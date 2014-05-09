# BEGIN_COPYRIGHT
# 
# Copyright (C) 2013-2014 CRS4.
# 
# This file is part of vispa.
# 
# vispa is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
# 
# vispa is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
# 
# You should have received a copy of the GNU General Public License along with
# vispa.  If not, see <http://www.gnu.org/licenses/>.
# 
# END_COPYRIGHT

"""
Runs the BLAST step from the TIGET workflow.

The first argument is the local path of input sequences in FASTA format;
the second one is the local path of the blast db archive (see below).

The blast db archive must:

 * be in one of the Hadoop-supported formats (zip, tar, tgz/tar.gz)
 * contain BLAST db files AT THE TOP LEVEL

BLAST db files, in turn, must begin with the db name, set through the
--blast-db option (if this is not set, the program uses the archive's
basename with any extensions removed).

Options can be provided through the 'tiget_blast.cfg' file, which is
looked for in the current directory. Example:

[DEFAULT]
log_level: DEBUG
blast_mappers: 100
blast_db: mm9

Note that option names are command line long option names with dashes
replaced by underscores. A command line option overrides its
corresponding configuration file option.
"""

import sys, os, logging, optparse, ConfigParser, uuid, hashlib
import subprocess as sp

import pydoop.hdfs as hdfs
import bl.tiget.mr.blast.al_type as al_type


LOG_FORMAT = '%(asctime)s|%(levelname)-8s|%(message)s'
LOG_DATEFMT = '%Y-%m-%d %H:%M:%S'
LOG_LEVELS = ['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL']

CONFIG_FILE = "tiget_blast.cfg"
BUFSIZE = 1024 * os.sysconf("SC_PAGE_SIZE")

DEFAULTS = {
  "log_level": "WARNING",
  "out_prefix": "",
  "disable_guardian": False,
  #--
  "hadoop_home": os.getenv("HADOOP_HOME", "/opt/hadoop"),
  "f2t_mappers": 1,
  "blast_mappers": 1,
  #--
  "blastall": "/usr/bin/blastall",
  "formatdb": "/usr/bin/formatdb",
  "blast_prog": "blastn",
  "blast_db": None,  # None = get from archive name
  "blast_evalue": 1.0,
  "blast_gap_cost": 1,
  "blast_word_size": 20,
  "blast_filters": False,
  #--
  "tiget_max_hits": 10,
  "tiget_max_start": 4,
  "tiget_homology_percent": 95.0,
  "tiget_min_al2seq_percent": 15.0,
  "tiget_min_score_diff": 20.0,
}

F2T_BASE_MR_OPT = {
  "mapred.job.name": "fasta2tab",
  "hadoop.pipes.java.recordreader": "false",
  "hadoop.pipes.java.recordwriter": "true",
  "mapred.map.tasks": str(DEFAULTS["f2t_mappers"]),
  "mapred.reduce.tasks": "0",
  "bl.mr.log.level": DEFAULTS["log_level"],
  }

BLAST_BASE_MR_OPT = {
  "mapred.job.name": "tiget_blast",
  "hadoop.pipes.java.recordreader": "true",
  "hadoop.pipes.java.recordwriter": "true",
  "mapred.create.symlink": "yes",
  "mapred.map.tasks": str(DEFAULTS["blast_mappers"]),
  "mapred.reduce.tasks": "0",
  "bl.mr.log.level": DEFAULTS["log_level"],
  }


class RandomStringGenerator(object):

  def __init__(self, prefix="mr_blast"):
    self.prefix = prefix

  def generate(self):
    return "%s%s" % (self.prefix, uuid.uuid4().hex)

STR_GENERATOR = RandomStringGenerator()

def rnd_str():
  return STR_GENERATOR.generate()


def checksum(f):
  md5 = hashlib.md5()
  while 1:
    s = f.read(BUFSIZE)
    if not s:
      break
    md5.update(s)
  return md5.hexdigest()


def write_launcher(outf, app_module):
  outf.write('#!/bin/bash\n')
  outf.write('""":"\n')
  for item in os.environ.iteritems():
    outf.write('export %s="%s"\n' % item)
  outf.write('exec python -u $0 $@\n')
  outf.write('":"""\n')
  outf.write('from %s import run_task\n' % app_module)
  outf.write('run_task()\n')


def build_d_options(opt_dict):
  d_options = []
  for name, value in opt_dict.iteritems():
    d_options.append("-D %s=%s" % (name, value))
  return " ".join(d_options)


def hadoop_pipes(pipes_opts, opt):
  cmd = "%s --config %s pipes %s" % (
    opt.hadoop, opt.hadoop_conf_dir, pipes_opts
    )
  logging.debug("cmd: %s" % cmd)
  p = sp.Popen(cmd, shell=True,
               stdout=opt.mr_dump_file, stderr=opt.mr_dump_file)
  return os.waitpid(p.pid, 0)[1]


class HelpFormatter(optparse.IndentedHelpFormatter):
  def format_description(self, description):
    return description + "\n" if description else ""


def make_parser():
  parser = optparse.OptionParser(
    usage="%prog [OPTIONS] INPUT DB_ARCHIVE",
    formatter=HelpFormatter(),
    )
  parser.set_description(__doc__.lstrip())
  add_hadoop_optgroup(parser)
  add_blast_optgroup(parser)
  add_tiget_optgroup(parser)
  parser.add_option("--log-file", metavar="FILE", help="log file [stderr]")
  parser.add_option("--log-level", metavar="STRING", choices=LOG_LEVELS,
                    help="log level ['INFO']")
  parser.add_option("--mr-dump-file", metavar="FILE",
                    help="MapReduce out/err dump file [stderr]")
  parser.add_option("--out-prefix", type="str", metavar="STRING",
                    help="prefix for output files ['%default']")
  parser.add_option("--disable-guardian", action="store_true",
                    help="disable guardian for BLAST processes [False]")
  config = ConfigParser.SafeConfigParser(DEFAULTS)
  config.read(CONFIG_FILE)
  defaults = config.defaults()
  try:
    defaults["blast_filters"] = config.getboolean("DEFAULT", "blast_filters")
  except ValueError:
    defaults["blast_filters"] = False
  try:
    defaults["disable_guardian"] = config.getboolean(
      "DEFAULT", "disable_guardian"
      )
  except ValueError:
    defaults["disable_guardian"] = False
  parser.set_defaults(**defaults)
  return parser


def add_hadoop_optgroup(parser):
  optgroup = optparse.OptionGroup(parser, "Hadoop Options")
  optgroup.add_option("--hadoop-home", type="str", metavar="STRING",
                      help="Hadoop home directory ['%default']")
  optgroup.add_option("--hadoop", type="str", metavar="STRING",
                      help="Hadoop executable ['%default']")
  optgroup.add_option("--hadoop-conf-dir", type="str", metavar="STRING",
                      help="Hadoop configuration directory ['%default']")
  optgroup.add_option("--f2t-mappers", type="int", metavar="INT",
                      help="n. mappers for fasta2tab [%default]")
  optgroup.add_option("--blast-mappers", type="int", metavar="INT",
                      help="n. mappers for blast [%default]")
  parser.add_option_group(optgroup)


def add_blast_optgroup(parser):
  optgroup = optparse.OptionGroup(parser, "BLAST Options")
  optgroup.add_option("--blastall", type="str", metavar="STRING",
                      help="Blastall executable ['%default']")
  optgroup.add_option("--formatdb", type="str", metavar="STRING",
                    help="Formatdb executable ['%default']")
  optgroup.add_option("-p", "--blast-prog", type="str", metavar="STRING",
                      help="BLAST program ['%default']")
  optgroup.add_option("-d", "--blast-db", type="str", metavar="STRING",
                      help="BLAST database [get from archive name]")
  optgroup.add_option("-e", "--blast-evalue", type="float", metavar="FLOAT",
                      help="BLAST expectation value [%default]")
  optgroup.add_option("-g", "--blast-gap-cost", type="int", metavar="INT",
                      help="BLAST gap opening cost [%default]")
  optgroup.add_option("-w", "--blast-word-size", type="int", metavar="INT",
                      help="BLAST word size [%default]")
  optgroup.add_option("-F", "--blast-filters", action="store_true",
                      help="BLAST filters [False]")
  parser.add_option_group(optgroup)


def add_tiget_optgroup(parser):
  optgroup = optparse.OptionGroup(parser, "TIGET Options")
  optgroup.add_option("-K", "--tiget-max-hits", type="int", metavar="INT",
                      help="number of BLAST hits to keep [%default]")
  optgroup.add_option("-S", "--tiget-max-start", type="int", metavar="INT",
                      help="upper threshold for query start [%default]")
  optgroup.add_option("-H", "--tiget-homology-percent", type="float",
                      metavar="FLOAT",
                      help="lower threshold for homology [%default]")
  optgroup.add_option("-v", "--tiget-min-al2seq-percent", type="float",
                      metavar="FLOAT",
                      help="min O1-O2 value for unique hits [%default]")
  optgroup.add_option("-z", "--tiget-min-score-diff", type="float",
                      metavar="FLOAT",
                      help="min S1-S2 value for unique hits [%default]")
  parser.add_option_group(optgroup)


def update_f2t_options(mr_opt, opt):
  mr_opt["mapred.map.tasks"] = opt.f2t_mappers
  mr_opt["bl.mr.log.level"] = opt.log_level_str


def update_blast_options(mr_opt, opt):
  mr_opt["mapred.map.tasks"] = opt.blast_mappers
  mr_opt["bl.mr.seq.blastall.exe"] = opt.blastall
  mr_opt["bl.mr.seq.formatdb.exe"] = opt.formatdb
  mr_opt["bl.mr.seq.blastall.program"] = opt.blast_prog
  mr_opt["bl.mr.seq.blastall.db.name"] = opt.blast_db
  mr_opt["bl.mr.seq.blastall.evalue"] = opt.blast_evalue
  mr_opt["bl.mr.seq.blastall.gap.cost"] = opt.blast_gap_cost
  mr_opt["bl.mr.seq.blastall.word.size"] = opt.blast_word_size
  mr_opt["bl.mr.seq.blastall.filter"] = 'true' if opt.blast_filters else 'false'
  mr_opt["bl.spawner.guardian"] = 'false' if opt.disable_guardian else 'true'
  mr_opt["bl.mr.log.level"] = opt.log_level_str


def update_tiget_options(mr_opt, opt):
  mr_opt["bl.mr.seq.tiget.max.hits"] = opt.tiget_max_hits
  mr_opt["bl.mr.seq.tiget.max.start"] = opt.tiget_max_start
  mr_opt["bl.mr.seq.tiget.min.identity.percent"] = opt.tiget_homology_percent
  mr_opt["bl.mr.seq.tiget.min.al2seq.percent"] = opt.tiget_min_al2seq_percent
  mr_opt["bl.mr.seq.tiget.min.score.diff"] = opt.tiget_min_score_diff


class Runner(object):

  def __init__(self, fs, lfs, logger):
    self.fs = fs
    self.lfs = lfs
    self.logger = logger
  
  def upload_archive(self, db_archive):
    db_archive_hdfs = os.path.basename(db_archive)  
    if self.fs.exists(db_archive_hdfs):
      with open(db_archive) as f:
        local_sum = checksum(f)
      with self.fs.open_file(db_archive_hdfs) as f:
        hdfs_sum = checksum(f)
      if hdfs_sum == local_sum:
        self.logger.warn(
          "using hdfs-cached db, remove it manually to replace it"
          )
        return db_archive_hdfs
    self.logger.info("uploading blast db archive")
    self.lfs.copy(db_archive, self.fs, db_archive_hdfs)
    return db_archive_hdfs
  
  def run_f2t(self, input_local, opt):
    input_hdfs, output_hdfs, f2t_launcher_hdfs = [rnd_str() for _ in xrange(3)]
    self.lfs.copy(input_local, self.fs, input_hdfs)
    with self.fs.open_file(f2t_launcher_hdfs, "w") as outf:
      write_launcher(outf, "bl.core.seq.mr.fasta2tab")
    mr_opt = {}
    mr_opt.update(F2T_BASE_MR_OPT)
    update_f2t_options(mr_opt, opt)
    d_options = build_d_options(mr_opt)
    self.logger.info("running fasta2tab, launcher='%s'" % f2t_launcher_hdfs)
    hadoop_pipes("%s -program %s -input %s -output %s" % (
      d_options, f2t_launcher_hdfs, input_hdfs, output_hdfs
      ), opt)
    return output_hdfs
  
  def run_blast(self, input_hdfs, db_archive_hdfs, opt):
    output_hdfs, blast_launcher_hdfs = [rnd_str() for _ in xrange(2)]
    with self.fs.open_file(blast_launcher_hdfs, "w") as outf:
      write_launcher(outf, "bl.tiget.mr.blast")
    mr_opt = {}
    mr_opt.update(BLAST_BASE_MR_OPT)
    mr_opt["mapred.cache.archives"] = "%s#%s" % (db_archive_hdfs, opt.blast_db)
    update_blast_options(mr_opt, opt)
    update_tiget_options(mr_opt, opt)
    d_options = build_d_options(mr_opt)
    self.logger.info("running blastall, launcher='%s'" % blast_launcher_hdfs)
    hadoop_pipes("%s -program %s -input %s -output %s" % (
      d_options, blast_launcher_hdfs, input_hdfs, output_hdfs
      ), opt)
    return output_hdfs
  
  def collect_output(self, output_hdfs, opt):
    ls = [r['name'] for r in self.fs.list_directory(output_hdfs)
          if r['name'].rsplit("/", 1)[1].startswith('part')]
    if os.path.sep in opt.out_prefix:
      os.makedirs(os.path.dirname(opt.out_prefix))
    output_filenames = {
      al_type.UNAMBIGUOUS: "unambiguous",
      al_type.REPEAT: "repeat",
      al_type.NO_HIT: "no_hit"
      }
    output_files = dict((k, open(opt.out_prefix+fn, "w"))
                        for k, fn in output_filenames.iteritems())
    blast_output_file = open(opt.out_prefix+"all_hits.tsv", "w")
    for i, path in enumerate(ls):
      self.logger.info("processing mapreduce output file %d/%d" %
                       (i+1, len(ls)))
      with self.fs.open_file(path) as f:
        old_seq_tag = None
        for line in f:
          line = line.strip()
          if not line:
            continue
          code, payload = line.split("\t", 1)
          code = int(code)
          outf = output_files[code]
          if code == al_type.NO_HIT:
            outf.write("%s\n" % payload)
            continue
          blast_output_file.write("%s\n" % payload)
          seq_tag, more_fields = payload.split("\t", 1)
          if seq_tag != old_seq_tag:
            outf.write("%s" % seq_tag)
            if code == al_type.UNAMBIGUOUS:
              outf.write("\t%s" % more_fields)
            outf.write("\n")
            old_seq_tag = seq_tag
    for f in output_files.itervalues():
      f.close()
    blast_output_file.close()


def main(argv):

  parser = make_parser()
  opt, args = parser.parse_args()
  try:
    input_fasta = args[0]
    db_archive = args[1]
  except IndexError:
    parser.print_help()
    sys.exit(2)

  STR_GENERATOR.prefix = os.path.basename(input_fasta)

  logger = logging.getLogger()
  for h in logger.handlers:
    logger.removeHandler(h)
  opt.log_level_str = opt.log_level
  opt.log_level = getattr(logging, opt.log_level)
  kwargs = {'format': LOG_FORMAT,
            'datefmt': LOG_DATEFMT,
            'level': opt.log_level}
  if opt.log_file:
    kwargs['filename'] = opt.log_file
  logging.basicConfig(**kwargs)

  logger.debug("cli args: %r" % (args,))
  logger.debug("cli opts: %s" % opt)

  if opt.mr_dump_file:
    opt.mr_dump_file = open(opt.mr_dump_file, "w")
  else:
    opt.mr_dump_file = sys.stderr
  
  if not opt.blast_db:
    opt.blast_db = os.path.basename(db_archive).split(".", 1)[0]
    logger.info("--blast-db not provided: setting to %r" % opt.blast_db)
  
  os.environ["HADOOP_HOME"] = opt.hadoop_home
  if not opt.hadoop:
    opt.hadoop = os.path.join(opt.hadoop_home, "bin/hadoop")
  if not opt.hadoop_conf_dir:
    opt.hadoop_conf_dir = os.path.join(opt.hadoop_home, "conf")
  os.environ["HADOOP_CONF_DIR"] = opt.hadoop_conf_dir
  hdfs.reset()

  fs = hdfs.hdfs()
  logger.debug("hdfs params: host=%s, port=%d" % (fs.host, fs.port))
  lfs = hdfs.hdfs("", 0)
  runner = Runner(fs, lfs, logger)

  try:
    db_archive_hdfs = runner.upload_archive(db_archive)
    blast_input_hdfs = runner.run_f2t(input_fasta, opt)
    blast_output_hdfs = runner.run_blast(blast_input_hdfs, db_archive_hdfs,
                                         opt)
    runner.collect_output(blast_output_hdfs, opt)
    logger.info("all done")
  finally:
    lfs.close()
    fs.close()
    if opt.mr_dump_file is not sys.stderr:
      opt.mr_dump_file.close()


if __name__ == "__main__":
  main(sys.argv)
