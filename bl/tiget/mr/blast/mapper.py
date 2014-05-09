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
import os, logging, math
logging.basicConfig(level=logging.DEBUG)

import pydoop.pipes as pp
import pydoop.utils as pu
from bl.core.seq.engines.blastall_2_2_21 import Engine
from bl.core.seq.stats.karlin_altschul import BlastallLKCalculator
from bl.tiget.repeats import is_repeat
import al_type


LN2 = 0.69314718055994529


## Tabular output fields: Query id, Subject id, % identity, alignment
## length, mismatches, gap openings, q. start, q. end, s. start,
## s. end, e-value, bit score

## This is for the next step (computing repeated sequences)
## --------------------------------------------------------
## Homology percent = 100 * LA / LT
##   LA = alignment length = |query start - query end|
##   LT = target (query) length


class Mapper(pp.Mapper):
  """
  Maps query sequences to blastall hits, discarding those that do not
  meet the following criteria:

    query start <= n
    % identity >= m

  The output key is al_type.UNAMBIGUOUS, al_type.REPEAT or al_type.NO_HIT.

  If there are no hits after filtering, the value is the sequence tag;
  otherwise, one k/v pair is emitted for each hit, where values are
  the tabular blast hits (the sequence tag is the first field).

  @input-record: C{key} does not matter (LineRecordReader), C{value} =
  whole sequence as output by fasta2tab (<HEADER>\t<SEQUENCE>)

  @output-record: tabular blastall hit against the specified db (first
  field is removed from the hit and emitted as key, so that the
  tab-separated k/v pair is the original tab-separated output).

  @jobconf-param: C{bl.mr.log.level} logging level, specified as a
  logging module literal; defaults to 'WARNING'.

  @jobconf-param: C{bl.mr.seq.blastall.db.name} The BLAST database
  name (REQUIRED). A BLAST db is typically obtained by running the
  C{formatdb} command on one or more fasta files. The archive provided
  through C{mapred.cache.archives} MUST contain BLAST db files
  beginning with the db name (DB_NAME.nin, etc.).

  @jobconf-param: C{bl.mr.seq.formatdb.exe} Full path to the formatdb
  executable.

  @jobconf-param: C{bl.mr.seq.blastall.exe} Full path to the blastall
  executable.

  @jobconf-param: C{bl.mr.seq.blastall.program} The BLAST program
  to use ('blastn', 'blastp', etc.); defaults to 'blastn'.

  @jobconf-param: C{bl.mr.seq.blastall.evalue} upper threshold for the
  expectation value.

  @jobconf-param: C{bl.mr.seq.blastall.gap.cost} gap opening cost.

  @jobconf-param: C{bl.mr.seq.blastall.word.size} length of best perfect match.

  @jobconf-param: C{bl.mr.seq.blastall.filter} filter options for DUST or SEG.

  @jobconf-param: C{mapred.cache.archives} distributed cache entry
  (HDFS_PATH#LINK_NAME) for an archive containing the pre-formatted db
  files at the top level, i.e., no directories.

  @jobconf_param: C{bl.mr.seq.tiget.max.hits} for each query seq, use
  only the first N blast hits.

  @jobconf_param: C{bl.mr.seq.tiget.max.start}: discard blast hits
  with sequence start higher than this value.

  @jobconf-param: C{bl.mr.seq.tiget.min.identity.percent} discard
  blast hits with identity lower than this value

  @jobconf-param: C{mapred.create.symlink} must be set to 'yes'.
  """
  COUNTER_CLASS = "BLASTALL"

  def __get_log_conf(self, jc):
    pu.jc_configure(self, jc, 'bl.mr.log.level', 'log_level', 'WARNING')
    try:
      self.log_level = getattr(logging, self.log_level)
    except AttributeError:
      raise ValueError("Unsupported log level: %r" % self.log_level)

  def __get_blastall_conf(self, jc):
    pu.jc_configure(self, jc, 'bl.mr.seq.blastall.exe',
                    'blastall_exe', '/usr/bin/blastall')
    pu.jc_configure(self, jc, 'bl.mr.seq.blastall.program', 'program', 'blastn')
    pu.jc_configure(self, jc, 'bl.mr.seq.blastall.db.name', 'db_name')
    pu.jc_configure_float(self, jc, 'bl.mr.seq.blastall.evalue', 'evalue', 1.0)
    pu.jc_configure_int(self, jc, 'bl.mr.seq.blastall.gap.cost', 'gap_cost', 1)
    pu.jc_configure_int(self, jc, 'bl.mr.seq.blastall.word.size',
                        'word_size', 20)
    pu.jc_configure_bool(self, jc, 'bl.mr.seq.blastall.filter',
                        'filter', False)

  def __get_tiget_conf(self, jc):
    pu.jc_configure_int(self, jc, 'bl.mr.seq.tiget.max.hits', 'max_hits', 10)
    pu.jc_configure_int(self, jc, 'bl.mr.seq.tiget.max.start', 'max_start', 4)
    pu.jc_configure_float(self, jc, 'bl.mr.seq.tiget.min.identity.percent',
                          'min_identity', 95.0)
    pu.jc_configure_float(self, jc, 'bl.mr.seq.tiget.min.al2seq.percent',
                          'min_al2seq', 15.0)
    self.min_al2seq /= 100
    pu.jc_configure_float(self, jc, 'bl.mr.seq.tiget.min.score.diff',
                          'min_score_diff', 20.0)

  def __get_conf(self, jc):
    self.__get_log_conf(jc)  # log always comes first
    self.__get_blastall_conf(jc)
    self.__get_tiget_conf(jc)
    pu.jc_configure(self, jc, 'bl.mr.seq.formatdb.exe',
                    'formatdb_exe', '/usr/bin/formatdb')
    pu.jc_configure_bool(self, jc, 'bl.spawner.guardian', 'guardian', True)

  def __get_counters(self):
    self.hit_counter = self.ctx.getCounter(self.COUNTER_CLASS, "TOTAL_HITS")
    self.hom_rej_hit_counter = self.ctx.getCounter(self.COUNTER_CLASS,
                                                   "IDENTITY_REJECTED_HITS")
    self.start_rej_hit_counter = self.ctx.getCounter(self.COUNTER_CLASS,
                                                     "START_REJECTED_HITS")

  def __init__(self, ctx):
    super(Mapper, self).__init__(ctx)
    self.ctx = ctx
    jc = self.ctx.getJobConf()
    self.__get_conf(jc)
    self.__get_counters()
    self.logger = logging.getLogger("mapper")
    self.logger.setLevel(self.log_level)
    self.input_file = "temp.in"
    self.output_file = "temp.out"
    engine_logger = logging.getLogger("blastall")
    engine_logger.setLevel(self.log_level)
    self.engine = Engine(exe_file=self.blastall_exe, logger=engine_logger,
                         create_guardian=self.guardian)
    try:
      self.db_dir = jc.get("mapred.cache.archives").split(",")[0].split("#")[1]
    except IndexError:
      raise ValueError('bad format for "mapred.cache.archives"')
    self.opts = {
      "blastall.program": self.program,
      "blastall.database": os.path.join(self.db_dir, self.db_name),
      "blastall.out.tabular": True,
      "blastall.input.file": self.input_file,
      "blastall.output.file": self.output_file,
      "blastall.evalue": self.evalue,
      "blastall.gap.cost": self.gap_cost,
      "blastall.word.size": self.word_size,
      "blastall.filter": self.filter,
      }
    c = BlastallLKCalculator(self.formatdb_exe,
                             self.blastall_exe,
                             log_level=self.log_level,
                             engine_opts=self.opts)
    self.lambda_, kappa = c.calculate()
    self.lnK = math.log(kappa)

  def map(self, ctx):
    header, query_seq = ctx.getInputValue().rstrip().split("\t", 1)
    # TODO: use stdin/stdout instead
    self.__write_input(header, query_seq)
    self.engine.blastall(opts=self.opts)
    results = list(self.__filter_results(self.__read_output()))
    if not results:
      ctx.emit(str(al_type.NO_HIT), header)
    else:
      repeat = is_repeat(len(query_seq), results,
                         self.min_al2seq, self.min_score_diff)
      key = str(al_type.REPEAT) if repeat else str(al_type.UNAMBIGUOUS)
      for r in results:
        ctx.emit(key, "\t".join(r))
    
  def __filter_results(self, results_stream):
    for i, r in enumerate(results_stream):
      if i > self.max_hits:
        break
      self.ctx.incrementCounter(self.hit_counter, 1)
      identity = float(r[2])  # percentage
      query_start, query_end = map(int, r[6:8])
      if query_start > self.max_start:
        self.ctx.incrementCounter(self.start_rej_hit_counter, 1)
        continue
      if identity < self.min_identity:
        self.ctx.incrementCounter(self.hom_rej_hit_counter, 1)
        continue
      r[-1] = str(self.__bit2raw(float(r[-1])))
      yield r
    
  def __write_input(self, header, query_seq):
    f = open(self.input_file, "w")
    f.write(">%s\n%s\n" % (header, query_seq))
    f.close()

  def __read_output(self):
    f = open(self.output_file)
    for line in f:
      yield line.rstrip().split()
    f.close()

  def __bit2raw(self, bit_score):
    return (LN2*bit_score+self.lnK)/self.lambda_
