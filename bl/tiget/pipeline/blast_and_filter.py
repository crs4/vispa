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
Reads sequences in tabular format, aligns them to a reference genome
with BLAST and filters results.
"""

import sys, os, optparse

import pydoop.hdfs as hdfs
import mr_common as mrc


BASE_MR_OPT = {
    "mapred.job.name": "blast_and_filter",
    "hadoop.pipes.java.recordreader": "true",
    "hadoop.pipes.java.recordwriter": "true",
    "mapred.reduce.tasks": "0",
    "mapred.create.symlink": "yes",
    }


def upload_db(db_path, logger):
    rfs, lfs = None, None
    upload = True
    try:
        rfs, lfs = hdfs.hdfs(), hdfs.hdfs("", 0)
        db_path_hdfs = os.path.basename(db_path)
        if rfs.exists(db_path_hdfs):
            with open(db_path) as f:
                local_sum = mrc.checksum(f)
            with rfs.open_file(db_path_hdfs) as f:
                hdfs_sum = mrc.checksum(f)
            if hdfs_sum == local_sum:
                logger.warn("BLAST db archive: using %s" % db_path_hdfs)
                upload = False
        if upload:
            logger.info("uploading BLAST db archive %s" % db_path)
            lfs.copy(db_path, rfs, db_path_hdfs)
    finally:
        for fs in rfs, lfs:
            if fs is not None:
                fs.close()
    return db_path_hdfs


def add_blast_optgroup(parser):
    optgroup = optparse.OptionGroup(parser, "BLAST Options")
    optgroup.add_option("-d", "--blast-db", metavar="FILE",
                        help="BLAST db archive (REQUIRED)")
    optgroup.add_option("--blastall", metavar="EXE",
                        default="/usr/bin/blastall",
                        help="blastall executable ['%default']")
    optgroup.add_option("--formatdb", metavar="EXE",
                        default="/usr/bin/formatdb",
                        help="Formatdb executable ['%default']")
    optgroup.add_option("-p", "--blast-prog", metavar="STRING",
                        default="blastn", help="BLAST program ['%default']")
    optgroup.add_option("-e", "--blast-evalue", type="float", metavar="FLOAT",
                        default=1., help="BLAST expectation value [%default]")
    optgroup.add_option("-g", "--blast-gap-cost", type="int", metavar="INT",
                        default=1, help="BLAST gap opening cost [%default]")
    optgroup.add_option("-w", "--blast-word-size", type="int", metavar="INT",
                        default=20, help="BLAST word size [%default]")
    optgroup.add_option("-F", "--blast-filters", action="store_true",
                        help="BLAST filters [False]")
    optgroup.add_option("--disable-guardian", action="store_true",
                        help="disable guardian for BLAST processes [False]")
    parser.add_option_group(optgroup)


def add_tiget_optgroup(parser):
    optgroup = optparse.OptionGroup(parser, "TIGET Options")
    optgroup.add_option("-K", "--max-hits", type="int", metavar="INT",
                        default=10, help="n. of BLAST hits to keep [%default]")
    optgroup.add_option("-S", "--max-start", type="int", metavar="INT",
                        default=4, help="max query start [%default]")
    optgroup.add_option("-H", "--min-homology", type="float", metavar="FLOAT",
                        default=95., help="min homology percent [%default]")
    optgroup.add_option("-v", "--min-al2seq", type="float", metavar="FLOAT",
                        default=15., help="min align-to-seq percent [%default]")
    optgroup.add_option("-z", "--min-score-diff", type="float", metavar="FLOAT",
                        default=20., help="min score difference [%default]")
    parser.add_option_group(optgroup)


def main(argv):
    parser = mrc.make_parser()
    add_blast_optgroup(parser)
    add_tiget_optgroup(parser)
    try:
        input_, output, opt, logger = mrc.parse_cl(parser, argv)
    except IndexError:
        sys.exit(2)
    if not opt.blast_db:
        sys.exit("--blast-db is required")
    mrc.config_pydoop(opt)
    db_path_hdfs = upload_db(opt.blast_db, logger)
    blast_db_name = os.path.basename(opt.blast_db).split(".", 1)[0]
    runner = mrc.PipesRunner(logger=logger)
    runner.set_input(input_, put=False)
    runner.set_output(output)
    runner.set_exe(mrc.build_launcher("bl.tiget.mr.blast"))
    mr_opt = BASE_MR_OPT.copy()
    mr_opt.update({
        "mapred.map.tasks": str(opt.mappers),
        "bl.mr.log.level": opt.log_level,
        "bl.spawner.guardian": 'false' if opt.disable_guardian else 'true',
        #-- BLAST --
        "mapred.cache.archives": "%s#%s" % (db_path_hdfs, blast_db_name),
        "bl.mr.seq.blastall.db.name": blast_db_name,
        "bl.mr.seq.blastall.exe": opt.blastall,
        "bl.mr.seq.formatdb.exe": opt.formatdb,
        "bl.mr.seq.blastall.program": opt.blast_prog,
        "bl.mr.seq.blastall.evalue": opt.blast_evalue,
        "bl.mr.seq.blastall.gap.cost": opt.blast_gap_cost,
        "bl.mr.seq.blastall.word.size": opt.blast_word_size,
        "bl.mr.seq.blastall.filter": 'true' if opt.blast_filters else 'false',
        #-- FILTER --
        "bl.mr.seq.tiget.max.hits": opt.max_hits,
        "bl.mr.seq.tiget.max.start": opt.max_start,
        "bl.mr.seq.tiget.min.identity.percent": opt.min_homology,
        "bl.mr.seq.tiget.min.al2seq.percent": opt.min_al2seq,
        "bl.mr.seq.tiget.min.score.diff": opt.min_score_diff,
      })
    runner.run(properties=mr_opt, mr_dump_file=opt.mr_dump_file)
    logger.info("all done")
    if opt.mr_dump_file is not sys.stderr:
        opt.mr_dump_file.close()


if __name__ == "__main__":
    main(sys.argv)
