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
Convert FASTA data to a tabular format on Hadoop.

Input:
>s1
ACAC
GTGT
>s2
CACA

Output:
s1	ACACGTGT
s2	CACA
"""

import sys
import mr_common as mrc


BASE_MR_OPT = {
    "mapred.job.name": "fasta2tab",
    "hadoop.pipes.java.recordreader": "false",
    "hadoop.pipes.java.recordwriter": "true",
    "mapred.reduce.tasks": "0",
    }
PREFIX = "fasta2tab_"


def main(argv):
    parser = mrc.make_parser()
    try:
        input_, output, opt, logger = mrc.parse_cl(parser, argv)
    except IndexError:
        sys.exit(2)
    mrc.config_pydoop(opt)
    runner = mrc.PipesRunner(logger=logger)
    runner.set_input(input_, put=False)
    runner.set_output(output)
    runner.set_exe(mrc.build_launcher("bl.core.seq.mr.fasta2tab"))
    mr_opt = BASE_MR_OPT.copy()
    mr_opt.update({
        "mapred.map.tasks": str(opt.mappers),
        "bl.mr.log.level": opt.log_level,
        })
    runner.run(properties=mr_opt, mr_dump_file=opt.mr_dump_file)
    logger.info("all done")
    if opt.mr_dump_file is not sys.stderr:
        opt.mr_dump_file.close()


if __name__ == "__main__":
    main(sys.argv)
