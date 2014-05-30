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
Common utilities for MapReduce apps.
"""
import os, optparse
from bl.core.utils import LOG_LEVELS


class HelpFormatter(optparse.IndentedHelpFormatter):
    def format_description(self, description):
        return description + "\n" if description else ""


def build_launcher(app_module):
    lines = [
        '#!/bin/bash\n',
        '""":"\n',
        ]
    for item in os.environ.iteritems():
        lines.append('export %s="%s"\n' % item)
    lines.extend([
        'exec python -u $0 $@\n',
        '":"""\n',
        'from %s import run_task\n' % app_module,
        'run_task()\n',
        ])
    return ''.join(lines)


def add_hadoop_optgroup(parser):
    optgroup = optparse.OptionGroup(parser, "Hadoop Options")
    optgroup.add_option("--hadoop-home", metavar="DIR",
                        default=os.getenv("HADOOP_HOME", "/opt/hadoop"),
                        help="Hadoop home directory ['%default']")
    optgroup.add_option("--hadoop", metavar="FILE",
                        help="Hadoop executable [${HADOOP_HOME}/bin/hadoop]")
    optgroup.add_option("--hadoop-conf-dir", metavar="DIR",
                        help="Hadoop config directory [${HADOOP_HOME}/conf]")
    optgroup.add_option("--mappers", type="int", metavar="INT", default=1,
                        help="n. mappers [%default]")
    parser.add_option_group(optgroup)


def make_parser():
    parser = optparse.OptionParser(
        usage="%prog [OPTIONS] INPUT OUTPUT",
        formatter=HelpFormatter(),
        )
    parser.set_description(__doc__.lstrip())
    add_hadoop_optgroup(parser)
    parser.add_option("--log-file", metavar="FILE", help="log file [stderr]")
    parser.add_option("--log-level", metavar="STRING", choices=LOG_LEVELS,
                      default="INFO", help="log level ['%default']")
    parser.add_option("--mr-dump-file", metavar="FILE",
                      help="MapReduce out/err dump file [stderr]")
    return parser
