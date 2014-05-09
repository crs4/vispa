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
Galaxy wrapper for the demux tool.
"""

import sys, os, cStringIO
import bl.tiget.pipeline.demux as demux


def run_demux(argv):
    stdout = sys.stdout
    sink = cStringIO.StringIO()
    sys.stdout = sink
    demux.main(argv)
    sys.stdout = stdout
    return sink.getvalue()


def write_output(output):
    row_template = '<tr><td>%s</td><td>%s</td><td>%s</td></tr>\n'
    print '<table border="1">\n'
    print '<tr><th>Barcode</th><th>Match Count</th><th>Output File</th></tr>\n'
    for line in output.splitlines():
        bc, count, out_fn = line.split()
        basename = os.path.basename(out_fn)
        link = '<a href="%s">%s</a>' % (basename, basename)
        print row_template % (bc, count, link)
    print '</table>\n'


def main(argv):
    output = run_demux(argv)
    write_output(output)


if __name__ == '__main__':
    main(sys.argv)
