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
Collate blast_and_filter output into the unambiguous, no-hit and
repeat categories.
"""

import sys, os, optparse, fileinput
import bl.tiget.mr.blast.al_type as al_type


def collate_output(in_stream, out_prefix=""):
    if os.path.sep in out_prefix:
        os.makedirs(os.path.dirname(out_prefix))
    output_filenames = {
        al_type.UNAMBIGUOUS: "unambiguous",
        al_type.REPEAT: "repeat",
        al_type.NO_HIT: "no_hit"
        }
    output_files = dict((k, open(out_prefix+fn, "w"))
                        for k, fn in output_filenames.iteritems())
    blast_output_file = open(out_prefix+"all_hits.tsv", "w")
    old_seq_tag = None
    for line in in_stream:
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


class HelpFormatter(optparse.IndentedHelpFormatter):
    def format_description(self, description):
        return description + "\n" if description else ""


def make_parser():
    parser = optparse.OptionParser(
        usage="%prog [OPTIONS] FILE [FILE]...", formatter=HelpFormatter()
        )
    parser.set_description(__doc__.lstrip())
    parser.add_option("--out-prefix", metavar="STRING", default="",
                      help="prefix for output files ['%default']")
    return parser


def main(argv):
    parser = make_parser()
    opt, args = parser.parse_args(argv[1:])
    if len(args) < 1:
        parser.print_help()
        sys.exit(2)
    collate_output(fileinput.input(args), out_prefix=opt.out_prefix)
    fileinput.close()


if __name__ == "__main__":
    main(sys.argv)
