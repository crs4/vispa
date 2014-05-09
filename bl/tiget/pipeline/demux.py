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
Demultiplex sequence data according to a set of barcodes.

Input format: FASTA; mismatches allowed: 0.
"""
import os, argparse, errno
from collections import Counter

from bl.core.seq.io.fasta import SimpleFastaReader as FastaReader


DEFAULT_OUTPUT_DIR = os.getcwd()


def get_barcodes(fn):
    with open(fn) as f:
        return set(line.strip() for line in f)


def demux_fasta(f, barcodes, outdir=DEFAULT_OUTPUT_DIR):
    out_map = dict.fromkeys(barcodes)
    seq_count = Counter()
    try:
        for bc in out_map:
            out_map[bc] = open(os.path.join(outdir, bc) + '.fa', 'w')
        reader = FastaReader(f)
        for header, seq in reader:
            for bc in barcodes:
                if seq.startswith(bc):
                    out_map[bc].write(">%s\n%s\n" % (header, seq))
                    seq_count[bc] += 1
    finally:
        for outf in out_map.itervalues():
            if outf is not None:
                outf.close()
    return seq_count, out_map


def make_parser():
    parser = argparse.ArgumentParser(
      description=__doc__.strip(),
      formatter_class=argparse.ArgumentDefaultsHelpFormatter,
      )
    parser.add_argument("input", metavar="INPUT", help="FASTA file")
    parser.add_argument("barcodes", metavar="BARCODES",
                        help="barcode file (one barcode per line)")
    parser.add_argument("--outdir", metavar="DIR_PATH",
                        default=DEFAULT_OUTPUT_DIR, help="output dir")
    return parser


def main(argv):
    parser = make_parser()
    args = parser.parse_args(argv[1:])
    try:
        os.makedirs(args.outdir)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise
    barcodes = get_barcodes(args.barcodes)
    with open(args.input) as f:
        seq_count, out_map = demux_fasta(f, barcodes, outdir=args.outdir)
    for bc, count in seq_count.most_common():
        print "%s\t%d\t%s" % (bc, count, out_map[bc].name)
