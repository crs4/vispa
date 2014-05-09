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
Merge redundant integration loci.
"""

import sys, argparse, csv
from collections import Counter


CHR_TAG = "s_id"
LOCUS_TAG = "s_start"
HEADER = [
  "q_id",
  CHR_TAG,
  "identity",
  "alignment_length",
  "mismatches",
  "gap_openings",
  "q_start",
  "q_end",
  LOCUS_TAG,
  "s_end",
  "e_value",
  "score",
  ]

HEADER_HELP = """
Header fields for the input file, separated by commas. Iff this is set
to 'auto', header fields will be looked for in the first row.
"""

def chrom_sorted(chrom_list):
  numbers, letters = [], []
  for c in chrom_list:
    try:
      c = int(c)
    except ValueError:
      letters.append(c)
    else:
      numbers.append(c)
  numbers.sort()
  letters.sort()
  return map(str, numbers) + letters


def first_pass(fn, delimiter="\t", header=None,
               chr_tag=CHR_TAG, locus_tag=LOCUS_TAG):
  """
  Split data by chromosome and merge identical loci.
  """
  unique_locs = {}
  with open(fn) as f:
    print "reading input from %s" % fn
    reader = csv.DictReader(f, fieldnames=header, delimiter=delimiter)
    for r in reader:
      try:
        chrom, locus = r[chr_tag], int(r[locus_tag])
      except KeyError as e:
        raise RuntimeError(
          "field %s not found in header %r" % (e, reader.fieldnames)
          )
      unique_locs.setdefault(chrom, Counter())[locus] += 1
  return unique_locs


def merge_single_chrom(sorted_data, window_size):
  """
  Merge redundant integration loci for a single chromosome.

  sorted_data: list of (locus, count) pairs, sorted by locus.
  """
  merged = []
  current = None
  for i, (locus, count) in enumerate(sorted_data):
    if i == 0:
      current = [locus, count]
      continue
    if (locus - sorted_data[i-1][0] <= window_size and
        locus - current[0] <= window_size):
      current[1] += count
    else:
      merged.append(current)
      current = [locus, count]
  merged.append(current)
  return merged


def make_parser():
  parser = argparse.ArgumentParser(
    description=__doc__.strip(),
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
  parser.add_argument("input", metavar="INPUT", help="BLAST tabular results")
  parser.add_argument("output", metavar="OUTPUT", help="merged int. sites")
  parser.add_argument("-d", "--delimiter", metavar="STRING",
                      help="field delimiter", default="\t")
  parser.add_argument("-w", "--window-size", metavar="INT", type=int,
                      help="window size", default=3)
  parser.add_argument("--dump-intermediate", action="store_true",
                      help="dump intermediate data (DEBUG)")
  parser.add_argument("--header", metavar="STRING",
                      help=HEADER_HELP, default=",".join(HEADER))
  parser.add_argument("--chr-tag", metavar="STRING",
                      help="chromosome header tag", default=CHR_TAG)
  parser.add_argument("--locus-tag", metavar="STRING",
                      help="integration locus header tag", default=LOCUS_TAG)
  return parser


def main(argv):
  parser = make_parser()
  args = parser.parse_args(argv[1:])
  if args.header == 'auto':
    args.header = None
  else:
    args.header = args.header.split(",")
  unique_locs = first_pass(
    args.input,
    delimiter=args.delimiter,
    header=args.header,
    chr_tag=args.chr_tag,
    locus_tag=args.locus_tag
    )
  if args.dump_intermediate:
    fdump = open("%s.intermediate" % args.output, "w")
  with open(args.output, "w") as fo:
    for chrom in chrom_sorted(unique_locs.keys()):
      data = unique_locs[chrom].items()
      data.sort()
      if args.dump_intermediate:
        for locus, count in data:
          fdump.write(args.delimiter.join(
            [chrom, str(locus), str(count)]
            )+"\n")
      merged = merge_single_chrom(data, args.window_size)
      for locus, count in merged:
        fo.write(args.delimiter.join([chrom, str(locus), str(count)])+"\n")
    print "wrote output to %s" % args.output
    if args.dump_intermediate:
      print " * DEBUG: wrote intermediate data to %s" % fdump.name


if __name__ == "__main__":
  main(sys.argv)
