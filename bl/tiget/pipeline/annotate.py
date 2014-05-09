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
Annotate integration sites.
"""

from __future__ import division

import os, argparse, csv, random
from bisect import bisect

import bl.core.io.bed as bed


DEFAULT_OUTPUT_FN = "./annotate.out"
OUTPUT_FIELDS = [
    ("chrom", "%s"),          # query chromosome
    ("pos", "%d"),            # query position
    ("name", "%s"),           # feature name
    ("start", "%d"),          # feature start
    ("end", "%d"),            # feature end
    ("strand", "%s"),         # '+' or '-'
    ("tss_d", "%d"),          # distance from transcription start site
    ("rel_pos", "%d"),        # -1: upstream; 0: in_gene; 1: downstream
    ("integration", "%.2f"),  # integration % (0 if not in-gene)
    ]


def query_iterator(fn, delimiter="\t", skip_first=False):
    with open(fn) as f:
        reader = csv.reader(f, delimiter=delimiter)
        if skip_first:
            reader.next()
        for row in reader:
            yield row[0], int(row[1])


def split_disjoint(interval_seq):
    """
    Split a sequence of intervals into subseqs of disjoint intervals.

    Intervals are assumed to be **sorted** (begin, end) tuples.
    """
    interval_seq = sorted(interval_seq)
    subseqs = [[]]
    i = 0
    while interval_seq:
        interval = interval_seq.pop(i)
        subseqs[-1].append(interval)
        i = bisect(interval_seq, (interval[1],))
        if i >= len(interval_seq) and interval_seq:
            subseqs.append([])
            i = 0
    return subseqs


def find_closest_on_disjoint(seq, pos):
    """
    In a sequence of disjoint intervals, find the one(s) closest to ``pos``.

    Return a (d, list_of_closest_intervals) tuple where:

    * if pos lies between the endpoints (included), d is 0
    * otherwise, d is the distance between pos and the closest endpoint

    Note that, in most cases, the list will contain only one interval.
    The only case where it will contain two is the one in which pos is
    equidistant from the end of an interval and the beginning of the
    next one.
    """
    i = bisect(seq, (pos, float("inf")))
    if i == 0:
        return seq[0][0] - pos, [seq[0]]
    elif i == len(seq):
        return max(0, pos-seq[-1][1]), [seq[-1]]
    else:
        left, right = seq[i-1], seq[i]
        dl, dr = pos - left[1], right[0] - pos
        if dr < dl:
            return dr, [right]
        elif dr > dl:
            return max(0, dl), [left]
        else:
            return dl, [left, right]


class AnnotationInfo(object):

    def __init__(self, interval_info=None, add_info=None):
        """
        Store annotation info for genomic intervals.

        You should not instantiate this class directly. To get an
        AnnotationInfo object, call :fun:`get_annotation_info` on a
        bed file name.

        interval_info: a dict that maps chromosome tags to sequences of
          unique intervals. Before calling search methods, these sequences
          must be split into subsequences of disjoint intervals via
          :meth:`split_all`.  Example::

            interval_info = {
                'chrUn_gl000218': [
                    [(38785, 97453)],
                    [(46844, 55048), (62408, 97453)],
                    ],
                }

        add_info: a dict that maps chromosome tags to subdicts that,
          in turn, map intervals to sets of (name, strand) tuples.
          Example::

            add_info = {
                'chrUn_gl000218': {
                    (38785, 97453): {('LOC100233156', '-')},
                    (46844, 55048): {('LOC389834', '-')},
                    (62408, 97453): {('LOC100233156', '-')},
                    },
                }
        """
        self.interval_info = interval_info or {}
        self.add_info = add_info or {}
        self.__loaded = self.interval_info and self.add_info
        self.__split = False

    def load_from_bed(self, bed_fn):
        if self.__loaded:
            return
        with bed.open(bed_fn) as fi:
            for r in fi:
                interval = (r["chromStart"], r["chromEnd"]-1)
                subd = self.add_info.setdefault(r["chrom"], {})
                if interval not in subd:
                    self.interval_info.setdefault(r["chrom"], []).append(
                        interval
                        )
                subd.setdefault(interval, set()).add((r["name"], r["strand"]))
        self.__loaded = True

    def split_all(self):
        if self.__split:
            return
        for chrom, data in self.interval_info.iteritems():
            self.interval_info[chrom] = split_disjoint(data)
        self.__split = True

    def find_closest(self, chrom, pos):
        self.split_all()
        ranking = {}
        for seq in self.interval_info[chrom]:
            d, closest_list = find_closest_on_disjoint(seq, pos)
            ranking.setdefault(d, []).extend(closest_list)
        min_d = min(ranking)
        return min_d, ranking[min_d]

    def annotate(self, chrom, pos, multi=False):
        min_d, interval_list = self.find_closest(chrom, pos)
        results = []
        for interval in interval_list:
            left, right = interval
            length = right - left
            for name, strand in self.add_info[chrom][interval]:
                res = {"chrom": chrom, "pos": pos}
                res["start"], res["end"] = left, right
                res["name"], res["strand"] = name, strand
                tss = left if strand == '+' else right
                res["tss_d"] = abs(pos - tss)
                if min_d == 0:
                    res["rel_pos"] = 0  # in-gene
                    res["integration"] = 100 * res["tss_d"] / length
                else:
                    if ((strand == '+' and pos < left) or
                        (strand == '-' and pos > right)):
                        res["rel_pos"] = -1  # upstream
                    else:
                        res["rel_pos"] = 1  # downstream
                    res["integration"] = 0
                results.append(res)
        if not multi:
            results = random.sample(results, 1)
        return results


def get_annotation_info(bed_fn):
    annotation_info = AnnotationInfo()
    annotation_info.load_from_bed(bed_fn)
    annotation_info.split_all()
    return annotation_info


def make_parser():
    parser = argparse.ArgumentParser(
      description=__doc__.strip(),
      formatter_class=argparse.ArgumentDefaultsHelpFormatter,
      )
    parser.add_argument("input", metavar="INPUT",
                        help="merged integration sites in tabular format")
    parser.add_argument("annot", metavar="ANNOTATION",
                        help="annotation info in BED format")
    parser.add_argument("-o", "--output", metavar="STRING",
                        default=DEFAULT_OUTPUT_FN, help="output file")
    parser.add_argument("-d", "--delimiter", metavar="STRING",
                        help="field delimiter", default="\t")
    parser.add_argument("--skip-first", action="store_true",
                        help="skip first line (e.g., header)")
    parser.add_argument("--multi", action="store_true",
                        help="output multiple lines per IS")
    return parser


def main(argv):
    parser = make_parser()
    args = parser.parse_args(argv[1:])
    annotation_info = get_annotation_info(args.annot)
    unknown_tags = set()
    with open(args.output, "w") as fo:
        writer = csv.writer(
            fo, delimiter=args.delimiter, lineterminator=os.linesep
            )
        for chrom, pos in query_iterator(
            args.input, args.delimiter, args.skip_first
            ):
            if chrom in annotation_info.add_info:
                for r in annotation_info.annotate(chrom, pos, multi=args.multi):
                    writer.writerow([t % r[n] for n, t in OUTPUT_FIELDS])
            else:
                unknown_tags.add(chrom)
    for t in sorted(unknown_tags):
        print 'no annotation for "%s"' % t
