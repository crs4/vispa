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

from __future__ import division

import unittest, copy, tempfile, os, itertools as it

import bl.tiget.pipeline.annotate as annotate
import bl.core.io.bed as bed


def write_bed(out_fn, annotation_info):
    with bed.open(out_fn, "w") as fo:
        for chrom, seq in annotation_info.interval_info.iteritems():
            r = {"chrom": chrom}
            for interval in seq:
                r["chromStart"], r["chromEnd"] = interval
                r["chromEnd"] += 1
                for name, strand in annotation_info.add_info[chrom][interval]:
                    r["name"] = name
                    r["score"] = 0  # no holes allowed
                    r["strand"] = strand
                    fo.writerow(r)


class SplitDisjointBase(object):

    def _check_split_disjoint(self, seq, subseqs):
        all_data = sum(subseqs, [])
        self.assertEqual(sorted(all_data), sorted(seq))
        for s in subseqs:
            for c in it.combinations(s, 2):
                t1, t2 = sorted(c)
                self.assertTrue(t1[1] <= t2[0])


class TestSplitDisjoint(unittest.TestCase, SplitDisjointBase):

    def setUp(self):
        self.data = [
          (0, 4), (5, 6), (8, 9),
          (1, 2), (3, 7),
          (1, 10),
          ]

    def runTest(self):
        subseqs = annotate.split_disjoint(self.data)
        self._check_split_disjoint(self.data, subseqs)


class TestFindClosestOnDisjoint(unittest.TestCase):

    def setUp(self):
        self.data = [(5, 7), (15, 25)]

    def _check_closest(self, pos, exp_d, exp_c_list):
        d, c_list = annotate.find_closest_on_disjoint(self.data, pos)
        self.assertEqual(d, exp_d)
        self.assertEqual(len(c_list), len(exp_c_list))
        for c, exp_c in it.izip(sorted(c_list), sorted(exp_c_list)):
            self.assertEqual(c, exp_c)

    def runTest(self):
        for i in xrange(1, 30):
            if i < 5:
                self._check_closest(i, 5-i, [self.data[0]])
            elif 5 <= i < 8:
                self._check_closest(i, 0, [self.data[0]])
            elif 8 <= i < 11:
                self._check_closest(i, i-7, [self.data[0]])
            elif i == 11:
                self._check_closest(i, 4, [self.data[0], self.data[1]])
            elif 12 <= i < 15:
                self._check_closest(i, 15-i, [self.data[1]])
            elif 15 <= i < 26:
                self._check_closest(i, 0, [self.data[1]])
            else:
                self._check_closest(i, i-25, [self.data[1]])


class TestAnnotationInfo(unittest.TestCase, SplitDisjointBase):

    def setUp(self):
        self.interval_info = {
            "chr1": [
                (1, 5), (8, 10),
                (1, 2), (4, 6),
                ],
            "chr2": [
                (5, 10),
                ],
            }
        self.add_info = {
            "chr1": {
                (1, 5): [("F11", "+")],
                (8, 10): [("F12", "+"), ("F13", "-")],
                (1, 2): [("F14", "+")],
                (4, 6): [("F15", "-")],
                },
            "chr2": {
                (5, 10): [("F21", "-")],
                }
            }
        self.annotation_info = annotate.AnnotationInfo(
            interval_info=copy.deepcopy(self.interval_info),
            add_info=copy.deepcopy(self.add_info),
            )
        fd, self.temp_fn = tempfile.mkstemp(prefix="bl_tiget_")
        os.close(fd)

    def tearDown(self):
        os.remove(self.temp_fn)

    def test_split_all(self):
        self.annotation_info.split_all()
        for seq, subseqs in it.izip(
            self.interval_info.itervalues(),
            self.annotation_info.interval_info.itervalues(),
            ):
            self._check_split_disjoint(seq, subseqs)

    def test_load_from_bed(self):
        write_bed(self.temp_fn, self.annotation_info)
        info = annotate.AnnotationInfo()
        info.load_from_bed(self.temp_fn)
        for a in "interval_info", "add_info":
            self.assertEqual(
                sorted(getattr(self.annotation_info, a)),
                sorted(getattr(info, a)),
                )

    def _check_closest(self, chrom, pos, exp_d, exp_closest):
        d, closest = self.annotation_info.find_closest(chrom, pos)
        self.assertEqual(d, exp_d)
        self.assertEqual(len(closest), len(exp_closest))
        for c, exp_c in it.izip(sorted(closest), sorted(exp_closest)):
            self.assertEqual(c, exp_c)

    def test_find_closest(self):
        self._check_closest("chr1", 0, 1, [(1, 5), (1, 2)])
        for pos in 1, 2:
            self._check_closest("chr1", pos, 0, [(1, 5), (1, 2)])
        self._check_closest("chr1", 3, 0, [(1, 5)])
        for pos in 4, 5:
            self._check_closest("chr1", pos, 0, [(1, 5), (4, 6)])
        self._check_closest("chr1", 6, 0, [(4, 6)])
        self._check_closest("chr1", 7, 1, [(4, 6), (8, 10)])
        for pos in 8, 9, 10:
            self._check_closest("chr1", pos, 0, [(8, 10)])
        self._check_closest("chr1", 11, 1, [(8, 10)])

    def _check_annotate(self, chrom, pos, exp_records):
        records = self.annotation_info.annotate(chrom, pos)
        self.assertEqual(len(records), len(exp_records))
        for r, exp_r in it.izip(sorted(records), sorted(exp_records)):
            self.assertEqual(len(r), len(exp_r))
            for k, v in r.iteritems():
                self.assertTrue(k in exp_r)
                exp_v = exp_r[k]
                if k == "integration":
                    self.assertAlmostEqual(v, exp_v)
                else:
                    self.assertEqual(v, exp_v)

    def test_annotate(self):
        self._check_annotate("chr1", 0, [
            {"chrom": "chr1", "pos": 0, "name": "F11", "start": 1, "end": 5,
             "strand": '+', "tss_d": 1, "rel_pos": -1, "integration": 0},
            {"chrom": "chr1", "pos": 0, "name": "F14", "start": 1, "end": 2,
             "strand": '+', "tss_d": 1, "rel_pos": -1, "integration": 0},
            ])
        for pos in 1, 2:
            I1, I2 = 100*(pos-1)/4, 100*(pos-1)
            self._check_annotate("chr1", pos, [
            {"chrom": "chr1", "pos": pos, "name": "F11", "start": 1, "end": 5,
             "strand": '+', "tss_d": pos-1, "rel_pos": 0, "integration": I1},
            {"chrom": "chr1", "pos": pos, "name": "F14", "start": 1, "end": 2,
             "strand": '+', "tss_d": pos-1, "rel_pos": 0, "integration": I2},
                ])
        self._check_annotate("chr1", 3, [
            {"chrom": "chr1", "pos": 3, "name": "F11", "start": 1, "end": 5,
             "strand": '+', "tss_d": 2, "rel_pos": 0, "integration": 50},
            ])
        for pos in 4, 5:
            I1, I2 = 100*(pos-1)/4, 100*(6-pos)/2
            self._check_annotate("chr1", pos, [
            {"chrom": "chr1", "pos": pos, "name": "F11", "start": 1, "end": 5,
             "strand": '+', "tss_d": pos-1, "rel_pos": 0, "integration": I1},
            {"chrom": "chr1", "pos": pos, "name": "F15", "start": 4, "end": 6,
             "strand": '-', "tss_d": 6-pos, "rel_pos": 0, "integration": I2},
                ])
        self._check_annotate("chr1", 6, [
            {"chrom": "chr1", "pos": 6, "name": "F15", "start": 4, "end": 6,
             "strand": '-', "tss_d": 0, "rel_pos": 0, "integration": 0},
            ])
        self._check_annotate("chr1", 7, [
            {"chrom": "chr1", "pos": 7, "name": "F15", "start": 4, "end": 6,
             "strand": '-', "tss_d": 1, "rel_pos": -1, "integration": 0},
            {"chrom": "chr1", "pos": 7, "name": "F12", "start": 8, "end": 10,
             "strand": '+', "tss_d": 1, "rel_pos": -1, "integration": 0},
            {"chrom": "chr1", "pos": 7, "name": "F13", "start": 8, "end": 10,
             "strand": '-', "tss_d": 3, "rel_pos": 1, "integration": 0},
            ])
        for pos in 8, 9, 10:
            I1, I2 = 100*(pos-8)/2, 100*(10-pos)/2
            self._check_annotate("chr1", pos, [
            {"chrom": "chr1", "pos": pos, "name": "F12", "start": 8, "end": 10,
             "strand": '+', "tss_d": pos-8, "rel_pos": 0, "integration": I1},
            {"chrom": "chr1", "pos": pos, "name": "F13", "start": 8, "end": 10,
             "strand": '-', "tss_d": 10-pos, "rel_pos": 0, "integration": I2},
                ])
        self._check_annotate("chr1", 11, [
            {"chrom": "chr1", "pos": 11, "name": "F12", "start": 8, "end": 10,
             "strand": '+', "tss_d": 3, "rel_pos": 1, "integration": 0},
            {"chrom": "chr1", "pos": 11, "name": "F13", "start": 8, "end": 10,
             "strand": '-', "tss_d": 1, "rel_pos": -1, "integration": 0},
            ])


def load_tests(loader, tests, pattern):
    test_cases = (
        TestSplitDisjoint,
        TestFindClosestOnDisjoint,
        TestAnnotationInfo,
        )
    suite = unittest.TestSuite()
    for tc in test_cases:
        suite.addTests(loader.loadTestsFromTestCase(tc))
    return suite


if __name__ == '__main__':
    suite = load_tests(unittest.defaultTestLoader, None, None)
    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(suite)
