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

import sys, cStringIO, unittest, tempfile, os, shutil

from bl.core.seq.io.fasta import SimpleFastaReader as FastaReader
import bl.tiget.pipeline.demux as demux


def run_and_get_stdout(func, *args, **kwargs):
    stdout = sys.stdout
    sink = cStringIO.StringIO()
    sys.stdout = sink
    func(*args, **kwargs)
    sys.stdout = stdout
    return sink.getvalue()


class TestDemuxFasta(unittest.TestCase):

    def setUp(self):
        self.barcodes = ['TTTTAAAA', 'GGGGCCCC']
        self.data = [
            ('0', 'TTTTAAAACATCGTAGTGACGTGCATGTCGACATGTGCTATGTGTCAGA'),
            ('1', 'TTTTAAATCATCGTAGTGACGTGCATGTCGACATGTGCTATGTGTCAGA'),
            ('2', 'GGGGCCCCCATCGTGATGCTGATCACTGAAGCTAGCTAGTCGTAGTCGT'),
            ('3', 'TTTTAAAAACTGATCGATGCTAGTCGATGTGCTAGTGCTAGTGTGCTAG'),
            ('4', 'ACTGATCGATGCTAGTCGATGTGCTAGTGCTAGTGTGCTAGTTTTAAAA'),
            ('5', 'ACTGATCGATGCTAGTTTTAAAATCGATGTGCTAGTGCTAGTGTGCTAG'),
            ]
        self.workdir = tempfile.mkdtemp(prefix='bl_tiget_')
        self.infn = os.path.join(self.workdir, 'input')
        self.bcfn = os.path.join(self.workdir, 'barcodes')
        with open(self.infn, 'w') as f:
            for t in self.data:
                f.write('>%s\n%s\n' % t)
        with open(self.bcfn, 'w') as f:
            for b in self.barcodes:
                f.write('%s\n' % b)

    def tearDown(self):
        shutil.rmtree(self.workdir)

    def runTest(self):
        print
        expected_content = {
            self.barcodes[0]: [self.data[i] for i in (0, 3)],
            self.barcodes[1]: [self.data[2]],
            }
        stdout = run_and_get_stdout(
            demux.main,
            ['test', self.infn, self.bcfn, '--outdir', self.workdir]
            )
        for bc in self.barcodes:
            fn = os.path.join(self.workdir, bc) + '.fa'
            self.assertTrue(os.path.isfile(fn))
            with open(fn) as f:
                self.assertEqual(list(FastaReader(f)), expected_content[bc])
        stdout = [_.split() for _ in stdout.splitlines()]
        match_count = dict((_[0], int(_[1])) for _ in stdout)
        for i in 0, 1:
            bc = self.barcodes[i]
            self.assertEqual(match_count[bc], len(expected_content[bc]))


def load_tests(loader, tests, pattern):
    test_cases = (
        TestDemuxFasta,
        )
    suite = unittest.TestSuite()
    for tc in test_cases:
        suite.addTests(loader.loadTestsFromTestCase(tc))
    return suite


if __name__ == '__main__':
    SUITE = load_tests(unittest.defaultTestLoader, None, None)
    RUNNER = unittest.TextTestRunner(verbosity=2)
    RUNNER.run(SUITE)
