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

# pylint: disable=C0103

"""
Remove LTR and LC from LAM-PCR output sequences.
"""

import sys, os, argparse

from Bio.Seq import Seq
from BlastPython import fasta_stream_from_stream
from BlastPython import seq_factory_from_fasta
from BlastPython import seq_factory_from_str
from BlastPython import blast_seq_stream
from BlastPython import blast_result_stream
from BlastPython import blaster
from BlastPython import blast_filter
from ncbi_toolkit import *
import ncbi_toolkit

STRAND = strand.both  # NCBI bl2seq default
LTR = "TGTGTGACTCTGGTAACTAGAGATCCCTCAGACCCTTTTAGTCAGTGTGGAAAATCTCTAGCA"
LINKER = "AGTGGCACAGCAGTTAGG"


class PrepareInputFile:

    def __init__(self, project, database, user, tag, infileName, outfileName,
                 trimmed):
        self.infileName = infileName
        self.infile = open(infileName, 'r')
        self.project = project
        self.database = database
        self.user = user
        self.tag = tag
        self.seqCounter = 0
        self.outfileName = outfileName
        self.trimmedSeqFile = open(trimmed, 'r')
        self._acquireTrimmed()
        self._scrollFastaFile()
        self.infile.close()
        self.trimmedSeqFile.close()

    def _acquireTrimmed(self):
        self.trimmed = {}
        try:
            while True:
                header = self.trimmedSeqFile.next().strip().replace(">", "")
                sequence = self.trimmedSeqFile.next().strip()
                self.trimmed[header] = sequence
        except StopIteration:
            return

    def _scrollFastaFile(self):
        self._setOutLine(["project", "database", "user", "filename", "tag",
                          "title", "sequence", "trimmed_sequence"])
        tmpFasta = []
        try:
            while True:
                line = self.infile.next().strip()
                if ">" in line:
                    if tmpFasta == []:
                        tmpFasta.append(line)
                    else:
                        self._setSequenceInTab(tmpFasta)
                        tmpFasta = [line]
                else:
                    tmpFasta.append(line)
        except StopIteration:
            self._setSequenceInTab(tmpFasta)
            return

    def _setSequenceInTab(self, tmpFasta):
        title = tmpFasta[0].split(None, 1)[0].replace(">", "")
        sequence = "".join(tmpFasta[1:]).strip()
        outLine = [self.project, self.database, self.user,
                   os.path.basename(self.infileName), self.tag,
                   title, sequence]
        if title in self.trimmed:
            outLine.append(self.trimmed[title])
        else: outLine.append("not_passed")
        self._setOutLine(outLine)
        self.seqCounter += 1
        return

    def _setOutLine(self, outLine):
        outFile = open(self.outfileName, 'a')
        outFile.write('\t'.join(outLine) + '\n')
        outFile.close()
        return


class seq_factory_from_str(ncbi_toolkit.blast_sseq_loc_from_str):

    def __init__(self, strand) :
        ncbi_toolkit.blast_sseq_loc_from_str.__init__(self)
        self.strand   = strand
        self.seq_gi   = 0

    def make(self, s):
        self.seq_gi += 1
        return super(seq_factory_from_str, self).make(
            s, False, self.seq_gi, 'fake title', self.strand, 0, 0
            )


class similarity_filter(blast_filter):

    len_query=None

    def __init__(self, min_sim_fraction,len_query, in_stream):
        blast_filter.__init__(self, in_stream)
        self.min_sim_fraction = min_sim_fraction
        self.len_query=len_query

    def is_acceptable(self, r):
        blast_filter.is_acceptable(self, r)
        subject, r = r
        hit_list = r[0]
        if hit_list is None:
            return False
        hsp_list = hit_list[0]
        hsp = hsp_list[0]
        off=self.len_query-hsp.query[2]
        sim_frac = hsp.num_ident/(float(self.len_query)-off*(off==3))
        if sim_frac < self.min_sim_fraction:
            pass
        return sim_frac > self.min_sim_fraction


class mapped_stream(object):

    def __init__(self, mapper, in_stream):
        self.in_stream = in_stream
        self.mapper    = mapper

    def __iter__(self):
        return self

    def next(self):
        r = self.in_stream.next()
        r = self.mapper.convert(r)
        if r:
            return r
        else:
            return self.next()


class cut_ltr_mapper_output_fasta(object):

    def convert(self, r):
        subject, r  = r
        hit_list = r[0]
        if hit_list is not None:
            hsp_list = hit_list[0]
            for hsp in hsp_list:
                my_strand= hsp.query[0]
                start= hsp.subject[1]
                end = hsp.subject[2]
                my_seq=subject.get_sequence()
                if my_strand == 1:
                    new_seq = my_seq[end:]
                elif my_strand == -1:
                    rev_seq = my_seq[:start]
                    rev_seq = Seq(rev_seq)
                    rev_seq = rev_seq.reverse_complement()
                    new_seq = str(rev_seq)
                if len(new_seq) > 10:
                    return '>%s\n%s' % (subject.id, new_seq)
                else:
                    return None
                    break


class cut_linker_mapper_output_fasta(object):

    def convert(self, r):
        subject, r = r
        hit_list = r[0]
        if hit_list is not None:
            hsp_list = hit_list[0]
            for hsp in hsp_list:
                my_strand= hsp.query[0]
                start= hsp.subject[1]
                end = hsp.subject[2]
                my_seq=subject.get_sequence()
                new_seq = my_seq[:start]
                return subject.id, new_seq
        else:
            seq = subject.get_sequence()
            return subject.id, seq


class ExecuteTrimming:

    def __init__(self, ltr, linker):
        self.ltr_fasta = ">ltr\n%s\n" % ltr
        self.linker_fasta = ">linker\n%s\n" % linker

    def trimming(self, raw_file, outname):
        self.raw_file = raw_file
        self.outname = outname
        default_strand = STRAND
        seq_factory_fasta = seq_factory_from_fasta(default_strand)
        ltr_seq = seq_factory_fasta.make(self.ltr_fasta)
        kwds = {'Program' : ncbi_toolkit.EProgram.eBlastn}
        ltr_blaster = blaster(ltr_seq, **kwds)
        linker_seq = seq_factory_fasta.make(self.linker_fasta)
        kwds = {'Program' : ncbi_toolkit.EProgram.eBlastn}
        linker_blaster = blaster(linker_seq, **kwds)
        fin = open(self.raw_file)
        fs = fasta_stream_from_stream(fin)
        ss = blast_seq_stream(seq_factory_fasta, fs)
        bs = blast_result_stream(ltr_blaster, ss)
        fbs = similarity_filter(0.89, len(ltr_seq.get_sequence()), bs)
        ltr_strs = mapped_stream(cut_ltr_mapper_output_fasta(), fbs)
        ltr_seqs = blast_seq_stream(seq_factory_fasta, ltr_strs)
        lnk_res = blast_result_stream(linker_blaster, ltr_seqs)
        lnk_strs = mapped_stream(cut_linker_mapper_output_fasta(), lnk_res)
        for x in lnk_strs:
            if len(x[1]) > 19:
                line = "\n".join(
                    [">"+x[0].replace("|", "_").replace("lcl_", ""),x[1]]
                    )
                self._setOut(line)
            else:
                pass

    def _setOut(self, line):
        out = open(self.outname, 'a')
        out.write(line + '\n')
        out.close()
        return


class TrimmingWrapper:

    def __init__(self, infilename, ltr, linker, outname):
        self.infilename = infilename
        if ltr == "def":
            self.ltr = LTR
        else: self.ltr = ltr
        if linker == "def":
            self.linker = LINKER
        else: self.linker = linker
        self.outname = outname
        self._runTrimm()

    def _runTrimm(self):
        E = ExecuteTrimming(self.ltr, self.linker)
        E.trimming(self.infilename, self.outname)
        return


def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', dest='infileName', action='store',
                        required=True, help='Fasta input file name')
    parser.add_argument('-ltr', dest='ltrSeq', action='store', required=True,
                        help='LTR sequence - def for default')
    parser.add_argument('-lk', dest='linkerSeq', action='store', required=True,
                        help='Linker sequence - def for default')
    parser.add_argument('-tr', dest='trimmed', action='store', required=True,
                        help='trimmed sequences')
    parser.add_argument('-usedb', dest='usedb', action='store',
                        required=True, help='project name')
    parser.add_argument('-p', dest='project', action='store',
                        default='tiget_pj', required=False,
                        help='project name')
    parser.add_argument('-d', dest='database', action='store',
                        default='tiget_db', required=False,
                        help='database name')
    parser.add_argument('-u', dest='user', action='store',
                        default='tiget_user', required=False, help='user name')
    parser.add_argument('-t', dest='tag', action='store', default='tiget_tag',
                        required=False, help='tag name')
    parser.add_argument('-o', dest='outfileName', action='store',
                        required=False, help='False outfile name')
    args = parser.parse_args()
    TW = TrimmingWrapper(args.infileName, args.ltrSeq, args.linkerSeq,
                         args.trimmed)
    if args.usedb == "advanced":
        P = PrepareInputFile(args.project, args.database, args.user, args.tag,
                             args.infileName, args.outfileName, args.trimmed)
    return


if __name__ == "__main__":
    sys.exit(main(sys.argv))
