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

def al2seq(result, seq_len):
  query_start = int(result[6])
  query_end = int(result[7])
  return abs(query_start - query_end) / float(seq_len)


score = lambda result: float(result[-1])


def is_repeat(seq_len, blast_results, min_al2seq=.15, min_score_diff=20.):
  """
  Mark a repeat according to tiget rules (FIXME: summarize).

  NOTE: sorts ``blast_results`` by decreasing score.

  :type seq_len: int
  :param seq_len: length of the query sequence

  :type blast_results: list
  :param blast_results: list of tabular blast hits as lists of strings
  """
  blast_results.sort(key=score, reverse=True)
  try:
    r1 = blast_results[0]
  except IndexError:
    return None  # no match
  try:
    r2 = blast_results[1]
  except IndexError:
    return False
  if abs(al2seq(r1, seq_len) - al2seq(r2, seq_len)) <= min_al2seq:
    return True
  if score(r1) - score(r2) <= min_score_diff:
    return True
  return False
