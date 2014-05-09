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
Merge a collection of FASTA files into a single one.
"""

import sys, os, argparse, zipfile, tarfile, contextlib, logging

from bl.core.seq.io.fasta import SimpleFastaReader as FastaReader
from bl.core.utils import NullLogger


LOG_LEVELS = ['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL']
DEFAULT_OUTPUT = "output.fa"


def make_logger(level_str="INFO", filename=None):
  formatter = logging.Formatter(
    fmt='%(asctime)s|%(levelname)-8s|%(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    )
  logger = logging.getLogger(__name__)
  for h in logger.handlers:
    logger.removeHandler(h)
  if filename:
    handler = logging.FileHandler(filename, 'w')
  else:
    handler = logging.StreamHandler()
  handler.setFormatter(formatter)
  logger.addHandler(handler)
  logger.setLevel(getattr(logging, level_str))
  return logger


class Multiplexer(object):

  SEP = "/"

  def __critical(self, msg, exc_type=RuntimeError):
    self.logger.critical(msg)
    raise exc_type(msg)

  def __set_protocol_tar(self):
    self.open = tarfile.open
    self.getmembers = lambda a: a.getmembers()
    self.getname = lambda m: m.name
    self.isdir = lambda m: m.isdir()
    self.extractfile = lambda a, m: contextlib.closing(a.extractfile(m))

  def __set_protocol_zip(self):
    self.open = zipfile.ZipFile
    self.getmembers = lambda a: a.infolist()
    self.getname = lambda m: m.filename
    self.isdir = lambda m: m.file_size <= 0
    self.extractfile = lambda a, m: a.open(m, "rU")

  def __set_protocol(self, fn):
    if not isinstance(fn, basestring):  # TODO: support a list of filenames
      self.__critical("not a string: %r" % (fn,), TypeError)
    try:
      if tarfile.is_tarfile(fn):
        self.logger.debug("%r looks like a tar file" % (fn,))
        self.__set_protocol_tar()
        return
    except IOError as e:
      self.__critical(str(e), IOError)
    if zipfile.is_zipfile(fn):
      self.logger.debug("%r looks like a zip file" % (fn,))
      self.__set_protocol_zip()
      return
    self.__critical("unsupported archive format: %r" % (fn,), ValueError)

  def __init__(self, fasta_in, sep=SEP, logger=None):
    self.logger = logger or NullLogger()
    self.__set_protocol(fasta_in)
    self.input = fasta_in
    self.sep = sep

  def __iter__(self):
    sep = self.sep
    with self.open(self.input) as a:
      members = self.getmembers(a)
      for i, m in enumerate(members):
        basename = os.path.basename(self.getname(m))
        self.logger.info(
          "processing [%d/%d]: %r" % (i+1, len(members), basename)
          )
        if self.isdir(m):
          continue
        with self.extractfile(a, m) as f:
          tag = os.path.splitext(basename)[0]
          reader = FastaReader(f)
          for header, seq in reader:
            yield sep.join((tag, header)), seq


def make_parser():
  desc = "Merge a collection of FASTA files into a single one"
  parser = argparse.ArgumentParser(description=desc)
  parser.add_argument('input', metavar="ZIP_ARCHIVE",
                      help='zip archive of FASTA files')
  parser.add_argument('-o', '--output', metavar="FILE",
                      help='output FASTA file', default=DEFAULT_OUTPUT)
  parser.add_argument('--log-file', metavar="FILE", help='log file [stderr]')
  parser.add_argument('--log-level', type=str, choices=LOG_LEVELS,
                      help='logging level', default='INFO')
  return parser


def main(argv):
  parser = make_parser()
  args = parser.parse_args(argv)
  logger = make_logger(level_str=args.log_level, filename=args.log_file)
  multiplexer = Multiplexer(args.input, logger=logger)
  with open(args.output, "w") as fo:
    for header, seq in multiplexer:
      fo.write(">%s\n%s\n" % (header, seq))
  logger.info("wrote %r" % (args.output,))


if __name__ == "__main__":
  main(sys.argv[1:])
