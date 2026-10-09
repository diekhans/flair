"""Which reads support a thing, one row per thing.

Columns, in order, with no header line:

  name         the isoform, allele group or allele the reads support
  reads        comma separated read names

Written by flair transcriptome, quantify, alleles and isoalleles, each with its
own file name.  Headerless, since these files are joined and filtered by id with
shell tools as often as they are read by a program.
"""
from flair.pycbio.tsv import TsvReader, TsvWriter

COLUMNS = ('name', 'reads')

class ReadMapReader(TsvReader):
    def __init__(self, read_map_tsv):
        super().__init__(read_map_tsv, columns=COLUMNS, defaultColType=str)

class ReadMapWriter(TsvWriter):
    def __init__(self, read_map_tsv, *, outFh=None):
        super().__init__(read_map_tsv, columns=COLUMNS, defaultColType=str,
                         outFh=outFh, writeHeader=False)

    def writeReads(self, name, reads):
        "reads are written in the order given, so callers sort when it matters"
        self.writeRow((name, ','.join(reads)))

def read_names(reads):
    "the read names of one row"
    return [] if reads == '' else reads.split(',')
