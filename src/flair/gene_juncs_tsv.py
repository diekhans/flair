"""Reads grouped by the gene and junction chain they support, written per region
by flair spliceevents and read back when the region's samples are joined.

Columns, in order, with no header line:

  gene         the gene the reads were assigned to
  juncs        the junction chain, each junction start.end; comma separated in
               the file and a list on the row
  start        the read's start
  end          the read's end
  strand       the strand of the chain
  read         the read name

One row per read, so a chain appears once per supporting read.
"""
from flair.pycbio.tsv import TsvReader, TsvWriter
from flair.tsv_column_types import commaListType

COLUMNS = ('gene', 'juncs', 'start', 'end', 'strand', 'read')
TYPE_MAP = {'gene': str, 'juncs': commaListType, 'start': int, 'end': int,
            'strand': str, 'read': str}

class GeneJuncsReader(TsvReader):
    def __init__(self, gene_juncs_txt):
        super().__init__(gene_juncs_txt, columns=COLUMNS, typeMap=TYPE_MAP)

class GeneJuncsWriter(TsvWriter):
    def __init__(self, gene_juncs_txt, *, outFh=None):
        super().__init__(gene_juncs_txt, columns=COLUMNS, typeMap=TYPE_MAP,
                         outFh=outFh, writeHeader=False)

def parse_junc_coord(coord):
    "one junction of the chain, written start.end"
    start, end = coord.split('.')
    return int(start), int(end)
