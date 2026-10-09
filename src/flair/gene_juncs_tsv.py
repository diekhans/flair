"""Reads grouped by the gene and junction chain they support, written per region
by flair spliceevents and read back when the region's samples are joined.

Columns, in order, with no header line:

  gene         the gene the reads were assigned to
  juncs        the junction chain, junctions comma separated, each start.end
  start        the read's start
  end          the read's end
  strand       the strand of the chain
  read         the read name

One row per read, so a chain appears once per supporting read.
"""
from flair.pycbio.tsv import TsvReader, TsvWriter

COLUMNS = ('gene', 'juncs', 'start', 'end', 'strand', 'read')
TYPE_MAP = {'gene': str, 'juncs': str, 'start': int, 'end': int,
            'strand': str, 'read': str}

class GeneJuncsReader(TsvReader):
    def __init__(self, gene_juncs_txt):
        super().__init__(gene_juncs_txt, columns=COLUMNS, typeMap=TYPE_MAP)

class GeneJuncsWriter(TsvWriter):
    def __init__(self, gene_juncs_txt, *, outFh=None):
        super().__init__(gene_juncs_txt, columns=COLUMNS, typeMap=TYPE_MAP,
                         outFh=outFh, writeHeader=False)
