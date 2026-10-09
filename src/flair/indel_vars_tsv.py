"""Insertions and deletions called by identify_vars, one row per variant.

Columns, in order, with no header line:

  chrom        the chromosome
  pos          the position
  indel_type   insertion or deletion
  ref_seq      the reference sequence
  var_seq      the variant sequence
  indel_reads  reads carrying the variant
  tot_cov      reads covering the position
  filters      filter tags, dn marking a dense cluster of variants around this
               one; comma separated in the file and a list on the row
"""
from flair.pycbio.tsv import TsvWriter
from flair.tsv_column_types import commaListType

COLUMNS = ('chrom', 'pos', 'indel_type', 'ref_seq', 'var_seq', 'indel_reads',
           'tot_cov', 'filters')

class IndelVarsWriter(TsvWriter):
    def __init__(self, indel_vars_tsv, *, outFh=None):
        super().__init__(indel_vars_tsv, columns=COLUMNS, defaultColType=str,
                         typeMap={'filters': commaListType},
                         outFh=outFh, writeHeader=False)
