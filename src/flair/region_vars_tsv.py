"""VCF variants grouped by the region they fall in, written by
annotate_group_vcf_vars.

Columns, in order, with no header line:

  chrom        the chromosome
  region       the region the variant falls in
  pos          the position
  ref          the reference allele
  alts         comma separated alternate alleles
  name         the variant name from the VCF
"""
from flair.pycbio.tsv import TsvWriter

COLUMNS = ('chrom', 'region', 'pos', 'ref', 'alts', 'name')

class RegionVarsWriter(TsvWriter):
    def __init__(self, region_vars_tsv, *, outFh=None):
        super().__init__(region_vars_tsv, columns=COLUMNS, defaultColType=str,
                         outFh=outFh, writeHeader=False)
