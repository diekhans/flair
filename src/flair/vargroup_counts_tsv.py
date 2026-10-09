"""Reads grouped by the set of variants they carry, from flair variantquant.

Columns, in order, with no header line:

  chrom            the reference the reads aligned to
  gene             the gene the group was assigned to, the commonest among the
                   reads' variants
  num_reads        how many reads carry this set of positions
  num_positions    how many variant positions the set has
  mods_per_read    comma separated count of carried variants, one per read
  varpos           comma separated variant positions

Headerless, as written before the file had a format of its own; the column names
here are the documentation.
"""
from contextlib import contextmanager
from flair.pycbio.sys import fileOps
from flair.pycbio.tsv import TsvWriter

COLUMNS = ('chrom', 'gene', 'num_reads', 'num_positions', 'mods_per_read', 'varpos')

class VarGroupCountsWriter(TsvWriter):
    def __init__(self, vargroup_counts_tsv, *, outFh=None):
        super().__init__(vargroup_counts_tsv, columns=COLUMNS, defaultColType=str,
                         outFh=outFh, writeHeader=False)

@contextmanager
def vargroup_counts_writer(vargroup_counts_tsv):
    "write atomically; the path appears only once complete"
    fileOps.ensureFileDir(vargroup_counts_tsv)
    with fileOps.AtomicFileCreate(vargroup_counts_tsv) as tmp_tsv:
        with VarGroupCountsWriter(tmp_tsv) as writer:
            yield writer
