"""Reads carrying each variant, from flair variantquant.

Columns, in order:

  varpos       the variant, chrom:position
  gene         the gene it falls in, the first of the sorted list when it falls
               in more than one
  transcript   the transcript, left empty for now
  <sample>     one column per sample, holding unmodified;modified read counts

A variant is written when some sample has at least --threshold reads over it and
some sample carries it, or always with --output_all.
"""
from contextlib import contextmanager
from flair.pycbio.sys import fileOps
from flair.pycbio.tsv import TsvWriter

ID_COLUMNS = ('varpos', 'gene', 'transcript')

def columns(sample_columns):
    return list(ID_COLUMNS) + list(sample_columns)

class VarCountsWriter(TsvWriter):
    def __init__(self, var_counts_tsv, sample_columns, *, outFh=None):
        super().__init__(var_counts_tsv, columns=columns(sample_columns),
                         defaultColType=str, outFh=outFh)

@contextmanager
def var_counts_writer(var_counts_tsv, sample_columns):
    "write atomically; the path appears only once complete"
    fileOps.ensureFileDir(var_counts_tsv)
    with fileOps.AtomicFileCreate(var_counts_tsv) as tmp_tsv:
        with VarCountsWriter(tmp_tsv, sample_columns) as writer:
            yield writer
