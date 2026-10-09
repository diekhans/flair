"""Transcripts per million, one row per isoform.

Columns, in order:

  gene_id      the gene the isoform belongs to
  isoform_id   the isoform
  <sample>     one column per sample, named by its sample id, holding TPM

The same shape as a counts matrix, with the counts replaced by TPM, so it is read
and written through counts_matrix_tsv with float values.  flair quantify --tpm
writes it beside the counts matrix.
"""
from flair.counts_matrix_tsv import (CountsRow, counts_matrix_writer,  # noqa: F401
                                     read_counts_rows, read_sample_columns)

def tpm_writer(tpm_tsv, sample_columns):
    "write a TPM matrix atomically; values are floats rather than counts"
    return counts_matrix_writer(tpm_tsv, sample_columns, countType=float)
