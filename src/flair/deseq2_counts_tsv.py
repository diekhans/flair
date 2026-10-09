"""The counts matrix handed to DESeq2, one row per feature.

Columns, in order:

  <empty>      the feature id, a gene id in the gene table and an isoform id in
               the isoform table.  The column name is empty because R reads the
               file with row.names = 1 and takes this column as the row names
  <sample>     one integer column per sample, named as in the formula matrix

Written by diffexp into its working directory after filtering by expression; not
an output anyone reads directly.
"""
from contextlib import contextmanager
from flair.pycbio.sys import fileOps
from flair.pycbio.tsv import TsvWriter

# R takes column one as row names, so it carries no name
FEATURE_COLUMN = ''

class Deseq2CountsWriter(TsvWriter):
    def __init__(self, deseq2_counts_tsv, sample_columns):
        super().__init__(deseq2_counts_tsv,
                         columns=[FEATURE_COLUMN] + list(sample_columns),
                         defaultColType=str)

    def writeFeature(self, feature_id, counts):
        self.writeRow([feature_id] + list(counts))

@contextmanager
def deseq2_counts_writer(deseq2_counts_tsv, sample_columns):
    "write atomically; the path appears only once complete"
    fileOps.ensureFileDir(deseq2_counts_tsv)
    with fileOps.AtomicFileCreate(deseq2_counts_tsv) as tmp_tsv:
        with Deseq2CountsWriter(tmp_tsv, sample_columns) as writer:
            yield writer

def write_deseq2_counts(deseq2_counts_tsv, sample_columns, feature_ids, feature_counts):
    with deseq2_counts_writer(deseq2_counts_tsv, sample_columns) as writer:
        for feature_id, counts in zip(feature_ids, feature_counts):
            writer.writeFeature(feature_id, counts)
