"""Per-event Fisher's exact results from diffsplice_fishers_exact.

Columns, in order:

  the columns of the event quant TSV it was given; see event_quant_tsv
  <sample1>-<sample2>_pval   the p-value of the test between those two samples

One row per side of each event, the two sides of an event carrying the same
p-value, since the test is of the event's two by two table.
"""
from contextlib import contextmanager
from flair.pycbio.sys import fileOps
from flair.pycbio.tsv import TsvWriter

def pval_column(colname1, colname2):
    return f'{colname1}-{colname2}_pval'

class DiffspliceFishersWriter(TsvWriter):
    def __init__(self, fishers_tsv, event_quant_columns, colname1, colname2):
        super().__init__(fishers_tsv,
                         columns=list(event_quant_columns) + [pval_column(colname1, colname2)],
                         defaultColType=str)

@contextmanager
def diffsplice_fishers_writer(fishers_tsv, event_quant_columns, colname1, colname2):
    "write atomically; the path appears only once complete"
    fileOps.ensureFileDir(fishers_tsv)
    with fileOps.AtomicFileCreate(fishers_tsv) as tmp_tsv:
        with DiffspliceFishersWriter(tmp_tsv, event_quant_columns, colname1, colname2) as writer:
            yield writer
