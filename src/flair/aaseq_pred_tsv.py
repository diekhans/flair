"""Protein predictions for variant-carrying transcripts, written by
predict_aaseq_withvar.

Columns, in order, with no header line:

  transcript     the transcript, with its variants, as the modified transcript
                 FASTA names it
  productivity   the prediction: productive, or why it is not
  utr_vars       which untranslated regions carry variants, 5utr, 3utr, both
                 comma separated, or empty
  aaseq          the predicted protein sequence
"""
from flair.pycbio.tsv import TsvWriter

COLUMNS = ('transcript', 'productivity', 'utr_vars', 'aaseq')

class AaSeqPredWriter(TsvWriter):
    def __init__(self, aaseq_pred_tsv, *, outFh=None):
        super().__init__(aaseq_pred_tsv, columns=COLUMNS, defaultColType=str,
                         outFh=outFh, writeHeader=False)
