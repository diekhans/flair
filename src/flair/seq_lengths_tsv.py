"""The length of every sequence in a FASTA, and the histogram of those lengths.

Sequence lengths, no header line:

  name         the FASTA record name, the whole line after '>'
  length       the number of residues

Histogram, no header line:

  length       a length that at least one sequence has
  count        how many sequences have it, in increasing length order

Both written by fasta_seq_lengths.  Headerless because they feed plotting and
shell pipelines that expect bare columns; the column names here are the
documentation.
"""
from flair.pycbio.tsv import TsvReader, TsvWriter

LENGTH_COLUMNS = ('name', 'length')
HISTOGRAM_COLUMNS = ('length', 'count')

class SeqLengthsReader(TsvReader):
    def __init__(self, seq_lengths_tsv):
        super().__init__(seq_lengths_tsv, columns=LENGTH_COLUMNS,
                         typeMap={'name': str, 'length': int})

class SeqLengthHistogramReader(TsvReader):
    def __init__(self, histogram_tsv):
        super().__init__(histogram_tsv, columns=HISTOGRAM_COLUMNS, defaultColType=int)

class SeqLengthsWriter(TsvWriter):
    def __init__(self, seq_lengths_tsv):
        super().__init__(seq_lengths_tsv, columns=LENGTH_COLUMNS,
                         typeMap={'name': str, 'length': int}, writeHeader=False)

class SeqLengthHistogramWriter(TsvWriter):
    def __init__(self, histogram_tsv):
        super().__init__(histogram_tsv, columns=HISTOGRAM_COLUMNS,
                         defaultColType=int, writeHeader=False)
