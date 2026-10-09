"""What each combined id was made from, written by flair combine.

Columns, in order, with no header line:

  new_id       the FLT or FLG id the combined transcriptome uses
  source_ids   the ids it was built from, '; ' separated, each written
               <samples>:<id> with the samples comma separated

flair combine numbers afresh, so these files are the only link between an id in
a combined transcriptome and the per-sample ids behind it; see ids.rst.  One
file for isoforms, one for genes.
"""
from flair.pycbio.tsv import TsvReader, TsvWriter

COLUMNS = ('new_id', 'source_ids')
SOURCE_SEP = '; '

class CombinedMapReader(TsvReader):
    def __init__(self, combined_map_txt):
        super().__init__(combined_map_txt, columns=COLUMNS, defaultColType=str)

class CombinedMapWriter(TsvWriter):
    def __init__(self, combined_map_txt, *, outFh=None):
        super().__init__(combined_map_txt, columns=COLUMNS, defaultColType=str,
                         outFh=outFh, writeHeader=False)

    def writeSources(self, new_id, source_ids):
        self.writeRow((new_id, SOURCE_SEP.join(source_ids)))
