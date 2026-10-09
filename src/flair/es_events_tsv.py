"""Skipped exon events called from an isoform BED, one row per exon.

Columns, in order, with no header line:

  exon           the skipped exon, chrom:acceptor-donor
  strand         the strand of its gene
  num_inclusion  how many isoforms include the exon
  num_exclusion  how many skip it
  inclusion_isos comma separated ids of the including isoforms
  exclusion_isos comma separated ids of the skipping isoforms

Written by es_as and read by es_as_inc_excl_to_counts, which attaches counts and
writes an event quant TSV; see event_quant_tsv.  A row with no excluding isoform
is not an event and is dropped there.

Headerless because the two programs are a pair and the columns are fixed; the
column names here are the documentation.
"""
from flair.pycbio.tsv import TsvReader, TsvWriter

COLUMNS = ('exon', 'strand', 'num_inclusion', 'num_exclusion',
           'inclusion_isos', 'exclusion_isos')

TYPE_MAP = {'exon': str, 'strand': str, 'num_inclusion': int, 'num_exclusion': int,
            'inclusion_isos': str, 'exclusion_isos': str}

class EsEventsReader(TsvReader):
    def __init__(self, es_events_tsv):
        super().__init__(es_events_tsv, columns=COLUMNS, typeMap=TYPE_MAP)

class EsEventsWriter(TsvWriter):
    def __init__(self, es_events_tsv, *, outFh=None):
        super().__init__(es_events_tsv, columns=COLUMNS, typeMap=TYPE_MAP,
                         outFh=outFh, writeHeader=False)

def isoform_ids(isos):
    "the comma separated isoform ids of one side of an event"
    return [] if isos == '' else isos.split(',')
