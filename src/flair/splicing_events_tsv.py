"""Alternative splicing events called by flair spliceevents, one row per event.

Columns, in order:

  eventname            the event, <type>-(<strand>)-<gene>
  eventtype            the kind of event
  gene                 the gene it was called in
  junctions_included   junctions supporting inclusion
  junctions_excluded   junctions supporting exclusion
  outer_junctions      the junctions flanking the event
  exons                the exons it involves
  <sample>             one column per sample

Three files share these columns and differ only in what the sample columns hold:

  .diffsplice.counts.tsv    reads supporting the event
  .diffsplice.PSIjunc.tsv   its share of the reads over its own junctions
  .diffsplice.PSItot.tsv    its share of the reads over the whole locus

A PSI column is empty when the sample has too few reads to take a fraction from,
rather than 0, which would read as "never included".

spliceevents writes one headerless file per partition and concatenates them, so
the header is written once when the parts are joined.
"""
from flair.pycbio.tsv import TsvWriter

EVENT_COLUMNS = ('eventname', 'eventtype', 'gene', 'junctions_included',
                 'junctions_excluded', 'outer_junctions', 'exons')

def columns(sample_columns):
    return list(EVENT_COLUMNS) + list(sample_columns)

class SplicingEventsWriter(TsvWriter):
    def __init__(self, splicing_events_tsv, sample_columns, *, outFh=None, writeHeader=True):
        super().__init__(splicing_events_tsv, columns=columns(sample_columns),
                         defaultColType=str, outFh=outFh, writeHeader=writeHeader)

    def writeEvent(self, event_fields, sample_values):
        self.writeRow(list(event_fields) + list(sample_values))
