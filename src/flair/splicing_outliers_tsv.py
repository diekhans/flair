"""Samples whose PSI of an alternative of a splicing event is far from the median of
all samples', from spliceevents_outliers.  One row per outlying sample of an
alternative.

Columns, in order:

  event                  the event, as in the event files (splicing_events_tsv)
  alternative            the alternative of the event
  event_type             the kind of event
  gene                   the gene the event was found in
  gene_name              its name
  additional_genes       other genes sharing the event's junctions
  additional_gene_names  their names
  sample                 the sample that stands out
  median_psi             the median PSI across samples
  dev_iqr_half           half the interquartile range, the spread PSI is judged
                         against
  sample_psi             this sample's PSI
  samples_with_psi       how many samples had enough reads to take a PSI from
  alternative_reads      this sample's reads of the alternative
  psi_total              the total its PSI divides by: the event's junction reads
                         with --psi junc, the locus reads with --psi tot
  delta_psi_to_median    sample_psi minus median_psi, as an absolute value
  dev_from_median        that difference in units of dev_iqr_half, or inf when the
                         other samples' PSIs are identical, so the spread it
                         divides by is zero

The four gene columns contain several values, comma separated in the file and a list
on the row.  The numbers are written as computed, rounded to six places, and read
back as strings.

Written unfiltered as .spliceevents.outliers.tsv and, after the deviation cutoff and
leaving out outliers a simpler event explains, as
.spliceevents.outliers.filtered.tsv.
"""
from flair.pycbio.tsv import TsvReader, TsvWriter
from flair.tsv_column_types import commaListType
from flair.splicing_events_tsv import GENE_COLUMNS

COLUMNS = ('event', 'alternative', 'event_type', *GENE_COLUMNS, 'sample', 'median_psi', 'dev_iqr_half',
           'sample_psi', 'samples_with_psi', 'alternative_reads', 'psi_total',
           'delta_psi_to_median', 'dev_from_median')
TYPE_MAP = {'alternative_reads': int, 'psi_total': int,
            **{col: commaListType for col in GENE_COLUMNS}}

class SplicingOutliersReader(TsvReader):
    def __init__(self, splicing_outliers_tsv):
        super().__init__(splicing_outliers_tsv, defaultColType=str, typeMap=TYPE_MAP)

class SplicingOutliersWriter(TsvWriter):
    def __init__(self, splicing_outliers_tsv, *, outFh=None):
        super().__init__(splicing_outliers_tsv, columns=COLUMNS, defaultColType=str, typeMap=TYPE_MAP,
                         outFh=outFh)
