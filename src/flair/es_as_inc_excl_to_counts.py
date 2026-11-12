#!/usr/bin/env python3
"""Attach counts to the skipped exon events es_as.py calls, writing the
es events.quant.tsv."""

import argparse
import os
import sys
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import numpy as np  # noqa: E402 - openblas setting must be before numpy import
from flair.counts_matrix_tsv import read_sample_columns, read_counts_rows  # noqa: E402
from flair.es_events_tsv import EsEventsReader  # noqa: E402
from flair.event_quant_tsv import EventQuantWriter  # noqa: E402
from flair.cli_help import FlairHelpFormatter  # noqa: E402
from flair.pycbio.sys import cli  # noqa: E402

def side_counts(data, n_samps, isos):
    "the counts of one side of an event, summed over its isoforms"
    vals = np.asarray([data.get(iso, np.zeros(n_samps, dtype=np.int64)) for iso in isos])
    return np.sum(vals, axis=0)

def write_es_events(counts_matrix_tsv, es_events_tsv):
    sample_names = read_sample_columns(counts_matrix_tsv)
    n_samps = len(sample_names)
    # int64, so a large count is written plainly rather than in scientific notation
    # and is not rounded at the 24 bit mantissa of float32
    data = {row.isoform_id: np.asarray(row.counts, dtype=np.int64)
            for row in read_counts_rows(counts_matrix_tsv)}

    with EventQuantWriter(None, sample_names, outFh=sys.stdout) as writer:
        for row in EsEventsReader(es_events_tsv):
            # an exon needs an isoform on each side to be a skipping event
            if row.inclusion_isos and row.exclusion_isos:
                writer.writeSide('inclusion', row.exon, row.exon,
                                 side_counts(data, n_samps, row.inclusion_isos),
                                 row.inclusion_isos)
                writer.writeSide('exclusion', row.exon, row.exon,
                                 side_counts(data, n_samps, row.exclusion_isos),
                                 row.exclusion_isos)

def build_parser():
    parser = argparse.ArgumentParser(prog='es_as_inc_excl_to_counts', description=__doc__,
                                     formatter_class=FlairHelpFormatter)
    parser.add_argument('counts_matrix', help='the isoform counts matrix of flair quantify')
    parser.add_argument('es_events_tsv', help='the skipped exon events of es_as')
    return parser


def main():
    with cli.ErrorHandler():
        opts = build_parser().parse_args()
        write_es_events(opts.counts_matrix, opts.es_events_tsv)


if __name__ == '__main__':
    main()
