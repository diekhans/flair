#!/usr/bin/env python3

import os
import sys
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import numpy as np  # noqa: E402 - openblas setting must be before numpy import
from flair.counts_matrix_tsv import read_sample_columns, read_counts_rows  # noqa: E402

sample_names = read_sample_columns(sys.argv[1])
nSamps = len(sample_names)
data = {row.isoform_id: np.asarray(row.counts, dtype=np.float32)
        for row in read_counts_rows(sys.argv[1])}

with open(sys.argv[2]) as fin2:
    print('\t'.join(['feature_id', 'coordinate'] + sample_names + ['isoform_ids']))
    for line in fin2:
        cols = line.rstrip().split()
        if int(cols[3]) == 0:
            continue
        else:
            exon, strand, _, exc, incIsos, excIsos = cols
        incVals = np.asarray([data.get(x, np.zeros(nSamps)) for x in incIsos.split(",")])
        excVals = np.asarray([data.get(x, np.zeros(nSamps)) for x in excIsos.split(",")])

        incVals = np.sum(incVals, axis=0)
        excVals = np.sum(excVals, axis=0)
        # totVals = incVals + excVals

        # must have at least 1 count to support the inc and exc of this exon
        # if incVals.all() < 1 or excVals.all() < 1:
        #     continue

        print("inclusion_%s" % exon, exon, "\t".join(str(x) for x in incVals), incIsos, sep="\t")
        print("exclusion_%s" % exon, exon, "\t".join(str(x) for x in excVals), excIsos, sep="\t")

        # print(exon,"\t".join("%.2f" % x for x in incVals/totVals))
