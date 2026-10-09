#!/usr/bin/env python3

import argparse
import csv
import os
import scipy.stats as sps
from flair import FlairInputDataError
from flair.counts_matrix import read_header, read_counts_rows, ID_COLUMNS
from flair.pycbio.sys import cli


def build_parser():
    desc = """Calculates the usage of each isoform as a fraction of the total expression
    of the gene and compares this between samples."""

    parser = argparse.ArgumentParser(prog='diff_iso_usage', description=desc)
    parser.add_argument('counts_matrix_tsv',
                        help='counts matrix TSV from flair-quantify')
    parser.add_argument('colname1',
                        help='the name of the column of the first sample')
    parser.add_argument('colname2',
                        help='the name of the column of the second sample')
    parser.add_argument('outfile',
                        help='output filename containing the p-value associated with differential '
                        'isoform usage for each isoform')
    return parser

def sample_column_index(header, colname, counts_matrix_tsv):
    "index of a named sample column within a row's counts"
    sample_columns = header[len(ID_COLUMNS):]
    if colname not in sample_columns:
        raise FlairInputDataError(
            f"{counts_matrix_tsv} has no sample column named {colname}; it has: "
            f"{' '.join(sample_columns)}")
    return sample_columns.index(colname)

def diff_iso_usage(counts_matrix_tsv, colname1, colname2, outfilename):  # noqa: C901 - FIXME: reduce complexity
    header = read_header(counts_matrix_tsv)
    col1 = sample_column_index(header, colname1, counts_matrix_tsv)
    col2 = sample_column_index(header, colname2, counts_matrix_tsv)

    counts = {}
    for row in read_counts_rows(counts_matrix_tsv):
        if row.gene_id not in counts:
            counts[row.gene_id] = {}
        counts[row.gene_id][row.isoform_id] = [float(row.counts[col1]), float(row.counts[col2])]

    with open(outfilename, 'wt') as outfile:
        writer = csv.writer(outfile, delimiter='\t', lineterminator=os.linesep)
        # the column names say which sample each number came from, so the sign of
        # delta_PSI can be read from the file without knowing the argument order
        writer.writerow(['geneID', 'isoID', 'fisher_pval',
                         f'this_iso_{colname1}_count', f'this_iso_{colname2}_count',
                         f'other_isos_{colname1}_count', f'other_isos_{colname2}_count',
                         f'{colname1}_PSI', f'{colname2}_PSI', 'delta_PSI'])
        geneordered = sorted(counts.keys())
        for gene in geneordered:
            generes = []
            for iso in counts[gene]:
                thesecounts = counts[gene][iso]
                othercounts = [0, 0]
                for iso_ in counts[gene]:
                    if iso_ != iso:
                        othercounts[0] += counts[gene][iso_][0]
                        othercounts[1] += counts[gene][iso_][1]
                ctable = [thesecounts, othercounts]
                if thesecounts[0] + othercounts[0] == 0 or thesecounts[1] + othercounts[1] == 0 or sum(thesecounts) == 0 or sum(othercounts) == 0:  # do not test this isoform if no gene exp in one sample
                    generes.append([gene, iso, 'NA'] + ctable[0] + ctable[1] + ['NA', 'NA', 'NA'])
                else:
                    s1PSI, s2PSI, deltaPSI = 'NA', 'NA', 'NA'
                    if ctable[1][0] + ctable[0][0] > 0:
                        s1PSI = round(ctable[0][0] / (ctable[1][0] + ctable[0][0]), 3)
                    if ctable[1][1] + ctable[0][1] > 0:
                        s2PSI = round(ctable[0][1] / (ctable[1][1] + ctable[0][1]), 3)
                    if s1PSI != 'NA' and s2PSI != 'NA':
                        deltaPSI = round(s2PSI - s1PSI, 3)
                    psi_data = [s1PSI, s2PSI, deltaPSI]

                    generes.append([gene, iso, sps.fisher_exact(ctable)[1]] + ctable[0] + ctable[1] + psi_data)

            # if not generes:
            #     writer.writerow([gene, iso, 'NA'] + ctable[0] + ctable[1] + psi_data)
            #     continue

            for res in generes:
                writer.writerow(res)


def main():
    args = cli.parseArgsWithLogging(build_parser())
    with cli.ErrorHandler():
        diff_iso_usage(args.counts_matrix_tsv, args.colname1, args.colname2, args.outfile)


if __name__ == '__main__':
    main()
