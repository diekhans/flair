"""
Tests for flair_diffExp module.
"""
import pytest
from flair import FlairInputDataError
from flair.counts_matrix import MAX_REPORTED_IDS
from flair.flair_diffExp import check_iso_gene_ids

def write_counts_matrix(path, iso_gene_ids):
    with open(path, 'w') as fh:
        fh.write('\t'.join(['ids', 's1_ctl_b1', 's2_test_b1']) + '\n')
        for iso_gene in iso_gene_ids:
            fh.write('\t'.join([iso_gene, '1', '2']) + '\n')

def write_inputs(tmp_path, iso_gene_ids):
    counts = str(tmp_path / 'x.counts.tsv')
    write_counts_matrix(counts, iso_gene_ids)
    return counts

def test_check_iso_gene_ids_annotated(tmp_path):
    check_iso_gene_ids(write_inputs(tmp_path, ['ENST1.1_ENSG1.1', 'ENST2.1_ENSG1.1']))

def test_check_iso_gene_ids_refseq(tmp_path):
    "the transcript half carries an underscore of its own"
    check_iso_gene_ids(write_inputs(tmp_path, ['NM_000123.4_GENEA', 'NM_000124.5_GENEA']))

def test_check_iso_gene_ids_novel_locus(tmp_path):
    "a novel isoform names its locus rather than a gene accession, which is a gene here"
    check_iso_gene_ids(write_inputs(tmp_path, ['read1_chr12:25205000', 'read2_chr12:25205000']))

def test_check_iso_gene_ids_no_gene(tmp_path):
    counts = write_inputs(tmp_path, ['ENST1.1_ENSG1.1', 'ENST2.1_ENSG1.1', 'ENST3.1'])
    with pytest.raises(FlairInputDataError, match=r"1 row ids .* do not name a gene: ENST3\.1;"):
        check_iso_gene_ids(counts)

def test_check_iso_gene_ids_no_gene_list_elided(tmp_path):
    counts = write_inputs(tmp_path, ['ENST%d.1' % num for num in range(MAX_REPORTED_IDS + 5)])
    with pytest.raises(FlairInputDataError, match=r"15 row ids .*, \.\.\.;"):
        check_iso_gene_ids(counts)

def test_check_iso_gene_ids_bare_refseq(tmp_path):
    "bare RefSeq ids split into a gene of their own, caught only by the grouping check"
    counts = write_inputs(tmp_path, ['NM_000123.4', 'NM_000124.5'])
    with pytest.raises(FlairInputDataError, match="holds more than one isoform"):
        check_iso_gene_ids(counts)

def test_check_iso_gene_ids_no_rows(tmp_path):
    with pytest.raises(FlairInputDataError, match="has no isoform rows"):
        check_iso_gene_ids(write_inputs(tmp_path, []))
