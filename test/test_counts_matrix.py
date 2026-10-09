"""
Tests for counts_matrix module.
"""
import pytest
from flair import FlairInputDataError
from flair.counts_matrix import (SampleInfo, sample_info_path, write_sample_info,
                                 read_sample_info, select_condition_pair,
                                 condition_column_indexes, CountsRow, ID_COLUMNS,
                                 read_header, read_counts_rows, read_sample_columns,
                                 write_counts_matrix)

SAMPLE_INFOS = (SampleInfo('s1', 'ctl', 'b1'),
                SampleInfo('s2', 'ctl', 'b1'),
                SampleInfo('s3', 'test', 'b2'))

def write_test_matrix(path, sample_columns):
    write_counts_matrix(path, sample_columns,
                        [CountsRow('gene1', 'iso1', ['1'] * len(sample_columns))])

def test_read_header_rejects_old_format(tmp_path):
    "a matrix with one isoform_gene column cannot be split reliably, so it is refused"
    counts = str(tmp_path / 'x.counts.tsv')
    with open(counts, 'w') as fh:
        fh.write('ids\ts1\ts2\n')
        fh.write('iso1_gene1\t1\t2\n')
    with pytest.raises(FlairInputDataError, match="must start with the columns"):
        read_header(counts)

def test_counts_matrix_round_trip(tmp_path):
    counts = str(tmp_path / 'x.counts.tsv')
    rows = [CountsRow('gene1', 'iso1', ['1', '2']), CountsRow('gene1', 'iso2', ['3', '4'])]
    write_counts_matrix(counts, ['s1', 's2'], rows)
    assert read_header(counts) == list(ID_COLUMNS) + ['s1', 's2']
    assert read_sample_columns(counts) == ['s1', 's2']
    assert read_counts_rows(counts) == rows

def test_sample_info_path_counts_tsv():
    assert sample_info_path('out/flair.quantify.counts.tsv') == 'out/flair.quantify.sample_info.tsv'

def test_sample_info_path_other_names():
    assert sample_info_path('matrix.tsv') == 'matrix.sample_info.tsv'
    assert sample_info_path('matrix') == 'matrix.sample_info.tsv'

def test_read_sample_info_round_trip(tmp_path):
    counts = str(tmp_path / 'x.counts.tsv')
    write_test_matrix(counts, [si.sample_id for si in SAMPLE_INFOS])
    write_sample_info(sample_info_path(counts), SAMPLE_INFOS)
    assert read_sample_info(counts) == list(SAMPLE_INFOS)

def test_read_sample_info_allows_underscores(tmp_path):
    "the fields are columns of their own, so they may contain the joining character"
    infos = [SampleInfo('s_1', 'ctl_a', 'b_1'), SampleInfo('s_2', 'test_a', 'b_2')]
    counts = str(tmp_path / 'x.counts.tsv')
    write_test_matrix(counts, [si.sample_id for si in infos])
    write_sample_info(sample_info_path(counts), infos)
    assert read_sample_info(counts) == infos

def test_read_sample_info_falls_back_to_column_names(tmp_path):
    "a counts matrix written before the sample info file existed"
    counts = str(tmp_path / 'x.counts.tsv')
    write_test_matrix(counts, ['s1_ctl_b1', 's2_test_b1'])
    assert read_sample_info(counts) == [SampleInfo('s1_ctl_b1', 'ctl', 'b1'),
                                        SampleInfo('s2_test_b1', 'test', 'b1')]

def test_read_sample_info_must_match_counts_columns(tmp_path):
    counts = str(tmp_path / 'x.counts.tsv')
    write_test_matrix(counts, ['s1', 's2'])
    write_sample_info(sample_info_path(counts), SAMPLE_INFOS)
    with pytest.raises(FlairInputDataError, match="does not describe the columns"):
        read_sample_info(counts)

def test_select_condition_pair_sorted_default():
    assert select_condition_pair(['test', 'ctl', 'ctl'], '', '', 'x.tsv') == ('ctl', 'test')

def test_select_condition_pair_ignores_column_order():
    "the whole point: the same conditions in any column order give the same pair"
    assert (select_condition_pair(['ctl', 'test', 'test', 'ctl'], '', '', 'x.tsv') ==
            select_condition_pair(['test', 'test', 'ctl', 'ctl'], '', '', 'x.tsv'))

def test_select_condition_pair_named():
    assert select_condition_pair(['ctl', 'test'], 'test', 'ctl', 'x.tsv') == ('test', 'ctl')

def test_select_condition_pair_unknown_name():
    with pytest.raises(FlairInputDataError, match="is not a condition in"):
        select_condition_pair(['ctl', 'test'], 'nope', 'ctl', 'x.tsv')

def test_select_condition_pair_needs_both():
    with pytest.raises(FlairInputDataError, match="must both be given"):
        select_condition_pair(['ctl', 'test'], 'ctl', '', 'x.tsv')

def test_select_condition_pair_needs_exactly_two():
    with pytest.raises(FlairInputDataError, match="has 3 conditions"):
        select_condition_pair(['a', 'b', 'c'], '', '', 'x.tsv')

def test_condition_column_indexes():
    assert condition_column_indexes(['ctl', 'test', 'ctl'], 'ctl') == [0, 2]
    assert condition_column_indexes(['ctl', 'test', 'ctl'], 'test') == [1]
