"""Checks that a read id names one read, in the three places a repeat is caught:
the clipping of the reads being realigned, the reads labelled with their allele, and
a sample's read map, which is the only one that sees the whole sample."""
import pysam
import pytest
from flair import FlairInputDataError
from flair.count_sam_transcripts import parse_sam, read_clipping, read_isoforms_bed
from flair.flair_isoalleles import LABEL_TAG, _write_labelled_read
from flair.read_map_tsv import ReadMapWriter, check_read_ids_unique

_HEADER = pysam.AlignmentHeader.from_dict({'SQ': [{'SN': 'chr1', 'LN': 10000}]})


def write_read_map(path, rows):
    with ReadMapWriter(path) as writer:
        for name, reads in rows:
            writer.writeReads(name, reads)
    return path


def write_clipping(path, rows):
    with open(path, 'w') as fh:
        for name, left, right in rows:
            fh.write(f'{name}\t{left}\t{right}\n')
    return str(path)


def aligned_read(name, start, length=100):
    read = pysam.AlignedSegment(_HEADER)
    read.query_name, read.reference_id, read.reference_start = name, 0, start
    read.cigartuples = [(pysam.CMATCH, length)]
    read.query_sequence = 'A' * length
    return read


class Collector():
    "stands in for the labelled BAM being written"
    def __init__(self):
        self.reads = []

    def write(self, read):
        self.reads.append(read)


def test_read_map_ids_are_unique(tmp_path):
    path = write_read_map(str(tmp_path / 'x.read.map.tsv'),
                          [('iso1', ['readA', 'readB']), ('iso2', ['readC'])])
    check_read_ids_unique(path)


def test_read_map_id_in_two_rows(tmp_path):
    "two molecules of one name full-length for different isoforms, as merged runs give"
    path = write_read_map(str(tmp_path / 'x.read.map.tsv'),
                          [('iso1', ['readA', 'readB']), ('iso2', ['readA'])])
    with pytest.raises(FlairInputDataError, match=r'1 read ids name more than one read .*readA'):
        check_read_ids_unique(path)


def test_read_map_id_twice_in_one_row(tmp_path):
    path = write_read_map(str(tmp_path / 'x.read.map.tsv'), [('iso1', ['readA', 'readA'])])
    with pytest.raises(FlairInputDataError, match=r'1 read ids name more than one read .*readA'):
        check_read_ids_unique(path)


def test_read_map_empty(tmp_path):
    check_read_ids_unique(write_read_map(str(tmp_path / 'x.read.map.tsv'), []))
    check_read_ids_unique(write_read_map(str(tmp_path / 'y.read.map.tsv'), [('iso1', ['readA'])]))


def test_clipping_of_each_read(tmp_path):
    clipping = read_clipping(write_clipping(tmp_path / 'c.tsv',
                                            [('readA', 10, 20), ('readB', 0, 5)]))
    assert clipping == {'readA': [10, 20], 'readB': [0, 5]}


def test_clipping_of_a_repeated_read_id(tmp_path):
    path = write_clipping(tmp_path / 'c.tsv', [('readA', 10, 20), ('readA', 0, 5)])
    with pytest.raises(FlairInputDataError, match=r'read id readA names more than one read'):
        read_clipping(path)


def test_labelled_read_fetched_again_for_an_overlapping_gene():
    "the same alignment twice is the batch's genes overlapping, not two reads"
    written, out = {}, Collector()
    for _ in range(2):
        _write_labelled_read(aligned_read('readA', 1000), 'h1', written, out, 'sample.bam')
    assert [read.query_name for read in out.reads] == ['readA']
    assert out.reads[0].get_tag(LABEL_TAG) == 'h1'


def test_labelled_read_id_at_two_positions():
    written, out = {}, Collector()
    _write_labelled_read(aligned_read('readA', 1000), 'h1', written, out, 'sample.bam')
    with pytest.raises(FlairInputDataError, match=r'read id readA names more than one read of sample.bam'):
        _write_labelled_read(aligned_read('readA', 5000), 'h2', written, out, 'sample.bam')


def write_transcript_sam(path, alignments, tlen=300):
    "a SAM of whole-transcript alignments, as minimap2 writes them, grouped by read"
    with open(path, 'w') as fh:
        fh.write(f'@SQ\tSN:t1\tLN:{tlen}\n')
        for name in alignments:
            fh.write('\t'.join([name, '0', 't1', '1', '60', f'{tlen}M', '*', '0', '0',
                                'A' * tlen, '*', f'AS:i:{tlen}', f'MD:Z:{tlen}']) + '\n')
    return str(path)


def write_transcript_bed(path, tlen=300):
    with open(path, 'w') as fh:
        fh.write('\t'.join(['chr1', '1000', str(1000 + tlen), 't1', '0', '+', '1000',
                            str(1000 + tlen), '0', '1', f'{tlen},', '0,']) + '\n')
    return str(path)


def count_sam(tmp_path, alignments):
    info = read_isoforms_bed(isoforms=write_transcript_bed(tmp_path / 'isoforms.bed'))
    return parse_sam(write_transcript_sam(tmp_path / 'aln.sam', alignments), info, {},
                     stringent=False, check_splice=False, fusion_breakpoints=None,
                     allow_UTR_indels=False, trimmedreads=None, soft_clipping_buffer=50,
                     output_endpos=None, trust_ends=False)


def test_alignments_of_one_read_each(tmp_path):
    counted = count_sam(tmp_path, ['readA', 'readB', 'readC'])
    assert sorted(read for reads in counted.values() for read, _, _ in reads) == ['readA', 'readB', 'readC']


def test_alignments_of_a_read_id_that_comes_back(tmp_path):
    "an id whose alignments are in two places in the stream is two reads"
    with pytest.raises(FlairInputDataError, match=r'read id readA names more than one read'):
        count_sam(tmp_path, ['readA', 'readB', 'readA'])
