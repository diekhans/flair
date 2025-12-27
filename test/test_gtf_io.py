"""
tests of gff_io module
"""
import pytest
from io import StringIO
from flair.gtf_io import gtf_data_parser, GtfIdError, gtf_write_row
from flair import SeqRange

@pytest.fixture(scope="session")
def basic_data():
    return gtf_data_parser("input/basic.annotation.gtf")

def test_parse(basic_data):
    trans_ids = sorted(basic_data.iter_transcript_ids())
    assert len(trans_ids) == 59
    assert trans_ids[0:10] == ['ENST00000225792.10', 'ENST00000256078.9', 'ENST00000279052.10', 'ENST00000311936.8', 'ENST00000348547.6',
                               'ENST00000357394.8', 'ENST00000411577.5', 'ENST00000413587.5', 'ENST00000416206.5', 'ENST00000438317.5']

def test_trans_lookup(basic_data):
    transcript = basic_data.get_transcript('ENST00000556131.1')
    assert transcript is not None
    assert transcript.transcript_id == 'ENST00000556131.1'
    assert transcript.gene_id == "ENSG00000133703.12"
    assert transcript.gene_name == "KRAS"
    assert transcript.coords == SeqRange(name='chr12', start=25233818, end=25250929, strand='-')
    assert transcript.coords_no_strand == SeqRange(name='chr12', start=25233818, end=25250929)
    assert transcript.attrs['level'] == 2
    assert transcript.attrs['transcript_support_level'] == "1"

def test_trans_error(basic_data):
    transcript = basic_data.get_transcript('Fred')
    assert transcript is None
    with pytest.raises(GtfIdError, match=r"unknown transcript id `Barney'"):
        basic_data.fetch_transcript('Barney')

def test_chroms(basic_data):
    chroms = basic_data.get_chroms()
    assert chroms == ['chr12', 'chr17', 'chr20']


# overlaps KRAS
KRAS_OVER_RANGE = SeqRange("chr12", 25205000, 25252000)
KRAS_TRANS_IDS = ["ENST00000256078.9", "ENST00000311936.8", "ENST00000556131.1", "ENST00000557334.5"]

def test_iter_overlap(basic_data):
    trans_ids = sorted([t.transcript_id
                        for t in basic_data.iter_overlap_transcripts_sr(KRAS_OVER_RANGE)])
    assert trans_ids == KRAS_TRANS_IDS

def test_iter_overlap_strand(basic_data):
    seq_range = SeqRange(*KRAS_OVER_RANGE[0:3], '-')
    trans_ids = sorted([t.transcript_id
                        for t in basic_data.iter_overlap_transcripts_sr(seq_range)])
    assert trans_ids == KRAS_TRANS_IDS

def test_iter_overlap_other_strand(basic_data):
    seq_range = SeqRange(*KRAS_OVER_RANGE[0:3], '+')
    trans_ids = sorted([t.transcript_id
                        for t in basic_data.iter_overlap_transcripts_sr(seq_range)])
    assert trans_ids == []

def test_rec_str(basic_data):
    trans_str = str(basic_data.fetch_transcript("ENST00000256078.9"))
    assert trans_str == ('chr12\tHAVANA\ttranscript\t25205246\t25250929\t.\t-\t.\t'
                         'gene_id "ENSG00000133703.12"; transcript_id "ENST00000256078.9"; gene_type "protein_coding"; '
                         'gene_name "KRAS"; transcript_type "protein_coding"; transcript_name "KRAS-201"; level 2; '
                         'protein_id "ENSP00000256078.4"; transcript_support_level "1"; hgnc_id "HGNC:6407"; '
                         'tag "RNA_Seq_supported_only"; tag "basic"; tag "appris_principal_4"; tag "CCDS"; ccdsid "CCDS8703.1"; '
                         'havana_gene "OTTHUMG00000171193.4"; havana_transcript "OTTHUMT00000412232.4";')

def test_rec_str1(basic_data):
    rec_fh = StringIO()

    gtf_write_row(rec_fh, 'chr12', 'HAVANA', 'transcript', 25205246, 25250929, None, '-', None,
                  gene_id="ENSG00000133703.12", transcript_id="ENST00000256078.9", gene_name="KRAS")

    assert rec_fh.getvalue() == ('chr12\tHAVANA\ttranscript\t25205247\t25250929\t.\t-\t.\t'
                                 'gene_id "ENSG00000133703.12"; gene_name "KRAS"; transcript_id "ENST00000256078.9";\n')

def test_rec_str2(basic_data):
    rec_fh = StringIO()

    gtf_write_row(rec_fh, 'chr12', 'HAVANA', 'transcript', 25205246, 25250929, 10, '-', 1,
                  gene_id="ENSG00000133703.12", transcript_id="ENST00000256078.9", gene_name="KRAS")

    assert rec_fh.getvalue() == ('chr12\tHAVANA\ttranscript\t25205247\t25250929\t10\t-\t1\t'
                                 'gene_id "ENSG00000133703.12"; gene_name "KRAS"; transcript_id "ENST00000256078.9";\n')

def test_rec_str3(basic_data):
    rec_fh = StringIO()

    gtf_write_row(rec_fh, 'chr12', 'HAVANA', 'transcript', 25205246, 25250929, 100, '-', 2,
                  attrs={"gene_id": "ENSG00000133703.12",
                         "transcript_id": "ENST00000256078.9",
                         "gene_name": "KRAS"})

    assert rec_fh.getvalue() == ('chr12\tHAVANA\ttranscript\t25205247\t25250929\t100\t-\t2\t'
                                 'gene_id "ENSG00000133703.12"; transcript_id "ENST00000256078.9"; gene_name "KRAS";\n')
