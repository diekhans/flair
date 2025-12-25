"""
tests of gff_io module
"""
import pytest
from flair.gtf_io import gtf_data_parser, GtfIdError
from flair import SeqRange

@pytest.fixture(scope="session")
def gtf_basic_data():
    return gtf_data_parser("input/basic.annotation.gtf")

def test_basic_parse(gtf_basic_data):
    trans_ids = sorted(gtf_basic_data.iter_transcript_ids())
    assert len(trans_ids) == 59
    assert trans_ids[0:10] == ['ENST00000225792.10', 'ENST00000256078.9', 'ENST00000279052.10', 'ENST00000311936.8', 'ENST00000348547.6',
                               'ENST00000357394.8', 'ENST00000411577.5', 'ENST00000413587.5', 'ENST00000416206.5', 'ENST00000438317.5']

def test_basic_trans_lookup(gtf_basic_data):
    transcript = gtf_basic_data.get_transcript('ENST00000556131.1')
    assert transcript is not None
    assert transcript.transcript_id == 'ENST00000556131.1'
    assert transcript.gene_id == "ENSG00000133703.12"
    assert transcript.gene_name == "KRAS"
    assert transcript.coords == SeqRange(name='chr12', start=25233818, end=25250929, strand='-')
    assert transcript.coords_no_strand == SeqRange(name='chr12', start=25233818, end=25250929)

def test_basic_trans_error(gtf_basic_data):
    transcript = gtf_basic_data.get_transcript('Fred')
    assert transcript is None
    with pytest.raises(GtfIdError, match=r"unknown transcript id `Barney'"):
        gtf_basic_data.fetch_transcript('Barney')
