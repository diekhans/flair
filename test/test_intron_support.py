"""Reading junction support: the score column of an intron BED is a read count, not
a browser shade."""
from flair.intron_support import IntronSupport

# chr1:1000-5000, supported by more reads than the BED score convention's 1000
INTRON_BED = "chr1\t1000\t5000\tGT/AG\t5000\t+\n"


def _loaded(tmp_path, text):
    bed = tmp_path / 'junctions.bed'
    bed.write_text(text)
    support = IntronSupport()
    support.load_introns_bed(str(bed))
    return list(support.introns())


def test_support_above_the_bed_score_cap_is_kept(tmp_path):
    introns = _loaded(tmp_path, INTRON_BED)
    assert len(introns) == 1
    assert introns[0].read_support_cnt == 5000


def test_the_same_intron_twice_sums_its_support(tmp_path):
    introns = _loaded(tmp_path, INTRON_BED + INTRON_BED)
    assert introns[0].read_support_cnt == 10000


def test_support_below_the_cap_is_unchanged(tmp_path):
    introns = _loaded(tmp_path, "chr1\t1000\t5000\tGT/AG\t7\t+\n")
    assert introns[0].read_support_cnt == 7
