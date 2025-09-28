
from flair.intron_support import IntronSupport, load_read_bed_introns, load_annot_gtf_introns, load_read_star_introns

def _assert_introns(over, expect):
    assert [str(o) for o in over] == expect

def _basic_load_reads_annot_support_test(intron_support):
    "same data, loaded from BED or STAR SJ"
    # chr12:25209911-25215436(-)
    start = 25209911
    end = 25215436
    over = intron_support.overlap("chr12", start, start + 1)
    _assert_introns(over, [
        'Intron(chr12:25209911-25215436(-) annot=True read=True read_cnt=12',
        'Intron(chr12:25209911-25225613(-) annot=True read=False read_cnt=0',
        'Intron(chr12:25209911-25245273(-) annot=True read=False read_cnt=0',
        'Intron(chr12:25209911-25225613(-) annot=False read=True read_cnt=96',
    ])

    over = intron_support.overlap("chr12", end - 1, end)
    _assert_introns(over, [
        'Intron(chr12:25209911-25215436(-) annot=True read=True read_cnt=12',
    ])

    over = intron_support.overlap_introns("chr12", start, end, 5)
    _assert_introns(over, [
        'Intron(chr12:25209911-25215436(-) annot=True read=True read_cnt=12',
    ])

def test_load_bed_annot_support():
    intron_support = IntronSupport()
    load_read_bed_introns(intron_support, "input/basic.shortread_junctions.bed")
    load_annot_gtf_introns(intron_support, "input/basic.annotation.gtf")
    _basic_load_reads_annot_support_test(intron_support)

def test_load_star_annot_support():
    intron_support = IntronSupport()
    load_read_star_introns(intron_support, "input/basic.shortread_junctions.tab")
    load_annot_gtf_introns(intron_support, "input/basic.annotation.gtf")
    _basic_load_reads_annot_support_test(intron_support)
