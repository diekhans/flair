"""
Isoform names are hashes of their structure, and isoform ids are given by name, so
isoforms of the same coordinates on different chromosomes need different names: the
single-exon ERCC spike-ins starting at 1 all had one id.

Two isoforms taking one name is a defect, and stops the run rather than one of them
replacing the other in the dict it is stored in.  FLAIR 3.0.1 named an end variant
after the column of the BED that holds the exon count, so variants of unrelated
junction chains were named alike and all but the last were lost with nothing
reported: BrooksLabUCSC/flair#732.
"""
import os
import subprocess
import sys
import pytest
from flair import FlairError
from flair.flair_bed import FlairBed
from flair.flair_transcriptome import CandidateIsoforms, assign_final_ends, get_new_ids, store_by_name
from flair.isoform_data import Isoform, Junc
from flair.read_ends_tsv import COLUMNS as ENDS_COLUMNS

JUNCS = (Junc(1100, 1500), Junc(1700, 2000))


def test_single_exon_isoforms_on_different_chromosomes():
    a = Isoform('ERCC-00003', '+', (), 0, 1006)
    b = Isoform('ERCC-00042', '+', (), 0, 1006)
    assert a.name != b.name
    assert a.name == Isoform('ERCC-00003', '+', (), 0, 1006).name
    assert a.name != Isoform('ERCC-00003', '-', (), 0, 1006).name


def test_spliced_isoforms_on_different_chromosomes():
    assert Isoform('chr1', '+', JUNCS, 800, 2200).name != Isoform('chr2', '+', JUNCS, 800, 2200).name
    # before the ends are known, by the junctions
    assert Isoform('chr1', '+', JUNCS).name != Isoform('chr2', '+', JUNCS).name


def test_names_are_the_same_between_runs():
    "Python's string hash changes between runs, so names must not depend on it"
    code = ("from flair.isoform_data import Isoform, Junc; "
            "print(Isoform('chr1', '+', (Junc(1100, 1500),), 800, 2200).name, Isoform('ERCC-00003', '+', (), 0, 1006).name)")
    names = set()
    for seed in ('1', '2'):
        env = dict(os.environ, PYTHONHASHSEED=seed)
        names.add(subprocess.run([sys.executable, '-c', code], env=env, capture_output=True, text=True,
                                 check=True).stdout)
    assert len(names) == 1


####
# two isoforms of one name, BrooksLabUCSC/flair#732
####
ENDS_HEADER = '\t'.join(ENDS_COLUMNS) + '\n'

# the junction chain of the end variant tests, with the ends line below written
# against its first junction start and its last junction end
END_JUNCS = (Junc(1200, 1300), Junc(1400, 1500))


def _isoform(start, end, chrom='chr1', strand='+'):
    return Isoform(chrom, strand, END_JUNCS, start, end)


def _ends_line(name, start, end):
    "count_sam_transcripts' ends line for a read aligned from start to end"
    return f'{name}\ttx\t0\t{1200 - start}\t\t1\t{end - 1500}\t\n'


def _final_ends(tmp_path, read_ends, max_ends):
    "the isoforms assign_final_ends gives a junction chain with these reads, by name"
    isoform = _isoform(800, 2200)
    isoform.name = 'tx'
    ends = tmp_path / 'ends.tsv'
    ends.write_text(ENDS_HEADER + ''.join(_ends_line(f'r{i}', s, e) for i, (s, e) in enumerate(read_ends)))
    return assign_final_ends({'tx': isoform}, str(ends), str(tmp_path / 'map.txt'), normalize_ends=False,
                             max_ends=max_ends, end_window=100, sjc_support=2, frac_support=0.05)


def test_the_end_variants_of_one_chain_each_keep_their_own_name(tmp_path):
    "what BrooksLabUCSC/flair#732 reports lost: three end groups give three isoforms"
    read_ends = [(1000, 2000)] * 3 + [(1100, 1800)] * 3 + [(900, 2150)] * 3
    final = _final_ends(tmp_path, read_ends, max_ends=3)
    assert len(final) == 3
    assert sorted((iso.start, iso.end) for iso in final.values()) == [(900, 2150), (1000, 2000), (1100, 1800)]


def test_an_end_variant_is_named_for_its_own_ends(tmp_path):
    "the names come from the structure, so no variant keeps the first-pass isoform's name"
    final = _final_ends(tmp_path, [(1000, 2000)] * 3 + [(1100, 1800)] * 3, max_ends=2)
    assert 'tx' not in final
    assert all(iso.name == name for name, iso in final.items())


def test_isoforms_of_one_structure_have_one_name():
    assert _isoform(1000, 2000).name == _isoform(1000, 2000).name


def test_isoforms_differing_in_an_end_are_named_apart():
    name = _isoform(1000, 2000).name
    assert _isoform(1001, 2000).name != name
    assert _isoform(1000, 2001).name != name


def test_a_repeated_candidate_name_stops_the_run():
    candidates = CandidateIsoforms()
    candidates.add(_isoform(1000, 2000))
    with pytest.raises(FlairError, match='two candidate isoforms have the name'):
        candidates.add(_isoform(1000, 2000))


def test_the_repeated_name_message_names_both_isoforms():
    isoforms = {}
    store_by_name(isoforms, _isoform(1000, 2000), 'candidate')
    with pytest.raises(FlairError, match=r'chr1:1000-2000 with 2 junctions and chr1:1000-2000'):
        store_by_name(isoforms, _isoform(1000, 2000), 'candidate')


def test_a_candidate_name_that_is_free_is_stored():
    candidates = CandidateIsoforms()
    candidates.add(_isoform(1000, 2000))
    candidates.add(_isoform(1000, 2100))
    assert len(candidates.isoforms) == 2


def _write_isoform_bed(path, names):
    with open(path, 'w') as fh:
        for i, name in enumerate(names):
            bed = FlairBed('chr1', 1000 + i, 2000 + i, name, score=3, strand='+', thickStart=1000 + i,
                           thickEnd=2000 + i, itemRgb='0', gene_id='g', read_support=3, frac_support=0.5)
            bed.addBlock(1000 + i, 2000 + i)
            bed.write(fh)


def test_two_isoforms_of_one_name_in_the_bed_stop_the_run(tmp_path):
    "a repeat here would give two isoforms one FLT id, in the BED, the counts and the read map"
    output = str(tmp_path / 'out')
    _write_isoform_bed(output + '.isoforms.bed', ['12345', '12345'])
    with pytest.raises(FlairError, match='two isoforms of .* have the name 12345'):
        get_new_ids(output)


def test_each_isoform_of_the_bed_takes_the_next_id(tmp_path):
    output = str(tmp_path / 'out')
    _write_isoform_bed(output + '.isoforms.bed', ['12345', '67890'])
    assert get_new_ids(output) == {'12345': 'FLT00000001', '67890': 'FLT00000002'}
