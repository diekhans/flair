"""
Tests for partition_runner module.
"""
import os
import pytest
from flair import SeqRange
from flair.gtf_io import gtf_data_parser, GtfAttrsSet
from flair.intron_support import IntronSupport
from flair.partition_runner import Partition, PartitionRunner

GTF_FILE = "input/basic.annotation.gtf"       # relative to test/ directory (run via make lib-tests)
INTRONS_BED = "input/basic.shortread_junctions.bed"

# Regions matching chroms in the test GTF
CHR12_REGION = SeqRange("chr12", 0, 30000000)
CHR17_REGION = SeqRange("chr17", 0, 30000000)
CHR20_REGION = SeqRange("chr20", 0, 70000000)

# KRAS transcripts expected in chr12 region
KRAS_TRANS_IDS = sorted(["ENST00000256078.9", "ENST00000311936.8",
                         "ENST00000556131.1", "ENST00000557334.5"])


@pytest.fixture(scope="session")
def gtf_data():
    return gtf_data_parser(GTF_FILE, attrs=GtfAttrsSet.FLAIR)


@pytest.fixture(scope="session")
def intron_support():
    is_db = IntronSupport()
    is_db.load_introns_bed(INTRONS_BED)
    return is_db


@pytest.fixture()
def runner(tmp_path, gtf_data, intron_support):
    regions = [CHR12_REGION, CHR17_REGION, CHR20_REGION]
    return PartitionRunner(regions, str(tmp_path), gtf_data=gtf_data, intron_support=intron_support)


def test_partition_dirs_created(runner):
    for part in runner:
        assert os.path.isdir(part.temp_dir)


def test_partition_pickle_files_exist(runner):
    for part in runner:
        assert os.path.exists(part.temp_path("gtf_data.pkl"))
        assert os.path.exists(part.temp_path("intron_support.pkl"))


def test_partition_load_gtf_data(runner):
    for part in runner:
        loaded = part.load_gtf_data()
        assert loaded is not None
        chroms = loaded.get_chroms()
        assert all(c == part.region.name for c in chroms)


def test_partition_load_intron_support(runner):
    chr12_part = next(p for p in runner if p.region.name == "chr12")
    loaded = chr12_part.load_intron_support()
    assert loaded is not None
    introns = list(loaded.introns("chr12"))
    assert len(introns) > 0


def test_partition_gtf_subset_chr12(runner):
    chr12_part = next(p for p in runner if p.region.name == "chr12")
    loaded = chr12_part.load_gtf_data()
    trans_ids = sorted(loaded.iter_transcript_ids())
    assert KRAS_TRANS_IDS == [t for t in trans_ids if t in KRAS_TRANS_IDS]


def test_partition_gtf_no_cross_contamination(runner):
    chr17_part = next(p for p in runner if p.region.name == "chr17")
    loaded = chr17_part.load_gtf_data()
    assert "chr12" not in loaded.get_chroms()
    assert "chr20" not in loaded.get_chroms()


def test_partition_intron_support_empty_for_chr17(runner):
    chr17_part = next(p for p in runner if p.region.name == "chr17")
    loaded = chr17_part.load_intron_support()
    assert loaded is not None
    assert list(loaded.introns("chr17")) == []


def test_runner_len(runner):
    assert len(runner) == 3


def test_runner_iter(runner):
    parts = list(runner)
    assert len(parts) == 3
    assert all(isinstance(p, Partition) for p in parts)


def test_runner_regions(runner):
    regions = [p.region for p in runner]
    assert regions == [CHR12_REGION, CHR17_REGION, CHR20_REGION]


def test_run_calls_func_for_each_partition(runner):
    seen_regions = []

    def record_region(*, partition, gtf_data, intron_support):
        seen_regions.append(partition.region)

    runner.run(record_region)
    assert seen_regions == [CHR12_REGION, CHR17_REGION, CHR20_REGION]


def test_run_passes_correct_kwargs(runner):
    received = {}

    def capture(*, partition, gtf_data, intron_support, extra):
        received[partition.region.name] = (gtf_data, intron_support, extra)

    runner.run(capture, extra="value")
    assert set(received.keys()) == {"chr12", "chr17", "chr20"}
    assert all(v[2] == "value" for v in received.values())


def test_run_gtf_data_is_gtf_data_type(runner):
    from flair.gtf_io import GtfData
    received_types = []

    def check_type(*, partition, gtf_data, intron_support):
        received_types.append(type(gtf_data))

    runner.run(check_type)
    assert all(t is GtfData for t in received_types)


def test_run_intron_support_is_intron_support_type(runner):
    received_types = []

    def check_type(*, partition, gtf_data, intron_support):
        received_types.append(type(intron_support))

    runner.run(check_type)
    assert all(t is IntronSupport for t in received_types)


def test_partition_no_gtf(tmp_path):
    regions = [CHR12_REGION]
    runner = PartitionRunner(regions, str(tmp_path))
    part = list(runner)[0]
    assert not os.path.exists(part.temp_path("gtf_data.pkl"))
    assert part.load_gtf_data() is None


def test_partition_no_intron_support(tmp_path):
    regions = [CHR12_REGION]
    runner = PartitionRunner(regions, str(tmp_path))
    part = list(runner)[0]
    assert not os.path.exists(part.temp_path("intron_support.pkl"))
    assert part.load_intron_support() is None
