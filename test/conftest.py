import os
import pytest
from flair.gtf_io import gtf_data_parser, GtfAttrsSet

@pytest.fixture(scope="session")
def make_output_dir():
    "where the Makefile put the files these tests read, set by test/Makefile"
    return os.environ.get('FLAIR_TEST_OUTPUT_DIR', 'output')

@pytest.fixture(scope="session")
def basic_gtf_data():
    return gtf_data_parser("input/basic.annotation.gtf")

@pytest.fixture(scope="session")
def basic_gtf_data_all():
    return gtf_data_parser("input/basic.annotation.gtf", attrs=GtfAttrsSet.ALL)
