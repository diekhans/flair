import pytest
from flair.gtf_io import gtf_data_parser

@pytest.fixture(scope="session")
def basic_gtf_data():
    return gtf_data_parser("input/basic.annotation.gtf")
