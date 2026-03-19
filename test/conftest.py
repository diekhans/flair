import pytest
from flair.gtf_io import gtf_data_parser, GtfAttrsSet

@pytest.fixture(scope="session")
def basic_gtf_data():
    return gtf_data_parser("input/basic.annotation.gtf")

@pytest.fixture(scope="session")
def basic_gtf_data_all():
    return gtf_data_parser("input/basic.annotation.gtf", attrs=GtfAttrsSet.ALL)
