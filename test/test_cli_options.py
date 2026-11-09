"""Options that take a list of values, comma separated or by repeating the option."""
from flair.cli_options import option_values


def test_one_comma_separated_option():
    assert option_values(['a,b,c']) == ['a', 'b', 'c']


def test_repeated_option():
    assert option_values(['a', 'b,c', 'd']) == ['a', 'b', 'c', 'd']


def test_trailing_and_empty_values():
    assert option_values(['a,b,', ',c']) == ['a', 'b', 'c']


def test_option_not_given():
    assert option_values(None) == []
