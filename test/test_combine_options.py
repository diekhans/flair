"""The input options of flair combine: --prefixes and --isoform_beds take a comma
separated list, and --prefixes may also be repeated."""
from flair.flair_combine import option_values, parse_input


def test_one_comma_separated_option():
    assert option_values(['a,b,c']) == ['a', 'b', 'c']


def test_repeated_option():
    assert option_values(['a', 'b,c', 'd']) == ['a', 'b', 'c', 'd']


def test_trailing_and_empty_values():
    assert option_values(['a,b,', ',c']) == ['a', 'b', 'c']


def test_prefixes_name_the_isoform_beds(tmp_path):
    for prefix in ('one', 'sub/two'):
        bed = tmp_path / (prefix + '.isoforms.bed')
        bed.parent.mkdir(exist_ok=True)
        bed.write_text('')
    paths = parse_input(None, [f'{tmp_path}/one', f'{tmp_path}/sub/two'], None)
    assert paths == [f'{tmp_path}/one.isoforms.bed', f'{tmp_path}/sub/two.isoforms.bed']
