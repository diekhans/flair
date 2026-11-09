"""The input options of flair combine: --prefixes names the transcriptomes by prefix,
--isoform_beds by path, and both take a list (cli_options)."""
from flair.flair_combine import parse_input


def write_beds(tmp_path, prefixes):
    for prefix in prefixes:
        bed = tmp_path / (prefix + '.isoforms.bed')
        bed.parent.mkdir(parents=True, exist_ok=True)
        bed.write_text('')


def test_prefixes_name_the_isoform_beds(tmp_path):
    write_beds(tmp_path, ('one', 'sub/two'))
    paths = parse_input(None, [f'{tmp_path}/one', f'{tmp_path}/sub/two'], None)
    assert paths == [f'{tmp_path}/one.isoforms.bed', f'{tmp_path}/sub/two.isoforms.bed']


def test_isoform_beds_are_paths():
    assert parse_input(None, None, ['a.bed,b.bed', 'c.bed']) == ['a.bed', 'b.bed', 'c.bed']
