"""An input of a run that is one of the files the run writes is rejected before the
run starts: the subcommands' output suffixes and run_checks.check_inputs_not_outputs.
"""
import pytest
from flair import FlairInputDataError
from flair import flair_combine, flair_diffExp, flair_diffSplice, flair_fusion
from flair import flair_isoalleles, flair_quantify, flair_spliceevents, flair_transcriptome
from flair.output_compress import output_compression, output_paths, set_output_compression
from flair.run_checks import check_inputs_not_outputs

# the subcommand modules that name what they write, and one output of each that an
# input could plausibly be
SUBCOMMANDS = ((flair_transcriptome, '.isoforms.bed'),
               (flair_combine, '.combined.isoforms.bed'),
               (flair_quantify, '.isoforms.counts.tsv'),
               (flair_spliceevents, '.spliceevents.bed'),
               (flair_fusion, '.fusions.isoforms.bed'),
               (flair_isoalleles, '.isoalleles.bed'),
               (flair_diffExp, '.isoforms_sig_exp_change_norm_by_gene.tsv'),
               (flair_diffSplice, '.alt3.events.quant.tsv'))


@pytest.fixture
def compress(request):
    "a run writing its outputs compressed or plain, the setting restored after the test"
    was = output_compression()
    set_output_compression(request.param)
    yield request.param
    set_output_compression(was)


def module_outputs(module, output):
    return output_paths(output, module.OUTPUT_TEXT_SUFFIXES, getattr(module, 'OUTPUT_BINARY_SUFFIXES', ()))


@pytest.mark.parametrize('module,suffix', SUBCOMMANDS)
def test_every_output_of_a_subcommand_is_rejected_as_an_input(module, suffix):
    outputs = module_outputs(module, 'out/flair')
    for path in outputs:
        with pytest.raises(FlairInputDataError, match='file is one this run writes'):
            check_inputs_not_outputs(outputs, {'isoform BED': path})


@pytest.mark.parametrize('module,suffix', SUBCOMMANDS)
def test_the_main_output_of_a_subcommand_is_one_of_its_paths(module, suffix):
    assert 'out/flair' + suffix in module_outputs(module, 'out/flair')


@pytest.mark.parametrize('module,suffix', SUBCOMMANDS)
def test_an_input_beside_the_output_is_not_rejected(module, suffix):
    check_inputs_not_outputs(module_outputs(module, 'out/flair'),
                             {'annotation GTF': 'out/annotation.gtf',
                              'aligned reads BAM': ['out/sample.bam'],
                              'region BED': None})


@pytest.mark.parametrize('compress', (False, True), indirect=True)
def test_a_compressed_run_claims_the_compressed_names(compress):
    outputs = module_outputs(flair_transcriptome, 'out/flair')
    assert ('out/flair.isoforms.bed.gz' in outputs) == compress
    assert ('out/flair.isoforms.bed' in outputs) != compress
    # the bigBed is written plain whatever the run compresses
    assert 'out/flair.isoforms.bb' in outputs


def test_the_same_file_by_another_path_is_rejected(tmp_path):
    "the paths are compared as realpaths, so a link or a dotted path is the same file"
    output = f'{tmp_path}/flair'
    bed = output + '.combined.isoforms.bed'
    open(bed, 'w').close()
    by_dots = f'{tmp_path}/./sub/../flair.combined.isoforms.bed'
    with pytest.raises(FlairInputDataError, match='file is one this run writes'):
        check_inputs_not_outputs(module_outputs(flair_combine, output), {'isoform BED': by_dots})
