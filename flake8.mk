
# relative to src/flair
FLAKE8_SRC = \
	__init__.py \
	flair_cli.py \
	flair_partition.py \
        flair_transcriptome.py \
	intron_support.py \
	filter_transcriptome_align.py \
	gtf_io.py

FIXME_NOT_WORKING =\
        flair_variantmodels.py \
	flair_spliceevents.py \
        flair_variantquant.py


FLAKE8_TEST = \
	test_correct_lib.py \
	test_gtf_io.py \
	bin/gtf_io_perf

FLAKE8_CHECK = ${FLAKE8_SRC:%=src/flair/%} ${FLAKE8_TEST:%=test/%}
