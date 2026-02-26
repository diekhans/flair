
# relative to src/flair
FLAKE8_SRC = \
	__init__.py \
	flair_cli.py \
	flair_partition.py \
        flair_transcriptome.py \
        flair_variantmodels.py \
	intron_support.py \
	filter_transcriptome_align.py \
	gtf_io.py

FIXME_NOT_WORKING =\
	flair_spliceevents.py \
        flair_variantquant.py


FLAKE8_TEST = \
	test_correct_lib.py

FLAKE8_CHECK = ${FLAKE8_SRC:%=src/flair/%} ${FLAKE8_TEST:%=test/%}
