
# relative to src/flair
FLAKE8_SRC = \
	__init__.py \
	flair_cli.py \
	flair_partition.py \
	flair_spliceevents.py \
        flair_variantmodels.py \
        flair_variantquant.py \
        flair_transcriptome.py \
	intron_support.py

FLAKE8_TEST = \
	test_correct_lib.py

FLAKE8_CHECK = ${FLAKE8_SRC:%=src/flair/%} ${FLAKE8_TEST:%=test/%}
