
# relative to src/flair
FLAKE8_SRC = \
    __init__.py \
    flair_cli.py \
    flair_partition.py \
    flair_transcriptome.py \
    flair_variantmodels.py

# FIXME: 
FLAKE8_SRC_BROKEN =  \
   flair_spliceevents.py \
   flair_variantquant.py \

FLAKE8_CHECK = ${FLAKE8_SRC:%=src/flair/%}
