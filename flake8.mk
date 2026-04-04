
# relative to src/flair
FLAKE8_SRC = \
	__init__.py \
	flair_cli.py \
	flair_partition.py \
        flair_transcriptome.py \
	intron_support.py \
	filter_transcriptome_align.py \
	partition_runner.py \
	gtf_io.py \
        flair_variantmodels.py \
	flair_spliceevents.py \
        flair_variantquant.py \
	fasta_seq_lengths.py \
	diffsplice_fishers_exact.py \
	counts_to_tpm.py \
	mark_intron_retention.py \
	read_processing.py \
	annotate_group_vcf_vars.py \
	annotate_aaseq_with_uniprot.py \
	filter_transcriptome_chim.py \
	flair_align.py \
	flair_diffExp.py \
	junction_correct.py \
	remove_internal_priming.py

FLAKE8_TEST = \
	test_correct_lib.py \
	test_gtf_io.py \
	bin/gtf_io_perf

FLAKE8_CHECK = ${FLAKE8_SRC:%=src/flair/%} ${FLAKE8_TEST:%=test/%}
