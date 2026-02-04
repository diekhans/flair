root = ..
include ${root}/defs.mk

##
# the test-installed targets set
#   use_installed_flair=yes
# to not add local targets
##
ifneq (${use_installed_flair},yes)
    export PATH := ${root}/bin:${PATH}
    export PYTHONPATH := ${root}/src:${PYTHONPATH}
endif

ifeq (${PYTHONWARNINGS},)
   # set to error to find were they are coming from
   # FIXME: gets lots of warnings by not closing files, and these break some
   # expect results.  Disable for now.
   #   ResourceWarning: Enable tracemalloc to get the object allocation traceback
   #export PYTHONWARNINGS=always
endif


##
# all targets
##
default:
	@echo To run all tests, use: >&2
	@echo "   " make test >&2
	@echo To run tests without R depedencies, use: >&2
	@echo "   " make test-base >&2
	@echo To run tests using installed FLAIR rather than source tree, use: >&2
	@echo "   " make test-installed >&2
	@echo To run tests without R depedencies using installed FLAIR, use: >&2
	@echo "   " make test-base-installed >&2
	@echo otherwise run individual tests by specifying test target, such as: >&2
	@echo "   " make align-tests >&2
	@echo to installed FLAIR for individual tests use the use_installed_flair option, such as: >&2
	@echo "   " make align-tests use_installed_flair=yes >&2
	@echo remove output files: >&2
	@echo "   " make clean >&2
	@echo remove output and downloaded genomes files: >&2
	@echo "   " make real-clean >&2
	@echo run flake8: >&2
	@echo "   " make flake8 >&2
	@echo "   " make lint >&2
	@exit 1

test: test-base test-with-R
test-base: align-tests quantify-tests \
	combine-tests predict-productivity-tests diff-iso-usage-tests \
	plot-usage-tests transcriptome-tests fusion-tests variants-tests \
	partition-tests lib-tests

test-with-R: diffexp-tests diffsplice-tests


test-installed:
	${MAKE} test use_installed_flair=yes

test-base-installed:
	${MAKE} test-base use_installed_flair=yes

test-help: test-main-help test-align-help \
           test-quantify-help test-diffexp-help test-diffsplice-help \
           test-predict-productivity-help test-diff-iso-usage-help \
	   test-plot-usage-help test-transcriptome-help test-variants-help \
	   test-fusion-help

lint: flake8
flake8:
	cd .. && ${MAKE} flake8

##
# Note: expected and output are sorted in cases where the order is not stable
##

I = input
E = expected
O = output
TI = tmp-input

# reference data
GENOME_FA = ${TI}/genome.fa

##
# basic: general tests that programs seem to work
##
BASIC_READS_FA = $I/basic.reads.fa
BASIC_ANNOTATION = $I/basic.annotation.gtf
BASIC_ANNOTATION_INCOMPLETE = $I/basic.annotation.incomplete.gtf
BASIC_SHORTREAD_JUNCTIONS = $I/basic.shortread_junctions.tab
BASIC_PROMOTER_REGIONS = $I/basic.promoter_regions.bed
BASIC_READS_MANIFEST = $I/basic.reads_manifest.txt

# FIXME: test are based on old test data, however flair quantify
# does not produce sufficient pvalues to work diffexpr
BASIC_COUNTS_MATRIX = $I/basic.counts_matrix_diffexp.tsv
BASIC_DIFFSPLICE_MATRIX = $I/basic.counts_matrix_diffsplice.tsv
BASIC_DIFFSPLICE_BED = $I/test-collapse-annot.isoforms.bed

# FIXME: recreate these files with transcriptome, and others below 
BASIC_ISOFORMS_FA = $E/test-collapse-annot.isoforms.fa
BASIC_ISOFORMS_BED = $E/test-collapse-annot.isoforms.bed

##
# Files generate by tests that are used by other tests.
# These enfoce files being create in the right order.
## 
TEST_ALIGN_BAM = $O/test-align.bam
TEST_QUANTIFY_B1_BAM = $O/test-quantify.B1.B.flair.aligned.bam

###
# seg1: single exon tests
###
SEG1_READS_FQ = input/seg1.cdna-ont.fastq
SEG1_ANNOTATION = input/seg1.gencodeV47.gtf
SEG1_PROMOTER_REGIONS = input/seg1.promoter-regions.bed

SEG1_ALIGN_BED = $E/test-align-seg1.bed
SEG1_ALIGN_BAM = $O/test-align-seg1.bam
SEG1_ISOFORMS_FA = expected/test-collapse-seg1.isoforms.fa
SEG1_ISOFORMS_BED = expected/test-collapse-seg1.isoforms.bed


# The help test results vary on screen width.  This forces them to be consistent.
# It might seem weird to see if help works, but some it did break in the past.
export COLUMNS = 80

##
# compare output, expected with result
##
define cmp_expect
	diff $E/$(1) $O/$(1)
endef

##
# Expected and output are sorted in cases where the order is not stable.  In
# this case, we first test for the file as sorted input subsitution failure
# isn't propigate by shell.
##
define cmp_expect_sort
	@test -e $O/${1} || (echo Error: output file missing: $O/${1}; exit 1)
	diff <(sort $E/${1}) <(sort $O/${1})
endef

###
# MAIN COMMAND
###

test-main-help: mkdirs
	flair --help | fgrep -v Version> $O/$@.out
	$(call cmp_expect,$@.out)

###
# FLAIR ALIGN
###

align-tests: test-align test-align-intpriming test-align-seg1 test-align-separate test-align-help

test-align : mkdirs ${GENOME_FA}
	flair align -r $(BASIC_READS_FA) --genome ${GENOME_FA} -t 4 -o $O/$@
	$(call cmp_expect,$@.bed)

${TEST_ALIGN_BAM}: test-align
${SEG1_ALIGN_BAM}: test-align-seg1

test-align-separate : mkdirs ${GENOME_FA}
	flair align -r $(BASIC_READS_FA) --genome ${GENOME_FA} -t 4 -o $O/$@ --filtertype separate
	$(call cmp_expect,$@.bed)

test-align-intpriming : mkdirs ${GENOME_FA}
	flair align -r $(BASIC_READS_FA) --genome ${GENOME_FA} -t 4 -o $O/$@ --remove_internal_priming -f $(BASIC_ANNOTATION) --remove_singleexon
	$(call cmp_expect,$@.bed)

test-align-seg1 : mkdirs ${GENOME_FA}
	flair align -r ${SEG1_READS_FQ} --genome ${GENOME_FA} -t 4 -o $O/$@
	$(call cmp_expect,$@.bed)

test-align-help: mkdirs
	flair align --help > $O/$@.out
	$(call cmp_expect,$@.out)

# macros to do comparions with expected results
# FIXME: make sort version go away by making results stable
define diff_bed_gtf_sort
	$(call cmp_expect_sort,$@.isoforms.bed)
	$(call cmp_expect_sort,$@.isoforms.gtf)
endef
define diff_bed_gtf
	$(call cmp_expect,$@.isoforms.bed)
	$(call cmp_expect,$@.isoforms.gtf)
endef

define collapse_promoter_diff
	$(call cmp_expect_sort,$@.promoter_supported.bed)
endef

# resulting read map differs between run (see #443)
# just check if it exists for now
# should do: $(call cmp_expect_sort,$@.isoform.read.map.txt)
define diff_bed_gtf_sort_read_map
	@test -e $O/${1} || (echo Error: output file missing: $O/${1}; exit 1)
endef

###
# FLAIR COMBINE
###

define combine_diff
	$(call cmp_expect,$@/combined.bed)
	$(call cmp_expect,$@/combined.counts.tsv)
	$(call cmp_expect,$@/combined.isoform.map.txt)
endef
define combine_fa_diff
	$(call combine_diff)
	$(call cmp_expect,$@/combined.fa)
endef

combine-tests: test-combine-bed-only test-combine-bed-fa test-combine-bed-fa-map

update-combine-tests:
	for tname in "test-combine-bed-only" "test-combine-bed-fa" "test-combine-bed-fa-map"; \
	do \
	cp $O/$$tname/combined.bed $E/$$tname/combined.bed; \
	cp $O/$$tname/combined.counts.tsv $E/$$tname/combined.counts.tsv; \
	cp $O/$$tname/combined.isoform.map.txt $E/$$tname/combined.isoform.map.txt; \
	done
	for tname in "test-combine-bed-fa" "test-combine-bed-fa-map"; \
	do \
	cp $O/$$tname/combined.fa $E/$$tname/combined.fa; \
	done

test-combine-bed-only: mkdirs
	rm -rf $O/$@
	mkdir -p $O/$@
	flair combine --manifest=input/bed-only.combine_manifest.txt -o $O/$@/combined
	$(call combine_diff)

test-combine-bed-fa: mkdirs
	rm -rf $O/$@
	mkdir -p $O/$@
	flair combine --manifest=input/bed-fa.combine_manifest.txt -o $O/$@/combined
	$(call combine_fa_diff)

test-combine-bed-fa-map: mkdirs
	rm -rf $O/$@
	mkdir -p $O/$@
	flair combine --manifest=input/bed-fa-map.combine_manifest.txt -o $O/$@/combined
	$(call combine_fa_diff)

###
# FLAIR QUANTIFY
###
quantify-tests: test-quantify test-quantify-nobam test-quantify-help

test-quantify: mkdirs
	rm -rf $O/$@.tmp
	flair quantify -r $(BASIC_READS_MANIFEST) -i $(BASIC_ISOFORMS_FA) --isoform_bed $(BASIC_ISOFORMS_BED) --generate_map --temp_dir $O/$@.tmp --tpm --sample_id_only -o $O/$@ --output_bam
	samtools view $O/$@.A1.A.flair.aligned.bam | cut -f 1,3 > $O/$@.A1.A.flair.aligned.readtotranscript.txt
	$(call cmp_expect,$@.counts.tsv)
	$(call cmp_expect_sort,$@.A1.A.flair.aligned.readtotranscript.txt)
	$(call cmp_expect,$@.tpm.tsv)

${TEST_QUANTIFY_B1_BAM}: test-quantify

# without map or bam outputs
test-quantify-nobam : mkdirs
	rm -rf $O/$@.tmp
	flair quantify -r $(BASIC_READS_MANIFEST) -i $(BASIC_ISOFORMS_FA) --isoform_bed $(BASIC_ISOFORMS_BED) --temp_dir $O/$@.tmp --tpm --sample_id_only -o $O/$@
	$(call cmp_expect,$@.counts.tsv)
	$(call cmp_expect,$@.tpm.tsv)

test-quantify-help: mkdirs
	flair quantify --help > $O/$@.out
	$(call cmp_expect,$@.out)

###
# FLAIR DIFFEXP
###
diffexp-tests: test-diffexp test-diffexp-help

test-diffexp : mkdirs
	flair diffexp -q $(BASIC_COUNTS_MATRIX) -o $O/$@ -e 1 -of
	$(call cmp_expect_sort,$@/genes_deseq2_MCF7_v_A549.tsv)
	$(call cmp_expect_sort,$@/isoforms_deseq2_MCF7_v_A549.tsv)
	$(call cmp_expect_sort,$@/isoforms_drimseq_MCF7_v_A549.tsv)

test-diffexp-help: mkdirs
	flair diffexp --help > $O/$@.out
	$(call cmp_expect,$@.out)

###
# FLAIR DIFFSPLICE
###
diffsplice-tests: test-diffsplice test-diffsplice-help

# various output files not long produce due to changes in test cases
test-diffsplice: mkdirs
	flair diffsplice -i $(BASIC_DIFFSPLICE_BED) -q $(BASIC_DIFFSPLICE_MATRIX) --test -o $O/$@ -of
	$(call cmp_expect_sort,$@/diffsplice.alt3.events.quant.tsv)
	$(call cmp_expect_sort,$@/diffsplice.alt5.events.quant.tsv)
	$(call cmp_expect,$@/diffsplice.es.events.quant.tsv)
	$(call cmp_expect_sort,$@/diffsplice.ir.events.quant.tsv)
	$(call cmp_expect,$@/drimseq_alt3_A_v_B.tsv)
	$(call cmp_expect,$@/drimseq_es_A_v_B.tsv)

test-diffsplice-cond: mkdirs
	flair diffsplice -i $(BASIC_ISOFORMS_BED) -q $(BASIC_DIFFSPLICE_MATRIX) --test -o $O/$@ -of --conditionA A --conditionB B
	$(call cmp_expect_sort,test-diffsplice/diffsplice.alt3.events.quant.tsv)
	$(call cmp_expect_sort,test-diffsplice/diffsplice.alt5.events.quant.tsv)
	$(call cmp_expect,test-diffsplice/diffsplice.es.events.quant.tsv)
	$(call cmp_expect_sort,test-diffsplice/diffsplice.ir.events.quant.tsv)
	$(call cmp_expect,test-diffsplice/drimseq_alt3_A_v_B.tsv)
	$(call cmp_expect,test-diffsplice/drimseq_es_A_v_B.tsv)

test-diffsplice-help: mkdirs
	flair diffsplice --help > $O/$@.out
	$(call cmp_expect,$@.out)


##
# PREDICT PRODUCTIVITY
##
predict-productivity-tests: test-predict-productivity test-predict-productivity-help

test-predict-productivity : mkdirs ${GENOME_FA}
	predictProductivity -i $(BASIC_ISOFORMS_BED) -g $(BASIC_ANNOTATION) -f ${GENOME_FA} -o $O/$@ --longestORF
	$(call cmp_expect,$@.bed)
	$(call cmp_expect,$@.info.tsv)

test-predict-productivity-help: mkdirs
	predictProductivity --help > $O/$@.out
	$(call cmp_expect,$@.out)

##
# DIFF ISO USAGE
##
diff-iso-usage-tests: test-diff-iso-usage test-diff-iso-usage-help

test-diff-iso-usage: mkdirs
	diff_iso_usage $(BASIC_COUNTS_MATRIX) MCF7cdnar1r3_MCF7_batch1 MCF7cdnastrandr3r2_MCF7_batch1 $O/$@.tsv
	$(call cmp_expect,$@.tsv)

test-diff-iso-usage-help: mkdirs
	diff_iso_usage --help > $O/$@.out
	$(call cmp_expect,$@.out)

##
# plot isoform usage
##
plot-usage-tests: test-plot-usage test-plot-usage-help

test-plot-usage: mkdirs
	plot_isoform_usage expected/test-predict-productivity.bed expected/test-collapse-annot.annotated_transcripts.alignment.counts ENSG00000133703.12 -o $O/$@
	file $O/test-plot-usage_isoforms.png | grep 'PNG image data' >/dev/null
	file $O/test-plot-usage_usage.png | grep 'PNG image data' >/dev/null

test-plot-usage-help: mkdirs
	plot_isoform_usage --help > $O/$@.out
	$(call cmp_expect,$@.out)

###
# flair transcriptome
###

transcriptome-tests: test-transcriptome test-transcriptome-partition test-transcriptome-noannot test-transcriptome-predcds \
	test-transcriptome-seg1 test-transcriptome-ends test-transcriptome-help

update-transcriptome-tests:
	for tname in "test-transcriptome" "test-transcriptome-partition" "test-transcriptome-noannot" "test-transcriptome-predcds" "test-transcriptome-seg1" "test-transcriptome-ends"; \
	do \
	cp $O/$$tname.isoforms.bed $E/$$tname.isoforms.bed; \
	cp $O/$$tname.isoforms.gtf $E/$$tname.isoforms.gtf; \
	cp $O/$$tname.isoform.read.map.txt $E/$$tname.isoform.read.map.txt; \
	done


test-transcriptome: mkdirs ${TEST_ALIGN_BAM}
	flair transcriptome  -b ${TEST_ALIGN_BAM} -g ${GENOME_FA} --junction_tab ${BASIC_SHORTREAD_JUNCTIONS} -t 4 -f $(BASIC_ANNOTATION) -o $O/$@ --keep_intermediate
	$(diff_bed_gtf)

test-transcriptome-ends: mkdirs ${TEST_ALIGN_BAM}
	flair transcriptome  -b ${TEST_ALIGN_BAM} -g ${GENOME_FA} --junction_tab ${BASIC_SHORTREAD_JUNCTIONS} -t 4 -f $(BASIC_ANNOTATION) -o $O/$@ --end_norm_dist 500 --output_endpos
	$(diff_bed_gtf)

test-transcriptome-partition: mkdirs ${TEST_ALIGN_BAM}
	flair transcriptome  -b ${TEST_ALIGN_BAM} -g ${GENOME_FA} --junction_tab ${BASIC_SHORTREAD_JUNCTIONS} -t 4 -f $(BASIC_ANNOTATION) -o $O/$@ --parallel_mode byregion
	$(diff_bed_gtf)

test-transcriptome-noannot: mkdirs ${TEST_ALIGN_BAM}
	flair transcriptome  -b ${TEST_ALIGN_BAM} -g ${GENOME_FA} -t 4 -f $(BASIC_ANNOTATION) -o $O/$@ --no_align_to_annot
	$(diff_bed_gtf)

test-transcriptome-predcds: mkdirs ${TEST_ALIGN_BAM}
	flair transcriptome  -b ${TEST_ALIGN_BAM} -g ${GENOME_FA} -t 4 -f $(BASIC_ANNOTATION) -o $O/$@ --predict_cds
	$(diff_bed_gtf)
	$(call cmp_expect,$@.isoforms.CDS.bed)

# FIXME: doesn't work, returns nothing. Need to debug later
test-transcriptome-seg1: mkdirs  ${SEG1_ALIGN_BAM}
	flair transcriptome -b ${SEG1_ALIGN_BAM} -g ${GENOME_FA} -t 4 -f $(SEG1_ANNOTATION) -o $O/$@
	$(diff_bed_gtf)

test-transcriptome-help: mkdirs
	flair transcriptome --help > $O/$@.out
	$(call cmp_expect,$@.out)


###
# flair variants
# FIXME need to add better checks for AAseq results, however the inconsistent sorting of the variants leads to inconsistent AAseq IDs which are hard to check
###
variants-tests:  test-variants test-variants-help

# longshot is not currently include in the docker container
have_longshot := $(shell which longshot >/dev/null 2>&1 && echo yes || echo no)

ifeq (${have_longshot},yes)
test-variants: mkdirs ${GENOME_FA} ${TEST_QUANTIFY_B1_BAM}
	minimap2 -ax map-ont --secondary=no ${GENOME_FA} input/basic.reads.4.fa | samtools view -hb - | samtools sort - > $O/$@.genomealigned.bam
	samtools index $O/$@.genomealigned.bam
	samtools faidx $(BASIC_ISOFORMS_FA)
	longshot --force_overwrite --bam $O/$@.genomealigned.bam --ref ${GENOME_FA} --out $O/$@.genomealigned.longshot.vcf --min_cov 3 --min_alt_count 3 --strand_bias_pvalue_cutoff 0.000001
	echo 's1    ${TEST_QUANTIFY_B1_BAM}   $O/$@.genomealigned.longshot.vcf' > $O/$@.varmanifest.tsv
	flair variants -m $O/$@.varmanifest.tsv -i $(BASIC_ISOFORMS_FA) -b $(BASIC_ISOFORMS_BED) -g ${GENOME_FA} -f $(BASIC_ANNOTATION) -o $O/$@
	$(call cmp_expect,$@.isoswithvars.counts.tsv)
	$(call cmp_expect,$@.aaseq.counts.tsv)
else
test-variants:
	@echo 'WARNING: ***' >&2
	@echo 'WARNING: *** longshot command not found, test-variants skipped' >&2
	@echo 'WARNING: ***' >&2
endif

test-variants-help: mkdirs
	flair variants --help > $O/$@.out
	$(call cmp_expect,$@.out)


###
# flair fusion
###
ERGIC3_KRAS_SIM_FA_GZ = $O/ERGIC3_KRAS_sim_fusion.badread20x.stranded.fastq.gz
FUSION_TEST_SHARED_OUT = $O/fusion_tests
FUSION_TEST_BAM = ${FUSION_TEST_SHARED_OUT}_chimeric.bam

fusion-tests: test-fusion test-fusion-pary test-fusion-help

test-fusion: mkdirs ${GENOME_FA} ${ERGIC3_KRAS_SIM_FA_GZ} ${FUSION_TEST_BAM}
	flair fusion -g ${GENOME_FA} -f $(BASIC_ANNOTATION) -r ${ERGIC3_KRAS_SIM_FA_GZ} -b ${FUSION_TEST_BAM} --maxloci 4 -o $O/$@ -s 2
	diff <(sort <(cut -f 1,2,3,5,6,7,8,9,10,11,12 $E/$@.fusions.isoforms.bed)) <(sort <(cut -f 1,2,3,5,6,7,8,9,10,11,12 $O/$@.fusions.isoforms.bed))
	diff <(sort <(grep -v ">" $E/$@.fusions.isoforms.fa)) <(sort <(grep -v ">" $O/$@.fusions.isoforms.fa))

test-fusion-pary: mkdirs ${GENOME_FA} ${ERGIC3_KRAS_SIM_FA_GZ} ${FUSION_TEST_BAM}
	flair fusion -g ${GENOME_FA} -f $I/parY.annotation.gtf -r ${ERGIC3_KRAS_SIM_FA_GZ} -b ${FUSION_TEST_BAM} --maxloci 4 -o $O/$@ -s 2
	$(call cmp_expect_sort,$@.fusions.isoforms.bed)
	$(call cmp_expect_sort,$@.fusions.isoforms.fa)

${FUSION_TEST_BAM}: mkdirs  ${ERGIC3_KRAS_SIM_FA_GZ}
	flair align -r $I/ERGIC3_KRAS_sim_fusion.badread20x.stranded.fastq \
	   --genome ${GENOME_FA} -t 4 -o ${FUSION_TEST_SHARED_OUT} \
	   --filtertype separate --minfragmentsize 40 --maxintronlen 350k


test-fusion-help: mkdirs
	flair fusion  --help > $O/$@.out
	$(call cmp_expect,$@.out)

${ERGIC3_KRAS_SIM_FA_GZ}: $I/ERGIC3_KRAS_sim_fusion.badread20x.stranded.fastq
	@mkdir -p $(dir $@)
	gzip -c $< >$@


###
# flair spliceevents
###
test-spliceevents: mkdirs
	python3 ../src/flair/flair_spliceevents.py -m $I/spliceevents.manifest.tsv -g ${GENOME_FA} -t 4 -f $(BASIC_ANNOTATION) -o $O/$@


##
# flair_partition command
##
partition-tests: test-partition-basic test-partition-adjacent test-partition-min-reads \
	test-partition-sam

test-partition-basic: mkdirs
	flair_partition --threads=2 --log-level=info --bed=input/partition-test.bed output/$@.ranges.bed >& output/$@.out
	$(call cmp_expect,$@.ranges.bed)
	$(call cmp_expect,$@.out)

# merge adjacent partitions
test-partition-adjacent: mkdirs
	flair_partition --threads=2 -part_merge_dist=10000 --log-level=info --bed=input/partition-test.bed output/$@.ranges.bed >& output/$@.out
	$(call cmp_expect,$@.ranges.bed)
	$(call cmp_expect,$@.out)

# minimum reads per partition
test-partition-min-reads: mkdirs
	flair_partition --threads=2 --min_partition_items=3  --log-level=info --bed=input/partition-test.bed output/$@.ranges.bed >& output/$@.out
	$(call cmp_expect,$@.ranges.bed)
	$(call cmp_expect,$@.out)

# with bed, sam, and gtf
test-partition-sam: mkdirs
	flair_partition --threads=2 --log-level=info --bed=input/partition-test.bed --bam=input/tiny.sam --gtf=input/seg1.gencodeV47.gtf \
	    output/$@.ranges.bed >& output/$@.out
	$(call cmp_expect,$@.ranges.bed)
	$(call cmp_expect,$@.out)

##
# library tests
##
lib-tests: mkdirs
	${PYTEST} --capture=no -vv


##
# build test input files
##
${GENOME_FA}:
	mkdir -p $(dir $@)
	rm -f $@.fai
	wget -nv https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr12.fa.gz
	wget -nv https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr17.fa.gz
	wget -nv https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr20.fa.gz
	zcat chr12.fa.gz chr17.fa.gz chr20.fa.gz > $@
	rm chr12.fa.gz chr17.fa.gz chr20.fa.gz
	samtools faidx $@

mkdirs:
	@mkdir -p $O

clean:
	rm -rf $O $D __pycache__
	rm -rf [0-9a-f]*-[0-9a-f]*-[0-9a-f]*-[0-9a-f]*-[0-9a-f]*

real-clean : clean
	rm -rf ${TI}
