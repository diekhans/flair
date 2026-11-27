.. _migration-3.0-3.1-label:

Migrating from FLAIR 3.0 to 3.1
===============================

What a command line written for FLAIR 3.0 has to change to run under FLAIR 3.1.
Each table gives what 3.0 had, what 3.1 has, and what to know about the change.
An empty first column is new in 3.1.  Options and outputs that are not listed are
unchanged.

An old name still given on a command line is reported with the name that replaced
it, rather than as unrecognized: ``flair transcriptome --out_dir out`` reports
``--out_dir is now --output, a file name base rather than a directory``, and
``flair collapse`` reports that ``flair transcriptome`` builds transcript models
from the aligned reads.  The old names are not accepted, and are not listed in any
``--help``.

All commands
------------

.. list-table::
   :header-rows: 1
   :widths: 30 30 40

   * - 3.0
     - 3.1
     - Comment
   * - ``--genomealignedbam``, ``--conditionA`` and other run-together or camelCase
       names
     - ``--genome_aligned_bam``, ``--condition_a``
     - Option names are snake_case.  Each one is listed under its command below.
   * -
     - ``-z``, ``--compress``
     - gzips the text files a run writes, adding ``.gz`` to their names, and bgzips
       FASTA so it can be indexed.  A global option, so it comes before the
       subcommand: ``flair -z transcriptome ...``.  Compressed input is read whether
       or not it is given.
   * - ``flair.read.map.txt`` and the other ``.txt`` tables
     - ``flair.isoforms.read.map.tsv``
     - The TSVs FLAIR writes are named ``.tsv``.
   * - ``<isoform>_<gene>`` names
     - ``FLT`` isoform, ``FLG`` gene and ``FLP`` protein identifiers
     - The gene is a column of its own rather than part of the isoform name.  See
       :ref:`ids-label`.
   * - ``--quiet``
     - ``--log-level``, ``--log-stderr``, ``--log-conf``, ``--log-debug``
     - Logging is set by these global options.

Commands replaced or removed
----------------------------

.. list-table::
   :header-rows: 1
   :widths: 30 30 40

   * - 3.0
     - 3.1
     - Comment
   * - ``flair correct``, ``flair collapse``
     - ``flair transcriptome``
     - One command does both, starting from a genome aligned BAM, with the junction
       evidence given by ``--gtf``, ``--junction_bed`` or ``--junction_star``.
   * - ``flair align``
     - gone
     - It ran minimap2 and filtered its output, which nothing else used.  Align the
       reads with the minimap2 command under this table.
   * - ``flair variants``
     - ``flair isoalleles``
     - Allele calling from reads is replaced by haplotyping with an upstream tool,
       such as WhatsHap, whose haplotagged reads ``flair isoalleles`` takes.  See
       :ref:`isoalleles-label`.

``flair transcriptome`` and ``flair fusion`` take a BAM of all the reads aligned to
the genome:

.. code:: text

   minimap2 -ax splice --MD --secondary=no -t 4 genome.fa reads.fq | samtools sort -o sample.bam -
   samtools index sample.bam

:ref:`alignment-label` states what each option does and when to change it.

flair transcriptome
-------------------

.. list-table::
   :header-rows: 1
   :widths: 30 30 40

   * - 3.0
     - 3.1
     - Comment
   * - ``--genomealignedbam``
     - ``-b``, ``--genome_aligned_bam``
     -
   * - ``--keep_sup``
     - ``--keep_supplementary``
     -
   * - ``--noaligntoannot``
     - ``--no_align_to_annot``
     -
   * - ``--junction_tab``
     - ``--junction_star``
     - The old name still works and warns.
   * - ``-s``, ``--support``
     - ``--sjc_support``, ``--single_exon_support``, ``--frac_support``
     - One option became three: spliced isoforms (default 1), single exon isoforms
       (default 3), and the fraction of the gene locus an isoform must make up
       (default 0.05).
   * - ``--no_stringent``, ``--no_check_splice``
     - gone
     - ``--stringent`` and ``--check_splice`` always apply: turning them off counted
       reads that do not support an isoform toward it.
   * - ``-n``, ``--no_redundant``
     - gone
     -
   * - ``--predictCDS``, ``--longestORF``
     - gone
     - An isoform's CDS and productivity are reported whenever the annotation given
       with ``-f`` has start codons.
   * - ``--annotation_reliant``, ``--annotated_bed``, ``--transcriptomefasta``
     - gone
     - Reads are aligned to the transcripts of ``-f`` by default, with no FASTA or
       BED of them to build or pass; ``--no_align_to_annot`` turns that alignment
       off.
   * - ``--remove_internal_priming``, ``--intprimingthreshold``,
       ``--intprimingfracAs``
     - gone
     - Single-exon reads whose 3' end is in A-rich genome, 60% A in the 20 bases
       after it or 80% A in its last 20 aligned bases, are always removed as
       oligo-dT mispriming; ``--keep_internally_primed`` keeps them.
   * - ``--nvrna``
     - gone
     - Correcting a read's strand against the annotation is what correction does;
       ``--trust_strand`` and ``--directRNA`` keep the alignment's strand instead.
   * - ``--promoters``, ``--3prime_regions``, ``--gene_tss``,
       ``--no_gtf_end_adjustment``
     - gone
     - An isoform's ends come from the ends of the reads assigned to it, clustered
       within ``--end_window``, at most ``--max_ends`` of them, or
       ``--normalize_ends``.  No promoter or transcription end site file, no
       gene-level ends, and the annotation's ends are used to keep subset isoforms
       rather than to move an isoform's ends.
   * - ``--mm2_args``
     - gone
     - The minimap2 runs are internal.
   * - ``--longshot_bam``, ``--longshot_vcf``
     - gone
     - Haplotypes come from ``flair isoalleles``.
   * - ``--allow_paralogs``
     - gone
     -
   * - ``--parallelmode``, ``--parallel_mode``
     - gone
     - The input is always partitioned into clusters of overlapping reads and
       annotated transcripts.
   * - ``--trust_ends``
     - gone
     - Each isoform's ends come from the ends of the reads assigned to it, and
       passing it is an error that reports that.
   * -
     - ``--sample_name``
     - Required, and recorded in the output files.
   * -
     - ``--trust_strand``, ``--trust_junctions``, ``--normalize_ends``
     -
   * -
     - ``--directRNA``, ``--keep_internally_primed``
     - For reads not primed on poly(A), where a single-exon read ending in A-rich
       sequence is a real transcript end rather than oligo-dT mispriming.
   * -
     - ``--region``, ``--region_bed``
     - Build transcripts in the given regions only, each grown to the whole cluster
       of reads and transcripts that overlap it.
   * -
     - ``--temp_dir``
     - ``--keep_intermediate`` now requires it.
   * -
     - ``--generate_map``, ``--fusion_breakpoints``, ``--subset_backup_support``
     -

Outputs:

.. list-table::
   :header-rows: 1
   :widths: 30 30 40

   * - 3.0
     - 3.1
     - Comment
   * - ``flair.isoforms.bed``, ``flair.isoforms.gtf``, ``flair.isoforms.fa``
     - the same names
     -
   * - ``flair.read.map.txt``
     - ``flair.isoforms.read.map.tsv``
     - Written with ``--generate_map``.

flair quantify
--------------

.. list-table::
   :header-rows: 1
   :widths: 30 30 40

   * - 3.0
     - 3.1
     - Comment
   * - ``-r``, ``--reads_manifest``
     - ``--manifest``
     - The manifest names each sample's reads as a sorted, indexed BAM aligned to
       the genome, where it named a reads file before.  Its first three columns,
       sample id, condition and batch, are unchanged.
   * - ``--isoformbed``
     - ``--isoform_bed``
     - Given with ``-g``, ``--genome``.
   * - ``-i``, ``--isoforms``
     - gone
     - The isoforms are given as ``--isoform_bed``, not as a FASTA.
   * - ``--output_bam``, ``--sample_id_only``, ``--quality``
     - gone
     -
   * - ``--tpm``
     - gone
     - The table it asked for is always written, as ``<output>.isoforms.cpm.tsv``.
   * - ``--stringent``, ``--check_splice``
     - gone
     - Both requirements always apply to read assignment, so there is nothing to
       turn on.
   * - ``--with_gene``
     - gone
     - The counts matrix always names the gene, in its own column.
   * -
     - ``--normalize_ends``
     -
   * -
     - ``--temp_dir``
     - The directory a run makes its own temporary directory in, by default
       ``$TMPDIR``; ``--keep_intermediate`` keeps it.

Outputs:

.. list-table::
   :header-rows: 1
   :widths: 30 30 40

   * - 3.0
     - 3.1
     - Comment
   * - ``<output>.counts.tsv`` with one ``ids`` column
     - ``<output>.isoforms.counts.tsv`` with ``gene_id`` and ``isoform_id``
     - The one column joined the two.
   * - sample columns named ``sample_condition_batch``
     - ``<output>.isoforms.sample_info.tsv``
     - The column names contain the sample id alone.  ``flair diffexp`` and ``flair
       diffsplice`` read the sample info file.
   * - the per sample read maps
     - ``<output>.isoforms.<sample>.read.map.tsv``
     -
   * - ``<output>.tpm.tsv``
     - ``<output>.isoforms.cpm.tsv``
     - Written by every run.  It contained counts per million under the name TPM,
       and is now named for what it contains: a long read is not one transcript,
       since reads are fragmented, so the counts are not normalized by isoform
       length.

flair combine
-------------

.. list-table::
   :header-rows: 1
   :widths: 30 30 40

   * - 3.0
     - 3.1
     - Comment
   * - ``--output_prefix``
     - ``-o``, ``--output``
     -
   * - ``--endwindow``
     - ``--end_window``
     -
   * - ``--minpercentusage``
     - ``--frac_support``
     - A fraction of the gene's reads where it was a percentage of them:
       ``--minpercentusage 5`` is ``--frac_support 0.05``.
   * - ``--include_se``
     - ``--remove_single_exon``
     - Single exon isoforms are kept by default, and this drops them.
   * - ``-c``, ``--convert_gtf``
     - gone
     - The GTF is always written.
   * - ``-f``, ``--filter``
     - gone
     - ``-f`` is now ``--gtf``, below.
   * - the manifest's sample name, type, FASTA and read map columns
     - one isoform BED path per line
     -
   * -
     - ``--prefixes``, ``--isoform_beds``
     - Name the transcriptomes to combine without a manifest file; they and
       ``--manifest`` are mutually exclusive.
   * -
     - ``-g``, ``--genome``
     -
   * -
     - ``-f``, ``--gtf``
     - With it, each combined isoform takes its gene, reference transcript, CDS,
       productivity and protein from that annotation, so samples run with different
       annotations can be combined.
   * -
     - ``--max_ends``, ``--min_reads``
     -

flair fusion
------------

.. list-table::
   :header-rows: 1
   :widths: 30 30 40

   * - 3.0
     - 3.1
     - Comment
   * - ``--genomechimbam``
     - ``-b``, ``--genome_aligned_bam``
     -
   * - ``--maxloci``
     - ``--max_loci``
     -
   * - ``--minfragmentsize``
     - ``--min_fragment_size``
     -
   * - ``-s``, ``--support``
     - ``--min_support``
     -
   * - ``--stringent``, ``--check_splice``
     - gone
     - Both always apply; neither can be turned off.
   * - ``-r``, ``--reads``, ``--transcriptchimbam``
     - gone
     -
   * -
     - ``--sample_name``
     - Required.
   * -
     - ``--max_dist_to_tss``, ``--min_dist_between_breakpoints``,
       ``--max_partners``
     -
   * -
     - ``--frac_support``, ``--quality``, ``--no_transcriptome_alignment``
     -
   * -
     - ``--keep_intermediate``, ``--temp_dir``
     - The intermediate files are in a directory of the run's own rather than beside
       the output.

Outputs:

.. list-table::
   :header-rows: 1
   :widths: 30 30 40

   * - 3.0
     - 3.1
     - Comment
   * -
     - ``OUTPUT.fusions.tsv``
     - A table of the fusions called: their reads, isoforms and best supported
       breakpoint.
   * -
     - ``OUTPUT.fusions.isoforms.gtf``, ``OUTPUT.fusions.isoforms.bb``
     - Beside the BED, as flair transcriptome writes them.

flair diffexp
-------------

.. list-table::
   :header-rows: 1
   :widths: 30 30 40

   * - 3.0
     - 3.1
     - Comment
   * - ``--out_dir``
     - ``-o``, ``--output``
     -
   * - ``-of``, ``--out_dir_force``
     - ``--overwrite_output``
     -
   * - ``-e``, ``--exp_thresh``
     - ``--min_expression``
     -
   * -
     - ``--condition_a``, ``--condition_b``
     - Name the two conditions to compare.  Without them the two are taken in sorted
       order rather than from the column order of the counts matrix.

Outputs:

.. list-table::
   :header-rows: 1
   :widths: 30 30 40

   * - 3.0
     - 3.1
     - Comment
   * -
     - ``isoforms_sig_exp_change_norm_by_gene.tsv``
     - Isoforms whose usage changes, from a t-test on counts normalized by gene.

flair diffsplice
----------------

.. list-table::
   :header-rows: 1
   :widths: 30 30 40

   * - 3.0
     - 3.1
     - Comment
   * - ``-i``, ``--isoforms``
     - ``--isoform_bed``
     -
   * - ``--outDir``
     - ``-o``, ``--output``
     -
   * - ``-of``, ``--out_dir_force``
     - ``--overwrite_output``
     -
   * - ``--conditionA``, ``--conditionB``
     - ``--condition_a``, ``--condition_b``
     -
   * - ``--drim1``, ``--drim2``, ``--drim3``, ``--drim4``
     - ``--min_samps_gene_expr``, ``--min_samps_feature_expr``,
       ``--min_gene_expr``, ``--min_feature_expr``
     -
   * - ``inclusion_chr1:500-650``
     - ``inclusion_chr1:500-650(+)``
     - Event names and the ``coordinate`` column include the strand, as ``flair
       spliceevents`` names an event.  Without the strand, an event called on both
       strands at one locus appeared twice under one name, and DRIMSeq rejected the
       file.
   * - a ``drimseq_<event>_A_v_B.tsv`` for each event type tested
     - all four, always
     - An event type that called nothing, or whose events DRIMSeq found nothing
       significant in, leaves a file of only its header.

Standalone programs
-------------------

.. list-table::
   :header-rows: 1
   :widths: 30 30 40

   * - 3.0
     - 3.1
     - Comment
   * - ``bed_to_sequence``
     - gone
     - Use ``bedtools getfasta``.
   * - ``counts_to_tpm``
     - gone
     - ``flair quantify`` writes ``<output>.isoforms.cpm.tsv`` itself.
   * - ``normalize_counts_matrix``
     - gone
     -
   * - ``predictProductivity``
     - gone
     - ``flair transcriptome`` writes the productivity call in the isoform BED, and
       ``--predictCDS`` is gone with it.
   * - ``identify_annotated_gene``
     - gone
     -
   * - ``identify_gene_isoform``
     - gone
     - It named isoforms transcript_gene for ``bed_to_gtf`` to split, and FLAIR no
       longer splits names (:ref:`ids <ids-label>`).
   * - ``annotate_aaseq_with_uniprot``
     - gone
     -
   * - ``junctions_from_sam``
     - gone
     - It read only alignments without an indel or clipping, which is 4% of a
       minimap2 BAM, and took the strand from an ``XS:A:`` tag that minimap2 does
       not write, so its junctions had a fraction of the read support and no strand.
       Call junctions with `intron-prospector
       <https://github.com/diekhans/intron-prospector>`__ and pass them with
       ``--junction_bed``, or use STAR's ``SJ.out.tab`` with ``--junction_star``.
   * - ``bed_to_gtf --dontuseCDS``
     - ``bed_to_gtf --no_cds``
     - ``--add_reference_transcript_id`` is gone, and ``--is_flair_bed`` turns the
       extra BED columns into GTF attributes.  A BED not from FLAIR is rejected
       unless given ``--force``, which uses each name as its gene id too: names are
       no longer split into a transcript and gene at an underscore.  A thick range
       is a coding region whenever ``thickStart`` is below ``thickEnd``, as the BED
       convention has it, so a transcript with no UTR keeps its CDS.  FLAIR before
       3.1 wrote a thick range covering the whole transcript for a transcript with
       no coding region, which was read as having none; convert such a BED with
       ``--no_cds``.
   * - ``plot_isoform_usage --vcf``
     - gone
     -
   * - ``mark_intron_retention``, ``fasta_seq_lengths``,
       ``diffsplice_fishers_exact``
     - the same arguments
     - They parse them with argparse, so each has ``--help``.
   * - ``flair_partition``
     - the same arguments
     - It takes the ``--log-*`` options.
   * - ``gtf_to_bed``
     - unchanged
     -
   * - ``diff_iso_usage``
     - the same arguments
     - It writes three more columns, below.

File formats
------------

.. list-table::
   :header-rows: 1
   :widths: 30 30 40

   * - 3.0
     - 3.1
     - Comment
   * - one ``ID`` column of ``<isoform>_<gene>``
     - ``gene_id`` and ``isoform_id``
     - The counts matrix and the CPM table open with the two columns.  A counts
       matrix from an earlier FLAIR, whose sample columns are named
       ``sample_condition_batch``, is still read: the fields are taken back out of
       the column name.
   * - ``.isoform.read.map.txt``, no header
     - ``.isoforms.read.map.tsv``, header ``name`` and ``reads``
     -
   * - ``diff_iso_usage`` with ``geneID`` and ``isoID``
     - ``gene_id`` and ``isoform_id``, and three columns at the end
     - The three are ``<sample1>_PSI``, ``<sample2>_PSI`` and ``delta_PSI``.  A
       p-value or a PSI that cannot be taken is an empty field where it was ``NA``.
   * - the isoform BED of 3.0
     - BED12 plus 13 columns
     - In this order: ``gene_id``, ``ref_transcript_id``, ``ref_gene_mappings``,
       ``read_support``, ``frac_support``, ``productivity``, ``transcript_class``,
       ``fused_genes``, ``pos_in_fusion``, ``samples``, ``source_isoform``,
       ``allele_group``, ``aaseq_id``.  A column that does not apply to a row is
       empty.  Anything reading the BED by column number past column 12 has to be
       changed to the new column numbers.
   * - isoform name ``<isoform>_<gene>``
     - isoform name ``FLT00000001``
     - The isoform id alone, with the gene in the ``gene_id`` column beside it.
