#! /usr/bin/env python3

import argparse
import os
import pipettor
import shutil
import pysam
import logging
from collections import Counter, namedtuple
from flair import PosRange
from flair.flair_align import intron_chain_to_exon_starts
from flair.gtf_io import gtf_data_parser, gtf_write_row, GtfTranscript, GtfExon, GtfAttrsSet, TRANSCRIPT_EXON_FEATURES
from flair.intron_support import IntronSupport
from flair.junction_correct import JunctionCorrector
from flair.partition_runner import parallel_mode_parse, partition_runner_factory
from flair.pycbio.hgdata.bed import Bed, BedReader
from flair.bed_to_sequence import bed_to_sequence
from flair.bed_to_gtf import bed_to_gtf

# FIXME: temporarily disabled C901 (too complex) in .flake8
# FIXME: add object for all file names
# FIXME: use real TSVs
# FIXME: need to document all the files
# FIXME: it seems overkill to discard a single read based on one unsupported junction.
#        These can be recovered from other reads. Simple way is to add all reads valid
#        junctions, although maybe faster to track with introns are support and figure
#        this out in correct.  Another might have a polymorphic different with the reference
#        that changes splice junctions.  Should this be discarded if multiple long-reads
#        support it, but it isn't annotated.  Maybe these can be identified.

def get_args():
    parser = argparse.ArgumentParser(description='generates confident transcript models directly from a bam file '
                                                 'of aligned long rna-seq reads')
    parser.add_argument('-b', '--genome_aligned_bam', required=True,
                        help='Sorted and indexed bam file aligned to the genome')
    parser.add_argument('-g', '--genome', type=str, required=True,
                        help='FastA of reference genome, can be minimap2 indexed')
    parser.add_argument('-o', '--output', default='flair',
                        help='output file name base for FLAIR isoforms (default: flair)')
    parser.add_argument('-t', '--threads', type=int, default=12,
                        help='number of threads to run with - related to parallel_mode')
    parser.add_argument('-f', '--gtf', dest="annot_gtf", default=None,
                        help='GTF annotation file, used for identifying annotated isoforms')

    mutexc = parser.add_mutually_exclusive_group(required=False)
    mutexc.add_argument('--junction_tab', help='short-read junctions in SJ.out.tab format. '
                                               'Use this option if you aligned your short-reads with STAR, '
                                               'STAR will automatically output this file')
    mutexc.add_argument('--junction_bed', help='short-read junctions in bed format '
                                               '(can be generated from long-read alignment with intronProspector)')
    parser.add_argument('--junction_support', type=int, default=1,
                        help='if providing short-read junctions, minimum junction support required to keep junction. '
                             'If your junctions file is in bed format, the score field will be used for read support.')
    parser.add_argument('--ss_window', type=int, default=15,
                        help='window size for correcting splice sites (15)')
    parser.add_argument('-w', '--end_window', type=int, default=100,
                        help='window size for comparing TSS/TES (100)')

    parser.add_argument('--sjc_support', type=int, default=1,
                        help='''minimum number of supporting reads for a spliced isoform''')
    parser.add_argument('--se_support', type=int, default=3,
                        help='''minimum number of supporting reads for a single exon isoform''')
    parser.add_argument('--frac_support', type=float, default=0.05,
                        help='''minimum fraction of gene locus support for isoform to be called
                        default: 0.05, only isoforms that make up more than 5 percent of the gene
                        locus are reported. Set to 0 for max recall''')
    parser.add_argument('--no_stringent', default=False, action='store_true',
                        help='''specify if all supporting reads don't need to be full-length
                        (aligned to first and last exons of transcript).  Use this for fragmented libraries,
                        with an understanding that it will impact precision.''')
    parser.add_argument('--no_check_splice', default=False, action='store_true',
                        help='''don't enforce accurate alignment around splice site.
                        Specify this for libraries with high error rates, but it will reduce precision''')
    parser.add_argument('--no_align_to_annot', default=False, action='store_true',
                        help='''related to old annotation_reliant, now specify if you don't want
                        an initial alignment to the annotated sequences and only want transcript
                        detection from the genomic alignment.
                         Will be slightly faster but less accurate if the annotation is good''')
    parser.add_argument('-n', '--no_redundant', default='none',
                        help='''For each unique splice junction chain, report options include:
                        none-- multiple supported TSSs/TESs chosen for each set of splice junctions (modulated by max_ends);
                        longest--single TSS/TES chosen to maximize length;
                        best_only--single most supported TSS/TES used in conjunction chosen (none)''')
    parser.add_argument('--max_ends', type=int, default=1,
                        help='maximum number of TSS/TES picked per isoform (1) make higher for more precise end detection')
    parser.add_argument('--filter', default='nosubset',
                        help='''Report options include:
                        nosubset--any isoforms that are a proper set of another isoform are removed;
                        bysupport--subset isoforms are removed based on support;
                        comprehensive--default set + all subset isoforms;
                        ginormous--comprehensive set + single exon subset isoforms''')
    parser.add_argument('--quality', default=1, type=int,
                        help='minimum mapping quality threshold to consider genomic alignments for defining transcripts')
    parser.add_argument('--parallel_mode', default='auto:1GB',
                        help='''parallelization mode. Default: "auto:1GB" This indicates an automatic threshold where
                            if the file is less than 1GB, parallelization is done by chromosome, but if it's larger,
                            parallelization is done by region of non-overlapping reads. Other modes: bychrom, byregion,
                            auto:xGB - for setting the auto threshold, it must be in units of GB.''')
    parser.add_argument('--predict_cds', default=False, action='store_true',
                        help='specify if you want to predict the CDS of the final isoforms. '
                             'Will be output in the final bed file but not the gtf file. '
                             'Productivity annotation is also added in the name field, '
                             'which is detailed further in the predictProductivity documentation')
    parser.add_argument('--keep_intermediate', default=False, action='store_true',
                        help='''specify if intermediate and temporary files are to be kept for debugging.
                        Intermediate files include: promoter-supported reads file,
                        read assignments to firstpass isoforms''')
    parser.add_argument('--keep_sup', default=False, action='store_true',
                        help='''specify if you want to keep supplementary alignments to define isoforms''')
    parser.add_argument('--end_norm_dist', type=int,
                        help='specify the number of basepairs to extend transcript ends if you want to '
                             'normalize them across transcripts in a gene and extend them')
    parser.add_argument('--output_endpos', default=False, action='store_true',
                        help='specify if you want to output a separate file with corrected read end positions. '
                             'For development purposes')
    parser.add_argument('--output_bam', default=False, action='store_true',
                        help='output intermediate bams aligned to the transcriptome. '
                             'Only works with --keep_intermediate, for debugging')
    parser.add_argument('--fusion_breakpoints',
                        help='''for fusion detection only - bed file containing locations of fusion breakpoints on the synthetic genome''')
    parser.add_argument('--allow_paralogs', default=False, action='store_true',
                        help='specify if want to allow reads to be assigned to multiple paralogs with equivalent alignment')
    parser.add_argument('--generate_map', default=False, action='store_true',
                        help='''specify this argument to generate a txt file of read-isoform assignments''')
    args = parser.parse_args()
    args.parallel_mode = parallel_mode_parse(parser, args.parallel_mode)
    args.trust_ends = False
    args.remove_internal_priming = False

    if not os.path.exists(args.genome_aligned_bam):
        parser.error(f'Aligned reads file path does not exist: {args.genome_aligned_bam}')
    if not os.path.exists(args.genome):
        parser.error(f'Genome file path does not exist: {args.genome}')
    return args


####
# basic types
####
# FIXME: junction_corrector isn't using this, should these be more global?
class Junc(PosRange):
    """Stores start, end, just adds a type name to SeqRange for clearer code and error messages"""
    pass

class Exon(namedtuple("Exon",("start", "end", "name"))):
    def __new__(cls, start, end, name=None):
        assert start <= end
        if name == None:
            name = ''
        return super(Exon, cls).__new__(cls, start, end, name)
    
    def __len__(self):
        return self.end - self.start

def exons_to_juncs(exons):
    """Convert exon ranges to junctions"""
    return [Junc(exons[i].end, exons[i + 1].start)
            for i in range(len(exons) - 1)]


ISO_SRC_ANNOT = 'annot'
ISO_SRC_NOVEL = 'novel'

# tolerance for terminal exon boundary comparisons
TERMINAL_EXON_BOUNDARY_TOLERANCE = 20

# margin for single-exon isoform overlap comparisons
SINGLE_EXON_OVERLAP_MARGIN = 10

# expression ratio threshold for filtering overlapping single-exon isoforms
SINGLE_EXON_EXPRESSION_RATIO = 1.2

# overlap fraction thresholds for gene assignment
MIN_ISOFORM_OVERLAP_FRAC = 0.5
MIN_ANNOT_OVERLAP_FRAC = 0.8

# search window for binary search of single-exon annotations
ANNOT_SE_SEARCH_WINDOW = 2


# class IsoIdSrc(namedtuple("IsoIdSrc",
#                           ("id", "src"))):
#     """isoform identifier along with the source of the isoform"""
#     # FIXME: it is unclear if this is the best way to store the information,
#     # this was create as a transition from iso (id) or (iso_id) (marker, id)
#     pass

####
# misc
###
def make_temp_dir(out_prefix):
    # FIXME: use TMPDIR unless directory explicitly specified
    temp_dir = out_prefix + ".intermediate"
    try:
        os.makedirs(temp_dir, exist_ok=True)
    except OSError as exc:
        raise OSError(f"Creation of the directory `{temp_dir}' failed") from exc
    return temp_dir + '/'


def binary_search(query, data):
    """ Query is a coordinate interval. Binary search for the query in sorted data,
        which is a list of coordinates. Finishes when an overlapping value of query and
        data exists and returns the index in data. """
    # FIXME: uses python bisect module
    i = int(round(len(data) / 2))  # binary search prep
    lower, upper = 0, len(data)
    while True:
        if upper - lower < 2:  # stop condition but not necessarily found
            break
        if data[i][1] < query[0]:
            lower = i
            i = int(round((i + upper) / 2))
        elif data[i][0] > query[1]:
            upper = i
            i = int(round((lower + i) / 2))
        else:  # found
            break
    return i

def bed_to_junctions(bed):
    # FIXME: a junctions object might be good
    return [Junc(bed.blocks[i - 1].end, bed.blocks[i].start)
            for i in range(1, bed.blockCount)]

####
# splice junction correction
####
def read_correct_to_readrec(junction_corrector, read):
    # FIXME: remove unnecessary initial build of ReadRec and make junctions from
    # read, correct, and then make bed
    readrec = ReadRec.from_read(read)
    corrected_bed = junction_corrector.correct_read_bed(readrec.to_bed())
    if corrected_bed is None:
        return None
    readrec.juncs = tuple(Junc(corrected_bed.blocks[i].end, corrected_bed.blocks[i + 1].start)
                          for i in range(len(corrected_bed.blocks) - 1))
    return readrec

def build_intron_support(gtf_file, junction_tab, junction_bed):
    """Build IntronSupport from all available sources."""
    is_db = IntronSupport()
    if junction_bed is not None:
        is_db.load_introns_bed(junction_bed)
    if junction_tab is not None:
        is_db.load_star(junction_tab)
    if gtf_file is not None:
        is_db.load_gtf(gtf_data_parser(gtf_file, attrs=GtfAttrsSet.FLAIR))
    return is_db

def setup_junction_corrector(intron_support, ss_window, junction_support):
    """Create a JunctionCorrector from a pre-built IntronSupport."""
    return JunctionCorrector(intron_support, ss_window, junction_support)

def get_rgb(strand, junclen):
    if junclen == 0:
        return "99,99,99"
    elif strand == '+':
        return "27,158,119"
    else:
        return "217,95,2"


def get_bed_exons_from_juncs(juncs, start, end):
    if len(juncs) == 0:
        exon_starts = [0]
        exon_sizes = [end - start]
    else:
        exon_starts = [0] + [j.end - start for j in juncs]
        exon_sizes = ([juncs[0].start - start] + [juncs[i + 1].start - juncs[i].end
                                                  for i in range(len(juncs) - 1)] +
                      [end - juncs[-1].end])
    return exon_starts, exon_sizes

def get_bed_exons_from_exons(exons, start):
    exon_starts = [e.end - start for e in exons]
    exon_sizes = [e.end - e.start for e in exons]
    return exon_starts, exon_sizes

def get_sequence_for_exons(genome, chrom, strand, exons):
    trans_seq = ''.join([genome.fetch(chrom, e.start, e.end)
                         for e in exons])
    # FIXME: this upper cases only if reverse strand
    if strand == '-':
        trans_seq = get_reverse_complement(trans_seq)
    return trans_seq


class ReadRec:
    """Read alignment with location, junction, and metadata fields.

    Stores chrom, start, end, name, score, strand, and juncs directly.
    Exons are computed on the fly from start, end, and juncs.
    A Bed record can be produced on demand via to_bed().
    """
    def __init__(self, chrom, strand, juncs, start=None, end=None, name=None):
        self.chrom = chrom
        self.strand = strand
        self.juncs = juncs
        self.start = start
        self.end = end
        self.name = name
    
    @property
    def exons(self):
        """Return exons as list of Exon objects, computed from start, end, and juncs."""
        if self.start == None or self.end == None:
            return None
        if not self.juncs:
            return [Exon(self.start, self.end)]
        exons = [Exon(self.start, self.juncs[0].start)]
        for i in range(len(self.juncs) - 1):
            exons.append(Exon(self.juncs[i].end, self.juncs[i + 1].start))
        exons.append(Exon(self.juncs[-1].end, self.end))
        return exons

    def generate_name(self):
        if self.name == None:
            self.name = str(hash(tuple(self.exons)))

    def to_bed(self):
        """Create and return a Bed object."""
        if self.start == None or self.end == None or self.name == None:
            return None
        exon_starts, exon_sizes = get_bed_exons_from_juncs(self.juncs, self.start, self.end)
        bed = Bed(self.chrom, self.start, self.end, self.name,
                  score=self.score, strand=self.strand,
                  thickStart=self.start, thickEnd=self.end,
                  itemRgb=get_rgb(self.strand, len(self.juncs)))
        for i in range(len(exon_starts)):
            blk_start = self.start + exon_starts[i]
            bed.addBlock(blk_start, blk_start + exon_sizes[i])
        return bed

    def get_bed_line(self):
        """Return BED format row."""
        return self.to_bed().toRow()

    def get_sequence(self, genome):
        if self.start == None or self.end == None:
            return None
        return get_sequence_for_exons(genome, self.chrom, self.strand, self.exons)

    def update_from_juncs(self, new_juncs):
        """Update juncs, keeping chrom, start, end, name, score, strand."""
        self.juncs = tuple(new_juncs)
    
    @property
    def genomic_length(self):
        if self.start == None or self.end == None:
            return None
        return self.end - self.start
    
    def reset_from_exons(self, exons):
        """Update ReadRec from a list of Exon objects."""
        self.start = exons[0].start
        self.end = exons[-1].end
        self.juncs = tuple(exons_to_juncs(sorted(exons)))
    
    def score(self):
        return 0
    
    @classmethod
    def from_read(cls, read, junc_direction=None):
        """Create a ReadRec from a pysam aligned read."""
        # FIXME switch to pycbio.hgdata.cigar
        align_start = read.reference_start
        ref_pos = align_start
        intron_blocks = []
        has_match = False
        for block in read.cigartuples:
            if block[0] == 3:  # intron
                if has_match:
                    intron_blocks.append([ref_pos, ref_pos + block[1]])
                # this fixes weird bug if there's an intron, then an insertion, then another intron???
                elif len(intron_blocks) > 0:
                    intron_blocks[-1][1] += block[1]
                has_match = False
                ref_pos += block[1]
            elif block[0] in {0, 7, 8, 2}:  # consumes reference
                ref_pos += block[1]
                if block[0] in {0, 7, 8}:
                    has_match = True
        if junc_direction not in {'+', '-'}:
            junc_direction = "-" if read.is_reverse else "+"
        juncs = tuple(Junc(blk[0], blk[1]) for blk in intron_blocks)
        return cls(read.reference_name, junc_direction, juncs, align_start, ref_pos, read.query_name)

    @classmethod
    def from_junctions(cls, chrom, start, end, name, score, strand, juncs):
        """Create a ReadRec from junction coordinates."""
        return cls(chrom, strand, juncs if isinstance(juncs, tuple) else tuple(juncs), start, end, name)


class IsoWithReads(ReadRec):
    def __init__(self, chrom, strand, juncs, start=None, end=None, reads=None):
        super().__init__(chrom, strand, juncs, start, end)
        self.reads = reads if reads != None else []
        self.gene_id = None
        self.transcript_id = None
    
    @classmethod
    def from_readrec(cls, readrec):
        return cls(readrec.chrom, readrec.strand, readrec.juncs)
    
    @property
    def starts(self):
        return [x.start for x in reads]
    
    @property
    def ends(self):
        return [x.end for x in reads]

    @property
    def num_reads(self):
        return len(self.reads)

    @property
    def score(self):
        return len(self.reads)
    


class ReadEndInfo:
    def __init__(self, start, end, name, mapq):
        self.start = start
        self.end = end
        self.name = name
        self.mapq = mapq
        self.polyA = None
        self.intprim = None
    
    @classmethod
    def from_readrec(cls, readrec):
        return cls(readrec.start, readrec.end, readrec.name, readrec.score)

    


class AnnotData(object):
    def __init__(self):
        # FIXME: what the keys of these dicts()?
        # FIXME: update names

        # map of (transcript_id, gene_id) -> (start, end)
        self.transcript_to_exons = {}

        # list of (transcript_id, gene_id, strand)
        self.transcripts = []

        # map of ((start0, end0), ...) -> (transcript_id, gene_id)
        self.juncchain_to_transcript = {}

        # map of (start, ent) -> set of (transcript_id, gene_id)
        self.junc_to_gene = {}

        # list of (start, end, strand, gene_id):
        # FIXME: rename once it is figured out how this works in get_single_exon_gene_overlaps
        # FIXME: make set
        self.all_annot_SE = []

        # map of strand to map of gene_id  to set of (start, end)
        # FIXME: why is strand needed here
        self.spliced_exons = {'+': {}, '-': {}}

        # map of gene_id to set of (start, end)
        self.gene_to_annot_juncs = {}

        # map of gene_id to strand
        self.gene_to_strand = {}


def annot_data_from_gtf(gtf_data, region):
    """Build AnnotData for a region from a pre-partitioned GtfData."""
    annots = AnnotData()
    if gtf_data is None:
        return annots
    region_map = {region: annots}
    for trans in gtf_data.transcripts:
        exons = [Exon(exon.start, exon.end) for exon in trans.exons]
        if not exons:
            continue
        sorted_exons = sorted(exons)
        t_start = sorted_exons[0].start
        t_end = sorted_exons[-1].end
        save_transcript_annot_to_region(trans.transcript_id, trans.gene_id, region,
                                        region_map, t_start, t_end, trans.strand, sorted_exons)
    return annots

def save_transcript_annot_to_region(transcript_id, gene_id, region, regions_to_annot_data, t_start, t_end, strand, t_exons):
    # regions is tuple of ('chr20', 0, 64444167)
    # FIXME: t_exons are list of (32186476, 32190360)
    assert isinstance(t_exons[0], Exon)  # FIXME tmp debugging
    sorted_exons = sorted(t_exons)
    annots = regions_to_annot_data[region]
    annots.transcript_to_exons[(transcript_id, gene_id)] = tuple(sorted_exons)
    juncs = exons_to_juncs(sorted_exons)
    annots.transcripts.append((transcript_id, gene_id, strand))
    if gene_id not in annots.gene_to_strand:
        annots.gene_to_strand[gene_id] = strand
    if len(juncs) == 0:
        annots.all_annot_SE.append((t_start, t_end, strand, gene_id))
    else:
        if gene_id not in annots.spliced_exons[strand]:
            annots.spliced_exons[strand][gene_id] = set()
        annots.spliced_exons[strand][gene_id].update(set(sorted_exons))
        annots.juncchain_to_transcript[tuple(juncs)] = (transcript_id, gene_id)
        if gene_id not in annots.gene_to_annot_juncs:
            annots.gene_to_annot_juncs[gene_id] = set()
        for j in juncs:
            if j not in annots.junc_to_gene:
                annots.junc_to_gene[j] = set()
            annots.junc_to_gene[j].add((transcript_id, gene_id))
            annots.gene_to_annot_juncs[gene_id].add(j)
    annots.all_annot_SE = sorted(annots.all_annot_SE)  # FIXME: make set?


def get_filter_tome_align_cmd(args, ref_bed, output_name, map_file, is_annot, clipping_file, unique_bound):
    # FIXME: convert filter_transcriptome_align.py to a library, however
    # minimap output needs to be piped through filter_transcriptome_align
    # without saving the bam file.

    # count sam transcripts ; the dash at the end means STDIN
    # use 1 thread in because this is already multithreaded here
    count_cmd = ['filter_transcriptome_align.py', '--sam', '-',
    # count_cmd = ['count_sam_transcripts.py', '--sam', '-',
                 '-o', output_name, '-t', 1,]
    if clipping_file:
        count_cmd.extend(['--trimmedreads', clipping_file])
    if map_file:
        count_cmd.extend(['--generate_map', map_file])
    if args.output_endpos or is_annot:
        count_cmd.extend(['--output_endpos', output_name.split('.counts.txt')[0] + '.ends.tsv'])
    if args.end_norm_dist:
        count_cmd.extend(['--end_norm_dist', args.end_norm_dist])
    if not is_annot and not args.no_stringent:
        count_cmd.extend(['--stringent'])
    if is_annot:
        count_cmd.extend(['--allow_UTR_indels'])
    if args.output_bam:
        count_cmd.extend(['--output_bam', output_name.split('.counts.txt')[0] + '.bam'])
    if not args.no_check_splice:
        count_cmd.append('--check_splice')
    if not args.no_check_splice or not args.no_stringent or is_annot or args.fusion_breakpoints:
        count_cmd.extend(['-i', ref_bed])  # annotated isoform bed file
    if args.trust_ends:
        count_cmd.append('--trust_ends')
    if unique_bound and (not args.no_stringent or is_annot):
        count_cmd.extend(['--unique_bound', unique_bound])
    if args.remove_internal_priming:
        count_cmd.extend(['--remove_internal_priming',
                          '--intprimingthreshold', str(args.intprimingthreshold),
                          '--intprimingfracAs', str(args.intprimingfracAs),
                          '--transcriptomefasta', args.transcriptfasta])
    if args.remove_internal_priming and is_annot:
        count_cmd.append('--permissive_last_exons')
    if args.fusion_breakpoints:
        count_cmd += ['--fusion_breakpoints', args.fusion_breakpoints]
    if args.allow_paralogs:
        count_cmd += ['--allow_paralogs']
    # print(' '.join([str(x) for x in count_cmd]))
    return count_cmd


def transcriptome_align_and_count(args, input_reads, align_ref_fasta, ref_bed, output_name, map_file, is_annot, clipping_file, unique_bound):
    # minimap (results are piped into count_sam_transcripts.py)
    # '--split-prefix', 'minimap2transcriptomeindex', doesn't work with MD tag
    if isinstance(input_reads, str):
        input_reads = [input_reads]
    mm2_cmd = ['minimap2', '-a', '-N', '4', '--MD'] + [align_ref_fasta] + input_reads

    # FIXME add in step to filter out chimeric reads here
    # FIXME really need to go in and check on how count_sam_transcripts is working
    count_cmd = get_filter_tome_align_cmd(args, ref_bed, output_name, map_file, is_annot, clipping_file, unique_bound)
    pipettor.run([mm2_cmd, count_cmd])

##
# Transcript end assignment
##

# class IsoformEndInfo:
#     """Represents a group of reads sharing similar transcript ends"""
#     def __init__(self, start, end, strand, read_id, weighted_score=0.0, supporting_reads=None):
#         self.start = start
#         self.end = end
#         self.strand = strand
#         self.read_id = read_id  # Representative read for this group
#         self.weighted_score = weighted_score
#         self.supporting_reads = supporting_reads or []

#     @property
#     def num_reads(self):
#         return len(self.supporting_reads)

#     @property
#     def length(self):
#         return self.end - self.start

#     @property
#     def score(self):
#         """Alias for num_reads for compatibility"""
#         return self.num_reads

# class IsoformInfo:
#     """Represents metadata about a detected isoform"""
#     def __init__(self, transcript_id, strand, exons):
#         self.transcript_id = transcript_id
#         self.strand = strand
#         self.exons = exons
#         self.gene_id = None  # Assigned later

#     @property
#     def start(self):
#         return self.exons[0][0] if self.exons else None

#     @property
#     def end(self):
#         return self.exons[-1][1] if self.exons else None

#     def set_gene_id(self, gene_id):
#         self.gene_id = gene_id

# class EndInfo:
#     """Represents transcript end information for final isoforms"""
#     def __init__(self, start, end, iso_id_src, read_names):
#         self.start = start
#         self.end = end
#         self.iso_id_src = iso_id_src
#         self.read_names = read_names

#     @property
#     def score(self):
#         return min(len(self.read_names), 1000)

#     @property
#     def length(self):
#         return self.end - self.start

JunctionChain = namedtuple('JunctionChain', ('chrom', 'strand', 'juncs'))

# class GeneIsoformData:
#     """Organizes all isoforms for a single gene"""
#     def __init__(self, gene_id):
#         self.gene_id = gene_id
#         self._by_junction = {}

#     def add_isoform(self, chrom, strand, juncs, end_info):
#         """Add an isoform to the gene's isoform data"""
#         junc_chain = JunctionChain(chrom, strand, juncs)
#         if junc_chain not in self._by_junction:
#             self._by_junction[junc_chain] = []
#         self._by_junction[junc_chain].append(end_info)

#     def get_isoforms(self, chrom, strand, juncs):
#         """Get all isoforms for a specific junction chain"""
#         junc_chain = JunctionChain(chrom, strand, juncs)
#         return self._by_junction.get(junc_chain, [])

#     def set_isoforms(self, chrom, strand, juncs, isoforms):
#         """Set the isoforms list for a specific junction chain"""
#         junc_chain = JunctionChain(chrom, strand, juncs)
#         self._by_junction[junc_chain] = isoforms

#     def junction_chains(self):
#         """Iterate over all junction chains in this gene"""
#         return iter(self._by_junction.keys())

#     @property
#     def total_read_support(self):
#         """Sum of all reads supporting all isoforms"""
#         return sum(len(iso.read_names)
#                    for isos in self._by_junction.values()
#                    for iso in isos)

def get_best_ends(curr_group, end_window):
    best_ends = []
    if len(curr_group) > int(end_window):
        all_starts = Counter([x.start for x in curr_group])
        all_ends = Counter([x.end for x in curr_group])
        for read_info in curr_group:
            weighted_score = all_starts[read_info.start] + all_ends[read_info.end]
            best_ends.append((weighted_score, read_info.start, read_info.end))
    else:
        for read_info1 in curr_group:
            score, weighted_score = 0, 0
            for read_info2 in curr_group:
                if abs(read_info1.start - read_info2.start) <= end_window and abs(read_info1.end - read_info2.end) <= end_window:
                    score += 2
                    weighted_score += (((end_window - abs(read_info1.start - read_info2.start)) / end_window) +
                                       ((end_window - abs(read_info1.end - read_info2.end)) / end_window))
            best_ends.append((weighted_score, read_info1.start, read_info1.end))
    best_ends.sort(reverse=True)
    # FIXME: DO I WANT TO ADD CORRECTION TO NEARBY ANNOTATED TSS/TTS????
    return best_ends[0]

def group_reads_by_ends(read_info_list, sort_index, end_window):
    sorted_ends = sorted(read_info_list, key=lambda x: x.start if sort_index == 0 else x.end)
    new_groups, group = [], []
    last_edge = 0
    for iso_info in sorted_ends:
        edge = iso_info.start if sort_index == 0 else iso_info.end
        if edge - last_edge <= end_window:
            group.append(iso_info)
        else:
            if len(group) > 0:
                new_groups.append(group)
            group = [iso_info]
        last_edge = edge
    if len(group) > 0:
        new_groups.append(group)
    return new_groups

# MAIN METHOD
# read_ends is a list containing elements with: (read.start, read.end, read.strand, read.name)
# If the reads are spliced, the group will contain only the info for reads with a shared splice junction
# if the reads are unspliced, the group will contain info for all unspliced reads in a given chromosome/region,
# The output is a list of IsoWithReads objects containing:
#    - weighted_score (represents how many reads have ends similar to this exact position)
#    - start, end, strand, read_id (representative read id)
#    - supporting_reads (list of all read names in group)

def collapse_end_groups(end_window, isoform):
    start_groups = group_reads_by_ends(isoform.reads, 0, end_window)
    all_end_groups, iso_end_groups = [], []
    for start_group in start_groups:
        all_end_groups.extend(group_reads_by_ends(start_group, 1, end_window))
    for end_group in all_end_groups:
        weighted_score, start, end = get_best_ends(end_group, end_window)
        read_end_info = IsoWithReads(isoform.chrom, isoform.strand, isoform.juncs, start, end, end_group)
        read_end_info.generate_name()
        iso_end_groups.append(read_end_info)
    return iso_end_groups


def get_isos_with_similar_juncs(juncs, firstpass_junc_to_name, junc_to_gene):
    """Find isoforms sharing junctions with the given junction set.
    Returns separate sets for novel (string UUIDs) and annotated ((transcript_id, gene_id) tuples)."""
    novel_isos = set()
    annot_isos = set()
    for j in juncs:
        if firstpass_junc_to_name and j in firstpass_junc_to_name:
            novel_isos.update(firstpass_junc_to_name[j])
        if j in junc_to_gene:
            annot_isos.update(junc_to_gene[j])
    return novel_isos, annot_isos

def _is_junction_subset(juncs, otheriso_juncs):
    """Check if juncs is a proper subset of otheriso_juncs using string matching."""
    if len(juncs) >= len(otheriso_juncs):
        return False
    iso_juncs_str = str(juncs)[1:-1].rstrip(',')
    otheriso_juncs_str = str(otheriso_juncs)[1:-1]
    return iso_juncs_str in otheriso_juncs_str


def _check_terminal_exon_overlap(first_exon, last_exon, other_exon, otheriso_score,
                                 terminal_exon_is_subset, superset_support):
    """Check overlap with terminal exon of other transcript (first or last).
    Only requires sharing the same terminal splice site."""
    if first_exon.end == other_exon.end:
        terminal_exon_is_subset[0] = 1
        superset_support.append(otheriso_score)
    elif last_exon.start == other_exon.start:
        terminal_exon_is_subset[1] = 1
        superset_support.append(otheriso_score)


def _check_internal_exon_overlap(first_exon, last_exon, other_exon, otheriso_score,
                                 terminal_exon_is_subset, superset_support, unique_seq_bound):
    """Check overlap with internal exon of other transcript.
    Records unique sequence boundaries and checks containment within tolerance."""
    if first_exon.end == other_exon.end:
        unique_seq_bound.append((0, first_exon.end - other_exon.start))
        if first_exon.start >= (other_exon.start - TERMINAL_EXON_BOUNDARY_TOLERANCE):
            terminal_exon_is_subset[0] = 1
            superset_support.append(otheriso_score)
    if last_exon.start == other_exon.start:
        unique_seq_bound.append((1, other_exon.end - last_exon.start))
        if last_exon.end <= (other_exon.end + TERMINAL_EXON_BOUNDARY_TOLERANCE):
            terminal_exon_is_subset[1] = 1
            superset_support.append(otheriso_score)


def _check_junction_subset(juncs, first_exon, last_exon, otheriso_score, otheriso_juncs, otheriso_exons,
                           terminal_exon_is_subset, superset_support, unique_seq_bound):
    """Check if juncs is a subset of otheriso_juncs and update tracking lists."""
    if not _is_junction_subset(juncs, otheriso_juncs):
        return
    for i, other_exon in enumerate(otheriso_exons):
        is_terminal = (i == 0 or i == len(otheriso_exons) - 1)
        if is_terminal:
            _check_terminal_exon_overlap(first_exon, last_exon, other_exon, otheriso_score,
                                         terminal_exon_is_subset, superset_support)
        else:
            _check_internal_exon_overlap(first_exon, last_exon, other_exon, otheriso_score,
                                         terminal_exon_is_subset, superset_support, unique_seq_bound)


def identify_spliced_iso_subset_annot(annot_iso_id, sup_annot_transcript_to_juncs, annots,
                                      juncs, first_exon, last_exon, terminal_exon_is_subset,
                                      superset_support, unique_seq_bound):
    """Check if query isoform is subset of an annotated isoform.
    annot_iso_id is (transcript_id, gene_id) tuple."""
    if sup_annot_transcript_to_juncs is not None:
        # using only supported annotated
        if annot_iso_id not in sup_annot_transcript_to_juncs:
            return
        otheriso_score, otheriso_juncs = sup_annot_transcript_to_juncs[annot_iso_id]
        otheriso_exons = annots.transcript_to_exons[annot_iso_id]
    elif annot_iso_id in annots.transcript_to_exons:
        # using all annotated
        otheriso_exons = annots.transcript_to_exons[annot_iso_id]
        otheriso_juncs = exons_to_juncs(otheriso_exons)
        otheriso_score = 0
    else:
        return
    _check_junction_subset(juncs, first_exon, last_exon, otheriso_score, otheriso_juncs, otheriso_exons,
                           terminal_exon_is_subset, superset_support, unique_seq_bound)


def identify_spliced_iso_subset_novel(novel_iso_id, firstpass_unfiltered,
                                      juncs, first_exon, last_exon, terminal_exon_is_subset,
                                      superset_support, unique_seq_bound):
    """Check if query isoform is subset of a novel (firstpass) isoform.
    novel_iso_id is a string UUID."""
    otheriso = firstpass_unfiltered[novel_iso_id]
    _check_junction_subset(juncs, first_exon, last_exon, otheriso.score, otheriso.juncs, otheriso.exons,
                           terminal_exon_is_subset, superset_support, unique_seq_bound)

def filter_spliced_iso(filter_type, support, juncs, exons, name, score, annots,
                       firstpass_junc_to_name, firstpass_unfiltered,
                       sup_annot_transcript_to_juncs, strand):
    assert isinstance(exons[0], Exon)  # FIXME: debugging
    novel_isos, annot_isos = get_isos_with_similar_juncs(juncs, firstpass_junc_to_name, annots.junc_to_gene)
    terminal_exon_is_subset = [0, 0]  # first exon is a subset, last exon is a subset
    first_exon, last_exon = exons[0], exons[-1]
    superset_support = []
    unique_seq_bound = []
    for novel_iso_id in novel_isos:
        if novel_iso_id != name:
            identify_spliced_iso_subset_novel(novel_iso_id, firstpass_unfiltered,
                                              juncs, first_exon, last_exon, terminal_exon_is_subset,
                                              superset_support, unique_seq_bound)
    unique_seq_bound = list(set(unique_seq_bound))
    if strand == '-':
        # just invert the indexes
        for i in range(len(unique_seq_bound)):
            unique_seq_bound[i] = f'{abs(unique_seq_bound[i][0] - 1)}_{unique_seq_bound[i][1]}'
    else:
        for i in range(len(unique_seq_bound)):
            unique_seq_bound[i] = f'{unique_seq_bound[i][0]}_{unique_seq_bound[i][1]}'

    if sum(terminal_exon_is_subset) < 2:  # both first and last exon have to overlap
        return True, unique_seq_bound
    elif filter_type != 'nosubset':
        if score >= support and score > max(superset_support) * 1.2:
            return True, unique_seq_bound
    return False, None

####
# terminal exon normalization
####
class GeneMaxTerminalExonsEnds:
    """Class to collect the maximal terminal exons ends for a gene.
    Exons are groups based on the location of the internal splice junction"""
    def __init__(self, gene_id):
        self.gene_id = gene_id
        self.left_ends = {}
        self.right_ends = {}

    def add_left_end(self, exon):
        if exon.end not in self.left_ends:
            self.left_ends[exon.end] = exon.start
        else:
            self.left_ends[exon.end] = min(self.left_ends[exon.end], exon.start)

    def add_right_end(self, exon):
        if exon.start not in self.right_ends:
            self.right_ends[exon.start] = exon.end
        else:
            self.right_ends[exon.start] = max(self.right_ends[exon.start], exon.end)

    def get_left_end(self, exon):
        return self.left_ends[exon.end]

    def get_right_end(self, exon):
        return self.right_ends[exon.start]

class MaxTerminalExonsEnds:
    """Collection of maximal terminal exons ends by gene.
    A genes terminal exons are grouped by the interior exon splice junction
    location.
    """
    def __init__(self):
        # FIXME: this is temporary.  The code groups by (gene_id, strand)
        # for reasons that are suspected to be bugs in stranding.  We keep
        # this but generate a warning until we are sure it is fixed
        self._by_gene_id = {}  # (gene_id, strand) -> GeneMaxTerminalExonsEnds
        self._gene_id_to_strand = {}
        self._genes_warned = set()

    def _obtain(self, gene_id, strand):
        "get current entry or create a new one"
        gene_key = (gene_id, strand)
        gene_entry = self._by_gene_id.get(gene_key)
        if gene_entry is None:
            gene_entry = GeneMaxTerminalExonsEnds(gene_id)
            self._by_gene_id[gene_key] = gene_entry

        # FIXME: tmp generate strand warning
        existing_strand = self._gene_id_to_strand.get(gene_id)
        if existing_strand is None:
            self._gene_id_to_strand[gene_id] = strand
        elif (strand != existing_strand) and (gene_id not in self._genes_warned):
            self._genes_warned.add(gene_id)
            logging.warning("BUG: gene id '%s' has transcripts on both strands", gene_id)

        return gene_entry

    def add_transcript(self, gene_id, strand, transcript_id, exons):
        # don't normalize ends for single exon transcripts, but still record gene
        # FIXME: do we actually want to add single-exon genes?
        gene_entry = self._obtain(gene_id, strand)
        if len(exons) > 1:
            gene_entry.add_left_end(exons[0])
            gene_entry.add_right_end(exons[-1])

    def fetch(self, gene_id, strand) -> GeneMaxTerminalExonsEnds:
        """return entry or error"""
        return self._by_gene_id[(gene_id, strand)]

def max_terminal_exons_ends_from_annots(annots):
    max_terminal_exons_ends = MaxTerminalExonsEnds()
    for transcript_id, gene_id, strand in annots.transcripts:
        exons = annots.transcript_to_exons[(transcript_id, gene_id)]
        max_terminal_exons_ends.add_transcript(gene_id, strand, transcript_id, exons)
    return max_terminal_exons_ends

def max_terminal_exons_ends_from_iso_infos(iso_to_info):
    max_terminal_exons_ends = MaxTerminalExonsEnds()
    for iso_name in iso_to_info:
        isoform = iso_to_info[iso_name]
        max_terminal_exons_ends.add_transcript(isoform.gene_id, isoform.strand, isoform.transcript_id, isoform.exons)
    return max_terminal_exons_ends

####
# transcriptome reference
####
def normalize_gene_terminal_exons(max_terminal_exons_ends, gene_id, strand, exons,
                                  *, add_length_at_ends=0):
    "updates terminal exons ends"
    gene_terminal_exons = max_terminal_exons_ends.fetch(gene_id, strand)
    exons[0] = Exon(gene_terminal_exons.get_left_end(exons[0]) - add_length_at_ends,
                    exons[0].end)
    exons[-1] = Exon(exons[-1].start,
                     gene_terminal_exons.get_right_end(exons[-1]) + add_length_at_ends)

def generate_transcriptome_reference_transcript(strand, transcript_to_strand, transcript_id, gene_id, annots, normalize_ends, max_terminal_exons_ends,
                                                add_length_at_ends, transcript_to_new_exons, chrom, genome, annot_bed_fh, annot_fa_fh, annot_uniqueseq_fh):
    transcript_to_strand[(transcript_id, gene_id)] = strand
    exons = list(annots.transcript_to_exons[(transcript_id, gene_id)])
    assert isinstance(exons[0], Exon)  # FIXME tmp debugging
    juncs = exons_to_juncs(exons)
    is_not_subset, unique_seq = filter_spliced_iso('nosubset', 0, juncs, exons, (transcript_id, gene_id),
                                                   0, annots, None, None, None, strand)
    if is_not_subset:
        if normalize_ends:
            normalize_gene_terminal_exons(max_terminal_exons_ends, gene_id, strand, exons,
                                          add_length_at_ends=add_length_at_ends)
            transcript_to_new_exons[(transcript_id, gene_id)] = tuple(exons)
        exons = tuple(exons)
        start, end = exons[0].start, exons[-1].end

        # FIXME: duplicated code
        exon_starts, exon_sizes = get_bed_exons_from_exons(exons, start)
        # FIXME: duplicated use BED class,
        bed_line = [chrom, start, end, transcript_id + '_' + gene_id, '.', strand, start, end, '0', len(exons),
                    ','.join([str(x) for x in exon_sizes]), ','.join([str(x) for x in exon_starts])]
        trans_seq = get_sequence_for_exons(genome, chrom, strand, exons)
        annot_bed_fh.write('\t'.join([str(x) for x in bed_line]) + '\n')
        annot_fa_fh.write('>' + transcript_id + '_' + gene_id + '\n')
        annot_fa_fh.write(''.join(trans_seq) + '\n')
        if len(unique_seq) > 0:
            annot_uniqueseq_fh.write(transcript_id + '_' + gene_id + '\t' + ','.join(unique_seq) + '\n')

def generate_transcriptome_reference_guts(normalize_ends, annots, add_length_at_ends, chrom, genome, annot_bed_fh, annot_fa_fh, annot_uniqueseq_fh):
    transcript_to_strand = {}
    transcript_to_new_exons = {}
    max_terminal_exons_ends = None
    if normalize_ends:
        max_terminal_exons_ends = max_terminal_exons_ends_from_annots(annots)

    for transcript_id, gene_id, strand in annots.transcripts:
        generate_transcriptome_reference_transcript(strand, transcript_to_strand, transcript_id, gene_id, annots, normalize_ends, max_terminal_exons_ends, add_length_at_ends,
                                                    transcript_to_new_exons, chrom, genome, annot_bed_fh, annot_fa_fh, annot_uniqueseq_fh)
    return transcript_to_strand, transcript_to_new_exons

def generate_transcriptome_reference(temp_prefix, annots, chrom, genome,
                                     normalize_ends=False, add_length_at_ends=0):
    with (open(temp_prefix + '.annotated_transcripts.bed', 'w') as annot_bed_fh,
          open(temp_prefix + '.annotated_transcripts.fa', 'w') as annot_fa_fh,
          open(temp_prefix + '.annotated_transcripts_uniquebound.txt', 'w') as annot_uniqueseq_fh):
        return generate_transcriptome_reference_guts(normalize_ends, annots, add_length_at_ends, chrom, genome, annot_bed_fh, annot_fa_fh, annot_uniqueseq_fh)


def identify_good_match_to_annot(args, temp_prefix, chrom, annots, genome):
    # FIXME: refactor
    # good_align_to_annot, firstpass_SE, sup_annot_transcript_to_juncs = [], set(), {}
    read_to_transcript = {}
    if not args.no_align_to_annot and len(annots.transcripts) > 0:
        # logging.info('generating transcriptome reference')
        # this part generates the fasta file for the annotation
        if args.end_norm_dist is not None:
            transcript_to_strand, transcript_to_new_exons = \
                generate_transcriptome_reference(temp_prefix, annots, chrom, genome,
                                                 normalize_ends=True,
                                                 add_length_at_ends=args.end_norm_dist)
        else:
            transcript_to_strand, transcript_to_new_exons = \
                generate_transcriptome_reference(temp_prefix, annots, chrom, genome)
        # FIXME: make a TSV
        clipping_file = temp_prefix + '.reads.genomicclipping.txt'
        logging.info('aligning to transcriptome reference')
        transcriptome_align_and_count(args, temp_prefix + '.reads.fasta',
                                      temp_prefix + '.annotated_transcripts.fa',
                                      temp_prefix + '.annotated_transcripts.bed',
                                      temp_prefix + '.matchannot.counts.txt',
                                      temp_prefix + '.matchannot.read.map.txt', True,
                                      clipping_file,
                                      temp_prefix + '.annotated_transcripts_uniquebound.txt')
        logging.info('processing good matches')
        for line in open(temp_prefix + '.matchannot.ends.tsv'):
            line = line.rstrip().split('\t')
            read, transcript = line[:2]
            startindex, startdist, endindex, enddist = [int(x) for x in line[2:]]
            read_to_transcript[read] = (transcript, startindex, startdist, endindex, enddist)
    else:
        # create empty output files
        # FIXME: why doesn't this create all of them?
        # FIXME: change above logic so file open all happens in one place
        with (open(temp_prefix + '.matchannot.counts.txt', 'w'),
              open(temp_prefix + '.matchannot.read.map.txt', 'w')):
            pass
        if args.output_endpos:
            with open(temp_prefix + '.ends.tsv', 'w') as _:
                pass
    # good_align_to_annot = set(good_align_to_annot)
    # return good_align_to_annot, firstpass_SE, sup_annot_transcript_to_juncs
    return read_to_transcript


def _should_process_read(read, region, min_quality, keep_sup, allow_secondary, allow_outside_range=False):
    """Check if read passes filtering criteria for processing"""
    if read.mapping_quality < min_quality:
        return False
    if read.is_secondary and not allow_secondary:
        return False
    if read.is_supplementary and not keep_sup:
        return False
    if read.reference_name != region.name:
        return False
    if not (region.start <= read.reference_start and read.reference_end <= region.end) and not allow_outside_range:
        return False
    return True


def _add_corrected_read_to_groups(corrected_read, sj_to_ends):
    """Add a corrected read to the junction-to-ends mapping"""
    junc_key = tuple(sorted(corrected_read.juncs))
    if junc_key not in sj_to_ends:
        sj_to_ends[junc_key] = IsoWithReads.from_readrec(corrected_read)
    sj_to_ends[junc_key].reads.append(ReadEndInfo.from_readrec(corrected_read))


def _correct_and_group_read(read, read_to_annot_transcript, annots, junction_corrector, sj_to_ends):
    """Correct a single read's splice junctions and add it to sj_to_ends groups."""
    if read.query_name in read_to_annot_transcript:
        transcript, startindex, startdist, endindex, enddist = read_to_annot_transcript[read.query_name]
        transcript, gene = transcript.split('_')
        exons = annots.transcript_to_exons[(transcript, gene)]
        juncs = [(exons[x].end, exons[x + 1].start) for x in range(len(exons) - 1)]
        # ignore unspliced transcripts
        if len(juncs) > 0:
            newstart = juncs[startindex][0] - startdist
            newend = juncs[endindex][1] + enddist
            juncs = tuple([Junc(x[0], x[1]) for x in juncs[startindex:endindex + 1]])
            strand = '-' if read.is_reverse else '+'
            corrected_read = ReadRec.from_junctions(read.reference_name, newstart, newend,
                                                   read.query_name, read.mapping_quality,
                                                   strand, juncs)
        else:
            corrected_read = read_correct_to_readrec(junction_corrector, read)
    else:
        corrected_read = read_correct_to_readrec(junction_corrector, read)
    if corrected_read:
        _add_corrected_read_to_groups(corrected_read, sj_to_ends)

def filter_correct_group_reads(args, temp_prefix, region, bam_file, read_to_annot_transcript, annots,
                               junction_corrector, generate_fasta=True, sj_to_ends=None,
                               allow_secondary=False):
    """Filter reads, correct splice junctions, and group by junction chain."""
    if sj_to_ends is None:
        sj_to_ends = {}
    for read in bam_file.fetch(region.name, region.start, region.end):
        if _should_process_read(read, region, args.quality, args.keep_sup, allow_secondary):
            _correct_and_group_read(read, read_to_annot_transcript, annots, junction_corrector, sj_to_ends)
    return sj_to_ends


def filter_ends_allow_multiple(isoforms, sjc_support, max_ends):
    """Allow multiple ends per junction chain.
    Returns list of IsoWithReads objects that meet support threshold."""
    if isoforms[0].num_reads < sjc_support:
        # If top candidate doesn't meet threshold, merge all reads into it
        best = isoforms[0]
        for iso in isoforms[1:]:
            best.reads.extend(iso.reads)
        return [best]
    else:
        # Filter to those meeting support threshold and limit to max_ends
        filtered = [x for x in isoforms if x.num_reads >= sjc_support]
        filtered = filtered[:max_ends]  # select only top most supported ends
        return filtered

def filter_ends_single_best(isoforms, no_redundant_mode):
    """Pick single best end from junction chain.
    Returns list with single IsoWithReads object."""
    # best_only uses the default sorting, doesn't require additional action
    if no_redundant_mode == 'longest':
        isoforms.sort(reverse=True, key=lambda x: x.length)

    # Pick single best end and merge all reads into it
    best = isoforms[0]
    for iso in isoforms[1:]:
        best.reads.extend(iso.reads)
    return [best]

def filter_ends_by_redundant_and_support(args, isoforms):
    """Sort ends, then select best ones based on support and value of args.no_redundant."""
    # First by weighted score, then by length
    isoforms.sort(key=lambda x: [x.num_reads, x.genomic_length],
                                  reverse=True)

    junc_support = sum([x.num_reads for x in isoforms])
    if junc_support < args.sjc_support:
        return []

    if args.no_redundant == 'none':
        # Allow multiple ends per junction chain
        return filter_ends_allow_multiple(isoforms, args.sjc_support, args.max_ends)
    else:
        # Pick single best end
        return filter_ends_single_best(isoforms, args.no_redundant)


def _write_unfiltered_ends(isoforms, fh):
    for iso_readrec in isoforms:
        fh.write('\t'.join(iso_readrec.get_bed_line()) + '\n')

def _add_end_to_firstpass(isoform, firstpass_unfiltered, firstpass_junc_to_name, firstpass_exons, iso_fh):
    iso_fh.write('\t'.join(isoform.get_bed_line()) + '\n')
    if isoform.juncs == ():
        firstpass_exons.add(Exon(isoform.start, isoform.end, isoform.name))
    else:
        for j in isoform.juncs:
            if j not in firstpass_junc_to_name:
                firstpass_junc_to_name[j] = set()
            firstpass_junc_to_name[j].add(isoform.name)
        firstpass_exons.update(set(isoform.exons))

def _process_junc_ends(args, isoforms, firstpass_unfiltered, firstpass_junc_to_name, firstpass_exons, iso_fh):
    if isoforms[0].juncs == ():
        filtered_isoforms = [x for x in isoforms if len(x.reads) >= args.sjc_support]
    else:
        filtered_isoforms = filter_ends_by_redundant_and_support(args, isoforms)
    for isoform in filtered_isoforms:
        firstpass_unfiltered[isoform.name] = isoform
        _add_end_to_firstpass(isoform, firstpass_unfiltered, firstpass_junc_to_name, firstpass_exons, iso_fh)

def _process_junc(args, isoform, firstpass_unfiltered, firstpass_junc_to_name, firstpass_exons, iso_fh, iso_unfilt_fh):
    these_firstpass = collapse_end_groups(args.end_window, isoform)
    _write_unfiltered_ends(these_firstpass, iso_unfilt_fh)
    _process_junc_ends(args, these_firstpass, firstpass_unfiltered, firstpass_junc_to_name, firstpass_exons, iso_fh)

def process_juncs_to_firstpass_isos(args, temp_prefix, sj_to_ends):
    firstpass_exons = set()
    firstpass_unfiltered, firstpass_junc_to_name = {}, {}
    with open(temp_prefix + '.firstpass.unfiltered.bed', 'w') as iso_fh, \
            open(temp_prefix + '.firstpass.reallyunfiltered.bed', 'w') as iso_unfilt_fh:
        for juncs, isoform in sj_to_ends.items():
            _process_junc(args, isoform, firstpass_unfiltered, firstpass_junc_to_name, firstpass_exons, iso_fh, iso_unfilt_fh)
    firstpass_exons = sorted(list(firstpass_exons))
    return firstpass_unfiltered, firstpass_junc_to_name, firstpass_exons

####
# single-exon transcript processing
####
def filter_single_exon_iso(args, single_exon, curr_group, firstpass_unfiltered):
    # FIXME: make object: grouped_iso (32186479, 32188247, '99bfe5c4-0f3a-4f4d-b5c9-bac459c45e5c')
    # FIXME: curr_group is a list of these
    iso_readrec = firstpass_unfiltered[single_exon.name]
    # FIXME: what does 'comp_' mean?
    expression_comp_with_superset = []
    is_contained = False
    for exon in curr_group:
        if exon != single_exon:
            if ((exon.start - SINGLE_EXON_OVERLAP_MARGIN) <= single_exon.start and
                    single_exon.end <= (exon.end + SINGLE_EXON_OVERLAP_MARGIN)):
                if exon.name != '' or args.filter == 'nosubset':  # is exon from spliced transcript
                    is_contained = True
                    break  # filter out
                else:  # is other single exon - check relative expression
                    other_score = firstpass_unfiltered[exon.name].score
                    score = iso_readrec.score
                    if score >= args.sjc_support and other_score * SINGLE_EXON_EXPRESSION_RATIO < score:
                        expression_comp_with_superset.append(True)
                    else:
                        expression_comp_with_superset.append(False)
    return not is_contained and all(expression_comp_with_superset)


def filter_single_exon_group(args, curr_group, firstpass_unfiltered, firstpass):
    #grouping single exon isoforms with exons from spliced isoform, using spliced isoform exons to filter single exon isoform
    for exon in curr_group:
        if exon.name != '':  # is single exon with name
            if filter_single_exon_iso(args, exon,  curr_group, firstpass_unfiltered):
                firstpass[exon.name] = firstpass_unfiltered[exon.name]
    return firstpass


def filter_all_single_exon(args, firstpass_exons, firstpass_unfiltered, firstpass):
    # group_start = 0
    last_end = 0
    curr_group = []

    for exon in firstpass_exons:
        if exon.start < last_end:
            curr_group.append(exon)
        else:
            if len(curr_group) > 0:
                firstpass = filter_single_exon_group(args, curr_group, firstpass_unfiltered, firstpass)
            curr_group = [exon]
            # group_start = start
        if exon.end > last_end:
            last_end = exon.end
    if len(curr_group) > 0:
        firstpass = filter_single_exon_group(args, curr_group, firstpass_unfiltered, firstpass)

    return firstpass


def filter_firstpass_isos(args, firstpass_unfiltered, firstpass_junc_to_name, firstpass_exons, annots,
                          sup_annot_transcript_to_juncs):
    # FIXME: firstpass_unfiltered is a dict of uuid to ReadRec
    iso_to_unique_bound = {}
    if args.filter == 'ginormous':
        firstpass = firstpass_unfiltered
    else:
        firstpass = {}
        for iso_name in firstpass_unfiltered:
            isoform = firstpass_unfiltered[iso_name]
            if isoform.juncs != ():
                if args.filter == 'comprehensive':
                    firstpass[iso_name] = isoform
                else:
                    assert isinstance(isoform.exons[0], Exon)  # FIXME tmp debugging
                    is_not_subset, unique_seq = filter_spliced_iso(args.filter, args.sjc_support, isoform.juncs, isoform.exons,
                                                                   iso_name, isoform.num_reads, annots,
                                                                   firstpass_junc_to_name, firstpass_unfiltered,
                                                                   sup_annot_transcript_to_juncs, isoform.strand)
                    if is_not_subset:
                        firstpass[iso_name] = isoform
                        if len(unique_seq) > 0:
                            iso_to_unique_bound[iso_name] = ','.join(unique_seq)
        # HANDLE SINGLE EXONS SEPARATELY - group first - one traversal of list
        firstpass = filter_all_single_exon(args, firstpass_exons, firstpass_unfiltered, firstpass)

    return firstpass, iso_to_unique_bound


def combine_temp_files_by_suffix(output, temp_prefixes, suffixes):
    for filesuffix in suffixes:
        with open(output + filesuffix, 'wb') as combined_fh:
            for temp_prefix in temp_prefixes:
                with open(temp_prefix + filesuffix, 'rb') as in_fh:
                    shutil.copyfileobj(in_fh, combined_fh, 1024 * 1024 * 10)


def get_genes_with_shared_juncs(juncs, annots):
    # FIXME: what does this actually return?
    gene_hits = {}
    if juncs != ():
        for j in juncs:
            if j in annots.junc_to_gene:
                for transcript_id, gene_id in annots.junc_to_gene[j]:
                    if gene_id not in gene_hits:
                        gene_hits[gene_id] = [0, -1 * len(annots.gene_to_annot_juncs[gene_id])]
                    gene_hits[gene_id][0] += 1
    return gene_hits


def get_single_exon_gene_overlaps(iso_readrec, annots):
    gene_hits = {}
    exon = iso_readrec.exons[0]
    index = binary_search(exon, annots.all_annot_SE)
    # FIXME: how does this ever work? all_annot_SE is [(start, end, strand, gene_id), ...]
    for annot_exon_info in annots.all_annot_SE[index - ANNOT_SE_SEARCH_WINDOW:index + ANNOT_SE_SEARCH_WINDOW]:
        # FIXME: make overlap a function
        overlap = min(exon.end, annot_exon_info[1]) - max(exon.start, annot_exon_info[0])
        if overlap > 0:
            # base coverage of long-read isoform by the annotated isoform
            frac_of_iso = float(overlap) / (exon.end - exon.start)
            # base coverage of the annotated isoform by the long-read isoform
            frac_of_annot = float(overlap) / (annot_exon_info[1] - annot_exon_info[0])
            if frac_of_iso > MIN_ISOFORM_OVERLAP_FRAC and frac_of_annot > MIN_ANNOT_OVERLAP_FRAC:
                if annot_exon_info[3] not in gene_hits or frac_of_iso > gene_hits[annot_exon_info[3]][0]:
                    gene_hits[annot_exon_info[3]] = [frac_of_iso, frac_of_annot]
    return gene_hits

def get_spliced_exon_overlaps(strand, exons, annots):
    gene_hits = []
    for annot_gene in annots.spliced_exons[strand]:
        annot_exons = sorted(list(annots.spliced_exons[strand][annot_gene]))
        # check if there is overlap in the genes
        # FIXME: not clear how this checks for overlap
        if (min((annot_exons[-1].end, exons[-1].end)) > max((annot_exons[0].start, exons[0].start))):
            covered_pos = set()
            for ex in exons:
                for aex in annot_exons:
                    for p in range(max((aex.start, ex.start)), min((aex.end, ex.end))):
                        covered_pos.add(p)
            if len(covered_pos) > sum([x.end - x.start for x in exons]) * 0.5:
                gene_hits.append([len(covered_pos), annot_gene, strand])
    return gene_hits

def _get_transcript_gene_from_annot(iso_readrec, annots, annot_name_to_used_counts):
    """Return (transcript_id, gene_id) if iso matches an annotated junction chain, else (None, None)."""
    if iso_readrec.juncs != () and iso_readrec.juncs in annots.juncchain_to_transcript:
        transcript_id, gene_id = annots.juncchain_to_transcript[iso_readrec.juncs]
        if transcript_id in annot_name_to_used_counts:
            annot_name_to_used_counts[transcript_id] += 1
            transcript_id = transcript_id + '-endvar' + str(annot_name_to_used_counts[transcript_id])
        else:
            annot_name_to_used_counts[transcript_id] = 1
        return transcript_id, gene_id
    else:
        return None, None


def _find_gene_id_by_overlap(iso_readrec, annots):
    """Find gene_id for an isoform without a matching junction chain, using junction or exon overlap."""
    if iso_readrec.juncs != ():
        gene_hits = get_genes_with_shared_juncs(iso_readrec.juncs, annots)
    else:
        gene_hits = get_single_exon_gene_overlaps(iso_readrec, annots)
    if gene_hits:
        return sorted(gene_hits.items(), key=lambda x: x[1], reverse=True)[0][0]
    else:
        # look for exon overlap
        if iso_readrec.strand != 'ambig':
            gene_hits = get_spliced_exon_overlaps(iso_readrec.strand, iso_readrec.exons, annots)
        else:
            gene_hits = (get_spliced_exon_overlaps('+', iso_readrec.exons, annots) +
                         get_spliced_exon_overlaps('-', iso_readrec.exons, annots))
        if gene_hits:
            gene_hits.sort(reverse=True)
            if iso_readrec.strand == 'ambig':
                iso_readrec.strand = gene_hits[0][2]
            return gene_hits[0][1]
        else:
            return None


def get_gene_name_firstpass(iso_name, isoform, annots, annot_name_to_used_counts, novel_gene_isos_to_group):
    transcript_id, gene_id = _get_transcript_gene_from_annot(isoform, annots, annot_name_to_used_counts)
    if transcript_id is None:
        transcript_id = isoform.name
        gene_id = _find_gene_id_by_overlap(isoform, annots)
    return gene_id, transcript_id


def get_gene_names_firstpass(firstpass, annots):
    annot_name_to_used_counts = {}
    novel_gene_isos_to_group = {'+': [], '-': []}
    for iso_name in firstpass:
        isoform = firstpass[iso_name]
        gene_id, isoform_id = get_gene_name_firstpass(iso_name, isoform, annots, annot_name_to_used_counts,
                                novel_gene_isos_to_group)
        firstpass[iso_name].transcript_id = isoform_id
        if gene_id is not None:
            # removing this strand correction breaks the unusual junction (due to underlying variant?) test
            firstpass[iso_name].strand = annots.gene_to_strand[gene_id]
            
            firstpass[iso_name].gene_id = gene_id
        else: 
            novel_gene_isos_to_group[isoform.strand].append((isoform.start, isoform.end, iso_name))

    return novel_gene_isos_to_group


def generate_non_gene_iso_groups_strand(novel_gene_isos_to_group, strand, chrom, firstpass):
    transcripts_to_group = sorted(novel_gene_isos_to_group[strand])
    last_end = 0
    group_start = 0
    curr_group = []
    for start, end, iso_name in transcripts_to_group:
        if start < last_end:
            curr_group.append((start, end, iso_name))
        else:
            if len(curr_group) > 0:
                group_name = f'{chrom}:{group_start}-{last_end}:{strand}'
                for s, e, n in curr_group:
                    firstpass[n].gene_id = group_name
            curr_group = [(start, end, iso_name)]
            group_start = start
        if end > last_end:
            last_end = end
    if len(curr_group) > 0:
        group_name = f'{chrom}:{group_start}-{last_end}:{strand}'
        for s, e, n in curr_group:
            firstpass[n].gene_id = group_name

def write_first_pass_isoforms(iso_name, normalize_ends, isoform, max_terminal_exons_ends, add_length_at_ends, unique_bound, unique_fh, iso_fh, seq_fh, genome):
    # FIXME: do normalization outside of write function
    if normalize_ends and len(iso_info.exons) > 1:  # don't normalize ends for single exon transcripts
        normalize_gene_terminal_exons(max_terminal_exons_ends, isoform.gene_id, isoform.strand, isoform.exons,
                                      add_length_at_ends=add_length_at_ends)
        isoform.reset_from_exons(iso_info.exons)
    isoform.name = isoform.transcript_id + '_' + isoform.gene_id

    if unique_bound and iso_name in unique_bound:
        unique_fh.write(isoform.name + '\t' + unique_bound[iso_name] + '\n')

    iso_fh.write('\t'.join(isoform.get_bed_line()) + '\n')
    seq_fh.write('>' + isoform.name + '\n')
    seq_fh.write(isoform.get_sequence(genome) + '\n')

def get_gene_names_and_write_firstpass(temp_prefix, chrom, firstpass, annots, genome, *,
                                       normalize_ends=False, add_length_at_ends=0, unique_bound=None):
    # THIS IS WHERE WE CAN GET GENES AND ADJUST NAMES

    novel_gene_isos_to_group = get_gene_names_firstpass(firstpass, annots)

    # generating non-gene iso groups
    for strand in novel_gene_isos_to_group:
        generate_non_gene_iso_groups_strand(novel_gene_isos_to_group, strand, chrom, firstpass)

    # generating standardized set of ends for gene
    if normalize_ends:
        max_terminal_exons_ends = max_terminal_exons_ends_from_iso_infos(firstpass)
    else:
        # FIXME: passing None is move obvious to flow control,
        # although making write_first_pass_isoforms less monolithic
        # it does more than writing
        max_terminal_exons_ends = {}

    with (open(temp_prefix + '.firstpass.bed', 'w') as iso_fh,
          open(temp_prefix + '.firstpass.fa', 'w') as seq_fh,
          open(temp_prefix + '.firstpass.uniquebound.txt', 'w') as unique_fh):
        for iso_name in firstpass:
            write_first_pass_isoforms(iso_name, normalize_ends, firstpass[iso_name], max_terminal_exons_ends,
                                      add_length_at_ends, unique_bound, unique_fh, iso_fh, seq_fh, genome)


def get_reverse_complement(seq):
    compbase = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N',
                'R': 'Y', 'Y': 'R', 'K': 'M', 'M': 'K', 'S': 'S',
                'W': 'W', 'B': 'V', 'V': 'B', 'D': 'H', 'H': 'D'}
    seq = seq.upper()
    new_seq = []
    for base in seq:
        new_seq.append(compbase[base])
    return ''.join(new_seq[::-1])

####
# results output
####
# def get_bed_gtf_from_info(end_info, chrom, strand, juncs, gene_id, genome):
#     # FIXME: what is end info??
#     # build gtf, bed, fasta data
#     start = end_info.start
#     end = end_info.end
#     iso_id_src = end_info.iso_id_src
#     read_names = end_info.read_names
#     score = min(len(read_names), 1000)
#     exon_starts, exon_sizes = get_bed_exons_from_juncs(juncs, start, end)
#     # FIXME: duplicate code to build BED record
#     bed_line = [chrom, start, end, iso_id_src.id + '_' + gene_id, score, strand, start, end,
#                 get_rgb(strand, len(juncs)), len(exon_starts), ','.join([str(x) for x in exon_sizes]),
#                 ','.join([str(x) for x in exon_starts])]
#     exons = [Exon(start + exon_starts[i], start + exon_starts[i] + exon_sizes[i])
#              for i in range(len(exon_starts))]
#     trans_seq = get_sequence_for_exons(genome, chrom, strand, exons)
#     if strand == '-':
#         exons.reverse()

#     # Create GtfTranscript with exons
#     gtf_transcript = GtfTranscript(chrom, 'FLAIR', 'transcript', start, end, score, strand, '.',
#                                    gene_id=gene_id, transcript_id=iso_id_src.id)
#     for i, exon in enumerate(exons, start=1):
#         gtf_exon = GtfExon(chrom, 'FLAIR', 'exon', exon.start, exon.end, score, strand, '.',
#                            gene_id=gene_id, transcript_id=iso_id_src.id, exon_number=i)
#         gtf_transcript.add_exon(gtf_exon)

#     return '\t'.join([str(x) for x in bed_line]) + '\n', gtf_transcript, trans_seq

# def write_iso_seq_map(iso_info, name_to_used_counts, chrom, strand, juncs, gene_id, genome, iso_fh, t_starts, t_ends, gtf_transcripts, map_fh,
#                       read_to_final_transcript, counts_fh, seq_fh):
#     iso_id_src = iso_info.iso_id_src
#     iso_id = iso_id_src.id.split('-endvar')[0]  # FIXME what is this all about?
#     if iso_id in name_to_used_counts:
#         name_to_used_counts[iso_id] += 1
#         iso_id = iso_id + '-endvar' + str(name_to_used_counts[iso_id])
#     else:
#         name_to_used_counts[iso_id] = 1
#     iso_info.iso_id_src = IsoIdSrc(iso_id, iso_id_src.src)
#     bed_line, gtf_data, tseq = get_bed_gtf_from_info(iso_info, chrom, strand, juncs, gene_id, genome)
#     iso_fh.write(bed_line)
#     t_starts.append(iso_info.start)
#     t_ends.append(iso_info.end)
#     gtf_transcripts.append(gtf_data)
#     map_fh.write(iso_id + '_' + gene_id + '\t' + ','.join(iso_info.read_names) + '\n')
#     for r in iso_info.read_names:
#         read_to_final_transcript[r] = (iso_id + '_' + gene_id, chrom, strand)
#     counts_fh.write(iso_id + '_' + gene_id + '\t' + str(len(iso_info.read_names)) + '\n')
#     seq_fh.write('>' + iso_id + '_' + gene_id + '\n')
#     seq_fh.write(tseq + '\n')

# def calculate_gene_total_reads(gene_isoform_data):
#     """Calculate total read support across all isoforms in a gene."""
#     gene_tot = 0
#     for junc_chain in gene_isoform_data.junction_chains():
#         for iso_info in gene_isoform_data.get_isoforms(junc_chain.chrom, junc_chain.strand, junc_chain.juncs):
#             gene_tot += len(iso_info.read_names)
#     return gene_tot

# def write_gene_isoforms(gene_isoform_data, gene_id, gene_tot, args, genome, iso_fh, map_fh, read_to_final_transcript, counts_fh, seq_fh):
#     """Write isoform data (BED, sequences, maps) for isoforms meeting the threshold.
#     Returns tuple of (gtf_transcripts, t_starts, t_ends) for GTF writing."""
#     gtf_transcripts, t_starts, t_ends = [], [], []

#     for junc_chain in gene_isoform_data.junction_chains():
#         ends_list = gene_isoform_data.get_isoforms(junc_chain.chrom, junc_chain.strand, junc_chain.juncs)

#         name_to_used_counts = {}
#         for iso_info in ends_list:
#             if len(iso_info.read_names) / gene_tot >= args.frac_support:
#                 write_iso_seq_map(iso_info, name_to_used_counts, junc_chain.chrom, junc_chain.strand, junc_chain.juncs,
#                                   gene_id, genome, iso_fh, t_starts, t_ends, gtf_transcripts, map_fh,
#                                   read_to_final_transcript, counts_fh, seq_fh)

#     return gtf_transcripts, t_starts, t_ends

# def write_gene_gtf(gtf_transcripts, t_starts, t_ends, gene_id, gtf_fh):
#     """Write GTF records for a gene (gene record, transcript records, and exon records)."""
#     first_transcript = gtf_transcripts[0]
#     gtf_write_row(gtf_fh, first_transcript.chrom, 'FLAIR', 'gene', min(t_starts), max(t_ends), '.',
#                   first_transcript.strand, '.', gene_id=gene_id)

#     # Write transcript and exon records
#     for gtf_transcript in gtf_transcripts:
#         # Write transcript record
#         print(gtf_transcript, file=gtf_fh)

#         # Write exon records
#         for gtf_exon in gtf_transcript.exons:
#             print(gtf_exon, file=gtf_fh)

# def write_gene_output(gene_to_juncs_to_ends, gene_id, args, genome, iso_fh, map_fh, read_to_final_transcript, counts_fh, seq_fh, gtf_fh):
#     """Write all output files for a single gene (isoforms, sequences, maps, and GTF)."""
#     gene_isoform_data = gene_to_juncs_to_ends[gene_id]

#     # Calculate total read support for the gene
#     gene_tot = calculate_gene_total_reads(gene_isoform_data)

#     # Write isoform data and collect GTF information
#     gtf_transcripts, t_starts, t_ends = write_gene_isoforms(
#         gene_isoform_data, gene_id, gene_tot, args, genome,
#         iso_fh, map_fh, read_to_final_transcript, counts_fh, seq_fh)

#     # Write GTF records
#     if gtf_transcripts:
#         write_gene_gtf(gtf_transcripts, t_starts, t_ends, gene_id, gtf_fh)

def get_transcirpts_to_reads(temp_prefix, suffix):
    transcript_to_reads = {}
    for line in open(temp_prefix + suffix):
        read_name, transcript_id, start, end = line.rstrip().split('\t')
        if transcript_id not in transcript_to_reads:
            transcript_to_reads[transcript_id] = []
        transcript_to_reads[transcript_id].append((read_name, start, end))
    return transcript_to_reads

def write_transcript_ends_bed(args, temp_prefix, suffix, read_to_final_transcript, ends_fh):
    transcript_to_reads = get_transcirpts_to_reads(temp_prefix, suffix)
    for t in transcript_to_reads:
        if len(transcript_to_reads[t]) >= args.sjc_support:  # FIXME this needs to be adjusted to consider single exons vs junction chains, also frac_support
            for r, start, end in transcript_to_reads[t]:
                if r in read_to_final_transcript:
                    t_name, chrom, strand = read_to_final_transcript[r]
                    ends_fh.write('\t'.join([chrom, start, end, t_name + '|' + r, '.', strand]) + '\n')

def write_transcript_ends_beds(args, temp_prefix, read_to_final_transcript, ends_fh):
    # FIXME: these are not real TSVs
    for suffix in ['.matchannot.ends.tsv']:
        write_transcript_ends_bed(args, temp_prefix, suffix, read_to_final_transcript, ends_fh)

def generate_genomic_alignment_read_to_clipping_file(temp_prefix, bam_file, region):
    c = 0
    with open(temp_prefix + '.reads.genomicclipping.txt', 'w') as clipping_fh:
        for read in bam_file.fetch(region.name, region.start, region.end):
            if not read.is_secondary and not read.is_supplementary:
                c += 1
                name = read.query_name
                cigar = read.cigartuples
                tot_clipped = 0
                if cigar[0][0] in {4, 5}:
                    tot_clipped += cigar[0][1]
                if cigar[-1][0] in {4, 5}:
                    tot_clipped += cigar[-1][1]
                clipping_fh.write(name + '\t' + str(tot_clipped) + '\n')
    return c, temp_prefix + '.reads.genomicclipping.txt'


def predict_productivity(out_prefix, genome_fasta, gtf):
    cmd = ('predictProductivity',
           '-i', out_prefix + '.isoforms.bed',
           '-o', out_prefix + '.isoforms.CDS',
           '--gtf', gtf,
           '--genome_fasta', genome_fasta,
           '--longestORF')
    pipettor.run(cmd)

def _iso_passes_support_filter(args, iso, num_exons, iso_to_counts, gene_to_tot):
    if iso not in iso_to_counts:
        return False
    else:
        count = iso_to_counts[iso]
        min_support = args.sjc_support if num_exons > 1 else args.se_support
        return (count >= min_support) and (count / gene_to_tot[iso.split('_')[-1]]) >= args.frac_support


def _run_region(*, partition, gtf_data, intron_support, args):
    region = partition.region
    annots = annot_data_from_gtf(gtf_data, region)

    # first extract reads for region as fasta
    pipettor.run([('samtools', 'view', '-h', args.genome_aligned_bam, region.name + ':' + str(region.start) + '-' + str(region.end)),
                  ('samtools', 'fasta', '-')],
                 stdout=partition.output_path('reads.fasta'))

    # then align reads to transcriptome and run count_sam_transcripts
    genome = pysam.FastaFile(args.genome)
    bam_file = pysam.AlignmentFile(args.genome_aligned_bam, 'rb')

    # genomic clipping: amount of clipping (from cigar) at ends of reads when aligned to genome
    # generates file with [read{\t}clipping amount] on each line
    # For comparing with amount of clipping after alignment to transcriptome
    # in order to check whether transcriptome alignment is comparable to or better than genomic alignment,
    # which can be considered to support isoform.
    # used in filter_transcriptome_align

    # logging.info('generating genomic clipping reference')
    num_reads, clipping_file = generate_genomic_alignment_read_to_clipping_file(partition.file_prefix, bam_file, region)

    if num_reads == 0:
        suffixes = ['.firstpass.reallyunfiltered.bed', '.firstpass.unfiltered.bed', '.firstpass.bed',
                        '.isoforms.bed',
                        '.isoform.read.map.txt', '.isoforms.gtf', '.isoforms.fa', '.isoform.counts.txt']
        if args.end_norm_dist:
            suffixes.extend(['.read_ends.bed', '.matchannot.ends.tsv'])
        if not args.no_align_to_annot:
            suffixes.extend(['.matchannot.counts.txt', '.matchannot.read.map.txt'])
        for s in suffixes:
            out = open(partition.file_prefix + s, 'w')
            out.close()
        return

    # aligning to reference transcriptome, then identifying reads that match well to reference transcripts
    # with filter_transcriptome_align
    # logging.info('identifying good match to annot')
    read_to_annot_transcript = \
        identify_good_match_to_annot(args, partition.file_prefix, region.name, annots, genome)

    # load splice junctions for chrom
    # logging.info('correcting splice junctions')
    junction_corrector = setup_junction_corrector(intron_support, args.ss_window, args.junction_support)

    # takes in bam file, for each read attempts to correct splice junctions (removes unsupported ones), then groups reads by junction chains
    # this also handles read strandedness if necessary
    sj_to_ends = filter_correct_group_reads(args, partition.file_prefix, region, bam_file, read_to_annot_transcript, annots,
                                            junction_corrector)
    bam_file.close()

    # for each junction chain, clusters ends - generates junction chain x ends
    # firstpass objects then does initial filtering by read support and
    # redundant ends also separates single exon isoforms from spliced isoforms
    # (because they're handled differently in future step for identifying
    # annotated gene/isoform names)
    firstpass_unfiltered, firstpass_junc_to_name, firstpass_exons = process_juncs_to_firstpass_isos(args, partition.file_prefix, sj_to_ends)

    # - filter isoforms - remove any that represent a subset of another
    # - identified isoform - based on what args.filter is set to also generate
    # - iso_to_unique_bound - a mapping of each isoform to the unique sequence
    #   at its ends (this is to better handle isoforms that represent junction
    #   subsets with additional sequence at the ends)
    firstpass, iso_to_unique_bound = filter_firstpass_isos(args, firstpass_unfiltered, firstpass_junc_to_name, firstpass_exons,
                                                           annots, {})
    if len(firstpass.keys()) > 0:
        # logging.info('getting gene names and writing firstpass')
        # this section identifies annotated gene and isoform names (primarily based on splice junction matching, secondarily by exon overlap)
        # also adjusts isoform strand, determines novel isoform and gene names
        # also normalizes transcript ends (temporarily extends ends so that transcript end alignment does not drive transcript assignment during transcriptome alignment)
        # writes out bed and fa files
        if args.end_norm_dist is not None:
            get_gene_names_and_write_firstpass(partition.file_prefix, region.name, firstpass, annots, genome,
                                               normalize_ends=True, add_length_at_ends=args.end_norm_dist, unique_bound=iso_to_unique_bound)
        else:
            get_gene_names_and_write_firstpass(partition.file_prefix, region.name, firstpass, annots, genome, unique_bound=iso_to_unique_bound)
            # logging.info('identifying good match to firstpass')

        # aligns to firstpass transcriptome, identifies best read -> isoform alignment for each read, then gets read counts per isoform
        transcriptome_align_and_count(args, partition.output_path('reads.fasta'),
                                      partition.output_path('firstpass.fa'),
                                      partition.output_path('firstpass.bed'),
                                      partition.output_path('isoform.counts.txt'),
                                      partition.output_path('isoform.read.map.txt'), False,
                                      partition.output_path('reads.genomicclipping.txt'),
                                      partition.output_path('firstpass.uniquebound.txt'))
    else:
        # create empty files
        with open(partition.output_path('firstpass.fa'), 'w') as _, \
             open(partition.output_path('firstpass.bed'), 'w') as _, \
             open(partition.output_path('isoform.counts.txt'), 'w') as _, \
             open(partition.output_path('isoform.read.map.txt'), 'w') as _:
            pass
        if args.output_endpos:
            with open(partition.output_path('ends.tsv'), 'w') as _:
                pass

    iso_to_counts = {}
    gene_to_tot = {}
    for line in open(partition.output_path('isoform.counts.txt')):
        iso, counts = line.rstrip().split('\t')
        gene = iso.split('_')[-1]
        if gene not in gene_to_tot:
            gene_to_tot[gene] = 0
        gene_to_tot[gene] += int(counts)
        iso_to_counts[iso] = int(counts)
    with open(partition.output_path('isoforms.bed'), 'w') as out:
        for line in open(partition.output_path('firstpass.bed')):
            temp = line.split('\t')
            iso = temp[3]
            num_exons = int(temp[-3])
            if _iso_passes_support_filter(args, iso, num_exons, iso_to_counts, gene_to_tot):
                out.write(line)
    bed_to_gtf(partition.output_path('isoforms.bed'), partition.output_path('isoforms.gtf'))
    bed_to_sequence(partition.output_path('isoforms.bed'), args.genome, partition.output_path('isoforms.fa'))

    genome.close()

def combine_chunks(args, output, partitions):
    files_to_combine = ['.firstpass.reallyunfiltered.bed', '.firstpass.unfiltered.bed', '.firstpass.bed',
                        '.isoforms.bed',
                        '.isoform.read.map.txt', '.isoforms.gtf', '.isoforms.fa', '.isoform.counts.txt']
    if args.end_norm_dist:
        files_to_combine.extend(('.read_ends.bed', '.matchannot.ends.tsv'))
    # if args.predict_cds:
    #     files_to_combine.append('.isoforms.CDS.bed')
    if not args.no_align_to_annot:
        files_to_combine.extend(['.matchannot.counts.txt', '.matchannot.read.map.txt'])
    combine_temp_files_by_suffix(output, [p.file_prefix for p in partitions], files_to_combine)

####
# main
####

def flair_transcriptome():
    # FIXME: split options out that are flags to indicate what to do
    # so args doesn't get passes but we don't have to pass so many options

    args = get_args()

    logging.info('loading genome')
    genome = pysam.FastaFile(args.genome)

    temp_dir = make_temp_dir(args.output)

    logging.info('building intron support database')
    is_db = build_intron_support(args.annot_gtf, args.junction_tab, args.junction_bed)

    annot_gtf_data = None
    if args.annot_gtf:
        logging.info('loading annotation GTF')
        annot_gtf_data = gtf_data_parser(args.annot_gtf, attrs=GtfAttrsSet.FLAIR, include_features=TRANSCRIPT_EXON_FEATURES)

    logging.info('partitioning genome')
    runner = partition_runner_factory(args.parallel_mode, genome, args.genome_aligned_bam,
                                      temp_dir, args.annot_gtf, args.threads,
                                      gtf_data=annot_gtf_data, intron_support=is_db)
    logging.info(f'number of partitions: {len(runner)}')

    logging.info('running partitions')
    runner.run(_run_region, args=args)
    combine_chunks(args, args.output, runner.partitions)

    # this needs to be done here outside of chunking because needs to load whole annotation gtf, which should only be done once
    if args.predict_cds:
        logging.info('predicting CDS')
        # FIXME: why is this passing in annotation GTF?
        predict_productivity(args.output, args.genome, args.annot_gtf)

    if not args.keep_intermediate:
        shutil.rmtree(temp_dir)

    genome.close()


if __name__ == "__main__":
    flair_transcriptome()
