"""Shared read-processing logic for FLAIR modules."""

import pysam
from flair.isoform_data import ReadRec, Junc


def should_process_read(read, region, min_quality, keep_sup, allow_secondary):
    """Check if read passes filtering criteria for processing"""
    if read.mapping_quality < min_quality:
        return False
    if read.is_secondary and not allow_secondary:
        return False
    if read.is_supplementary and not keep_sup:
        return False
    if read.reference_name != region.name:
        return False
    if not (region.start <= read.reference_start and read.reference_end <= region.end):
        return False
    return True


def add_corrected_read_to_groups(corrected_read, sj_to_ends):
    """Add a corrected read to the junction-to-ends mapping"""
    junc_key = tuple(sorted(corrected_read.juncs))
    if junc_key not in sj_to_ends:
        sj_to_ends[junc_key] = []
    sj_to_ends[junc_key].append(corrected_read)


def read_correct_to_readrec(junction_corrector, read):
    # FIXME: remove unnecessary initial build of ReadRec and make junctions from
    # read, correct, and then make bed
    readrec = ReadRec.from_read(read)
    corrected_bed = junction_corrector.correct_read_bed(readrec.to_bed())
    if corrected_bed is None:
        return None
    readrec.juncs = ReadRec._intern_juncs(tuple(Junc(corrected_bed.blocks[i].end, corrected_bed.blocks[i + 1].start)
                                                for i in range(len(corrected_bed.blocks) - 1)))
    return readrec


def generate_genomic_alignment_read_to_clipping_file(temp_prefix, bam_file, region):
    with open(temp_prefix + '.reads.genomicclipping.txt', 'w') as clipping_fh:
        for read in bam_file.fetch(region.name, region.start, region.end):
            if not read.is_secondary and not read.is_supplementary:
                name = read.query_name
                cigar = read.cigartuples
                tot_clipped = 0
                if cigar[0][0] in {4, 5}:
                    tot_clipped += cigar[0][1]
                if cigar[-1][0] in {4, 5}:
                    tot_clipped += cigar[-1][1]
                clipping_fh.write(name + '\t' + str(tot_clipped) + '\n')
    return temp_prefix + '.reads.genomicclipping.txt'
