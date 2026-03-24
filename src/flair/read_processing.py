"""Shared read-processing logic for FLAIR modules."""

import pysam
import pipettor
from flair.isoform_data import ReadRec, Junc, IsoWithReads, convert_to_bed


def should_process_read(read, region, min_quality, keep_sup, allow_secondary, allow_outside_range=False):
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


def get_sequence_from_bed(genome, input_bed, output_fa):
    bed_cmd = ('bedtools','getfasta','-nameOnly', '-s', '-split',
                '-fi', genome,
                '-bed',input_bed, 
                '-fo', output_fa)
    pipettor.run([bed_cmd])
    out = open(output_fa.split('.fa')[0] + '.fixed.fa', 'w')
    for line in open(output_fa):
        if line[0] == '>':
            line = line.split('(')[0] + '\n'
        out.write(line)
    out.close()
    pipettor.run([('mv', output_fa.split('.fa')[0] + '.fixed.fa', output_fa)])


def add_corrected_read_to_groups(corrected_read, sj_to_ends):
    """Add a corrected read to the junction-to-ends mapping"""
    junc_key = tuple(sorted(corrected_read.juncs)) # FIXME add chromosome and strand to key
    # FIXME check to see if a given junction chain is on both strands, throw error
    if junc_key not in sj_to_ends:
        sj_to_ends[junc_key] = IsoWithReads.from_readrec(corrected_read)
    sj_to_ends[junc_key].reads.append(corrected_read)


def read_correct_to_readrec(junction_corrector, readrec):
    # FIXME: remove unnecessary initial build of ReadRec and make junctions from
    # read, correct, and then make bed
    corrected_bed = junction_corrector.correct_read_bed(convert_to_bed(readrec))
    if corrected_bed is None:
        return None
    readrec.juncs = ReadRec._intern_juncs(tuple(Junc(corrected_bed.blocks[i].end, corrected_bed.blocks[i + 1].start)
                                                for i in range(len(corrected_bed.blocks) - 1)))
    return readrec


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
