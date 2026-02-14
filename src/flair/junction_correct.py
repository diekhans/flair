"""
Correction of read splice junctions from external evidence.
"""
import logging
import copy
from math import inf
from flair import PosRange
from flair.pycbio.hgdata.bed import Bed

##
# somewhat arbitrary sizes to keep from going off ends
# or overlapping other introns.
##
MIN_INTERNAL_EXON_SIZE = 3   # there is one this small!
MIN_TERMINAL_EXON_SIZE = 32

class JunctionCorrector:
    """Correction of read splice sites from orthogonal evidence

       * flank_window - the number of based +/- a read junction to search for
         an supporting intron match.
       * min_read_support - minimum of reads to support a junction.
    """

    def __init__(self, intron_support, flank_window, min_read_support):
        self.intron_support = intron_support
        self.flank_window = flank_window
        self.min_read_support = min_read_support

    @property
    def chroms(self):
        return self.intron_support.chroms

    def overlap_introns(self, chrom, start, end, strand):
        def _filter_intron(intron):
            return ((intron.strand == strand) and
                    (intron.annot_supported or (intron.read_support_cnt > self.min_read_support)))
        return list(filter(_filter_intron,
                           self.intron_support.overlap_introns(chrom, start, end, self.flank_window)))

    def correct_read_junctions(self, read_bed):
        """correct a read based on support from introns.  Return None if there
        is no support for an intron."""
        assert len(read_bed.blocks) > 0
        return _correct_junctions(self, read_bed)

    def correct_read_bed(self, read_bed):
        """correct a read based on support from introns.  Return None if there
        is no support for an intron. Return a new BED"""
        new_junctions = self.correct_read_junctions(read_bed)
        if new_junctions is None:
            return None
        return _build_corrected_read_bed(read_bed, new_junctions)

###
# intron support search
###
def _calc_possible_junction_range(read_bed, new_junctions):
    """prevent going off ends of read or overlapping small exons"""
    min_start = read_bed.chromStart + MIN_TERMINAL_EXON_SIZE
    if len(new_junctions) > 0:
        # adjust for previous intron
        min_start = max(min_start, new_junctions[-1].end + MIN_INTERNAL_EXON_SIZE)
    max_end = read_bed.chromEnd - MIN_INTERNAL_EXON_SIZE
    return (min_start, max_end)

def _filter_too_close(read_bed, new_junctions, intron_hits):
    """drop introns overlapping the previous intron or making a too short
    an exon at ends"""
    min_start, max_end = _calc_possible_junction_range(read_bed, new_junctions)
    return list(filter(lambda ih: (ih.start >= min_start) and (ih.end <= max_end),
                       intron_hits))

def _find_best_hit(strand, start, end, intron_hits):
    closest_introns = _collect_closest_hits(start, end, intron_hits)
    if len(closest_introns) > 1:
        # sort to prefer annotated or then more reads
        closest_introns.sort(key=lambda ci: (ci.annot_supported, ci.read_support_cnt), reverse=True)
    return closest_introns[0]

def _collect_closest_hits(start, end, intron_hits):
    """Return list by introns with closest total distance from ends.
    Multiple are return"""
    # don't want to prefer annotated due to NAGNAG junctions
    min_dist = inf
    closest_introns = []
    for intron in intron_hits:
        dist = abs(start - intron.start) + abs(end - intron.end)
        if dist < min_dist:
            # new minimum
            closest_introns.clear()
            min_dist = dist
        if dist <= min_dist:
            closest_introns.append(intron)
    return closest_introns

def _correct_junction(corrector, read_bed, start, end, new_junctions):
    """add and update an intron junctions. Return False if any are not supported."""

    intron_hits = corrector.overlap_introns(read_bed.chrom, start, end, read_bed.strand)
    if intron_hits is None:
        logging.debug(f"No intron support for '{read_bed.name}' {read_bed.chrom}:{start}-{end}")
        return False
    intron_hits = _filter_too_close(read_bed, new_junctions, intron_hits)
    if len(intron_hits) == 0:
        logging.debug("Supporting introns too close to ends or another intron for "
                      f"'{read_bed.name}' {read_bed.chrom}:{start}-{end}")
        return False
    best_intron = _find_best_hit(read_bed.strand, start, end, intron_hits)
    new_junctions.append(PosRange(best_intron.start, best_intron.end))
    return True

def _correct_junctions(corrector, read_bed):
    """create a list of new introns for a read, or None if can't be correct"""
    new_junctions = []
    prev_end = read_bed.blocks[0].end
    for blk in read_bed.blocks[1:]:
        if not _correct_junction(corrector, read_bed, prev_end, blk.start,
                                 new_junctions):
            return None
        prev_end = blk.end
    return new_junctions

###
# create correct bed
###
def _build_corrected_read_bed(read_bed, new_junctions):
    new_bed = Bed(read_bed.chrom, read_bed.chromStart, read_bed.chromEnd,
                  read_bed.name, score=read_bed.score, strand=read_bed.strand,
                  thickStart=read_bed.thickStart, thickEnd=read_bed.thickEnd,
                  itemRgb=read_bed.itemRgb,
                  extraCols=copy.deepcopy(read_bed.extraCols),
                  numStdCols=read_bed.numStdCols)
    prev_end = read_bed.chromStart
    for junction in new_junctions:
        new_bed.addBlock(prev_end, junction.start)
        prev_end = junction.end
    new_bed.addBlock(prev_end, read_bed.chromEnd)
    return new_bed
