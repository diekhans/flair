"""Read/exon structures and junction utilities, shared across FLAIR modules."""

from collections import namedtuple
from flair import PosRange
from flair.pycbio.hgdata.bed import Bed


####
# basic types
####
class Junc(PosRange):
    """Stores start, end, just adds a type name to PosRange for clearer code and error messages"""
    pass


class Exon(PosRange):
    """Stores start, end, just adds a type name to PosRange for clearer code and error messages"""
    pass


ISO_SRC_ANNOT = 'annot'
ISO_SRC_NOVEL = 'novel'


class IsoIdSrc(namedtuple("IsoIdSrc",
                          ("id", "src"))):
    """isoform identifier along with the source of the isoform"""
    # FIXME: it is unclear if this is the best way to store the information,
    # this was create as a transition from iso (id) or (iso_id) (marker, id)
    pass


def exons_to_juncs(exons):
    """Convert exon ranges to junctions"""
    return [Junc(exons[i].end, exons[i + 1].start)
            for i in range(len(exons) - 1)]


def bed_to_junctions(bed):
    # FIXME: a junctions object might be good
    return [Junc(bed.blocks[i - 1].end, bed.blocks[i].start)
            for i in range(1, bed.blockCount)]


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


def get_reverse_complement(seq):
    compbase = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N',
                'R': 'Y', 'Y': 'R', 'K': 'M', 'M': 'K', 'S': 'S',
                'W': 'W', 'B': 'V', 'V': 'B', 'D': 'H', 'H': 'D'}
    seq = seq.upper()
    new_seq = []
    for base in seq:
        new_seq.append(compbase[base])
    return ''.join(new_seq[::-1])


def get_sequence_for_exons(genome, chrom, strand, exons):
    trans_seq = ''.join([genome.fetch(chrom, e.start, e.end)
                         for e in exons])
    # FIXME: this upper cases only if reverse strand
    if strand == '-':
        trans_seq = get_reverse_complement(trans_seq)
    return trans_seq


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
        if data[i].end < query.start:
            lower = i
            i = int(round((i + upper) / 2))
        elif data[i].start > query.end:
            upper = i
            i = int(round((lower + i) / 2))
        else:  # found
            break
    return i


class ReadRec:
    """Read alignment with location, junction, and metadata fields.

    Stores chrom, start, end, name, score, strand, and juncs directly.
    Exons are computed on the fly from start, end, and juncs.
    A Bed record can be produced on demand via to_bed().

    Junction tuples are interned via a class-level cache so identical junction
    chains share a single tuple object.
    """

    # keep on one copy of junction chain
    _juncs_cache = {}

    @classmethod
    def _intern_juncs(cls, juncs):
        return cls._juncs_cache.setdefault(juncs, juncs)

    def __init__(self, chrom, start, end, name, score, strand, juncs):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.name = name
        self.score = score
        self.strand = strand
        self.juncs = self._intern_juncs(juncs)

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
        return cls(read.reference_name, align_start, ref_pos, read.query_name,
                   read.mapping_quality, junc_direction, juncs)

    @classmethod
    def from_junctions(cls, chrom, start, end, name, score, strand, juncs):
        """Create a ReadRec from junction coordinates."""
        return cls(chrom, start, end, name, score, strand, tuple(juncs))

    @property
    def exons(self):
        """Return exons as list of Exon objects, computed from start, end, and juncs."""
        if not self.juncs:
            return [Exon(self.start, self.end)]
        exons = [Exon(self.start, self.juncs[0].start)]
        for i in range(len(self.juncs) - 1):
            exons.append(Exon(self.juncs[i].end, self.juncs[i + 1].start))
        exons.append(Exon(self.juncs[-1].end, self.end))
        return exons

    def to_bed(self):
        """Create and return a Bed object."""
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
        return get_sequence_for_exons(genome, self.chrom, self.strand, self.exons)

    def reset_from_exons(self, exons):
        """Update ReadRec from a list of Exon objects."""
        self.start = exons[0].start
        self.end = exons[-1].end
        self.juncs = self._intern_juncs(tuple(exons_to_juncs(sorted(exons))))

    def update_from_juncs(self, new_juncs):
        """Update juncs, keeping chrom, start, end, name, score, strand."""
        self.juncs = self._intern_juncs(tuple(new_juncs))
