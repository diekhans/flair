"""
Annotation and orthogonal intron support.
"""
import sys
from collections import defaultdict
from intervaltree import IntervalTree
from flair import MIN_INTRON_SIZE, MAX_INTRON_SIZE, FlairInputDataError
from flair.ssUtils import gtfToSSBed
from flair.pycbio.hgdata.bed import BedReader
from flair.pycbio.tsv import TsvReader

class Intron:
    """Coordinates and support.  Strand maybe None if not available"""
    __slots__ = ("chrom", "start", "end", "strand",
                 "annot_supported", "read_supported", "read_support_cnt")

    def __init__(self, chrom, start, end, strand):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand
        self.annot_supported = self.read_supported = False
        self.read_support_cnt = 0

    def __str__(self):
        return (f"Intron({self.chrom}:{self.start}-{self.end}({self.strand}) "
                f"annot={self.annot_supported} read={self.read_supported} read_cnt={self.read_support_cnt}")

class IntronSupport:
    """
    Table of intron support index by both start and end positions
    """
    def __init__(self, *, min_intron_size=MIN_INTRON_SIZE, max_intron_size=MAX_INTRON_SIZE):
        # dict index by chrom of interval trees, keyed on first base of donor and last base of the acceptor sites
        self.coords_maps = defaultdict(IntervalTree)
        self.min_intron_size = min_intron_size
        self.max_intron_size = max_intron_size

    def _find_point_strand(self, chrom, point, strand):
        for entry in self.coords_maps[chrom][point]:
            if entry.data.strand == strand:
                return entry.data
        return None

    def _find_intron(self, chrom, start, end, strand):
        donor = self._find_point_strand(chrom, start, strand)
        if donor is not None:
            accept = self._find_point_strand(chrom, end - 1, strand)
            if accept is donor:
                return donor  # same intron
        return None

    def _add_intron(self, chrom, start, end, strand, read_count):
        assert start < end
        intron = Intron(chrom, start, end, strand)
        self.coords_maps[chrom].addi(start, start + 1, intron)
        self.coords_maps[chrom].addi(end - 1, end, intron)
        return intron

    def _add_support(self, chrom, start, end, strand, read_count):
        intron = self._find_intron(chrom, start, end, strand)
        if intron is None:
            intron = self._add_intron(chrom, start, end, strand, read_count)
        if read_count is None:
            intron.annot_supported = True
        else:
            intron.read_supported = True
            intron.read_support_cnt += max(read_count, 1)

    def add_support(self, chrom, start, end, strand, read_count=None):
        """A read_count None indicates annot_support.  Drop intron outside of
        configured size range"""
        if self.min_intron_size <= (end - start) <= self.max_intron_size:
            self._add_support(chrom, start, end, strand, read_count)
            return True
        return False

    def overlap(self, chrom, start, end, margin=0):
        """Get list of overlapping introns where either ends overlaps this range with
        margin padding."""
        return [entry.data for entry in self.coords_maps[chrom].overlap(start - margin, end + margin)]

    def overlap_introns(self, chrom, start, end, margin=0):
        """get introns were splice junctions overlap each end of this range,
        with margin bases on either of ends of the range"""
        # only return introns that hit both ends
        starts = self.overlap(chrom, start, start + 1, margin)
        ends = self.overlap(chrom, end - 1, end, margin)
        ends_ids = set([id(e) for e in ends])
        return [intron for intron in starts
                if id(intron) in ends_ids]

    def dump(self, fh=sys.stderr):
        print("IntronSupport:", file=fh)
        for chrom, coords_map in self.coords_maps.items():
            for entry in coords_map.items():
                print(f"\t{chrom}:{entry.begin}:{entry.end}", entry.data, file=fh)


def _load_read_bed_intron(intron_support, bed):
    if not (6 <= bed.numStdCols <= 9):
        raise FlairInputDataError(f"intron BED must have 6 to 9 columns, found {bed.numStdCols}")
    if bed.strand not in ('+', '-'):
        raise FlairInputDataError(f"Invalid strand '{bed.strand}' in intron BED must be '+', '-', or '.'")
    strand = bed.strand if bed.strand != '.' else None
    return intron_support.add_support(bed.chrom, bed.chromStart, bed.chromEnd, strand, bed.score)

def load_read_bed_introns(intron_support, bed_file):
    """load introns from a BED6, with score being the counts of reads.
    Return number of introns loaded."""
    cnt = 0
    line_num = 0
    try:
        for bed in BedReader(bed_file):
            line_num += 1
            if _load_read_bed_intron(intron_support, bed):
                cnt += 1
    except Exception as exc:
        raise FlairInputDataError(f"parsing intron BED failed: {bed_file} line {line_num}") from exc
    if cnt == 0:
        raise FlairInputDataError(f"No introns loaded from BED: {bed_file}")
    return cnt

def _load_annot_gtf_introns(intron_support, juncs):
    cnt = 0
    for chrom in juncs.keys():
        for start, end, strand in juncs[chrom]:
            if intron_support.add_support(chrom, start, end, strand):
                cnt += 1
    return cnt

def load_annot_gtf_introns(intron_support, gtf_file):
    """Load introns from a annotation GTF. Return number of introns loaded."""
    try:
        knownSS = {}  # flotsam
        juncs, _, _ = gtfToSSBed(gtf_file, knownSS, False, None, False)
        cnt = _load_annot_gtf_introns(intron_support, juncs)
    except Exception as exc:
        raise FlairInputDataError(f"parsing annotation GTF failed: {gtf_file}") from exc
    if cnt == 0:
        raise FlairInputDataError(f"No introns loaded from GTF: {gtf_file}")
    return cnt

def _load_star_intron(intron_support, rec):
    strand = (None, '+', '-')[rec.strand]
    return intron_support.add_support(rec.chrom, rec.start - 1, rec.end, strand, rec.uniq_map_cnt)

def load_read_star_introns(intron_support, sj_file):
    """load introns from STAR junction file"""
    columns = ("chrom", "start", "end", "strand", "motif", "annot", "uniq_map_cnt",
               "multi_map_cnt", "max_overhang")
    cnt = 0
    line_num = 1  # include header
    try:
        for rec in TsvReader(sj_file, columns=columns, typeMap={"chrom": str}, defaultColType=int):
            line_num += 1
            if _load_star_intron(intron_support, rec):
                cnt += 1
    except Exception as exc:
        raise FlairInputDataError(f"parsing STAR SJ failed: {sj_file} line {line_num}") from exc
    if cnt == 0:
        raise FlairInputDataError(f"No introns loaded from STAR SJ file: {sj_file}")
    return cnt
