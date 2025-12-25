"""
Simplistic, non-validating GTF parser.
"""
import re
from typing import Optional
from collections import defaultdict
from flair.pycbio.sys import fileOps
from flair import SeqRange

StrNone = Optional[str]
Attrs = dict[str, str]

class GtfParseError(Exception):
    """Error parsing GTF record."""
    pass

class GtfIdError(KeyError):
    """GTF id lookup error"""
    pass

class GtfRecord:
    """Base GTF record."""
    def __init__(self, chrom: str, source: str, feature: str,
                 start: int, end: int, score: str, strand: str, frame: str,
                 attrs: Attrs):
        self.chrom = chrom
        self.source = source
        self.feature = feature
        self.start = start  # 0-based
        self.end = end
        self.score = score
        self.strand = strand
        self.frame = frame
        self.gene_id = attrs.get("gene_id")
        self.gene_name = attrs.get("gene_name")
        self.transcript_id = attrs.get("transcript_id")
        self.attrs = attrs

    def __len__(self) -> int:
        return self.end - self.start

    @property
    def coords(self) -> SeqRange:
        """coordinates of sequence"""
        return SeqRange(self.chrom, self.start, self.end, self.strand)

    @property
    def coords_no_strand(self) -> SeqRange:
        """coordinates of sequence, without strand"""
        return SeqRange(self.chrom, self.start, self.end)

def gtf_record_sort_key(rec):
    return (rec.chrom, rec.start, rec.end)

class GtfExon(GtfRecord):
    """GTF exon."""

    def __init__(self, chrom: str, source: str, feature: str, start: int, end: int,
                 score: str, strand: str, frame: str,
                 attrs: Attrs):
        super().__init__(chrom, source, feature, start, end, score, strand, frame, attrs)

class GtfCDS(GtfRecord):
    """GTF CDS (coding sequence)."""

    def __init__(self, chrom: str, source: str, feature: str, start: int, end: int,
                 score: str, strand: str, frame: str,
                 attrs: Attrs):
        super().__init__(chrom, source, feature, start, end, score, strand, frame, attrs)

class GtfTranscript(GtfRecord):
    """GTF transcript with exons."""
    def __init__(self, chrom: str, source: str, feature: str, start: int, end: int,
                 score: str, strand: str, frame: str, attrs: Attrs):
        super().__init__(chrom, source, feature, start, end, score, strand, frame, attrs)

        self.exons: list[GtfExon] = []
        self.cds_recs: list[GtfCDS] = []

    def add_exon(self, exon: GtfExon) -> None:
        """Add exon to transcript."""
        self.exons.append(exon)

    def add_cds(self, cds: GtfCDS) -> None:
        """Add CDS to transcript."""
        self.cds_recs.append(cds)

    def sort_children(self) -> None:
        self.exons.sort(key=gtf_record_sort_key)
        self.cds_recs.sort(key=gtf_record_sort_key)

class GtfData:
    """data from a GTF file"""
    def __init__(self):
        self.transcripts = []
        self.transcripts_by_id: dict[str, GtfTranscript] = {}

    def add_transcript(self, gtf_transcript):
        self.transcripts.append(gtf_transcript)
        self.transcripts_by_id[gtf_transcript.transcript_id] = gtf_transcript

    def get_transcript(self, transcript_id):
        """return transcript for id or None if not found"""
        return self.transcripts_by_id.get(transcript_id)

    def fetch_transcript(self, transcript_id):
        """return transcript for id or error if not found"""
        try:
            return self.transcripts_by_id[transcript_id]
        except KeyError:
            raise GtfIdError(f"unknown transcript id `{transcript_id}'")

    def iter_transcript_ids(self):
        return self.transcripts_by_id.keys()


def _parse_attributes(attrs_str: str) -> Attrs:
    """Parse GTF attributes string into dict."""
    return dict(re.findall(r'(\w+)\s+"([^"]+)"', attrs_str))

def _parse_coordinates(start_str: str, end_str: str) -> tuple[int, int]:
    try:
        start = int(start_str) - 1  # Convert to 0-based
        end = int(end_str)
    except ValueError as exc:
        raise GtfParseError(f"Invalid coordinates `{start_str}' `{end_str}'") from exc
    if start >= end:
        raise GtfParseError(f"Coordinates `{start_str}' > `{end_str}'")
    return start, end

def _parse_strand(strand: str) -> str:
    if strand not in {'+', '-', '.'}:
        raise GtfParseError(f"Invalid strand '{strand}'")
    return strand

def _parse_score(score_str: str) -> str:
    try:
        if score_str == '.':
            return None
        else:
            return float(score_str)
    except Exception as exc:
        raise GtfParseError(f"Invalid score '{score_str}'") from exc

def _parse_frame(frame_str: str) -> str:
    try:
        if frame_str == '.':
            return None
        else:
            frame = int(frame_str)
    except Exception as exc:
        raise GtfParseError(f"Invalid frame '{frame_str}'") from exc

    if not (0 <= frame <= 2):
        raise GtfParseError(f"Frame must be in the range 0..2 or `.', got'{frame}'")
    return frame


TRANSCRIPT_FEATURES = frozenset([
    'transcript', 'mRNA', 'lncRNA', 'miRNA',
    'ncRNA', 'rRNA', 'tRNA', 'snRNA', 'snoRNA',
    'processed_transcript', 'pseudogenic_transcript'
])

def _gtf_record_class(feature):
    if feature in TRANSCRIPT_FEATURES:
        return GtfTranscript
    elif feature == "exon":
        return GtfExon
    elif feature == "CDS":
        return GtfCDS
    else:
        return GtfRecord

def _parse_gtf_line(line: str) -> GtfRecord:
    """Parse a single GTF line into a GtfRecord or derived class."""
    fields = line.split('\t')
    if len(fields) != 9:
        raise GtfParseError(f"Expected 9 fields, got {len(fields)}")

    start, end = _parse_coordinates(fields[3], fields[4])
    attrs = _parse_attributes(fields[8])

    cls = _gtf_record_class(fields[2])
    return cls(chrom=fields[0],
               source=fields[1],
               feature=fields[2],
               start=start,
               end=end,
               score=_parse_score(fields[5]),
               strand=_parse_strand(fields[6]),
               frame=_parse_frame(fields[7]),
               attrs=attrs)

def gtf_record_parser(gtf_file: str):
    """Parse GTF file, yields GtfRecord GtfTranscript, GtfExon, or GtfCDS objects.
    File maybe compressed"""
    with fileOps.opengz(gtf_file) as fh:
        for line_num, line in enumerate(fh, start=1):
            try:
                line = line.rstrip()
                if not ((len(line) == 0) or line.startswith('#')):
                    rec = _parse_gtf_line(line)
                    if rec is not None:
                        yield rec
            except GtfParseError as exc:
                raise GtfParseError(f"{gtf_file}:{line_num}: invalid GTF record") from exc

def _load_gtf_record(gtf_file, gtf_data, transcript_id_to_exons, transcript_id_to_cds_recs):
    for rec in gtf_record_parser(gtf_file):
        if isinstance(rec, GtfTranscript):
            gtf_data.add_transcript(rec)
        elif isinstance(rec, GtfExon):
            transcript_id_to_exons[rec.transcript_id].append(rec)
        elif isinstance(rec, GtfCDS):
            transcript_id_to_cds_recs[rec.transcript_id].append(rec)

def _add_children(transcript, exons, cds_recs):
    if exons is not None:
        for exon in exons:
            transcript.add_exon(exon)
    if cds_recs is not None:
        for cds_rec in cds_recs:
            transcript.add_cds(cds_rec)
    transcript.sort_children()

def _resolve_gtf_records(gtf_data, transcript_id_to_exons, transcript_id_to_cds_recs):
    """add exon and CDS records to transcripts and sort"""
    for transcript in gtf_data.transcripts:
        _add_children(transcript,
                      transcript_id_to_exons.get(transcript.transcript_id),
                      transcript_id_to_cds_recs.get(transcript.transcript_id))

def gtf_data_parser(gtf_file):
    """parse a GTF file into a GtfData object"""
    # must save up exons and CDS, as sorting of GTF files is not required
    gtf_data = GtfData()
    transcript_id_to_exons = defaultdict(list)
    transcript_id_to_cds_recs = defaultdict(list)
    _load_gtf_record(gtf_file, gtf_data, transcript_id_to_exons, transcript_id_to_cds_recs)
    _resolve_gtf_records(gtf_data, transcript_id_to_exons, transcript_id_to_cds_recs)
    return gtf_data
