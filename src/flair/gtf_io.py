"""
Simplistic, non-validating GTF parser.
"""
from typing import Optional
from flair.pycbio.sys import fileOps

StrNone = Optional[str]

class GtfParseError(Exception):
    """Error parsing GTF record."""
    pass

class GtfRecord:
    """Base GTF record."""
    def __init__(self, chrom: str, source: str, feature: str,
                 start: int, end: int, score: str, strand: str, frame: str,
                 gene_id: StrNone = None,
                 gene_name: StrNone = None,
                 transcript_id: StrNone = None):
        self.chrom = chrom
        self.source = source
        self.feature = feature
        self.start = start  # 0-based
        self.end = end
        self.score = score
        self.strand = strand
        self.frame = frame
        self.gene_id = gene_id
        self.gene_name = gene_name
        self.transcript_id = transcript_id

    def __len__(self) -> int:
        return self.end - self.start


class GtfExon(GtfRecord):
    """GTF exon."""

    def __init__(self, chrom: str, source: str, feature: str, start: int, end: int,
                 score: str, strand: str, frame: str,
                 gene_id: StrNone = None,
                 gene_name: StrNone = None,
                 transcript_id: StrNone = None):
        super().__init__(chrom, source, feature, start, end,
                         score, strand, frame, gene_id, gene_name, transcript_id)

class GtfCDS(GtfRecord):
    """GTF CDS (coding sequence)."""

    def __init__(self, chrom: str, source: str, feature: str, start: int, end: int,
                 score: str, strand: str, frame: str,
                 gene_id: StrNone = None,
                 gene_name: StrNone = None,
                 transcript_id: StrNone = None):
        super().__init__(chrom, source, feature, start, end,
                         score, strand, frame, gene_id, gene_name, transcript_id)

class GtfTranscript(GtfRecord):
    """GTF transcript with exons."""
    def __init__(self, chrom: str, source: str, feature: str, start: int, end: int,
                 score: str, strand: str, frame: str,
                 gene_id: StrNone = None,
                 gene_name: StrNone = None,
                 transcript_id: StrNone = None):
        super().__init__(chrom, source, feature, start, end,
                         score, strand, frame, gene_id, gene_name, transcript_id)

        self.exons: list[GtfExon] = []
        self.cds_recs: list[GtfCDS] = []

    def add_exon(self, exon: GtfExon) -> None:
        """Add exon to transcript."""
        self.exons.append(exon)

    def add_cds(self, cds: GtfCDS) -> None:
        """Add CDS to transcript."""
        self.cds_recs.append(cds)


class GtfData:
    """data from a GTF file"""
    def __init__(self):
        self.transcripts_by_id: dict[GtfTranscript] = {}

    def add_transcript(self, gxf_transcript):
        self.transcripts_by_id[gxf_transcript.transcripts_id] = gxf_transcript


def _parse_attr_value(attr: str) -> tuple[str, str]:
    """Parse a single GTF attribute."""
    if ' ' not in attr:
        raise GtfParseError(f"Invalid attribute format '{attr}'")
    key, _, value = attr.partition(' ')
    value = value.strip('"')
    return key, value

def _parse_attribute(attr_str: str, attrs: dict) -> dict[str, str]:
    try:
        key, value = _parse_attr_value(attr_str)
        attrs[key] = value
    except Exception as exc:
        raise GtfParseError(f"Failed to parse attribute `{attr_str}'") from exc

def _parse_attributes(attrs_str: str) -> dict[str, str]:
    """Parse GTF attributes string into dict."""
    attrs = {}
    for attr_str in attrs_str.split(';'):
        _parse_attribute(attr_str.strip(), attrs)
    return attrs

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

def _gxf_record_class(feature):
    if feature in TRANSCRIPT_FEATURES:
        return GtfTranscript
    elif feature == "exon":
        return GtfExon
    elif feature == "CDS":
        return GtfCDS
    else:
        return GtfRecord

def _parse_gtf_line(line: str) -> GtfRecord:
    """Parse a single GTF line"""
    fields = line.split('\t')
    if len(fields) != 9:
        raise GtfParseError(f"Expected 9 fields, got {len(fields)}")

    start, end = _parse_coordinates(fields[3], fields[4])
    attrs = _parse_attributes(fields[8])

    cls = _gxf_record_class(fields[2])
    return cls(chrom=fields[0],
               source=fields[1],
               feature=fields[2],
               start=start,
               end=end,
               score=_parse_score(fields[5]),
               strand=_parse_strand(fields[6]),
               frame=_parse_frame[7],
               gene_id=attrs.get("gene_id"),
               gene_name=attrs.get("gene_name"),
               transcript_id=attrs.get("transcript_id"))

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
            except GtfParseError as e:
                raise GtfParseError(f"{gtf_file}:{line_num}: invalid GTF record") from e

def _load_gtf_record(gxf_file, gtf_data, transcript_id_to_exons, transcript_id_to_cdses):
    for rec in gtf_record_parser(gtf_file):
        if isinstance(rec, GtfTranscript):
            gtf_data.add_transcript(rec)
        elif isinstance(rec, GtfExon):
            transcript_id_to_exons[rec.transcript_id].append(rec)
        elif isinstance(rec, GtfCDS):
            transcript_id_to_cdses[rec.transcript_id].append(rec)

def _add_exons(transcript, exons):


def _resolve_gtf_records(gtf_data, transcript_id_to_exons, transcript_id_to_cdses):
    """add exon and CDS records to transcripts and sort"""
    for transcript_id in transcript_id_to_exons.keys():
        _add_exons(
        transcript = gtf_data.fetch_transcripts(transcript_id)



def gtf_data_parser(gtf_file):
    """parse a GTF file into a GtfData object"""

    # must save up exons and CDS, as sorting of GTF files is not required
    gtf_data = GtfData()
    transcript_id_to_exons = defaultlist(list)
    transcript_id_to_cdses = defaultlist(list)
    _load_gtf_record(gxf_file, gtf_data, transcript_id_to_exons, transcript_id_to_cdses)
