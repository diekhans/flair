#!/usr/bin/env python3
"""Call alternative 3' splice site, alternative 5' splice site and intron
retention events from an isoform BED."""
import argparse
from collections import namedtuple
from flair.counts_matrix_tsv import read_sample_columns, read_counts_rows
from flair import PosRange
from flair.diffsplice_event_ids import junction_id, site_id
from flair.isoform_data import Junc
from flair.event_quant_tsv import event_quant_writer
from flair.pycbio.hgdata.bed import BedReader
from flair.pycbio.sys import cli

# minimum distance apart for alt SS to be tested
wiggle = 10


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('isoforms_bed', help="isoform BED to call events from")
    parser.add_argument('out_prefix', help="prefix for the three events.quant.tsv outputs")
    parser.add_argument('counts_tsv', nargs='?', help="[optional] isoform counts matrix")
    parser.add_argument('--compress', action='store_true',
                        help="gzip the event files, as flair --compress requires")
    return parser.parse_args()


class Isoform(namedtuple("Isoform", ("name", "span", "exons"))):
    """One isoform of the BED: the span it covers and its exons, both PosRange."""
    __slots__ = ()

class FlankedJunc(namedtuple("FlankedJunc", ("junc", "prev_exon_start", "next_exon_end"))):
    """A junction with the outer bounds of the exons on either side, which is what
    distinguishes an alternative splice site from a skipped exon."""
    __slots__ = ()

def event_quant_path(args, event):
    "the file one event type is written to, gzipped when asked for"
    path = f'{args.out_prefix}.{event}.events.quant.tsv'
    return path + '.gz' if args.compress else path

def get_junctions_bed(exons):
    "the junctions between an isoform's exons, each with its flanking exon bounds"
    return [FlankedJunc(Junc(exons[b].end, exons[b + 1].start),
                        exons[b].start, exons[b + 1].end)
            for b in range(len(exons) - 1)]

def isoform_retains_junction(isoform, junc):
    "true when one exon of this isoform spans the junction, which is a retention"
    if (isoform.span.start > junc.end) or (isoform.span.end < junc.start):
        return False
    # every exon, not exons[1:]: a junction retained inside the first exon of
    # another isoform was not called, so this and mark_intron_retention differed
    # about the same event
    return any((exon.start < junc.start) and (exon.end > junc.end)
               for exon in isoform.exons)


def update_altsplice_dict(jdict, chrom, strand, fiveprime, threeprime, exon_start, exon_end,
                          sample_names, iso_counts, name, search_threeprime=True):
    if fiveprime not in jdict[chrom]:
        jdict[chrom][fiveprime] = {}  # 5' end anchor if search_threeprime
    if threeprime not in jdict[chrom][fiveprime]:
        jdict[chrom][fiveprime][threeprime] = {}
        jdict[chrom][fiveprime][threeprime]['counts'] = [0] * len(sample_names)
        jdict[chrom][fiveprime][threeprime]['isos'] = []  # isoform list for this junction
        jdict[chrom][fiveprime][threeprime]['exon_end'] = exon_end  # for detecting exon skipping
    elif (search_threeprime and strand == '+') or (not search_threeprime and strand == '-'):
        if exon_end < jdict[chrom][fiveprime][threeprime]['exon_end']:  # select shorter exon end
            jdict[chrom][fiveprime][threeprime]['exon_end'] = exon_end
    else:
        if exon_end > jdict[chrom][fiveprime][threeprime]['exon_end']:
            jdict[chrom][fiveprime][threeprime]['exon_end'] = exon_end
    jdict[chrom][fiveprime][threeprime]['isos'] += [name]
    for c in range(len(sample_names)):
        jdict[chrom][fiveprime][threeprime]['counts'][c] += iso_counts[name][c]
    return jdict


def _is_alternative_site(fiveprime, tp1, tp2, exon_end, search_threeprime, strand):
    """is tp2 an alternative site to tp1, both spliced to fiveprime: a site of its own,
    wiggle or more away, not already tested in the other order, and overlapping tp1's
    exon, since a site past its end is an exon skip"""
    if tp1 == tp2 or abs(tp2 - tp1) < wiggle or abs(fiveprime - tp1) > abs(fiveprime - tp2):
        return False
    if (search_threeprime and strand == '+') or (not search_threeprime and strand == '-'):
        return tp2 <= exon_end
    else:
        return tp2 >= exon_end


def _write_altss_events(juncs_of_site, fiveprime, chrom, writer, search_threeprime):
    "the alternative site events of one anchoring splice site"
    all_tp = list(juncs_of_site.keys())
    strand, chrom_clean = chrom[0], chrom[1:]
    n = 0
    for tp1 in all_tp:  # tp1 = three prime SS number 1
        exon_end = juncs_of_site[tp1]['exon_end']
        for tp2 in all_tp:  # tp2 is also a 3' SS with the same 5' anchor as tp1
            if _is_alternative_site(fiveprime, tp1, tp2, exon_end, search_threeprime, strand):
                inclusion, exclusion = tp1, tp2
                anchor = site_id(chrom_clean, strand, fiveprime)
                feature_suffix = anchor if n == 0 else anchor + '-' + str(n)
                event = (junction_id(chrom_clean, strand, fiveprime, inclusion) + '_' +
                         junction_id(chrom_clean, strand, fiveprime, exclusion))

                writer.writeSide('inclusion', feature_suffix, event,
                                 juncs_of_site[inclusion]['counts'],
                                 sorted(juncs_of_site[inclusion]['isos']))
                writer.writeSide('exclusion', feature_suffix, event,
                                 juncs_of_site[exclusion]['counts'],
                                 sorted(juncs_of_site[exclusion]['isos']))
                n += 1


def find_altss(alljuncs, writer, search_threeprime=True):
    """ If fiveprimeon is True, then alternative 5' SS will be reported instead. """
    for chrom in alljuncs:
        for fiveprime in alljuncs[chrom]:
            # one 3' end is not an alt 3' junction
            if len(alljuncs[chrom][fiveprime]) > 1:
                _write_altss_events(alljuncs[chrom][fiveprime], fiveprime, chrom, writer, search_threeprime)


def _ir_event_kept(ir_junction, sample_names):
    """An intron retention event needs an isoform retaining the intron, and, when
    counts were given, reads on that side."""
    if len(ir_junction['inclusion']['isos']) == 0:
        return False
    else:
        return (not sample_names) or (sum(ir_junction['inclusion']['counts']) > 0)


def _write_ir_event(writer, chrom, junc, ir_junction, sample_names):
    strand, chrom_clean = chrom[0], chrom[1:]
    event = junction_id(chrom_clean, strand, junc.start, junc.end)
    for side in ('inclusion', 'exclusion'):
        counts = ir_junction[side]['counts'] if sample_names else []
        writer.writeSide(side, event, event, counts, sorted(ir_junction[side]['isos']))


def main():
    args = parse_args()
    with cli.ErrorHandler():
        call_events(args)


def _read_iso_counts(counts_tsv):
    "isoform -> its counts by sample, and the sample names; both empty with no counts matrix"
    if not counts_tsv:
        return {}, []
    sample_names = read_sample_columns(counts_tsv)
    return ({row.isoform_id: [float(x) for x in row.counts] for row in read_counts_rows(counts_tsv)},
            sample_names)


def _new_ir_junction(nsamples):
    "the counts and isoforms of a junction's two sides, before any isoform is added"
    return {side: {'counts': [0] * nsamples, 'isos': []} for side in ('exclusion', 'inclusion')}


def _add_isoform_junctions(bed, chrom, exons, tables, sample_names, iso_counts):
    "one isoform's junctions, in the alternative splice site tables and the IR table"
    a3_junctions, a5_junctions, ir_chrom = tables
    for flanked in get_junctions_bed(exons):
        j = flanked.junc
        fiveprime, threeprime = j.start, j.end
        exon_end, exon_start = flanked.next_exon_end, flanked.prev_exon_start

        if bed.strand == '-':
            fiveprime, threeprime = threeprime, fiveprime
            exon_end, exon_start = exon_start, exon_end

        update_altsplice_dict(a3_junctions, chrom, bed.strand, fiveprime, threeprime,
                              exon_start, exon_end, sample_names, iso_counts, bed.name)
        update_altsplice_dict(a5_junctions, chrom, bed.strand, threeprime, fiveprime,
                              exon_end, exon_start, sample_names, iso_counts, bed.name,
                              search_threeprime=False)

        # IR detection needs the junction alone, without the flanking exons
        exclusion = ir_chrom.setdefault(j, _new_ir_junction(len(sample_names)))['exclusion']
        exclusion['isos'] += [bed.name]
        for c in range(len(sample_names)):
            exclusion['counts'][c] += iso_counts[bed.name][c]


def _collect_isoform_junctions(isoforms_bed, sample_names, iso_counts):
    """the isoforms and their junctions by stranded chromosome: the isoforms, the
    junctions for intron retention, and those for alternative 3' and 5' splice sites"""
    isoforms, ir_junctions, a3_junctions, a5_junctions = {}, {}, {}, {}
    for bed in BedReader(isoforms_bed, fixScores=True):
        # an isoform of the BED that the counts matrix has no row for is left out
        if (not iso_counts) or (bed.name in iso_counts):
            chrom = bed.strand + bed.chrom  # stranded comparisons
            if chrom not in isoforms:
                isoforms[chrom], ir_junctions[chrom] = {}, {}
                a3_junctions[chrom], a5_junctions[chrom] = {}, {}
            exons = tuple(PosRange(blk.start, blk.end) for blk in bed.blocks)
            isoforms[chrom][bed.name] = Isoform(bed.name, PosRange(bed.chromStart, bed.chromEnd), exons)
            _add_isoform_junctions(bed, chrom, exons, (a3_junctions, a5_junctions, ir_junctions[chrom]),
                                   sample_names, iso_counts)
    return isoforms, ir_junctions, a3_junctions, a5_junctions


def _add_ir_inclusions(chrom_junctions, chrom_isoforms, sample_names, iso_counts):
    "the isoforms retaining each junction of one chromosome, with their counts"
    for j, junction in chrom_junctions.items():
        for iname, isoform in chrom_isoforms.items():  # every other isoform, to find IR
            if (iname not in junction['exclusion']['isos']) and isoform_retains_junction(isoform, j):
                junction['inclusion']['isos'] += [iname]
                for c in range(len(sample_names)):
                    junction['inclusion']['counts'][c] += iso_counts[iname][c]


def _write_ir_events(writer, chrom, chrom_junctions, sample_names):
    for j, junction in chrom_junctions.items():
        if _ir_event_kept(junction, sample_names):
            _write_ir_event(writer, chrom, j, junction, sample_names)


def call_events(args):
    "the alternative 3' and 5' splice site and intron retention events of the isoforms"
    iso_counts, sample_names = _read_iso_counts(args.counts_tsv)
    isoforms, ir_junctions, a3_junctions, a5_junctions = _collect_isoform_junctions(args.isoforms_bed,
                                                                                    sample_names, iso_counts)
    with event_quant_writer(event_quant_path(args, 'alt3'), sample_names) as writer:
        find_altss(a3_junctions, writer)

    with event_quant_writer(event_quant_path(args, 'alt5'), sample_names) as writer:
        find_altss(a5_junctions, writer, search_threeprime=False)

    with event_quant_writer(event_quant_path(args, 'ir'), sample_names) as writer:
        for chrom in ir_junctions:
            _add_ir_inclusions(ir_junctions[chrom], isoforms[chrom], sample_names, iso_counts)
            _write_ir_events(writer, chrom, ir_junctions[chrom], sample_names)
            # a chromosome's junctions are not read again, and they are large
            ir_junctions[chrom] = None


if __name__ == '__main__':
    main()
