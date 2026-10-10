#!/usr/bin/env python3
"""Call alternative 3' splice site, alternative 5' splice site and intron
retention events from an isoform BED."""
import argparse
from collections import namedtuple
from flair.counts_matrix_tsv import read_sample_columns, read_counts_rows
from flair.isoform_data import Junc
from flair.event_quant_tsv import event_quant_writer
from flair.pycbio.hgdata.bed import BedReader

# minimum distance apart for alt SS to be tested
wiggle = 10


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('isoforms_bed', help="isoform BED to call events from")
    parser.add_argument('out_prefix', help="prefix for the three events.quant.tsv outputs")
    parser.add_argument('counts_tsv', nargs='?', help="[optional] isoform counts matrix")
    return parser.parse_args()


class FlankedJunc(namedtuple("FlankedJunc", ("junc", "prev_exon_start", "next_exon_end"))):
    """A junction with the outer bounds of the exons on either side, which is what
    tells an alternative splice site from a skipped exon."""
    __slots__ = ()

def get_junctions_bed(starts, sizes):
    "the junctions of one isoform, each with its flanking exon bounds"
    return [FlankedJunc(Junc(starts[b] + sizes[b], starts[b + 1]),
                        starts[b], starts[b + 1] + sizes[b + 1])
            for b in range(len(starts) - 1)]


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
        if exon_end < jdict[chrom][fiveprime][threeprime]['exon_end']:  # pick shorter exon end
            jdict[chrom][fiveprime][threeprime]['exon_end'] = exon_end
    else:
        if exon_end > jdict[chrom][fiveprime][threeprime]['exon_end']:
            jdict[chrom][fiveprime][threeprime]['exon_end'] = exon_end
    jdict[chrom][fiveprime][threeprime]['isos'] += [name]
    for c in range(len(sample_names)):
        jdict[chrom][fiveprime][threeprime]['counts'][c] += iso_counts[name][c]
    return jdict


def find_altss(alljuncs, writer, search_threeprime=True):
    """ If fiveprimeon is True, then alternative 5' SS will be reported instead. """
    for chrom in alljuncs:
        for fiveprime in alljuncs[chrom]:
            if len(alljuncs[chrom][fiveprime]) == 1:  # if there is only one 3' end, not an alt 3' junction
                continue

            all_tp = list(alljuncs[chrom][fiveprime].keys())

            n = 0
            for tp1 in all_tp:  # tp1 = three prime SS number 1
                exon_end = alljuncs[chrom][fiveprime][tp1]['exon_end']
                for tp2 in all_tp:  # tp2 is also a 3' SS with the same 5' anchor as tp1
                    if tp1 == tp2 or abs(tp2 - tp1) < wiggle or abs(fiveprime - tp1) > abs(fiveprime - tp2):
                        # two sites are the same, too close together, or have already been tested in another order
                        continue
                    inclusion = tp1
                    exclusion = tp2
                    strand, chrom_clean = chrom[0], chrom[1:]

                    if (search_threeprime and strand == '+') or (not search_threeprime and strand == '-'):
                        if tp2 > exon_end:  # exon skipping. tp2 does not overlap tp1's exon
                            continue
                    elif tp2 < exon_end:  # exon skipping for alt SS upstream of anchor
                        continue

                    feature_suffix = chrom_clean + ':' + str(fiveprime) if n == 0 else chrom_clean + ':' + str(fiveprime) + '-' + str(n)
                    event = chrom_clean + ':' + str(fiveprime) + '-' + str(inclusion) + '_' + chrom_clean + ':' + str(fiveprime) + '-' + str(exclusion)

                    writer.writeSide('inclusion', feature_suffix, event,
                                     alljuncs[chrom][fiveprime][inclusion]['counts'],
                                     sorted(alljuncs[chrom][fiveprime][inclusion]['isos']))
                    writer.writeSide('exclusion', feature_suffix, event,
                                     alljuncs[chrom][fiveprime][exclusion]['counts'],
                                     sorted(alljuncs[chrom][fiveprime][exclusion]['isos']))
                    n += 1


def main():  # noqa: C901 - FIXME: reduce complexity
    args = parse_args()
    bedfh = open(args.isoforms_bed)
    outfilenamebase = args.out_prefix

    iso_counts = {}
    sample_names = []
    if args.counts_tsv:
        sample_names = read_sample_columns(args.counts_tsv)
        for row in read_counts_rows(args.counts_tsv):
            iso_counts[row.isoform_id] = [float(x) for x in row.counts]

    isoforms = {}  # ir detection
    ir_junctions = {}  # ir detection
    a3_junctions = {}  # alt 3' ss detection
    a5_junctions = {}  # alt 5' ss detection
    for bed in BedReader(bedfh, fixScores=True):
        chrom, name, start, end, strand = bed.chrom, bed.name, bed.chromStart, bed.chromEnd, bed.strand

        if iso_counts and name not in iso_counts:
            continue

        blockstarts = [blk.start for blk in bed.blocks]
        blocksizes = [len(blk) for blk in bed.blocks]

        chrom = strand + chrom  # stranded comparisons
        if chrom not in isoforms:
            isoforms[chrom] = {}
            ir_junctions[chrom] = {}
            a3_junctions[chrom] = {}
            a5_junctions[chrom] = {}

        isoforms[chrom][name] = {}
        isoforms[chrom][name]['sizes'] = blocksizes
        isoforms[chrom][name]['starts'] = blockstarts
        isoforms[chrom][name]['range'] = start, end

        for flanked in get_junctions_bed(blockstarts, blocksizes):
            j = flanked.junc
            fiveprime, threeprime = j.start, j.end
            exon_end, exon_start = flanked.next_exon_end, flanked.prev_exon_start

            if strand == '-':
                fiveprime, threeprime = threeprime, fiveprime
                exon_end, exon_start = exon_start, exon_end

            a3_junctions = update_altsplice_dict(a3_junctions, chrom, strand, fiveprime, threeprime,
                                                 exon_start, exon_end, sample_names, iso_counts, name)
            a5_junctions = update_altsplice_dict(a5_junctions, chrom, strand, threeprime, fiveprime,
                                                 exon_end, exon_start, sample_names, iso_counts, name,
                                                 search_threeprime=False)

            # IR detection needs the junction alone, without the flanking exons
            if j not in ir_junctions[chrom]:  # ir detection
                ir_junctions[chrom][j] = {}
                ir_junctions[chrom][j]['exclusion'] = {}
                ir_junctions[chrom][j]['inclusion'] = {}
                ir_junctions[chrom][j]['exclusion']['counts'] = [0] * len(sample_names)
                ir_junctions[chrom][j]['inclusion']['counts'] = [0] * len(sample_names)
                ir_junctions[chrom][j]['exclusion']['isos'] = []
                ir_junctions[chrom][j]['inclusion']['isos'] = []
            ir_junctions[chrom][j]['exclusion']['isos'] += [name]
            for c in range(len(sample_names)):
                ir_junctions[chrom][j]['exclusion']['counts'][c] += iso_counts[name][c]

    with event_quant_writer(outfilenamebase + '.alt3.events.quant.tsv', sample_names) as writer:
        find_altss(a3_junctions, writer)

    with event_quant_writer(outfilenamebase + '.alt5.events.quant.tsv', sample_names) as writer:
        find_altss(a5_junctions, writer, search_threeprime=False)

    with event_quant_writer(outfilenamebase + '.ir.events.quant.tsv', sample_names) as writer:
        for chrom in ir_junctions:  # noqa: C901 - FIXME: reduce complexity
            for j in ir_junctions[chrom]:
                for iname in isoforms[chrom]:  # compare with all other isoforms to find IR
                    if iname in ir_junctions[chrom][j]['exclusion']['isos']:  # is an exclusion isoform
                        continue
                    start, end = isoforms[chrom][iname]['range']
                    if start > j.end or end < j.start:  # isoform boundaries do not overlap junction
                        continue
                    starts, sizes = isoforms[chrom][iname]['starts'], isoforms[chrom][iname]['sizes']
                    # every block, not starts[1:]: a junction retained inside the first
                    # exon of another isoform was not called, so this and
                    # mark_intron_retention disagreed about the same event
                    for start, size in zip(starts, sizes):
                        estart, eend = start, start + size  # exon start, exon end
                        if estart < j.start and eend > j.end:  # retention
                            ir_junctions[chrom][j]['inclusion']['isos'] += [iname]
                            for c in range(len(sample_names)):
                                ir_junctions[chrom][j]['inclusion']['counts'][c] += iso_counts[iname][c]

            for j in ir_junctions[chrom]:
                incounts = ir_junctions[chrom][j]['inclusion']['counts']
                if sum(incounts) == 0:
                    continue
                if not sample_names:
                    ir_junctions[chrom][j]['exclusion']['counts'] = ir_junctions[chrom][j]['inclusion']['counts'] = []

                chrom_clean = chrom[1:]
                event = chrom_clean + ':' + str(j.start) + '-' + str(j.end)
                writer.writeSide('inclusion', event, event,
                                 ir_junctions[chrom][j]['inclusion']['counts'],
                                 sorted(ir_junctions[chrom][j]['inclusion']['isos']))
                writer.writeSide('exclusion', event, event,
                                 ir_junctions[chrom][j]['exclusion']['counts'],
                                 sorted(ir_junctions[chrom][j]['exclusion']['isos']))
            ir_junctions[chrom] = None


if __name__ == '__main__':
    main()
