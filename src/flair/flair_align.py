#! /usr/bin/env python3

import sys
import argparse
import os
import pipettor
import pysam
import logging
from flair.pycbio.sys import cli
from flair import FlairInputDataError

def parse_args():
    desc = "FLAIR align outputs an unfiltered bam file and a filtered bed file for use in the downstream pipeline"
    parser = argparse.ArgumentParser(description=desc)

    reads = parser.add_argument_group('required named arguments')
    reads.add_argument('-r', '--reads', nargs='+', type=str, required=True,
                          help='FASTA/FASTQ file(s) of raw reads, either space or comma separated')
    genome = parser.add_argument_group('Either one of the following arguments is required')
    genome.add_argument('-g', '--genome', type=str,
                        help='FASTA of reference genome, can be minimap2 indexed')
    genome.add_argument('--mm_index', type=str, default='',
                        help='minimap2 index .mmi file')
    parser.add_argument('-o', '--output', default='flair.aligned',
                        help='output file name base (default: flair.aligned)')
    parser.add_argument('-t', '--threads', type=int, default=4,
                        help='minimap2 number of threads (4)')
    parser.add_argument('--junction_bed', default='',
                        help='annotated isoforms/junctions bed file for splice site-guided minimap2 genomic alignment')
    parser.add_argument('--nvrna', action='store_true', default=False,
                        help='specify this flag to use native-RNA specific alignment parameters for minimap2')
    parser.add_argument('--minfragmentsize', type=int, default=80,
                        help='minimum size of alignment kept, used in minimap -s. More important when doing downstream fusion detection')
    parser.add_argument('--maxintronlen', default='200k',
                        help='maximum intron length in genomic alignment. Longer can help recover more novel isoforms with long introns')
    parser.add_argument('--quiet', default=False, action='store_true', dest='quiet',
                        help='''Suppress minimap progress statements from being printed''')
    args = cli.parseArgsWithLogging(parser)

    reads = []
    for rfiles in args.reads:
        for rfile in rfiles.split(','):
            if not os.path.exists(rfile):
                raise FlairInputDataError(f'Error: read file does not exist: {rfile}')
            reads.append(rfile)
    args.reads = reads
    return args

def intron_chain_to_exon_starts(ichain, start, end):
    esizes, estarts = [], [0,]
    for i in ichain:
        esizes.append(i[0] - (start + estarts[-1]))
        estarts.append(i[1] - start)
    esizes.append(end - (start + estarts[-1]))
    return esizes, estarts


def inferMM2JuncStrand(read):
    # minimap gives junction strand denoted as 'ts'
    # the sign corresponds to the alignment orientation, where + agrees and - disagrees
    orientation = read.flag
    try:
        juncDir = read.get_tag('ts')
    except:
        juncDir = None

    # Try to resolve strand by looking for polyA
    if not juncDir:
        left, right = read.cigartuples[0], read.cigartuples[-1]
        s1, s2 = read.query_sequence[:100], read.query_sequence[-100:]
        if ("T" * 10 in s1 and left[0] == 4 and left[1] >= 10) and \
            ("A" * 10 in s2 and right[0] == 4 and right[1] >= 10):
            # probably internal priming
            juncDir = "ambig"
        elif ("T" * 10 in s1 and left[0] == 4 and left[1] >= 10):
            # maps to positive strand but has a rev comp polyA
            juncDir = '-'
        elif ("A" * 10 in s2 and right[0] == 4 and right[1] >= 10):
            # maps to positive strand but has a sense polyA
            juncDir = "+"
        else:
            juncDir = "ambig"

    else: # only executes for processing ts tag strand
        if orientation == 0 and juncDir == "+":
            juncDir = "+"
        elif orientation == 0 and juncDir == "-":
            juncDir = "-"
        elif orientation == 16 and juncDir == "+":
            juncDir = "-"
        elif orientation == 16 and juncDir == "-":
            juncDir = "+"
    # print(read.query_name, orientation, ogjuncdir, juncDir)
    return juncDir

def bed_from_cigar(alignstart, is_reverse, cigartuples, readname, referencename, qualscore, juncDirection):
    positiveTxn = "27,158,119"  # green
    negativeTxn = "217,95,2"  # orange
    unknownTxn = "99,99,99"
    refpos = alignstart
    intronblocks = []
    hasmatch = False
    for block in cigartuples:
        if block[0] == pysam.CIGAR_OPS.CREF_SKIP and hasmatch:
            # intron
            intronblocks.append([refpos, refpos + block[1]])
            refpos += block[1]
        elif block[0] in frozenset((pysam.CIGAR_OPS.CMATCH,
                                    pysam.CIGAR_OPS.CEQUAL,
                                    pysam.CIGAR_OPS.CDIFF,
                                    pysam.CIGAR_OPS.CDEL)):
            # consumes reference
            refpos += block[1]
            if block[0] in frozenset((pysam.CIGAR_OPS.CMATCH,
                                      pysam.CIGAR_OPS.CEQUAL,
                                      pysam.CIGAR_OPS.CDIFF)):
                hasmatch = True
    esizes, estarts = intron_chain_to_exon_starts(intronblocks,alignstart, refpos)
    rgbcolor = unknownTxn
    if juncDirection == "+":
        rgbcolor = positiveTxn
    elif juncDirection == "-":
        rgbcolor = negativeTxn
    else:
        juncDirection = "-" if is_reverse else "+"
    outline = [referencename, str(alignstart), str(refpos), readname, str(qualscore), juncDirection, str(alignstart), str(refpos), rgbcolor, str(len(intronblocks) + 1), ','.join([str(x) for x in esizes]) + ',', ','.join([str(x) for x in estarts]) + ',']
    return outline

def doalignment(args):
    # minimap
    mm2_cmd = ['minimap2', '-ax', 'splice', '-s', str(args.minfragmentsize),
               '-G', args.maxintronlen, '--MD', '-t', str(args.threads)]
    if args.nvrna:
        mm2_cmd += ['-uf', '-k14']
    if args.junction_bed:
        mm2_cmd += ['--junc-bed', args.junction_bed]
    mm2_cmd += ['-secondary=no']
    if args.mm_index:
        mm2_cmd += [args.mm_index]
    else:
        mm2_cmd += [args.genome]
    mm2_cmd += args.reads

    samtools_sort_cmd = ('samtools', 'sort', '-@', str(args.threads), '-o', args.output + '.bam', '-')
    samtools_index_cmd = ('samtools', 'index', args.output + '.bam')
    pipettor.run([mm2_cmd, samtools_sort_cmd])
    pipettor.run([samtools_index_cmd])

def align():
    args = parse_args()
    doalignment(args)

def main():
    align()

if __name__ == "__main__":
    main()
