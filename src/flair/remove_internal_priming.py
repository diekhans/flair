#!/usr/bin/env python3

import sys
import argparse
from bisect import bisect_left
import pysam


def checkIsNearAnnotEnd(read3endpos, annotends):
    """
    Implements binary search for nearest transcript end, return true if pos is <=200bp from nearest end
    """
    pos1 = bisect_left(annotends, read3endpos)
    if pos1 == len(annotends): return abs(annotends[pos1 - 1] - read3endpos) <= 200
    disttoend = min(abs(annotends[pos1 - 1] - read3endpos), abs(annotends[pos1] - read3endpos))
    return disttoend <= 200


###add annotation-reliant check for transcript end, implement binary search
def checkInternalPriming(read3endpos, thischr, genome, reqfreq, threshold):
    """
    Checks the genomic sequence adjacent to the read end position for a stretch of As with
    a frequency >= reqfreq and a length >= threshold
    """
    genomeseqnearend = genome.fetch(thischr, max(read3endpos - 30, 0), min(read3endpos + 30, genome.get_reference_length(thischr))).upper()
    # FIXME: maxfreq never used
    maxlen, maxfreq = 0, 0
    if len(genomeseqnearend) > threshold*2:
        halfseqlen = int(len(genomeseqnearend)/2)
        for i in list(range(-1 * halfseqlen, -1*threshold)) + list(range(threshold, halfseqlen)):
            thisseq = genomeseqnearend[min(i + halfseqlen, halfseqlen): max(i + halfseqlen, halfseqlen)]
            thiscount = max(thisseq.count('A'), thisseq.count('T'))
            thisfreq = thiscount / len(thisseq) if len(thisseq) > 0 else 0
            if thisfreq >= reqfreq and len(thisseq) > maxlen:
                maxlen, maxfreq = len(thisseq), thisfreq
    return maxlen >= threshold

def removeinternalpriming(refname, refstart, refend, isrev, genome, annottranscriptends, annotexons, threshold, fracAs):
    """
    Given info from a an aligned bam read, check whether it has internal priming
    refname, refstart, refend, isrev - all info about read alignment
    genome: pysam.FastaFile object
    annottranscriptends: [genomic alignment only] dictionary of chrom to sorted list of transcript end pos
    annotexons: [transcriptomic alignment only] dictionary of transcript name to list of exon lengths
    threshold: max length of stretch of As before something is internal priming
    fracAs: minimum frequency of As in sequence to qualify as polyA (6/8 bp=A -> threshold=0.75)
    """
    read3endpos = refend if not isrev else refstart
    # if aligned to transcriptome, check distance to transcript end
    # if read end is close enough to transcript end, return True (no internal priming)
    if not annottranscriptends:
        if annotexons and refname in annotexons:
            theseexons = annotexons[refname]
            # multi exon transcript
            if len(theseexons) > 1 and read3endpos > sum(theseexons) - theseexons[-1]:
                return True
            # single exon transcript
            elif len(theseexons) == 1 and read3endpos >= theseexons[0]-200:
                return True
    # if read doesn't have stretch of As beyond threshold, doesn't have internal priming,
    # The refname check is from #629 when this is called on a transcriptome alignment
    # FIXME: this function should not be called on a transcriptome alignment
    if ((refname not in genome.references) or
        (not checkInternalPriming(read3endpos, refname, genome, fracAs, threshold))):
        return True
    elif annottranscriptends and refname in annottranscriptends:
        isnearannotend = checkIsNearAnnotEnd(read3endpos, annottranscriptends[refname])
        if isnearannotend:
            return True
    return False
