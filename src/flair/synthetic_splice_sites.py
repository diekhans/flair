import sys
from collections import Counter
import pysam

alignedbedfile = sys.argv[1]
referencegtffile = sys.argv[2]
outfilename = sys.argv[3]
refbpfile = sys.argv[4]
sjwiggle = int(sys.argv[5]) #15
readcov = int(sys.argv[6]) #2
refgenomefile = sys.argv[7]

def grouper(iterable):
    prev = None
    group = []
    for item in iterable:
        if prev is None or item - prev <= sjwiggle:
            group.append(item)
        else:
            yield group
            group = [item]
        prev = item
    if group:
        yield group

fusiontoannotsj = {}
for line in open(referencegtffile):
    line = line.rstrip().split('\t', 5)
    if line[2] == 'exon':
        if line[0] not in fusiontoannotsj: fusiontoannotsj[line[0]] = set()
        fusiontoannotsj[line[0]].add((int(line[3]), int(line[4])))

fusiontobp = {}
for line in open(refbpfile):
    line = line.rstrip().split('\t')
    fusiontobp[line[0]] = int(line[1])


genome = pysam.FastaFile(refgenomefile)


splicejunctosupport = {}
chrtonovelss = {}
reconsideredss = {}
c = 0
introns_to_reads = {}
for line in open(alignedbedfile):
    line = line.rstrip().split('\t')
    thischr, iso, strand, start, esizes, estarts = line[0], line[3], line[5], int(line[1]), \
                            [int(x) for x in line[10].rstrip(',').split(',')], [ int(x) for x in line[11].rstrip(',').split(',')]
    for i in range(len(esizes) - 1):
        thisintron = tuple([thischr, start + estarts[i] + esizes[i] , start + estarts[i + 1] + 1])
        if thisintron not in introns_to_reads:
            introns_to_reads[thisintron] = 0
        introns_to_reads[thisintron] += 1

for thisintron in introns_to_reads:
    readcount = introns_to_reads[thisintron]
    thisintron = list(thisintron)
    
    foundSJ = False
    if genome.fetch(thischr, thisintron[1], thisintron[1]+2) == 'GT' and \
            genome.fetch(thischr, thisintron[2]-3, thisintron[2]-1) == 'AG':
        foundSJ = True
        strand = '+'
    elif genome.fetch(thischr, thisintron[1], thisintron[1]+2) == 'CT' and \
            genome.fetch(thischr, thisintron[2]-3, thisintron[2]-1) == 'AC':
        foundSJ = True
        strand = '-'
    thisintron.append(strand)

    # print(thisintron, readcount, foundSJ, thisintron[1] <= fusiontobp[thischr] <= thisintron[2])

    
    if thisintron[1] <= fusiontobp[thischr] <= thisintron[2] or foundSJ: ##only process introns that have correct motifs OR cross the fusion breakpoint
        doreconsider = True
        if thischr in fusiontoannotsj:
            # closestdist, closestpos = 1000, None
            close_ref = False
            for sj in fusiontoannotsj[thischr]:
                d1 = abs(thisintron[1] - sj[0])
                d2 = abs(thisintron[2] - sj[1])
                if d1 <= sjwiggle and d2 <= sjwiggle:#== 0:
                    # print(iso, thisintron[i], genome.fetch(thischr, thisintron[i], thisintron[i]+2))
                    # if i == 1: print(iso, thisintron[i], genome.fetch(thischr, thisintron[i], thisintron[i]+2))
                    # else: print(iso, thisintron[i], genome.fetch(thischr, thisintron[i]-3, thisintron[i]-1))
                    # if thisdist < closestdist: closestdist, closestpos = thisdist, sj
                    close_ref = True
                    break
            # if closestpos: doreconsider = False
        # print(close_ref)
        if not close_ref:
            if thischr not in chrtonovelss: chrtonovelss[thischr] = []
            chrtonovelss[thischr].extend([thisintron[1]*readcount])
            chrtonovelss[thischr].extend([thisintron[2]*readcount])
            reconsideredss[tuple(thisintron)] = readcount


chrtogoodss = {}
for chr in chrtonovelss:
    chrtogoodss[chr] = set()
    chrtonovelss[chr] = sorted(chrtonovelss[chr])
    for group in grouper(chrtonovelss[chr]):
        if len(group) >= readcov:
            groupsize = len(group)
            finalpos = []
            for pos, count in Counter(group).most_common():
                if count >= readcov and count >= groupsize/10: ###needs to also be 1/10 of locus
                    isdifferent = True
                    for fp in finalpos:
                        if fp-sjwiggle <= pos <= fp+sjwiggle: isdifferent = False
                    if isdifferent: finalpos.append(pos)
            for fp in finalpos:
                chrtogoodss[chr].add(fp)

# print(chrtogoodss)

for thisintron in reconsideredss:
    readcount = reconsideredss[thisintron]
    thisintron = list(thisintron)
    thischr = thisintron[0]
    for i in range(1, 3):  ##checking each splice site
        ###FIXME currently processing each splice site individually
        closestdist, closestpos = 1000, None
        if thischr in chrtogoodss:
            for sj in chrtogoodss[thischr]:
                thisdist = abs(thisintron[i] - sj)
                if thisdist <= sjwiggle:
                    if thisdist < closestdist: closestdist, closestpos = thisdist, sj
        if closestpos: thisintron[i] = closestpos
    thisintron[2] -= 1
    thisintron = tuple(thisintron)
    if thisintron not in splicejunctosupport: splicejunctosupport[thisintron] = 0
    splicejunctosupport[thisintron] += readcount


goodsj, lowsupsj = {}, {}

out = open(outfilename, 'w')
for j in splicejunctosupport:
    # print(j, splicejunctosupport[j])
    if splicejunctosupport[j] >= 2:
        # goodsj[j] = splicejunctosupport[j]
        out.write('\t'.join([str(x) for x in j[:3]] + ['.', str(splicejunctosupport[j]), j[3]]) + '\n')

out.close()
