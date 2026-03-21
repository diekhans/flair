#! /usr/bin/env python3

import sys
import argparse
import os, glob
import pipettor
import pysam
import logging
import gzip
from flair.gtf_to_bed import gtf_to_bed
from flair.convert_synthetic_to_genome_bed import convert_synthetic_isos, get_paralog_ref
from flair.identify_prelim_fusions import id_chimeras
from collections import defaultdict
from flair import FlairInputDataError
from flair.read_processing import get_sequence_from_bed



def def_value():
    return set()

def report_nofusions(outputprefix):
    logging.info('no fusions detected. Exiting')
    for file in [outputprefix + '.fusions.isoforms.bed', outputprefix + '.fusions.isoforms.fa']:
        f = open(file, 'w')
        f.close()

def align_to_synth_genome(genome, reads, output, additional_options):
    mm2_cmd = ['minimap2', '-ax', 'splice'] + additional_options + [genome, reads]
    samtools_filter_cmd = ('samtools', 'view', '-F', '2048', '-hb', '-')
    samtools_sort_cmd = ('samtools', 'sort', '-o', output, '-')
    samtools_index_cmd = ('samtools', 'index', output)
    pipettor.run([mm2_cmd, samtools_filter_cmd, samtools_sort_cmd])
    pipettor.run([samtools_index_cmd])


def detectfusions():
    parser = argparse.ArgumentParser()
    required = parser.add_argument_group('required named arguments')
    required.add_argument('-g', '--genome',
                          type=str, required=True, help='FastA of reference genome')
    parser.add_argument('-f', '--gtf',
                        type=str, required=True, help='GTF annotation file, used for renaming FLAIR isoforms to annotated isoforms and adjusting TSS/TESs')
    parser.add_argument('-r', '--reads', nargs='+',
                          type=str, help='FastA/FastQ files of raw reads, can specify multiple files')
    required.add_argument('-b', '--genomechimbam',
                          type=str, required=True, help='bam file of chimeric reads from genomic alignment from flair align')
    parser.add_argument('--transcriptchimbam',
                        help='Optional: bam file of chimeric reads from transcriptomic alignment. If not provided, this will be made for you')
    parser.add_argument('-o', '--output', default='flair.fusion',
                        help='output file name base for FLAIR isoforms (default: flair.collapse)')
    parser.add_argument('-t', '--threads', type=int, default=4,
                        help='minimap2 number of threads (4)')
    parser.add_argument('--minfragmentsize', type=int, default=40,
                        help='minimum size of alignment kept, used in minimap -s. More important when doing downstream fusion detection')
    parser.add_argument('-s', '--support', type=float, default=3.0,
                        help='''minimum number of supporting reads for a fusion (3)''')
    parser.add_argument('--maxloci', type=int, default=2,
                        help='''max loci detected in fusion. Set higher for detection of 3-gene+ fusions''')
    parser.add_argument('--max_dist_to_TSS', type=int, default=15000,
                        help='''maximum allowed distance of 5' alignment to TSS of annotated transcript. To not check this, set to -1''')
    parser.add_argument('--keep_intermediate', default=False, action='store_true',
                        help='''keep intermediate and temporary files for debugging purposes''')
    path = os.path.dirname(os.path.realpath(__file__)) + '/'

    # FIXME: incorrect way to check for missing arguments
    no_arguments_passed = len(sys.argv) == 1
    if no_arguments_passed:
        parser.print_help()
        parser.error("No arguments passed. Please provide a bam file, reads, genome, and annotation file")

    args = parser.parse_args()
    if args.max_dist_to_TSS == -1:
        args.max_dist_to_TSS = None

    if not (args.reads or args.transcriptchimbam):
        raise FlairInputDataError(f'please provide either a fastq/fasta reads file or a bam file of chimeric alignments to the transcriptome')
    if args.reads:
        if ',' in args.reads[0]:
            args.reads = args.reads[0].split(',')
        for rfile in args.reads:
            if not os.path.exists(rfile):
                raise FlairInputDataError(f'Read file path does not exist: {rfile}')

    if not os.path.exists(args.genome):
        raise FlairInputDataError(f'Genome file path does not exist: {args.genome}')
    if not os.path.exists(args.gtf):
        if not args.gtf:
            raise FlairInputDataError('Please specify annotated gtf with -f ')
        else:
            raise FlairInputDataError('GTF file path does not exist')


    # if args.annotated_fa == 'generate':
    # get transcript sequences
    
    

    ####NEED TO REMEMBER THAT FUSION DETECTION RELIES ON HAVING PROPERLY STRANDED READS - need to add stranding step and/or better documentation on this

    ###Processing the gtf file so many times is really inefficient, how can we resolve this??

    ###align to transcriptome with --secondary=no

    if not args.transcriptchimbam:
        args.annotated_bed = args.output + '.annotated_transcripts.bed'
        gtf_to_bed(args.annotated_bed, args.gtf, include_gene=True)
        args.annotated_fa = args.output + '.annotated_transcripts.fa'
        get_sequence_from_bed(args.genome, args.output + '.annotated_transcripts.bed', args.annotated_fa)
        
        mm2_cmd = ['minimap2', '-a', '-s', str(args.minfragmentsize), '-t', str(args.threads), '--secondary=no',
                   args.annotated_fa] + args.reads
        mm2_cmd = tuple(mm2_cmd)

        # samtools; the dash at the end means STDIN
        samtools_sort_cmd = ('samtools', 'sort', '-@', str(args.threads), '-o', args.output + '_unfilteredtranscriptome.bam', '-')
        samtools_index_cmd = ('samtools', 'index', args.output + '_unfilteredtranscriptome.bam')
        pipettor.run([mm2_cmd, samtools_sort_cmd])
        pipettor.run([samtools_index_cmd])

        ##filter transcriptome alignment to chimeric only and remove the rest

        ##run filtering
        samfile = pysam.AlignmentFile(args.output + '_unfilteredtranscriptome.bam', "rb")
        withsup = pysam.AlignmentFile(args.output + '_transcriptomechimeric.bam', "wb", template=samfile)
        for read in samfile.fetch():
            if read.is_mapped and not read.is_secondary:
                if read.has_tag('SA'):
                    withsup.write(read)
        samfile.close()
        withsup.close()
        pysam.index(args.output + '_transcriptomechimeric.bam')

        pipettor.run([('rm', args.output + '_unfilteredtranscriptome.bam', args.output + '_unfilteredtranscriptome.bam.bai')])
        args.transcriptchimbam = args.output + '_transcriptomechimeric.bam'
        print('aligned to transcriptome')


    print('reading gtf')
    genetoinfo, genetoexons, genetoname = {}, {}, {}
    chrom_to_gene_pos = {}
    gene_to_all_exons, juncs_to_gene = {}, {}
    for line in open(args.gtf):
        if line.startswith('#'):
            continue
        line = line.rstrip().split('\t', 8)
        chrom, ty, start, end, strand = line[0], line[2], int(line[3]) - 1, int(line[4]), line[6]
        if ty in {'gene', 'exon', 'transcript'}:
            gene_id = line[8].split('gene_id "')[1].split('"')[0]
            gene_id = gene_id.replace('_', '-')
            gene_id = gene_id.split('.')[0]
            if ty == 'gene':
                genetoinfo[gene_id] = [chrom, start, end, strand, []]
                genename = line[-1].split('gene_name "')[1].split('"')[0]
                genetoname[gene_id] = genename
                if chrom not in chrom_to_gene_pos:
                    chrom_to_gene_pos[chrom] = []
                    chrom_to_gene_pos[chrom].append((start, end, strand, gene_id))
                    juncs_to_gene[chrom] = {}

            elif ty == 'exon':
                transcript_id = line[8].split('transcript_id "')[1].split('"')[0]
                if gene_id not in genetoexons: 
                    genetoexons[gene_id] = {}
                if transcript_id not in genetoexons[gene_id]: genetoexons[gene_id][transcript_id] = []
                genetoexons[gene_id][transcript_id].append((start, end))

                
            elif ty == 'transcript':
                end5 = start if strand == '+' else end
                genetoinfo[gene_id][-1].append(end5)
    
    print('continuing to parse annot')
    ## FOR JUNCS TO GENE, DO BY CHROM AS WELL
    for chrom in chrom_to_gene_pos:
        chrom_to_gene_pos[chrom] = sorted(chrom_to_gene_pos[chrom])

    for gene in genetoexons:
        gene_to_all_exons[gene] = set()
        chrom = genetoinfo[gene][0]
        for t in genetoexons[gene]:
            exons = sorted(list(genetoexons[gene][t]))
            gene_to_all_exons[gene].update(set(exons))
            if len(exons) > 1:
                juncs = [(exons[x][1], exons[x+1][0]) for x in range(len(exons)-1)]
                for j in juncs:
                    if j not in juncs_to_gene[chrom]:
                        juncs_to_gene[chrom][j] = set()
                    juncs_to_gene[chrom][j].add(gene)
        allexons = sorted(list(gene_to_all_exons[gene]))
        newexons = []
        laststart, lastend = -1, -1
        for s, e in allexons:
            if s > lastend:
                if laststart != -1:
                    newexons.append((laststart, lastend))
                laststart, lastend = s, e
            else:
                lastend = max(lastend, e)
        if laststart != -1:
            newexons.append((laststart, lastend))
        gene_to_all_exons[gene] = newexons

    
    
    intronLocs, intronToGenome = {}, {}
    for g in genetoexons:
        chrom, start, end, strand, end5s = genetoinfo[g]
        for t in genetoexons[g]:
            myexons = sorted(genetoexons[g][t])
            first, last = myexons[0], myexons[-1]
            mylocs = [[0, myexons[0][0] - 500, myexons[0][0]]]  ##add start of transcript
            runningtot = 0
            for i in range(len(myexons) - 1):
                runningtot += myexons[i][1] - myexons[i][0]  ##add size of last exon
                mylocs.append([runningtot, myexons[i][1], myexons[i + 1][0]])  ##add intron
            runningtot += myexons[-1][1] - myexons[-1][0]
            mylocs.append([runningtot, myexons[-1][1], myexons[-1][1] + 500])
            if strand == '-':
                mylocs = [[runningtot - mylocs[x][0], mylocs[x][1], mylocs[x][2]] for x in range(len(mylocs))]
            intronLocs[t] = sorted([x[0] for x in mylocs])
            intronToGenome[t] = {x[0]: (x[1], x[2]) for x in mylocs}
    

    maxpromiscuity = 4

    gene_to_paralogs = get_paralog_ref(os.path.realpath(__file__).split('flair_fusion')[0] + 'dgd_Hsa_all_v71.tsv')
    
    print('loading transcriptomic chimeras')
    tchim = id_chimeras('transcriptomic', args.transcriptchimbam, genetoinfo, chrom_to_gene_pos, gene_to_all_exons, juncs_to_gene, gene_to_paralogs, genetoname, args.support, maxloci=args.maxloci, reqdisttostart=args.max_dist_to_TSS, maxpromiscuity=4, intronLocs=intronLocs, intronToGenome=intronToGenome)
    print('loading genomic chimeras')
    
    combchim = id_chimeras('genomic', args.genomechimbam, genetoinfo, chrom_to_gene_pos, gene_to_all_exons, juncs_to_gene, gene_to_paralogs, genetoname, args.support, maxloci=args.maxloci, reqdisttostart=args.max_dist_to_TSS, maxpromiscuity=4)
    
    print('combining genomic and transcriptomic')
    #
    for f in tchim:
        if f not in combchim:
            combchim[f] = tchim[f]
        else:
            combchim[f]['reads'] = combchim[f]['reads'] | tchim[f]['reads']
            combchim[f]['disttostart'].extend(tchim[f]['disttostart'])
            combchim[f]['qdist'].extend(tchim[f]['qdist'])
            genes = f.split('__')
            for g in genes:
                gends, tends = combchim[f][g][1:], tchim[f][g][1:]
                if gends[0] < gends[1]:
                    combchim[f][g][1] = min(gends[0], tends[0])
                    combchim[f][g][2] = max(gends[1], tends[1])
                else:
                    combchim[f][g][1] = max(gends[0], tends[0])
                    combchim[f][g][2] = min(gends[1], tends[1])

    bedout = open(args.output + '.prelimfusions.bed', 'w')
    fusionreads = set()
    for f in combchim:
        fusionreads.update(combchim[f]['reads'])
        genes = f.split('__')
        for g in genes:
            startdiststr = ','.join([str(x) for x in combchim[f]['disttostart']])
            qdiststr = ','.join([str(x) for x in combchim[f]['qdist']])
            bedname = g + '__' + '--'.join(genes) #+ '__' + startdiststr + '__' + qdiststr
            strand = '+' if combchim[f][g][1] < combchim[f][g][2] else '-'
            bedpos = [str(x) for x in sorted(combchim[f][g][1:])]

            bedline = [combchim[f][g][0]] + bedpos + [bedname, str(len(combchim[f]['reads'])), strand, startdiststr]#, qdiststr] #+ [str(x) for x in combchim[f][g][1:]] + ['255,0,0']
            bedout.write('\t'.join(bedline) + '\n')

    bedout.close()

    if os.path.getsize(args.output + '.prelimfusions.bed') == 0:
        report_nofusions(args.output)
        return
    
    print('obtaining fusion reads')

    seenreads = set()
    freadsname = args.output + '.fusionreads.prelim.fa'
    out_fa = open(freadsname, 'w')
    bamfile = pysam.AlignmentFile(args.genomechimbam, 'rb')
    for a in bamfile:
        if not a.is_secondary and not a.is_supplementary and a.query_name in fusionreads:
            out_fa.write('>' + a.query_name + '\n' + a.get_forward_sequence() + '\n')
            seenreads.add(a.query_name)
    bamfile.close()
    bamfile = pysam.AlignmentFile(args.genomechimbam, 'rb')
    for a in bamfile:
        if not a.is_secondary and not a.is_supplementary and a.query_name in fusionreads and not a.query_name in seenreads:
            out_fa.write('>' + a.query_name + '\n' + a.get_forward_sequence() + '\n')
    bamfile.close()
    out_fa.close()


    print('generating synthetic reference')

    makesynthcommand = ['python3', path + 'make_synthetic_fusion_reference.py', '-a', args.gtf, '-g', args.genome,
                        '-o', args.output, '-c', args.output + '.prelimfusions.bed']
    pipettor.run([makesynthcommand])
    if os.path.getsize(args.output + '-syntheticFusionGenome.fa') == 0:
        report_nofusions(args.output)
        return

    # FIXME: pipettor by default captures stderr to include in an error message, this hides logging from
    # lower level.  Maybe don't capture when running flair subtools
    faidxcommand = ['samtools', 'faidx', args.output + '-syntheticFusionGenome.fa']
    pipettor.run([faidxcommand])

    print('aligning to synthetic fusion genome')
    align_to_synth_genome(args.output + '-syntheticFusionGenome.fa', freadsname, args.output + '.syntheticAligned.nosplice.bam', 
                          ['-s', str(args.minfragmentsize), '-t', str(args.threads), '-un', '--secondary=no', '-G', '1000k'])
    align_to_synth_genome(args.output + '-syntheticFusionGenome.fa', freadsname, args.output + '.syntheticAligned.withsplice.bam', 
                          ['-s', str(args.minfragmentsize), '-t', str(args.threads), '--secondary=no', '-G', '1000k'])
    
    rname_to_read = {}
    with pysam.AlignmentFile(args.output + '.syntheticAligned.withsplice.bam', 'rb') as bamfile:
        for a in bamfile:
            if a.is_mapped and not a.is_secondary and not a.is_supplementary:
                rname_to_read[a.query_name] = a
    with pysam.AlignmentFile(args.output + '.syntheticAligned.nosplice.bam', 'rb') as bamfile:
        for a in bamfile:
            if a.is_mapped and not a.is_secondary and not a.is_supplementary:
                if a.query_name not in rname_to_read:
                    rname_to_read[a.query_name] = a
                else:
                    nosplice_match = a.get_cigar_stats()[0][0]
                    nosplice_insert = a.get_cigar_stats()[0][1]
                    withsplice_match = rname_to_read[a.query_name].get_cigar_stats()[0][0]
                    withsplice_insert = rname_to_read[a.query_name].get_cigar_stats()[0][1]
                    if nosplice_match > withsplice_match and nosplice_insert < withsplice_insert:
                        rname_to_read[a.query_name] = a
    template = pysam.AlignmentFile(args.output + '.syntheticAligned.withsplice.bam', 'rb')
    outbam = pysam.AlignmentFile(args.output + '.syntheticAligned.unsorted.bam', 'wb', template=template)
    for a in rname_to_read:
        outbam.write(rname_to_read[a])
    template.close()
    outbam.close()

    pipettor.run([('samtools', 'sort', '-o', args.output + '.syntheticAligned.bam', args.output + '.syntheticAligned.unsorted.bam')])
    pipettor.run([('samtools', 'index', args.output + '.syntheticAligned.bam')])
    pipettor.run([('rm', args.output + '.syntheticAligned.withsplice.bam', args.output + '.syntheticAligned.withsplice.bam.bai', 
                   args.output + '.syntheticAligned.nosplice.bam', args.output + '.syntheticAligned.nosplice.bam.bai',
                   args.output + '.syntheticAligned.unsorted.bam')])
    
    
    print('getting ss')
    # bamtobedcmd = ('bedtools', 'bamtobed', '-bed12', '-i', args.output + '.syntheticAligned.bam')
    # getsscommand = ['python3', path + 'synthetic_splice_sites.py', args.output + '.syntheticAligned.bed',
                        # args.output + '-syntheticReferenceAnno.gtf', args.output + '.syntheticAligned.SJ.bed', args.output + '-syntheticBreakpointLoc.bed', '8', '2', args.output + '-syntheticFusionGenome.fa']#'15', '2']
    # pipettor.run([bamtobedcmd], stdout=args.output + '.syntheticAligned.bed')
    # pipettor.run([getsscommand])
    
    ipcmd = ('intronProspector', f'--genome-fasta={args.output}-syntheticFusionGenome.fa', f'--intron-bed6={args.output}.syntheticAligned.IPSJ.bed', '-C', '0.0', '--sj-filter=all', f'{args.output}.syntheticAligned.bam')
    pipettor.run([ipcmd])
    
    fusiontobp = {}
    for line in open(f'{args.output}-syntheticBreakpointLoc.bed'):
        line = line.rstrip().split('\t')
        fusiontobp[line[0]] = int(line[1])
    

    fusion_to_bp_sj = {f:False for f in fusiontobp}
    good_sj = []
    
    for line in open(f'{args.output}.syntheticAligned.IPSJ.bed'):
        line = line.rstrip().split('\t')
        fusion = line[0]
        start, end = int(line[1]), int(line[2])
        readsup = int(line[4])
        sjmotif = line[3].split('_')[-1]
        strand = line[5]
        if readsup >= 2:
            if sjmotif in {"GT/AG", "GC/AG", "AT/AC"} and strand == '+': #for synthetic alignment, all junctions should be '+'
                good_sj.append(line)
                if start < fusiontobp[fusion] < end:
                    fusion_to_bp_sj[fusion] = True
    
    for line in open(f'{args.output}.syntheticAligned.IPSJ.bed'):
        line = line.rstrip().split('\t')
        fusion = line[0]
        start, end = int(line[1]), int(line[2])
        readsup = int(line[4])
        sjmotif = line[3].split('_')[-1]
        strand = line[5]
        if readsup >= 2 and fusion_to_bp_sj[fusion] == False and start < fusiontobp[fusion] < end: ##no good breakpoint junctions yet
            good_sj.append(line)


    out = open(f'{args.output}.syntheticAligned.SJ.bed', 'w')
    for line in good_sj:
        out.write('\t'.join(line) + '\n')
    out.close()


    transcriptome_command = ['flair', 'transcriptome',
                             '--genome_aligned_bam', args.output + '.syntheticAligned.bam',
                             '--genome', args.output + '-syntheticFusionGenome.fa',
                             '--gtf', args.output + '-syntheticReferenceAnno.gtf',
                             '--ss_window', '8',
                             '--generate_map',
                             '--quality', '0',
                             '--sjc_support', '2',
                             '--allow_paralogs',
                             '--end_window', '300',
                            #  '--no_check_splice',
                            #  '--no_stringent',
                             '--no_align_to_annot',
                             '--fusion_breakpoints', args.output + '-syntheticBreakpointLoc.bed',
                             '--output', args.output + '.syntheticAligned.flair',]
    # only include junction if any where found
    junc_bed = args.output + '.syntheticAligned.SJ.bed'
    if os.path.exists(junc_bed) and (os.path.getsize(junc_bed) > 0):
        transcriptome_command.extend(['--junction_bed', junc_bed])
    
    print('generating fusion transcriptome')
    print(' '.join(transcriptome_command))
    pipettor.run(transcriptome_command)

    ##clean up isoform/gene names for args.output + '.combined.isoform.read.map.txt', '.isoforms.bed', '.isoforms.fa'
    oldnametonewname = {}
    out = open(args.output + '.syntheticAligned.isoforms.bed', 'w')
    c = 0
    for line in open(args.output + '.syntheticAligned.flair.isoforms.bed'):
        c +=1
        line = line.rstrip().split('\t')
        newname = 'fusioniso' + str(c) + '_' + line[0]
        oldnametonewname[line[3]] = newname
        line[3] = newname
        out.write('\t'.join(line) + '\n')
    out.close()
    out = open(args.output + '.syntheticAligned.isoform.read.map.txt', 'w')
    for line in open(args.output + '.syntheticAligned.flair.isoform.read.map.txt'):
        line = line.split('\t', 1)
        if line[0] in oldnametonewname:
            line[0] = oldnametonewname[line[0]]
            out.write('\t'.join(line))
    out.close()
    out = open(args.output + '.syntheticAligned.isoform.counts.txt', 'w')
    for line in open(args.output + '.syntheticAligned.flair.isoform.counts.txt'):
        line = line.split('\t', 1)
        if line[0] in oldnametonewname:
            line[0] = oldnametonewname[line[0]]
            out.write('\t'.join(line))
    out.close()

    print('converting coordinates from synthetic to genomic')
    convert_synthetic_isos(args.output + '.syntheticAligned.isoforms.bed',
                      args.output + '.syntheticAligned.isoform.read.map.txt', freadsname,
                      args.output + '-syntheticBreakpointLoc.bed', args.output + '.fusions.isoforms.bed')
    goodisos = set()
    for line in open(args.output + '.fusions.isoforms.bed'):
        line = line.rstrip().split('\t')
        goodisos.add('_'.join(line[3].split('_')[1:]))

    out = open(args.output + '.fusions.isoforms.fa', 'w')
    good = False
    for line in open(args.output + '.syntheticAligned.flair.isoforms.fa'):
        if line[0] == '>':
            oldname = line[1:].rstrip()
            if oldnametonewname[oldname] in goodisos: good = True
            else: good = False
            if good: out.write('>' + oldnametonewname[oldname] + '\n')
        elif good: out.write(line)
    out.close()

    os.rename(args.output + '.syntheticAligned.isoform.read.map.txt', args.output + '.fusion.isoform.read.map.txt')

    # removing extra FLAIR files
    if not args.keep_intermediate:
        for filename in glob.glob(args.output + '.syntheticAligned.flair*'):
            os.remove(filename)






if __name__ == '__main__':
    detectfusions()
