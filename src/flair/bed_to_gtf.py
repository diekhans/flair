#!/usr/bin/env python3
import sys
import argparse
from flair import FlairInputDataError

def main():
    parser = argparse.ArgumentParser(description='options')
    parser.add_argument('inputfile', type=str,
            action='store', help='isoforms in bed format')
    parser.add_argument('--force', action='store_true', dest='force',
            help='specify to not split isoform name by underscore into isoform and gene ids')
    parser.add_argument('--add_reference_transcript_id', action='store_true', dest='reference_transcript_id',
            help='specify to add reference_transcript_id attribute')
    parser.add_argument('--noCDS', action='store_true',
                        help='do not carry forward CDS from bed file (thickstart and thickend) to gtf file')
    args = parser.parse_args()
    bed_to_gtf(query=args.inputfile, force=args.force, outputfile='/dev/stdout',
               reference_transcript_id=args.reference_transcript_id, useCDS= not args.noCDS)


def split_iso_gene(iso_gene):
    if '_chr' in iso_gene:
        iso = iso_gene[:iso_gene.rfind('_chr')]
        gene = iso_gene[iso_gene.rfind('_chr')+1:]
    elif '_XM' in iso_gene:
        iso = iso_gene[:iso_gene.rfind('_XM')]
        gene = iso_gene[iso_gene.rfind('_XM')+1:]
    elif '_XR' in iso_gene:
        iso = iso_gene[:iso_gene.rfind('_XR')]
        gene = iso_gene[iso_gene.rfind('_XR')+1:]
    elif '_NM' in iso_gene:
        iso = iso_gene[:iso_gene.rfind('_NM')]
        gene = iso_gene[iso_gene.rfind('_NM')+1:]
    elif '_NR' in iso_gene:
        iso = iso_gene[:iso_gene.rfind('_NR')]
        gene = iso_gene[iso_gene.rfind('_NR')+1:]
    elif '_R2_' in iso_gene:
        iso = iso_gene[:iso_gene.rfind('_R2_')]
        gene = iso_gene[iso_gene.rfind('_R2_')+1:]
    else:
        iso = iso_gene[:iso_gene.rfind('_')]
        gene = iso_gene[iso_gene.rfind('_')+1:]
    return iso, gene

def bed_to_gtf(query, outputfile, force=False, reference_transcript_id=False, useCDS=True):
    outfile = open(outputfile, 'w')
    gene_to_transcript_lines = {}
    gene_to_chrom_strand = {}
    for line in open(query):
        line = line.rstrip().split('\t')
        start = int(line[1])
        chrom, strand, score, name, start = line[0], line[5], line[4], line[3], int(line[1])
        tstarts = [int(n) + start for n in line[11].rstrip(',').split(',')]
        bsizes = [int(n) for n in line[10].rstrip(',').split(',')]
        end, thick_start, thick_end = int(line[2]), int(line[6]), int(line[7])

        if '_' not in name and not force:
            raise FlairInputDataError('Entry name should contain underscore-delimited transcriptid and geneid like so: \n'
                             'ENST00000318842.11_ENSG00000156313.12 or a4bab8a3-1d28_chr8:232000\n'
                             'So no GTF conversion was done. Please run identify_gene_isoform first\n'
                             'for best results, or run with --force')

        if ';' in name:
            name = name.replace(';', ':')

        if force == True:
            transcript_id, gene_id = name, name
        else:
            transcript_id, gene_id = split_iso_gene(name)

        if gene_id not in gene_to_transcript_lines:
            gene_to_transcript_lines[gene_id] = []
            gene_to_chrom_strand[gene_id] = (chrom, strand)

        attributes = f'gene_id "{gene_id}"; transcript_id "{transcript_id}";'
        if reference_transcript_id and '-referencetranscript' in transcript_id:
            trimmed_transcript_id = transcript_id[:transcript_id.find('-referencetranscript')]
            attributes = f'gene_id "{gene_id}"; transcript_id "{trimmed_transcript_id}"; reference_transcript_id "{trimmed_transcript_id}";'
        gene_to_transcript_lines[gene_id].append([chrom, 'FLAIR', 'transcript', start+1, tstarts[-1]+bsizes[-1], '.', strand, '.', attributes])
        if thick_start != thick_end and (thick_start != start or thick_end != end) and useCDS:
            gene_to_transcript_lines[gene_id].append([chrom, 'FLAIR', 'CDS', thick_start+1, thick_end, '.', strand, '.', attributes])
            if strand == '+':
                gene_to_transcript_lines[gene_id].append([chrom, 'FLAIR', 'start_codon', thick_start+1, thick_start+3, '.', strand, '.', attributes])
                gene_to_transcript_lines[gene_id].append([chrom, 'FLAIR', '5UTR', start+1, thick_start+1, '.', strand, '.', attributes])
                gene_to_transcript_lines[gene_id].append([chrom, 'FLAIR', '3UTR', thick_end, tstarts[-1]+bsizes[-1], '.', strand, '.', attributes])
            elif strand == '-':
                gene_to_transcript_lines[gene_id].append([chrom, 'FLAIR', 'start_codon', thick_end-2, thick_end, '.', strand, '.', attributes])
                gene_to_transcript_lines[gene_id].append([chrom, 'FLAIR', '3UTR', start+1, thick_start+1, '.', strand, '.', attributes])
                gene_to_transcript_lines[gene_id].append([chrom, 'FLAIR', '5UTR', thick_end, tstarts[-1]+bsizes[-1], '.', strand, '.', attributes])
        # if strand == '-':  # to list exons in 5'->3'
        #       for b in range(len(tstarts)):  # exon number
        #               bi = len(tstarts) - 1 - b  # block index
        #               attributes = 'gene_id \"{}\"; transcript_id \"{}\"; exon_number \"{}\";'\
        #                                               .format(gene_id, transcript_id, b)
        #               print('\t'.join([chrom, 'FLAIR', 'exon', str(tstarts[bi]+1), \
        #                       str(tstarts[bi]+bsizes[bi]), '.', strand, '.', attributes]))
        # else:
        for b in range(len(tstarts)):
            attributes = f'gene_id "{gene_id}"; transcript_id "{transcript_id}"; exon_number "{b}";'
            if reference_transcript_id and '-referencetranscript' in transcript_id:
                attributes = f'gene_id "{gene_id}"; transcript_id "{trimmed_transcript_id}"; exon_number "{b}"; reference_transcript_id "{trimmed_transcript_id}";'
            gene_to_transcript_lines[gene_id].append([chrom, 'FLAIR', 'exon', tstarts[b]+1, tstarts[b]+bsizes[b], '.', strand, '.', attributes])
    for gene_id in gene_to_transcript_lines:
        attributes = f'gene_id "{gene_id}";'
        chrom, strand = gene_to_chrom_strand[gene_id]
        tlines = gene_to_transcript_lines[gene_id]
        gene_line = [chrom, 'FLAIR', 'gene', min([x[3] for x in tlines]), max([x[4] for x in tlines]), '.', strand, '.', attributes]
        outfile.write('\t'.join([str(x) for x in gene_line]) + '\n')
        for line in tlines:
            outfile.write('\t'.join([str(x) for x in line]) + '\n')
    outfile.close()
        

if __name__ == "__main__":
    main()
