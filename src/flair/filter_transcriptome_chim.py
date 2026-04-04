import sys
import pysam

# FIXME: delete this and replace it with
#  samtools view -F 0x104 -e '[SA] != ""' in.bam

infile = pysam.AlignmentFile(sys.argv[1], 'r') ##sam file input
outname = sys.argv[2] ##bam suffix
outfile = pysam.AlignmentFile(outname, 'wb', template=infile)
for align in infile:
    if align.is_mapped and not align.is_secondary:
        if align.has_tag('SA'):
            outfile.write(align)
infile.close()
outfile.close()
