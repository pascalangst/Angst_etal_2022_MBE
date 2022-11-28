# The basic pipeline we'll follow here is to put all the de novo transcriptomes together and then use isoseq type analyses to
# simplify the hint dataset. The necessary scripts can be pulled from the cDNA_cupcake repo, https://github.com/Magdoll/cDNA_Cupcake
# first thing we cat the individual transcript files, then we move to mapping
minimap2 -t 55 -ax splice -uf --secondary=no -C5 -O6,24 -B4 WJBH01.fasta \ 
all_transcripts.fasta > hq_isoforms.fasta.sam 2> hq_isoforms.fasta.sam.log

# now we use cDNA-Cupcake scripts to collapse isoforms
sort -k 3,3 -k 4,4n hq_isoforms.fastq.sam > hq_isoforms.fastq.sorted.sam
collapse_isoforms_by_sam.py --input all_transcripts.fasta \
-s hq_isoforms.fasta.sorted.sam --dun-merge-5-shorter -o d_sinensis.transcripts

# we now have our target transcripts in the form of d_sinensis.transcripts.collapsed.rep.fa
