# Here we'll start by pulling the D. sinensis genome
wget ftp://ftp.ebi.ac.uk/pub/databases/ena/wgs/suppressed/wjb/WJBH01.fasta.gz

# Let's create a STAR database from the reference
STAR --runThreadN 12 --runMode genomeGenerate --genomeDir WJBH01.1_star_index --genomeFastaFiles WJBH01.fasta --genomeSAindexNbases 12

# Now we can align the read datasets to the reference
STAR --genomeDir WJBH01.1_star_index/ --runThreadN 12 --readFilesIn rnaseq/trimmed/SRR10389290.R1.fq.gz rnaseq/trimmed/SRR10389290.R2.fq.gz \ 
--readFilesCommand zcat --outFileNamePrefix SRR10389290 --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard
STAR --genomeDir WJBH01.1_star_index/ --runThreadN 12 --readFilesIn rnaseq/trimmed/SRR10389293.R1.fq.gz rnaseq/trimmed/SRR10389293.R2.fq.gz \ 
--readFilesCommand zcat --outFileNamePrefix SRR10389293--outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard
STAR --genomeDir WJBH01.1_star_index/ --runThreadN 12 --readFilesIn rnaseq/trimmed/SRR10389294.R1.fq.gz rnaseq/trimmed/SRR10389294.R2.fq.gz \ 
--readFilesCommand zcat --outFileNamePrefix SRR10389294 --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard

# We can use these bam files to generate reference assisted transcriptomes with Trinity
Trinity --version
# Trinity version: Trinity-v2.12.0
Trinity --genome_guided_bam SRR10389290Aligned.sortedByCoord.out.bam --genome_guided_max_intron 10000 --max_memory 150G --CPU 16 \ 
  --output SRR10389290_trinity_out_dir
Trinity --genome_guided_bam SRR10389293Aligned.sortedByCoord.out.bam --genome_guided_max_intron 10000 --max_memory 150G --CPU 16 \
  --output SRR10389293_trinity_out_dir
Trinity --genome_guided_bam SRR10389294Aligned.sortedByCoord.out.bam --genome_guided_max_intron 10000 --max_memory 150G --CPU 16 \
  --output SRR10389294_trinity_out_dir
  
# Let's have a look at the BUSCO scores of the assemblies
run_BUSCO.py  -i SRR10389290.transcripts.fasta -o SRR10389290 -m tran -l ~/bioinformatics/busco_db/arthropoda_odb9 -c 55
run_BUSCO.py  -i SRR10389293.transcripts.fasta -o SRR10389293 -m tran -l ~/bioinformatics/busco_db/arthropoda_odb9 -c 55
run_BUSCO.py  -i SRR10389294.transcripts.fasta -o SRR10389294 -m tran -l ~/bioinformatics/busco_db/arthropoda_odb9 -c 55
