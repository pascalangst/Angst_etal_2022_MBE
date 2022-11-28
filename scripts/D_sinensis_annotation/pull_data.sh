# pull the raw data from EMBL
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/090/SRR10389290/SRR10389290_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/090/SRR10389290/SRR10389290_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/093/SRR10389293/SRR10389293_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/093/SRR10389293/SRR10389293_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/094/SRR10389294/SRR10389294_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/094/SRR10389294/SRR10389294_2.fastq.gz

# run fastqc and multiqc
fastqc -t 6 *.fastq.gz
multiqc SRR*

# There is definite adapter contamination so let's run fastp
mkdir trim; cd trim; ln -s ../*.fastq.gz .
fastp -i SRR10389290_1.fastq.gz -I SRR10389290_2.fastq.gz -o SRR10389290_out.R1.fq.gz -O SRR10389290_out.R2.fq.gz -w 16
fastp -i SRR10389293_1.fastq.gz -I SRR10389293_2.fastq.gz -o SRR10389293_out.R1.fq.gz -O SRR10389293_out.R2.fq.gz -w 16
fastp -i SRR10389294_1.fastq.gz -I SRR10389294_2.fastq.gz -o SRR10389294_out.R1.fq.gz -O SRR10389294_out.R2.fq.gz -w 16

# re-run fastqc and multiqc
fastqc -t 6 *.fq.gz
multiqc SRR*

# trimming worked fine, we can now move to assemlby
