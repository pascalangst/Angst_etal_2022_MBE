#!/bin/bash

### SNPGenie (for piS and piN)
## SNPGenie runs on individual contigs and in one direction only (split contigs and then revert them)

# split reference into contigs
mkdir ref_split
awk -F "|" '/^>/ {close(F);F = "ref_split/"$1".fasta"} {print >> F}' 29082016_Xinb3_ref.fasta

# split vcf into contigs (masked / missing sites should have entries for each field)
sed -e 's~\./\.~./.:.,.:.:.:.,.,.~g' magna_spr2014_210930_filter.vcf > magna_spr2014_210930_filter.masked.vcf
Rscript vcfsplit.r magna_spr2014_210930_filter.masked.vcf

# split annotation file into contigs
Rscript gffsplit.r FI-XINB3_draft_annotation_23112017.gff

# revert contigs
for contig in ref_split/*.fasta
do
	contigname=$(echo ${contig:11:${#contig}} | sed 's/.fasta//' )
	echo $contigname
	mv vcf_split/$contigname\|quiver vcf_split/$contigname.vcf
	mv gff_split/falcon_$contigname gff_split/$contigname.gff
	gffread -E gff_split/$contigname.gff -T -o gff_split/$contigname.gtf
	mv ref_split/\>$contigname.fasta ref_split/$contigname.fasta
	contiglength=$(bioawk -c fastx '{print length($seq)}' < ref_split/$contigname.fasta)
	echo $contiglength
	perl ~/programs/SNPGenie/vcf2revcom.pl vcf_split/$contigname.vcf $contiglength
	perl ~/programs/SNPGenie/gtf2revcom.pl gff_split/$contigname.gtf $contiglength
	perl ~/programs/SNPGenie/fasta2revcom.pl ref_split/$contigname.fasta
	mv *_re* ref_split/.
done

# create input directory
mkdir input/
cp gff_split/*.gtf input/.
cp vcf_split/*.vcf input/.
cp ref_split/*.fasta input/.

# check if each contig has at least one annotation and one SNP, else remove contig (see below removed contigs)
for contig in input/*[^m].fasta
do
	contigname=$(echo ${contig:6:${#contig}} | sed 's/.fasta//' )
	if find input/*$contigname* -maxdepth 1 -empty | grep -q .; then
		echo "One file is empty, continue with next contig"
		echo $contigname
		rm input/$contigname*
		continue
  	fi
  	sed -i 's/transcript_id/gene_id/g' input/$contigname.gtf
  	sed -i 's/transcript_id/gene_id/g' input/${contigname}_revcom.gtf
done

# create output directory and run software
mkdir output
for contig in input/*[^m].fasta
do 
	contigname=$(echo ${contig:6:${#contig}} | sed 's/.fasta//' )
	echo $contigname
	perl ~/programs/SNPGenie/snpgenie.pl --vcfformat=4 --snpreport=input/$contigname.vcf --fastafile=input/$contigname.fasta --gtffile=input/$contigname.gtf --outdir output/SNPGenie_Results_${contigname}
	perl ~/programs/SNPGenie/snpgenie.pl --vcfformat=4 --snpreport=input/${contigname}_revcom.vcf --fastafile=input/${contigname}_revcom.fasta --gtffile=input/${contigname}_revcom.gtf --outdir output/SNPGenie_Results_${contigname}_rev
done





















# Contigs without annotation and/or SNP
#One file is empty, continue with next contig
#000041F
#One file is empty, continue with next contig
#000100F
#One file is empty, continue with next contig
#000130F
#One file is empty, continue with next contig
#000133F
#One file is empty, continue with next contig
#000144F
#One file is empty, continue with next contig
#000145F
#One file is empty, continue with next contig
#000147F
#One file is empty, continue with next contig
#000165F
#One file is empty, continue with next contig
#000174F
#One file is empty, continue with next contig
#000176F
#One file is empty, continue with next contig
#000182F
#One file is empty, continue with next contig
#000187F
#One file is empty, continue with next contig
#000190F
#One file is empty, continue with next contig
#000201F
#One file is empty, continue with next contig
#000208F
#One file is empty, continue with next contig
#000209F
#One file is empty, continue with next contig
#000211F
#One file is empty, continue with next contig
#000212F
#One file is empty, continue with next contig
#000213F
#One file is empty, continue with next contig
#000221F
#One file is empty, continue with next contig
#000222F
#One file is empty, continue with next contig
#000227F
#One file is empty, continue with next contig
#000229F
#One file is empty, continue with next contig
#000233F
#One file is empty, continue with next contig
#000235F
#One file is empty, continue with next contig
#000236F
#One file is empty, continue with next contig
#000237F
#One file is empty, continue with next contig
#000238F
#One file is empty, continue with next contig
#000239F
#One file is empty, continue with next contig
#000240F
#One file is empty, continue with next contig
#000241F
#One file is empty, continue with next contig
#000242F
#One file is empty, continue with next contig
#000244F
#One file is empty, continue with next contig
#000247F
#One file is empty, continue with next contig
#000249F
#One file is empty, continue with next contig
#000254F
#One file is empty, continue with next contig
#000255F
#One file is empty, continue with next contig
#000258F
#One file is empty, continue with next contig
#000259F
#One file is empty, continue with next contig
#000262F
#One file is empty, continue with next contig
#000263F
#One file is empty, continue with next contig
#000265F
#One file is empty, continue with next contig
#000266F
#One file is empty, continue with next contig
#000267R
#One file is empty, continue with next contig
#000268F
#One file is empty, continue with next contig
#000271F
#One file is empty, continue with next contig
#000275F
#One file is empty, continue with next contig
#000280F
#One file is empty, continue with next contig
#000283F
#One file is empty, continue with next contig
#000284F
#One file is empty, continue with next contig
#000287F
#One file is empty, continue with next contig
#000291F
#One file is empty, continue with next contig
#000292F
#One file is empty, continue with next contig
#000295F
#One file is empty, continue with next contig
#000300F
#One file is empty, continue with next contig
#000301F
#One file is empty, continue with next contig
#000304F
#One file is empty, continue with next contig
#000310F
#One file is empty, continue with next contig
#000312F
#One file is empty, continue with next contig
#000318F
#One file is empty, continue with next contig
#000321F
#One file is empty, continue with next contig
#000322F
#One file is empty, continue with next contig
#000323F
#One file is empty, continue with next contig
#000329F
#One file is empty, continue with next contig
#000332F
#One file is empty, continue with next contig
#000336F
#One file is empty, continue with next contig
#000343F
#One file is empty, continue with next contig
#000348F
#One file is empty, continue with next contig
#000349F
#One file is empty, continue with next contig
#000350F
#One file is empty, continue with next contig
#000354F
#One file is empty, continue with next contig
#000355F
#One file is empty, continue with next contig
#000357F
#One file is empty, continue with next contig
#000358F
#One file is empty, continue with next contig
#000359F
#One file is empty, continue with next contig
#000362F
#One file is empty, continue with next contig
#000365F
#One file is empty, continue with next contig
#000367F
#One file is empty, continue with next contig
#000370F
#One file is empty, continue with next contig
#000371F
#One file is empty, continue with next contig
#000373F
#One file is empty, continue with next contig
#000375F
#One file is empty, continue with next contig
#000378F
#One file is empty, continue with next contig
#000379F
#One file is empty, continue with next contig
#000380F
#One file is empty, continue with next contig
#000381F
#One file is empty, continue with next contig
#000384F
#One file is empty, continue with next contig
#000386F
#One file is empty, continue with next contig
#000388F
#One file is empty, continue with next contig
#000390F
#One file is empty, continue with next contig
#000395F
#One file is empty, continue with next contig
#000397F
#One file is empty, continue with next contig
#000398F
#One file is empty, continue with next contig
#000404F
#One file is empty, continue with next contig
#000405F
#One file is empty, continue with next contig
#000408F
#One file is empty, continue with next contig
#000409F
#One file is empty, continue with next contig
#000410F
#One file is empty, continue with next contig
#000412F
#One file is empty, continue with next contig
#000417F
#One file is empty, continue with next contig
#000419F
#One file is empty, continue with next contig
#000422F
#One file is empty, continue with next contig
#000424F
#One file is empty, continue with next contig
#000426F
#One file is empty, continue with next contig
#000427F
#One file is empty, continue with next contig
#000429F
#One file is empty, continue with next contig
#000432F
#One file is empty, continue with next contig
#000433F
#One file is empty, continue with next contig
#000434F
#One file is empty, continue with next contig
#000440F
#One file is empty, continue with next contig
#000444F
#One file is empty, continue with next contig
#000446F
#One file is empty, continue with next contig
#000447F
#One file is empty, continue with next contig
#000450F
#One file is empty, continue with next contig
#000452F
#One file is empty, continue with next contig
#000454F
#One file is empty, continue with next contig
#000455F
#One file is empty, continue with next contig
#000457F
#One file is empty, continue with next contig
#000458F
#One file is empty, continue with next contig
#000461F
#One file is empty, continue with next contig
#000464F
#One file is empty, continue with next contig
#000465F
#One file is empty, continue with next contig
#000468F
#One file is empty, continue with next contig
#000470F
#One file is empty, continue with next contig
#000471F
#One file is empty, continue with next contig
#000475F
#One file is empty, continue with next contig
#000476F
#One file is empty, continue with next contig
#000478F
#One file is empty, continue with next contig
#000479F
#One file is empty, continue with next contig
#000480F
#One file is empty, continue with next contig
#000481F
#One file is empty, continue with next contig
#000482F
#One file is empty, continue with next contig
#000483F
#One file is empty, continue with next contig
#000484F
#One file is empty, continue with next contig
#000485F
#One file is empty, continue with next contig
#000488F
#One file is empty, continue with next contig
#000489F
#One file is empty, continue with next contig
#000491F
#One file is empty, continue with next contig
#000492F
#One file is empty, continue with next contig
#000495F
#One file is empty, continue with next contig
#000496F
#One file is empty, continue with next contig
#000497F
#One file is empty, continue with next contig
#000500F
#One file is empty, continue with next contig
#000501F
#One file is empty, continue with next contig
#000504F
#One file is empty, continue with next contig
#000506F
#One file is empty, continue with next contig
#000508F
#One file is empty, continue with next contig
#000509F
#One file is empty, continue with next contig
#000511F
#One file is empty, continue with next contig
#000512F
#One file is empty, continue with next contig
#000514F
#One file is empty, continue with next contig
#000515F
#One file is empty, continue with next contig
#000517F
#One file is empty, continue with next contig
#000520F
#One file is empty, continue with next contig
#000522F
#One file is empty, continue with next contig
#000524F
#One file is empty, continue with next contig
#000525F
#One file is empty, continue with next contig
#000527F
#One file is empty, continue with next contig
#000528F
#One file is empty, continue with next contig
#000529F
#One file is empty, continue with next contig
#000531F
#One file is empty, continue with next contig
#000533F
#One file is empty, continue with next contig
#000534F
#One file is empty, continue with next contig
#000535F
#One file is empty, continue with next contig
#000536F
#One file is empty, continue with next contig
#000537F
#One file is empty, continue with next contig
#000538F
#One file is empty, continue with next contig
#000539F
#One file is empty, continue with next contig
#000540F
#One file is empty, continue with next contig
#000542F
#One file is empty, continue with next contig
#000543F
#One file is empty, continue with next contig
#000545F
#One file is empty, continue with next contig
#000546F
#One file is empty, continue with next contig
#000547F
#One file is empty, continue with next contig
#000548F
#One file is empty, continue with next contig
#000550F
#One file is empty, continue with next contig
#000551F
#One file is empty, continue with next contig
#000552F
#One file is empty, continue with next contig
#000553F
#One file is empty, continue with next contig
#000554F
#One file is empty, continue with next contig
#000555F
#One file is empty, continue with next contig
#000556F
#One file is empty, continue with next contig
#000557F
#One file is empty, continue with next contig
#000558F
#One file is empty, continue with next contig
#000559F
#One file is empty, continue with next contig
#000560F
#One file is empty, continue with next contig
#000561F
#One file is empty, continue with next contig
#000563F
#One file is empty, continue with next contig
#000564F
#One file is empty, continue with next contig
#000565F
#One file is empty, continue with next contig
#000566F
#One file is empty, continue with next contig
#000567F
#One file is empty, continue with next contig
#000568F
#One file is empty, continue with next contig
#000569F
#One file is empty, continue with next contig
#000571F
#One file is empty, continue with next contig
#000572F
#One file is empty, continue with next contig
#000574F
#One file is empty, continue with next contig
#000575F
#One file is empty, continue with next contig
#000576F
#One file is empty, continue with next contig
#000577F
#One file is empty, continue with next contig
#000578F
#One file is empty, continue with next contig
#000579F
#One file is empty, continue with next contig
#000580F
#One file is empty, continue with next contig
#000581F
#One file is empty, continue with next contig
#000582F
#One file is empty, continue with next contig
#000583F
#One file is empty, continue with next contig
#000584F
#One file is empty, continue with next contig
#000585F
#One file is empty, continue with next contig
#000586F
#One file is empty, continue with next contig
#000587F
#One file is empty, continue with next contig
#000588F
#One file is empty, continue with next contig
#000589F
#One file is empty, continue with next contig
#000590F
#One file is empty, continue with next contig
#594
#One file is empty, continue with next contig
#597
#One file is empty, continue with next contig
#598
#One file is empty, continue with next contig
#600
#One file is empty, continue with next contig
#602
#One file is empty, continue with next contig
#611
