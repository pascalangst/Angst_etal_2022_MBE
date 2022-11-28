# We're going to generate protein hints from the de novo transcriptomes of D. sinensis and the 
# proteome of D. magna v.3.0
# Let's use the Transdecoder approach to get the proteins out of hte transcriptomes
# For each of the following the commands were run inside the Trinity output directories
TransDecoder.LongOrfs -t Trinity-GG.fasta

# Now that we have the ORFs we're going to use diamond and hmmer to find matches to swissprot and pfam
diamond blastp --query Trinity-GG.transdecoder_dir/longest_orfs.pep --db ~/bioinformatics/swissprot/uniprot_sprot.dmnd \
--outfmt 6 --ultra-sensitive --max-target-seqs 1 --evalue 1e-5 > diamond.Trinity-GG.outfmt6

# Now let's use hmmer to generate matches against pfam. This is super slow so I'm going to use parallels, see script hmmscan.sh, though the basic form is:
hmmscan --cpu 1 --domtblout ${sample}.pfam.domtblout ~/bioinformatics/pfam/Pfam-A.hmm ${sample}

# let's split up the longest_orfs.pep into 100 smaller fasta files; pyfasta can be pulled using bioconda
pyfasta split -n 100 longest_orfs.pep

# and we call it as follows:
ls longest_orfs.pep.* | parallel -j55 -k bash hmscan.sh {}

# because the resultant files have a header we need to cut this and then cat these files; see trim_header.sh and use the same as hmmscan.sh
ls *.pfam.domtblout | parallel -j12 -k bash trim_headers.sh {}
cat *.out > compile.pfam.domtblout
head -n 3 longest_orfs.pep.000.pfam.domtblout | cat - compile.pfam.domtblout > pfam.Trinity-GG.domtblout

# now we can finish up the Transdecoder pipeline ande get our target proteins
TransDecoder.Predict -t Trinity-GG.fasta --retain_pfam_hits pfam.Trinity-GG.domtblout \
 --retain_blastp_hits diamond.Trinity-GG.outfmt6

# Because of the original renaming we did above we need to make sure our different protein files have unique names, so we do:
awk '/^>/{print ">SRR10389290" ++i; next}{print}' Trinity-GG.fasta.transdecoder.pep \
> SRR10389290.fasta.transdecoder.rename.pep

# Once we move the individual protein files into the same directory we can run the following:
cat *transdecoder.rename.pep > proteins.fasta

# Once we can these with the D. magna proteins we have our protein hints for MAKER2.
