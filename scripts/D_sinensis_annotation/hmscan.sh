#!/bin/sh

sample=$1

hmmscan --cpu 1 --domtblout ${sample}.pfam.domtblout ~/bioinformatics/pfam/Pfam-A.hmm ${sample}
