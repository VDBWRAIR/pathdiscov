#!/bin/bash

# path to scripts
d=$1	

sample=VLP-4

r1=../data/VLP-4_S1_L001_R1_001.fastq
r2=../data/VLP-4_S1_L001_R2_001.fastq

$d/pathogen.pl --example > param.txt

$d/pathogen.pl --sample $sample --command step1 quality_filter host_map ray2_assembly iterative_blast_phylo --paramfile param.txt --outputdir ../results/$sample --R1 $r1 --R2 $r2

$d/pathogen.pl --sample $sample --command step1 quality_filter host_map iterative_blast_phylo_2 --paramfile param.txt --outputdir ../results/$sample --R1 $r1 --R2 $r2

