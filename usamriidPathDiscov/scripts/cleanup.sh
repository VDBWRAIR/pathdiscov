#!/bin/bash

# clean up

# to do - add file compression (gzip) ... everything should be compressed ! 

outputdir=$1	# output directory

echo "[echo] clean step1"

if cd ${outputdir}/step1; then

	rm -f R1.fastq R2.fastq step1.R1 step1.R2

fi 

echo "[echo] quality_filter"

if cd ${outputdir}/quality_filter; then

	rm -f R1.cut2.fastq R1.cut.fastq R1.prinseq.bad.fastq 
	rm -f R2.cut2.fastq R2.cut.fastq R2.prinseq.bad.fastq

fi 

echo "[echo] clean host_map"

if cd ${outputdir}/host_map_1; then

	rm -f map*/*paired.fastq
	rm -f map*/*single.fastq

fi 

echo "[echo] clean ray2_assembly"

if cd ${outputdir}/ray2_assembly_1; then

	rm -f R1.paired.fastq R1.single.fastq R2.paired.fastq R2.single.fastq
	rm -f head*
	
fi 

echo "[echo] clean iterative_blast_phylo"

if cd ${outputdir}/iterative_blast_phylo_1; then

	rm -rf tmp*

fi 

if cd ${outputdir}/iterative_blast_phylo_2; then

	rm -rf tmp*

fi 
