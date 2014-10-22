#!/bin/bash

# count counts

outputdir=$1

# scripts directory
d=$( dirname $( readlink -m $0 ) )

if cd ${outputdir}; then

	cat reads.txt | sort -k1,1n > reads.sort.txt
#	rm read.txt

	# get unblasted contigs 
	cat../iterative_blast_phylo_1/iterative_blast_phylo_1.R1 | awk '$1 ~ />/{print substr($1,2)}' > contig_unblast.txt
	
	# get number of reads (R1+R2) that comprise each contig
	ln -s ../ray2_assembly_1/contig_numreads.txt	
	
	r1_assem=$( cat *.R1.count.txt | awk '$1=="ray_1"{print $2}' )
	r2_assem=$( cat *.R2.count.txt | awk '$1=="ray_1"{print $2}' )

fi

