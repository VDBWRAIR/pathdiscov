#!/bin/bash

# bowtie2 align to a reference; then extract from the original fastq file(s) the reads that didn't map
# assume file out.sam exists (bad style)
# this script wont properly account for mates where one mate maps and the other doesnt

outputdir=$1			# outputdir
d=$2					# scripts dir
r1=$3					# reads 1 file
r2=$4					# reads 2 file

if ! cd ${outputdir}; then
 
	echo "[error] changing directories"; 

else
			
	samtools view -bS out.sam > out.bam
	
	echo "[stats] paired stats"
	samtools flagstat out.bam
	
	# get unmapped read IDs
	cat out.sam | awk '$3=="*"' | cut -f1 > tmp1.paired.unmap.id
	sort -u tmp1.paired.unmap.id > tmp2.paired.unmap.id
	cat tmp2.paired.unmap.id | awk '{print "@"$0}' > paired.unmap.id
	
	rm tmp1.paired.unmap.id tmp2.paired.unmap.id
	rm out.sam
	
	# extract the 4 line chunks from a fastq file, given in argument 1 file, for IDs give in argument 2 file		
	cmd="${d}/fastq_extract_id.pl ${r1} paired.unmap.id > R1.unmap.fastq"
	echo "[cmd] "$cmd
	eval ${cmd}
	
	if [ -s ${r2} ]; then
			
		cmd="${d}/fastq_extract_id.pl ${r2} paired.unmap.id > R2.unmap.fastq"
		echo "[cmd] "$cmd
		eval ${cmd}
	
	fi
	
fi # cd

