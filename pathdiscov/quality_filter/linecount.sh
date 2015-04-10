#!/bin/bash

# count the lines in either a fastq, a fasta, or a regular file and write to a file
# NOTE: if fasta ASSUME A SINGLE LINE PER SEQUENCE ROW

# example
# linecount.sh input name output 0 1

myfile=$1		# input file
filter=$2		# filtering_program_name
outfile=$3		# output file
format=$4		# 1 if fastq, 2 if fasta, 0 otherwise
concat=$5		# 1 if concat to output file, 0 otherwise

# if fastq, divide by 4 to get read counts; if fasta, divide by 4 to get read counts; 
num=$( cat ${myfile} | wc -l | if [ ${format} == 1 ]; then awk '{print $0/4}'; elif [ ${format} == 2 ]; then awk '{print $0/2}'; else cat; fi; )

echo -e ${filter}"\t"${num} | if [ ${concat} == 1 ]; then cat >> ${outfile}; else cat > ${outfile}; fi
