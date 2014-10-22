#!/bin/bash

# convert sff file to fastq

infile=$1								# input file
outfile=$2								# output file

d=$( dirname $( readlink -m $0 ) )		# scripts directory

sffinfo -s ${infile} > all_reads.fna
sffinfo -q ${infile} > all_reads.qual

${d}/fastq_convert.pl all_reads.fna

rm all_reads.fna 
rm all_reads.qual

# output is called "all_reads.fastq" by default

mv all_reads.fastq ${outfile}