#!/bin/bash

# wrapper for the following scripts:

# annotate_blast.sh
# make_phylogeny_pipe.pl
# weighted_count.pl

# to run:
# phylogeny_wrapper.sh . tmp.blast outphy.txt /data/columbia/scratch/ref/taxonomy/taxdump/nodes.dmp /data/columbia/scratch/ref/taxonomy/taxdump/names.dmp

# takes a blast file and annotates it by using the gi id to get taxid + description
# IMPORTANT! assume the SECOND column is subject id of the form, e.g., gi|237900821|gb|FJ968794.1|
# assume BLAST file has NO header

outputdir=$( readlink -m $1 )			# output (or tmp) dir
inputfile=$( readlink -m $2 )			# inputfile (the output of BLAST), a file of the form: query_id, taxid, description
outfile_annotate=$( readlink -m $3 )	# outputfile annotate
outfile_t2q=$( readlink -m $4 )			# outputfile taxid 2 query id
outfile=$( readlink -m $5 )				# outputfile
nodes=$6								# nodes.dmp - e.g., /data/taxonomy/nodes.dmp
names=$7								# names.dmp - e.g., /data/taxonomy/names.dmp
ntdb=$8									# ncbi nt db prefix

# scripts directory
d=$( dirname $( readlink -m $0 ) )

if ! cd ${outputdir}; then
 
	echo "[error] changing directories"; 

else


	# check to make sure the input has an NCBI GI number in the second column ( e.g., should be of the form gi|237900821|gb|FJ968794.1| )
	num=$( head ${inputfile} | awk '$2 ~ /gi\|/' | wc -l )
	
	if [ $num == 0 ]; then
	 
		echo "NCBI GI numbers not found for "${inputfile}; 
		
	else	

	#	cat ${inputfile} | sed '1d' > tmp.nohead.blast
	
		# annotate the blast file	
		# args: outputdir, inputfile, outputfile, nt db
		${d}/annotate_blast.sh . ${inputfile} ${outfile_annotate} ${ntdb}
		
		# make a file that maps taxid to description
		cut -f2- ${outfile_annotate} > tmp.tableconcat

		# for each map each uniq taxid to the number of queries mapping to it
		# args: input, output
		${d}/taxid2queryid.pl ${outfile_annotate} tmp.t2q
	
		# a bit of formatting
		cat tmp.t2q | sort -k3,3nr | cut -f1,2 > ${outfile_t2q}

		# get counts (weighted) for each taxid	
		# assumption: col1=qid, col2=taxid
		cat ${outfile_annotate} | ${d}/weighted_count.pl > tmp.count	
	
		# for each taxid, get info about the kingdom, genus, class, etc
		# pipe in: taxid
		# remove trailing tab
		cat tmp.count | cut -f1 | ${d}/make_phylogeny_pipe.pl ${nodes} ${names} | awk '{if (NR==1) {print > "header_phy"} else {print}}' | sed 's|\t$||' > tmp.phy

		# get descriptions
		${d}/tableconcatlines tmp.count tmp.tableconcat | cut -f3 > tmp.descrip

		# paste counts to info about the kingdom, genus, class, etc	
		paste tmp.count tmp.phy tmp.descrip | sort -k2,2nr | awk 'BEGIN{print "taxid\tcount\tsuperkingdom\tkingdom\tclass\torder\tfamily\tgenus\tspecies\tdescrip"}{print}' > ${outfile}
		
		rm tmp.count tmp.phy tmp.tableconcat tmp.descrip tmp.t2q
	
	fi

fi
