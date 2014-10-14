#!/bin/bash

# wrapper for the following scripts:

# annotate_blast.sh
# make_phylogeny_pipe.pl
# weighted_count.pl

# to run:
# phylogeny_wrapper.sh . tmp.blast outphy.txt /data/taxonomy/nodes.dmp /data/taxonomy/names.dmp

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
	
		# args: outputdir, inputfile, outputfile, nt db
		${d}/annotate_blast.sh . ${inputfile} ${outfile_annotate} ${ntdb}
		
		# args: input, output
		${d}/taxid2queryid.pl ${outfile_annotate} tmp.t2q
	
		cat tmp.t2q | sort -k3,3nr | cut -f1,2 > ${outfile_t2q}
	
		# assumption: col1=qid, col2=taxid
		cat ${outfile_annotate} | ${d}/weighted_count.pl > tmp.count	
	
		# pipe in: taxid
		cat tmp.count | cut -f1 | ${d}/make_phylogeny_pipe.pl ${nodes} ${names} | awk '{if (NR==1) {print > "header_phy"} else {print}}' > tmp.phy
	
		paste tmp.count tmp.phy | sort -k2,2nr | awk 'BEGIN{print "taxid\tcount\tsuperkingdom\tkingdom\tclass\torder\tfamily\tgenus\tspecies"}{print}' > ${outfile}
		
	#	rm tmp.nohead.blast 
		rm tmp.count tmp.phy
		rm tmp.t2q
	
	fi

fi
