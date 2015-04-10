#!/bin/bash

# takes a blast file and annotates it by using the gi id to get taxid + description
# IMPORTANT! assume the SECOND column is subject id of the form, e.g., gi|237900821|gb|FJ968794.1|
# IMPORTANT! assume BLAST file has no header

outputdir=$( readlink -m $1 )		# output (or tmp) dir
inputfile=$( readlink -m $2 )		# inputfile (the output of BLAST)
outfile=$3							# a file of the form: query_id, taxid, description
ntdb=$4								# NCBI nt database prefix

# example: 
# args: outputdir, inputfile, outputfile, nt db
# annotate_blast.sh /my/output/dir ${inputfile} tmp.annotate.blast /data/db/nt

# scripts directory
d=$( dirname $( readlink -m $0 ) )

if ! cd ${outputdir}; then
 
	echo "[error] changing directories"; 

else

	# IMPORTANT! assume the SECOND column is subject id of the form, e.g., gi|237900821|gb|FJ968794.1|
	# IMPORTANT! BLAST file has no header

	# get query_id and subject_id
#	cat ${inputfile} | awk '{print $1 > "tmp_qid"; print $2 > "tmp_gi";}'	
	# description necessary? 
#	blastdbcmd -db ${ntdb} -entry_batch tmp_gi -outfmt '%T	%t' > tmp_taxid

	# for some reason -target_only doesnt work with long form gi|237900821|gb|FJ968794.1|

	# get query_id and subject_id
	cat ${inputfile} | awk '{print $1 > "tmp_qid"; split($2,a,"|"); print a[2] > "tmp_gi";}'

	# description necessary? "-target_only" is crucial! o.w., the output wil be degenerate with respect to gi s 
	blastdbcmd -db ${ntdb} -entry_batch tmp_gi -outfmt '%T	%t' -target_only > tmp_taxid

	
	paste tmp_qid tmp_taxid > ${outfile}
	
	rm tmp_qid tmp_gi tmp_taxid

fi
