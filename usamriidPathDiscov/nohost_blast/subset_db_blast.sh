#!/bin/bash

# BLAST a subset of your input file to NCBI nt db and make a subset indexed db from the BLAST indexed NT db from tax id s 

# to run:
# subset_db_blast.sh /my/tmp/dir input.fasta subset_db /data/db/nt /data/columbia/scratch/scripts/my_scripts/pathogen_discov_prep/current/gi2taxid.txt

# ASSUME fasta file is in the form: one id per new line one seq per new line
# HARDWIRED take the first 50 seqs

# output BLAST prefix is subset_db

outputdir=$( readlink -m $1 )		# output (or tmp) dir
inputfile=$( readlink -m $2 )		# inputfile fasta file to BLAST
blastdbprefix=$3					# prefix of the blast indexed db you want to create
ntdb=$4								# ncbi nt db prefix
gi2taxid=$( readlink -m $5 )		# gi2taxid.txt file (i.e., a file that has GI in col1, taxid in col2
num=$((2*${6}))						# number of seqs to blast 

# scripts directory
d=$( dirname $( readlink -m $0 ) )

if ! cd ${outputdir}; then
 
	echo "[error] changing directories"; 

else	
	
	# take the first 50 sequences
	cat ${inputfile} | head -${num} > tmp.filehead.fa
	${d}/blast_wrapper.pl --type blastn --query tmp.filehead.fa --db ${ntdb} --task megablast --out tmp.filehead.blast --options '-evalue 1e-4 -word_size 28'
	
	if [ -s tmp.filehead.blast ]; then	

		# take GIs of top hits for 50 qid s
		cat tmp.filehead.blast | awk '{if (!x[$1]) {x[$1]=1; print $2}}' > tmp_gi
		
		# get taxid for these entries (doesnt need to be unique - next step takes care of that)
		blastdbcmd -db ${ntdb} -entry_batch tmp_gi -outfmt '%T' > tmp_taxid	
	
		# get all other GIs with same taxid
		cat ${gi2taxid} | ${d}/get_gi_from_taxid.pl tmp_taxid > tmp_gi_num
		
		blastdb_aliastool -db ${ntdb} -gilist tmp_gi_num -dbtype nucl -out ${blastdbprefix} -title "subset_db"	
	
	# 	rm tmp.filehead.fa tmp.filehead.blast tmp_gi tmp_gi_num tmp_taxid

	else
	
		echo "[error] no blast hits"
	
	fi


fi
