#!/bin/bash

# blast in chunks as separate processes

# notes: 
# IMPORTANT! inputfasta must be in format ID \n seq --- that is to say, the sequences must be on single lines (so it will split properly)
# This script can only handle an input file of with a maximum number of lines ${chunk}*1000
# ${chunk} should be even (so it will split properly)

# e.g., 
# block_blast.sh mydir /my/input/file /data/db/nt blastn megablast 1000 25 /my/dir/blastout /my/dir/headerout

tmpdir=$( readlink -m $1 )			# tmp dir (i.e., output dir)
inputfasta=$( readlink -m $2 )		# input file 
db=$3								# blast db
blast_type=$4						# type of blast (options are: blastn blastx)
task=$5								# type of blast algorithm (options are: megablast dc-megablast blastn)
blast_options=$6					# blast options
chunk=$7							# number of lines in files to blast (note - IMPORTANT - must be even for fasta file!)
ninst=$8							# number of instances of BLAST you want to run in parallel
outfile=$9							# output file
outheader=${10}						# output header

# suggested:
# task=megablast	
# chunk=1000
# ninst=25;

# scripts directory
d=$( dirname $( readlink -m $0 ) )

if ! cd ${tmpdir}; then
 
	echo "[error] changing directories"; 

else

	# make sure file is in format ID \n seq, so it will split properly
	# cat ${inputfasta} | fastajoinlines > tmp.fa
	
	# n is the number of files which will be made during the split (round up to get last file which many have less than ${chunk} number of lines)
	# note: awk's int function is actually a floor function. it doesnt round
	len=$( cat ${inputfasta} | wc -l )
	n=$( echo ${len} | awk -v c=${chunk} '{print int(0.999+$1/c)}' )
	
	echo "splitfiles="${n}	
	echo "wc_input_length="${len}	
	echo "wc_split_length="${chunk}
	echo "instances="${ninst}
	
	# make sure the header matches the BLAST output format you choose below
	echo "qseqid	sseqid	pident	length	mismatch	gapopen	qstart	qend	sstart	send	evalue	bitscore	qlen	slen" > $outheader
	
	# KEEP IN MIND, MAX ALLOWED LINES IN FILE LIMITED TO TO 1000*${chunk}
	if [ ${len} -gt $((1000*${chunk})) ]; then 
		
		echo "[error] file size "${len}" too big. The max allowed number of rows is 1000 x "${chunk}"."

	else

		# -a 3 => 1000 files max
		split -a 3 -d -l ${chunk} ${inputfasta} tmpsplit

		c=0; 							# counter for files
		
		# while counter < #files
		while [ ${c} -lt ${n} ]; do
		
			pids=""						# PIDs (this variable will hold a concatenated string of ${ninst} number of pids) 

			c2=$((${c}+${ninst})); 		# submit ${ninst} number of processes
			echo ${c}"-"$((${c2}-1));

			# while counter < #files && c < c+ninst, submit jobs
			while [ ${c} -lt ${n} -a ${c} -lt ${c2} ]; do 
			{
				# echo "**"${c}"**"; 
				# make 3 digits
				suffix=$( padtowidth=3; printf "%0*d" $padtowidth $c; )

				# blast
				${d}/blast_wrapper.pl --type ${blast_type} --query tmpsplit${suffix} --db ${db} --task ${task} --out blastout${suffix} --options "${blast_options}" &
				pids=${!}" "${pids}

				c=$((1+${c}));			# increment c
			}
			done; 
			
			echo "wait on PIDs: "${pids}
			wait ${pids}
						
		done 
	fi
	
	cat blastout* > $outfile

fi