#!/bin/bash

# e.g., 
# prepare_plot_phylo_percents.sh /my/output/dir 1

# makes a file that looks like this:

# iteration	1
# Viruses	1479.26
# Eukaryota	124.88
# Unannotated	325.84
# Bacteria	437.96
# iteration	2
# Viruses	1739.26
# Eukaryota	391.4
# Unannotated	18.24

outputdir=$( readlink -m $1 )		# output dir
j=$2								# R1 or R2

# scripts directory
d=$( dirname $( readlink -m $0 ) )

if ! cd ${outputdir}; then
 
	echo "[error] changing directories"; 

else

	# regular 
	counter=1; 	
	for i in *.${j}.blast.phylo; do
	 
		echo -e "iteration\t"$counter; 
		
		cat $i | sed '1d;' | awk '{x[$3]+=$2}END{for (y in x) print y"\t"x[y]}' | sed 's|\-|Unannotated|'; 
		
		counter=$(($counter+1)); 
	
	done > ${j}.count.superclass 
 
	# top hits 
	counter=1; 
	for i in *.${j}.top.blast.phylo; do
	 
		echo -e "iteration\t"$counter; 
		
		cat $i | sed '1d;' | awk '{x[$3]+=$2}END{for (y in x) print y"\t"x[y]}' | sed 's|\-|Unannotated|'; 
		
		counter=$(($counter+1)); 
	
	done > ${j}.top.count.superclass

fi