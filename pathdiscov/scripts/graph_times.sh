#!/bin/bash

d=$( dirname $( readlink -m $0 ) )

logfile=$1		# .o log file
outputdir=$2	# output directory
name=$3			# output file name
sample=$4		# sample

# example: 
# d=/ifs/scratch/c2b2/rr_lab/oe2118/projects/pathogen_discovery/scripts/pathogen_discov_pipeline/bin12/scripts
# ${d}/graph_times.sh /ifs/scratch/c2b2/rr_lab/oe2118/projects/201211_mds_ng/rna/input/run1.o1017722 . deltat 201211_mds_ng_rna
	
cat ${logfile} | awk '{if ($1=="[module]") {printf $2"\t"} else if ($1=="[deltat]") {print $2/(60*60)}}' > ${outputdir}/all.deltat

# make a graph of read times
echo "plot delta times";       
# args: inputfile outputfile(.png) name
cmd="R --slave --vanilla --args ${outputdir}/all.deltat ${outputdir}/${name} ${sample} < ${d}/plot_times.r";
echo $cmd;
eval ${cmd}
