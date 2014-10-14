#!/bin/bash

d=$( dirname $( readlink -m $0 ) )

projdir=$1		# project directory
outputdir=$2	# output directory
r1name=$3		# r1 name
r2name=$4		# r2 name
sample=$5		# sample
commands=$6		# list of commands

# example: 
# dir=/ifs/scratch/c2b2/rr_lab/oe2118/projects/pathogen_discovery/scripts/pathogen_discov_pipeline/bin11/scripts
# ${dir}/graph_counts.sh ../results_1/201211_mds_ng_rna . R1.count R2.count 201211_mds_ng_rna "step1 host_map_1 quality_filter host_map_2"

echo "commands: "${commands}

rm -rf ${outputdir}/${r1name}
rm -rf ${outputdir}/${r2name}

# modify commands according to your order: 
for i in ${commands}; do 

	echo "***"${i}"***";

	for j in ${projdir}/${i}/*.count; do 

		if ! [[ $j =~ R2 ]]; then 
			cat $j | grep -v input | sed 's|_||g; s|genome|gen|; s|transcriptome|tran|; s|transcript|tran|' >> ${outputdir}/${r1name}
		fi; 

		if ! [[ $j =~ R1 ]]; then 
			cat $j | grep -v input | sed 's|_||g; s|genome|gen|; s|transcriptome|tran|; s|transcript|tran|' >> ${outputdir}/${r2name}
		fi; 

	done

#	cat ${projdir}/${i}/*.count | grep -v input >> ${outputdir}/${r1name}
#	cat ${projdir}/${i}/*.count | grep -v input >> ${outputdir}/${r2name}

done

# make a graph of read counts
echo "plot sequence counts R1";       
# args: inputfile outputfile(.png) name
cmd="R --slave --vanilla --args ${outputdir}/${r1name} ${outputdir}/${r1name} ${sample}_R1 < ${d}/plot_counts.r";
echo $cmd;
eval ${cmd}

# make a graph of read counts
echo "plot sequence counts R2";       
# args: inputfile outputfile(.png) name
cmd="R --slave --vanilla --args ${outputdir}/${r2name} ${outputdir}/${r2name} ${sample}_R2 < ${d}/plot_counts.r";
echo $cmd;
eval ${cmd}
