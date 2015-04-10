#!/bin/bash

d=$( dirname $( readlink -m $0 ) )

projdir=$1		# iterative_blast_phylo directory
outputdir=$2	# output directory
sample=$3		# sample
blast_alg=$4	# sample

# example: 
# dir=/ifs/scratch/c2b2/rr_lab/oe2118/projects/pathogen_discovery/scripts/pathogen_discov_pipeline/bin12/scripts
# ${dir}/graph_supercounts.sh ../iterative_blast_phylo_1 . sample "step1 host_map_1 quality_filter host_map_2"

for j in ${projdir}/*.superclass; do

	echo "**"${j}"**"
	
	name=$( basename $j );

	cat $j | awk '
	{
		if ($1=="iteration") 
		{
			if (NR==1) 
			{
				print "Eukaryota\tBacteria\tViruses\tArchaea\tUnannotated";
			}
			else
			{
				print euk"\t"bac"\t"vir"\t"arch"\t"un;
			} 
			euk=0; bac=0; vir=0; arch=0; un=0;
		}
		else
		{
			if ($1=="Eukaryota") {euk+=$2}
			if ($1=="Bacteria") {bac+=$2}
			if ($1=="Viruses") {vir+=$2}
			if ($1=="Archaea") {arch+=$2}
			if ($1=="Unannotated") {un+=$2}		 		
		}
	}END{print euk"\t"bac"\t"vir"\t"arch"\t"un}' > ${outputdir}/${sample}.${name}.txt

done