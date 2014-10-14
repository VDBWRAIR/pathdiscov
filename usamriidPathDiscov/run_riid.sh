#!/bin/bash

THIS=$(cd $(dirname $0) && pwd)

###### TODO ###########
# Remove hard coded paths and put in settings.sh

# Setup Pipeline
source ${THIS}/../files/settings.sh

##### GLOBALS #####

# EIDRU Home
EIDRU=/media/VD_Research/People/Dereje.Jima/data
test -d $EIDRU || echo "EIDRU directory $EIDRU doesn't exist"

# Root of all the NGS data
NGSDATA=${EIDRU}/Rickettsia
test -d $NGSDATA || echo "NGSData directory $NGSData doesn't exist"
####### END GLOBALS ######

# Get samplename from parameter?
samplename=$1
if [ "$samplename" == "" ] 
then
	echo "Please supply a samplename"
	exit 1
fi
echo "Samplename: $samplename"

# Setup analysisdir for this analysis(current directory)
analysisdir=${samplename}
echo "Analysis Directory: $analysisdir"
# Create analysis directory
if [ -d $analysisdir ]
then
    rm -rf ${analysisdir}.bk 2>/dev/null
    mv ${analysisdir} ${analysisdir}.bk
fi
# Create analysis dir with input and results folder
mkdir -p ${analysisdir}/{input,results}
# Get the absolute path for analysis dir
analysisdir=$(cd $analysisdir && pwd)

# Where are the reads located
samplereadsdir=${NGSDATA}/combineRawFastq/${samplename} 
echo "Looking for reads in $samplereadsdir"

# Make sure reads dir exists and has reads
if [ ! -d $samplereadsdir ] || [ $(ls $samplereadsdir | wc -l) -eq 0 ]
then
    echo "No reads for ${samplename} in $samplereadsdir"
    rm -rf ${analysisdir}
    exit 1
fi

####### Run RIID #######
# Drop env to a file for logging
env > ${analysisdir}/env.txt

# Setup params and reads
cd ${analysisdir}/input
# Get all reads and combine into F.fastq and R.fastq
echo "Combining and converting all read files"
combine_reads.py ${NGSDATA}/combineRawFastq/${samplename}/
echo "Finished combining read files"
# Setup params.txt file
pathogen.pl --example > param.txt

# Run RIID pipeline
echo "RIID Pipeline started -- $(date)"
run_standard.pl --sample $samplename --outputdir ../results --R1 F.fastq --R2 R.fastq | tee ${analysisdir}/analysis.log
echo "RIID Pipeline finished -- $(date)"

# Make some convienience symlinks
cd ..
ln -s results/iterative_blast_phylo_2/reports unassembledread_reports 
ln -s results/iterative_blast_phylo_1/reports contig_reports

# Do some graphics
# Requires R
if [ "$(rpm -qa | grep R-3.0)" != "" ]
then
	# Graphs read counts at each stage
	${THIS}/graph_counts.sh results ./ R1.count R2.count ${samplename} "step1 quality_filter host_map_1 ray_assembly_1 iterative_blast_phylo_1 iterative_blast_phylo_2"

	# Graphs how long it took to run each stage
	${THIS}/graph_times.sh ${analysisdir}/analysis.log . times ${samplename}
fi
