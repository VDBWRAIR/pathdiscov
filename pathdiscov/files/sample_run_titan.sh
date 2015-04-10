#!/bin/bash

. ~/titan_setup.sh scratch

d=$1

# note: even though you're using $RealBin as the path to script, you still need export PERL5LIB. why? because SGE copies the script to a different location and, thus, using itself as a reference point, it wont find the modules dir

PERL5LIB=$PERL5LIB:${d}/Local_Module; 
export PERL5LIB

mkdir -p logs

sample=BCN2
r1=/ifs/scratch/c2b2/rr_lab/hk2524/Cancer/201202_ptcl_af/RNAseq/rawfiles/120911_SN828_0145_BC0WE9ACXX_lane3_BCN2_CGATGT_1.fastq.gz
r2=/ifs/scratch/c2b2/rr_lab/hk2524/Cancer/201202_ptcl_af/RNAseq/rawfiles/120911_SN828_0145_BC0WE9ACXX_lane3_BCN2_CGATGT_3.fastq.gz

message=$( qsub -V -N step1.${sample} -e ./logs -o ./logs -l mem=14G,time=4:: -S /bin/sh -cwd -b y $d/pathogen.pl --sample $sample --command step1 --paramfile param.txt --outputdir ../results/$sample --R1 $r1 --R2 $r2 --SGE )
echo $message
h1=$( echo $message | cut -f3 -d" " )

message=$( qsub -hold_jid $h1 -V -N qual.${sample} -e ./logs -o ./logs -l mem=14G,time=9:: -S /bin/sh -cwd -b y $d/pathogen.pl --sample $sample --command step1 quality_filter --paramfile param.txt --outputdir ../results/$sample --SGE )
echo $message
h2=$( echo $message | cut -f3 -d" " )

message=$( qsub -hold_jid $h2 -V -N map.${sample} -e ./logs -o ./logs -l mem=14G,time=12:: -S /bin/sh -cwd -b y $d/pathogen.pl --sample $sample --command step1 quality_filter host_map_2 --paramfile param.txt --outputdir ../results/$sample --SGE )
echo $message
h3=$( echo $message | cut -f3 -d" " )

# for Ray, use MPI parallelization 
message=$( qsub -hold_jid $h3 -V -R y -N assem.${sample} -e ./logs -o ./logs -l mem=2.5G,time=6:: -pe smp 4 -S /bin/sh -cwd -b y $d/pathogen.pl --sample $sample --command step1 quality_filter host_map_2 ray2_assembly --paramfile param.txt --outputdir ../results/$sample --SGE )
echo $message
h4=$( echo $message | cut -f3 -d" " )

message=$( qsub -hold_jid $h4 -V -N itblast.${sample} -e ./logs -o ./logs -l mem=12G,time=10:: -S /bin/sh -cwd -b y $d/pathogen.pl --sample $sample --command step1 quality_filter host_map_2 ray2_assembly iterative_blast_phylo --paramfile param.txt --outputdir ../results/$sample --SGE )
echo $message
h5=$( echo $message | cut -f3 -d" " )



