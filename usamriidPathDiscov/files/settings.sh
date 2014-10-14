######## This file contains the global variables which will be available to the script after this file is sourced. It should be changed to suit your particular setup. After it is sourced, all dependencies should be taken care of.

# Point this to the directory just below where you have usamriidPathDiscov  installed
APPS=/media/VD_Research/People/Dereje.Jima/2014_programming/python/usamriidPathDiscov
# Point this to directory just below where all of your databases are located
DBDIR=/home/AMED/dereje.jima/databases
# open mpi
PATH=/usr/lib64/openmpi/bin:$PATH
PATH=${APPS}/usamriidPathDiscov/download/blast-2.2.28/bin:$PATH
#module load openmpi-x86_64
export INNO_SCRIPTS_PATH=${APPS}/usamriidPathDiscov

######## variables pertaining general properties of the data

export INNO_PHRED_OFFSET=33
export INNO_SEQUENCE_PLATFORM=illumina	# choices are: illumina 454

######## variables pertaining to parallelization

# how many processes (non-SGE) or jobs (SGE) to run in parallel:
export INNO_NODE_NUM=10

######## variables pertaining to reference databases

export INNO_BOWTIE_HUMAN_GENOME_DB=${DBDIR}/humandna/human_dna
export INNO_BOWTIE_HUMAN_TRAN_DB=${DBDIR}/humanrna/h_sapiens_rna

export INNO_BLAST_NT_DB=${DBDIR}/ncbi/blast/nt/nt

######## variables pertaining to taxonomy

export INNO_TAX_NODES=${DBDIR}/ncbi/taxonomy/nodes.dmp
export INNO_TAX_NAMES=${DBDIR}/ncbi/taxonomy/names.dmp

######## variables pertaining to PATH:

#INNO_PATH_PRINSEQ
#INNO_PATH_CUTADAPT
#INNO_PATH_BWA
#INNO_PATH_BLAST
#INNO_PATH_BOWTIE
#INNO_PATH_RAY

# put in front of path and export:

PATH=$APPS:$INNO_SCRIPTS_PATH:$INNO_SCRIPTS_PATH/bin:$INNO_SCRIPTS_PATH/scripts:$INNO_SCRIPTS_PATH/step1:$PATH

######## variables pertaining to perl modules:

# put in front of path and export:

export PERL5LIB=$INNO_PATH_PERL_MODULES:$INNO_SCRIPTS_PATH/Local_Module:$PERL5LIB

export R_LIBS=$R_LIBS:$INNO_SCRIPTS_PATH/scripts
PATH=${APPS}/usamriidPathDiscov/bin/usamriidPathDiscov_cli:$PATH
export LD_LIBRARY_PATH=/usr/lib64/openmpi/lib:$LD_LIBRARY_PATH
