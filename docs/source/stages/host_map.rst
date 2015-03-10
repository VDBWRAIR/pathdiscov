========
host_map
========

Takes a series of indexed databases and runs a series of alignments on them. For example, you might run bowtie2 with db=genome, followed bowtie2 with db=transcriptome. At each stage, the input is what didn't align in the previous stage.

#. Count input files (R1 & R2)
#. Iterate through mapper_db_list in :ref:`sample-param-base` and map input reads against mapper_db_list entries
#. Count reads after each mapping
#. Generate a list of unmapped reads 

Configuration Options
=====================

command host_map

* mapper_program_list
    Comma separated list of mappers to use. By default bowtie2,bowtie2 are used

    Choices: ``bwa``, ``bowtie2``

    Example: ``bowtie2,bowtie2``
* mapper_db_list
    Paths to bowtie indexed genomes as a comma separated list
    These values will be replaced when the pipeline installs by setting human_dna and h_sapiens_rna in the :ref:`config-yaml-base`
    The pipeline looks for HUMAN_DNA and H_SAPIENS_RNA and replaces them during install

    Example: ``HUMAN_DNA,H_SAPIENS_RNA``
* mapper_name_list
    Simply names the mapping assemblies for graphic generation

    Example: ``bowtie2_genome_local,bowtie2_transcript_local``
* mapper_options_list
    Comma separated list of options to pass to the mapper

    Example: ``--local,--local``

Output
======

* host_map_1.R1, host_map_1.R2
    Reads remaining that did not map in any of the host genome mappings that can be used
    further downstream in the pipeline.
* R1.count, R2.count
    Count of unmapped reads at the end of all mappings

    These are the number of reads that did not map to host genomes in mapper_db_list
* R2.discard, R1.discard
    Reads that mapped in all mappings to the host genomes set in mapper_db_list

    This is the number of reads that will be discarded due to mapping to host DBs
* map_1, map_2, ..., map_n
    Each mapping to a host genome in mapping_db_list will generate a new map_X directory
    which will contain all files related to that mapping

   Each map dir has `out.bam` which contians alignment of mapped reads,
   `R1.map.id, R2.map.id`- ids of mapped reads, `R1.unmap.fastq,
   R2.unmap.fastq` -fastq of umapped reads, `R1.unmap.id, R2.unmap.id` -
   ids of unmapped reads.

