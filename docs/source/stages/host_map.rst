========
host_map
========

Takes a series of indexed databases and runs a series of alignments on them. For example, you might run bowtie2 with db=genome, followed bowtie2 with db=transcriptome. At each stage, the input is what didn't align in the previous stage.

#. Count input files (R1 and/or R2)
#. Iterate through mapper_db_list in :ref:`sample-param-base` and map input 
   reads against ``mapper_db_list`` entries using ``mapper_program_list``
#. Count reads after each mapping
#. Generate a list of unmapped reads 

Configuration Options
=====================

command host_map

* mapper_program_list
    Comma separated list of mappers to use. By default bowtie2,bowtie2 are used

    Choices: ``snap``, ``bowtie2``

    Example: ``bowtie2,snap``

    *NOTE*:
        Selecting Snap takes a lot more memory than bowtie2 to run. For the
        built hg38 genome, approximately 30GB of memory is needed for Snap
        to run. You can estimate memory requirements by checking the size of
        the snap database with:

        .. code-block:: bash

            du -hs /path/to/snap/database

* mapper_db_list
    Paths to bowtie indexed genomes as a comma separated list
    These values will be replaced when the pipeline installs by setting ``host_dna`` and ``host_rna`` in the :ref:`config-yaml-base`
    The pipeline looks for HOST_DNA and HOST_RNA and replaces them during install

    Example: ``HOST_DNA,HOST_RNA``
* mapper_name_list
    Simply names the mapping assemblies for graphic generation

    Example: ``bowtie2_genome_local,snap_transcript_local``
* mapper_options_list
    Comma separated list of options to pass to the mapper

    Example: ``--local -p NUMINST,-t NUMINST``

.. _host-map-output-interpretation:

Output Interpretation
=====================

The host_map stage is fairly simple. For each database and mapping program you list
in the param.txt, a ``map_X`` directory will exist in the ``host_map`` directory for
it.

By default there will be 2 mapping directories, ``map_1`` and ``map_2``. This is
because param.txt has by default ``mapper_program_list bowtie2,bowtie2`` and 
``mapper_db_list /path/to/hg38,/path/to/hg38_mrna``

Mapping Directories
-------------------

``map_1`` relates to the first entry listed in the comma separated list in the 
``mapper_program_list`` inside of param.txt, and ``map_2`` relates to the second
entry and so forth.

Inside of each of the ``map_X`` directories you will find a few files. The only
files of main concern are the ``R1.unmap.fastq`` and ``R2.unmap.fastq`` files which
represent the reads that did not map to the database for that step. The other files
are just supplemental files that were used to extract the unmapped reads from the
outputted bam file from bowtie(``out.bam``)

All output files
================

* host_map_1.R1, host_map_1.R2
    Reads remaining that did not map in any of the host genome mappings that can 
    be used further downstream in the pipeline. May just be a symlink to 
    map_N/R1.unmap.fastq
* R1.count, R2.count
    Count of unmapped reads at the end of all mappings
    Each line represents the count of unmapped reads after each 
    ``mapper_program_list`` has completed
* R2.discard, R1.discard
    These are a list of read ids that were mapped to any of the 
    ``mapper_db_list`` databases
* map_1, map_2, ..., map_n
    Each mapping to a host genome in ``mapping_db_list`` will generate a new 
    map_X directory which will contain all files related to that mapping

    * out.bam
        Contains alignment of mapped/unmapped reads
    * R1.map.id, R2.map.id
        ids of mapped reads
    * R1.unmap.fastq, R2.unmap.fastq
        fastq of umapped reads. These are input for the next mapping.
    * R1.unmap.id, R2.unmap.id
        ids of unmapped reads
