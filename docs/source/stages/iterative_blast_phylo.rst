=====================
iterative_blast_phylo
=====================

#. For each mate supplied
    #. Symlink input reads and/or convert fastq->fasta
    #. Count input lines
    #. For each blast_db_list
        #. If diamond or blastx then run orf_filter to filter
           input reads prior to blast
        #. Blast in chunks defined by ninst_list
        #. Generate x.mate.blast.phylo for blast report
        #. Pull out top blast results only
        #. Then annotate top blast result
        #. Get all reads that didn't blast
        #. Count reads that didn't blast and add to .count file with name of 
           blast_task_list
        #. Get superclass counts file
        #. Symlink x.contig.noblast.fasta as input to next blast_db_list
#. Create reports for all blast_task_list items

Configuration Options
=====================

command iterative_blast_phylo

* blast_db_list
    Database[s] that are used to blast/diamond against

    The following two values are available for replacement from 
    :ref:`config-yaml-base`

        * ``BLAST_NT`` will be replaced using nt_db
        * ``DIAMOND_NR`` will be replaced using nr_db

    Choices: ``BLAST_NT``, ``DIAMOND_NR``

    Example: ``BLAST_NT,BLAST_NT,DIAMOND_NR``
* blast_task_list
    Defines the blast tasks that will be run. Must match blast_db_list in size

    Choices: ``megablast``, ``dc-megablast``, ``blastn``, ``blastx``, ``diamond``

    Example: ``megablast,dc-megablast``
* blast_options_list
    Blast options except for: ``-task -query -db -out -outfmt -num_descriptions``
    
    Example: ``-evalue 1e-4 -word_size 28,-evalue 1e-4 -word_size 12,-evalue 1e-4``
* ninst_list
    The input file will be broken into chunks and blasted in parallel - this 
    parameter is the number of instances of BLAST you want to run in parallel
    If ``NUMINST`` is specified then ``NODE_NUM`` from 
    :ref:`config-yaml-base` will be used to replace it during install
    
    Must be the same size as blast_task_list

    Example: ``NUMINST,NUMINST``
* taxonomy_names
    NCBI taxonomy names dump file

    If ``TAX_NAMES`` is specified then ``tax_names`` from 
    :ref:`config-yaml-base` will be used to replace it during install

    Example: ``TAX_NAMES``
* taxonomy_nodes
    NCBI taxonomy nodes dump file

    If ``TAX_NODES`` is specified then ``tax_nodes`` from 
    :ref:`config-yaml-base` will be used to replace it during install

    Example: ``TAX_NODES``
* blast_pro_db
    NCBI protein database for taxonomy information lookup with blastdbcmd

    If ``BLASTNR`` is specified then ``nr_db`` from :ref:`config-yaml-base` will 
    be used to replace it during install

    Example: ``BLASTNR``

Output
======

Since iterative_blast_phylo can be run more than once and also it can be
(and usually is) configured to run multiple blast steps, there are quite a few
files generated.

The files will contain some combination of the following as a prefix or
sometimes somewhere in the middle of the file.

The following files and directories use \* where one of the following
will replace it in the actual output name

* ``contig``
* ``R1`` and/or ``R2``

Additionally, anywhere you see ``X`` it will be replaced with
a number representing the blast_task_list item it was generated from.


Output Files/Directories
------------------------

* \*.fasta
    Input reads for the stage. May be a symlink if the input files
    were originally fasta files.
* \*.count
    Counts for all stages with names from blast_task_list
* tmp\_\*_X
    Read file's are divided into ninst numbers and placed into this directory

    * orf_filter
        If diamond or blastx were in blast_task_list then :doc:`orf_filter` is run
        in this directory to filter out reads
* X.\*.blast
    Blast results from each blast task run
* X.\*.top.blast
    Top result from each X.\*.blast file for each read
* X.\*.blast.phylo
    Blast report, with counts for each taxid
* X.\*.top
    Top result from each X.\*.blast.phylo 
* X.\*.noblast.fasta
    Reads that did not blast. The final result will be in the highest
    number. So if you have 3 blast_task_list items, 4.\*.noblast.fasta will
    be the final result.
* X.\*
    Symlink to X.\*.noblast.fasta
* \*.count.superclass
    Superclass count file from blast.phylo
* \*.top.count.superclass
    Superclass count from top.blast.phylo
* iterative_blast_phylo_N.\*
    Symlink to final noblast.fasta file. Represents the resulting
    reads that had no blast results for the entire stage. Would be used
    by next stage in the pipeline as input.
* reports/
    Contains all reports from all blast_task_list joined together

    * \*.samplename.phylo.txt
    * \*.samplename.top.phylo.txt
    * \*.samplename.top.report.txt
    * \*.samplename.top.smallreport.txt
        sequence columns removed from report
