=====================
iterative_blast_phylo
=====================

#. For each mate supplied
    #. Symlink input reads and/or convert fastq->fasta
    #. Count input lines
    #. For each blast_db_list
        #. Blast in chunks defined by ninst_list
        #. Generate x.mate.blast.phylo for blast report
        #. Pull out top blast results only
        #. Then annotate top blast result
        #. Get all reads that didn't blast
        #. Get superclass counts file
#. Create reports for all blast_task_list items

Configuration Options
=====================

command iterative_blast_phylo

* blast_db_list
    blast db prefix

    Example: ``BLAST_NT,BLAST_NT``
* blast_task_list
    Defines the blast tasks that will be run. Must match blast_db_list in size

    Choices: ``megablast``, ``dc-megablast``, ``blastn``, ``blastx``

    Example: megablast,dc-megablast
* blast_options_list
    Blast options except for: ``-task -query -db -out -outfmt -num_descriptions``
    
    Example: ``-evalue 1e-4 -word_size 28,-evalue 1e-4 -word_size 12,-evalue 1e-4``
* ninst_list
    The input file will be broken into chunks and blasted in parallel - this parameter is the number of instances of BLAST you want to run in parallel
    If NUMINST is specified then NODE_NUM from :ref:`config-yaml-base` will be used to replace it during install
    
    Must be the same size as blast_task_list

    Example: ``NUMINST,NUMINST``
* taxonomy_names
    NCBI taxonomy names dump file

    If TAX_NAMES is specified then tax_names from :ref:`config-yaml-base` will be used to replace it during install

    Example: ``TAX_NAMES``
* taxonomy_nodes
    NCBI taxonomy nodes dump file

    If TAX_NODES is specified then tax_nodes from :ref:`config-yaml-base` will be used to replace it during install

    Example: ``TAX_NODES``

Output
======

* 1.contig.fasta or 1.R1.fasta,1.R2.fasta
    Symlink to input reads
* contig.count
    Counts for all stages
* tmp_contig_1, tmp_contig_2
    Contains files that were split during initial blasting
* x.contig.blast
    Blast results for each blast_task_list where x is the index of the blast_task in the list
* x.contig.top.blast
    Only the top result for each read
* x.contig.blast.phylo
    Blast report, with counts for each taxid
* x.contig.top.blast.phylo
    Top blast results, with counts for each taxid
* x.contig.noblast.fasta
    Reads that didn't blast for an iteration
* 2.contig.fasta, or 2.R1.fasta,2.R2.fasta
    Symlink to 1.mate.noblast.fasta
* iterative_blast_phylo_1.contig
    Symlink to final noblast.fasta
* contig.count.superclass
    Superclass count file from blast.phylo
* contig.top.count.superclass
    Superclass count from top.blast.phylo
* reports/
    Contains all reports from all blast_task_list joined together

    * x.contig.top.blast
    * contig.samplename.phylo.txt
        x.mate.blast.phylo joined
    * contig.samplename.top.phylo.txt
        x.mate.top.blast.phylo joined
    * contig.samplename.top.report.txt
        x.mate.top.blast joined
    * contig.samplename.top.smallreport.txt
        sequence columns removed from report
