=====================
iterative_blast_phylo
=====================

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

1.contig.blast
1.contig.blast.ann
1.contig.blast.phylo
1.contig.blast.t2q
1.contig.fasta
1.contig.noblast.fasta
1.contig.top.blast
1.contig.top.blast.ann
1.contig.top.blast.phylo
1.contig.top.blast.t2q
2.contig.blast
2.contig.fasta
2.contig.noblast.fasta
2.contig.top.blast
2.contig.top.blast.phylo
3.contig.fasta
contig.count
contig.count.superclass
contig.top.count.superclass
iterative_blast_phylo_1.contig
reports/
contig.mock_project.phylo.txt
contig.mock_project.top.phylo.txt
contig.mock_project.top.report.txt
contig.mock_project.top.smallreport.txt
reports 2/
1.contig.top.blast
2.contig.top.blast
contig.mock_project.phylo.txt
contig.mock_project.top.phylo.txt
contig.mock_project.top.report.txt
contig.mock_project.top.smallreport.txt
mock_project.joinR1R2.smallreport.txt
testy.joinR1R2.smallreport.txt
