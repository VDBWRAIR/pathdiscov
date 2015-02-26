=========
Databases
=========

The pipeline requires that you have blast databases and host genome indexes available for some of the stages such as :doc:`stages/host_map` and :doc:`stages/iterative_blast_phylo`

Directory Setup
===============

If you don't have enough space in your home directory for the genome you plan to use, you may download the `databases` anywhere in your network and make a symoblic link to `$HOME/databases`
as follows.

.. code-block:: bash
      
     ln -s /path/to/databases $HOME/databases

.. code-block:: bash
    
    mkdir -p ~/databases/{humandna,humanrna,ncbi}
    mkdir -p ~/databases/ncbi/blast/{nt,taxonomy}

Blast
=====

In general you just need to unpack the nt/nr databases from ncbi(or wherever) into ~/databases/ncbi/blast/nt,nr,taxdb
There is a shell script you should be able to use to do this for you as well.
This may take longer time depending on your network connection.

.. code-block:: bash

    usamriidPathDiscov/scripts/get_blast_dbs.sh ~/databases/ncbi/blast nt taxdb

Taxonomy
========

You need to download and extract the taxonomy databases as well so the pipeline
can extract taxonomy names for each of the blast results

.. code-block:: bash

    pushd ~/databases/ncbi/blast/taxonomy
    wget http://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz -O - | tar xzvf -
    popd

Host Genome
===========

The host genomes can be built manually (takes a very long time) or they can be downloaded from a few various sources such as UCSC, Ensemble or NCBI.

Download the version of the genome you like, the current version is GRCh38/hg38. Note the version of the genome may be updated in the future so you will need to update it periodically.

The page can be found at http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/

There you can download the host DNA and RNA. Just note that the download is ~700 mb and could take a long time depending on your connection. Once you have downloaded the file you will need to extract it and index it using 'bowtie2-build' and then set the correct path in the :ref:`config-yaml-base`

Example Setup
-------------

Configure pipeline to use NCBI's build38

Ensure you are in the usamriidPathDiscov git cloned directory then execute the following:

.. code-block:: bash

    _cwd=$(pwd)
    pushd ~/databases/humandna
    wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromFa.tar.gz
    tar -xzvf hg38.chromFa.tar.gz
    rm chroms/*_random.fa
    rm chroms/*alt.fa
    # NOTE: If you have multiple hosts, you may download the fasta files of all hosts to same folder ('chroms/') and concatinate as show below. You may also modify the names accordingly, exmaple instead of hg38, you may name 'allHost.fa'
    cat chroms/*.fa > hg38_all.fa
    #index the database using bowite2-build
    ${_cwd}/usamriidPathDiscov/download/bowtie2/bowtie2-build hg38_all.fa hg38
    popd
    # replace the location of indexed database in the template config file 'usamriidPathDiscov/files/config.yaml'
    sed -i 's%humandna/human_dna%humandna/hg38%' usamriidPathDiscov/files/config.yaml

Download human rna from the same URL, the version of the geome might be different.

.. code-block:: bash
   
   _cwd=$(pwd)
   pushd ~/databases/humanrna
   wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/mrna.fa.gz 
   gunzip mrna.fa.gz
   # index the database suing bowtie2-build
   ${_cwd}/usamriidPathDiscov/download/bowtie2/bowtie2-build mrna.fa hg38_mrna
   popd
   # replace the location of indexed database in the template config file 'usamriidPathDiscov/files/config.yaml'
   sed -i 's%humanrna/h_sapiens_rna%humanrna/hg38_mrna%' usamriidPathDiscov/files/config.yaml

Verify Databases
================

**Note**: This command is only available after you install. Unfortuneatly at this 
point you cannot use ``verifydatabases`` until after you have finished the entire
installation.

You will probably want to ensure that the pipeline can find all of your databases.
There is now a handy script that you can use to do this prior to installing.

:doc:`scripts/verifydatabases`
