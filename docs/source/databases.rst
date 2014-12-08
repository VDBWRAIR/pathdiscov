=========
Databases
=========

The pipeline requires that you have blast databases and host genome indexes available for some of the stages such as :doc:`stages/host_map` and :doc:`stages/iterative_blast_phylo`

Directory Setup
===============

.. code-block:: bash
    
    mkdir -p ~/databases/{humandna,humanrna,ncbi}
    mkdir -p ~/databases/ncbi/blast/{nt}

Blast
=====

In general you just need to unpack the nt/nr databases from ncbi(or wherever) into ~/databases/ncbi/blast/nt,nr,taxdb
There is a shell script you should be able to use to do this for you as well.
This may take longer time depending on your network connection.

.. code-block:: bash

    usamriidPathDiscov/scripts/get_blast_dbs.sh ~/databases/ncbi/blast nt taxdb

Host Genome
===========

The host genomes can be built manually(takes a very long time) or they can be downloaded from a few various sources such as ucsc, Ensemble or NCBI.

Download the version of the gonome you like, the current version is GRCh38/hg38. This may change when you download and you may modify the following path accordingly.

The page can be found at http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/

There you can download the host DNA and RNA. Just note that the linke is ~ 700 mb  and can take a few time depending on your connection. Once you downloaded you will need to extract, index it using 'bowtie2-build' and then set the correct path in the :ref:`config-yaml-base`

Example Setup
-------------

Configure pipeline to use NCBI's build38

Ensure you are in the usamriidPathDiscov git cloned directory then execute the following:

.. code-block:: bash

    pushd ~/databases/humandna
    wget  http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromFa.tar.gz
    tar -xzvf  hg38.chromFa.tar.gz
    rm  chroms/*_random.fa
    rm  chroms/*alt.fa
    # NOTE: If you have multiple hosts, you may download the fasta files of all hosts to same folder ('chroms/') and concatinate as show below. You may also modify the names accordingly, exmaple instead of hg38, you may name  'allHost.fa'
    cat chroms/*.fa > hg38_all.fa
    #index the database using bowite2-build
    bowtie2-build hg38_all.fa hg38
    popd
    # replace the location of indexed database in the template config  file  'usamriidPathDiscov/files/config.yaml.base'
    sed -i 's%GENOMEDIR/humandna/human_dna%GENOMEDIR/humandna/hg38%' usamriidPathDiscov/files/config.yaml.base

Download human rna from the same URL, the version of the geome might be different.

.. code-block:: bash
   
   pushd ~/databases/humanrna
   wget  http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/mrna.fa.gz 
   gunzip mrna.fa.gz
   # index the database suing bowtie2-build
   bowtie2-build mrna.fa hg38_mrna
   popd
   # replace the location of indexed database in the template config   file  'usamriidPathDiscov/files/config.yaml.base'
   sed -i 's%GENOMEDIR/humanrna/h_sapiens_rna%GENOMEDIR/humanrna/hg38_mrna%'  usamriidPathDiscov/files/config.yaml.base
   

