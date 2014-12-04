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

The host genomes can be built manually(takes a very long time) or they can be downloaded from a few various sources such as Illumina, Ensemble or NCBI.

Illumina seems to have the easiest option as they are pre-built for you where the others may require you to build them yourself.

Illumina's iGenomes page can be found at http://support.illumina.com/sequencing/sequencing_software/igenome.html

There you can download any hosts you want. Just note that they are likely very large downloads(20GB+) and can take a few hours. Once they are downloaded you will need to extract them and then set the correct path in the :ref:`config-yaml-base`

Right now the pipeline only allows you to set a single human_dna

Example Setup
-------------

Configure pipeline to use NCBI's build37.2

Ensure you are in the usamriidPathDiscov git cloned directory then execute the following:

.. code-block:: bash

    pushd ~/databases/humandna
    wget ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Homo_sapiens/NCBI/build37.2/Homo_sapiens_NCBI_build37.2.tar.gz
    tar xzvf Homo_sapiens_NCBI_build37.2.tar.gz
    popd
    sed -i 's%GENOMEDIR/humandna/human_dna%GENOMEDIR/humandna/Homo_sapiens/NCBI/build37.2/Sequence/Bowtie2Index/genome%' usamriidPathDiscov/files/config.yaml.base

The h_sapiens_rna needs more documentation, but can simply be removed from the :ref:`config-yaml-base` and then the :ref:`sample-param-base` needs to be modified as well to have the second bowtie mapping removed from the host_map section
