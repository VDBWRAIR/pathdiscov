=========
Databases
=========

The pipeline requires that you have blast databases and host genome indexes 
available for some of the stages such as :doc:`stages/host_map` 
and :doc:`stages/iterative_blast_phylo`

Prereqs
=======

Most of the commands below may not exist on your system until after you
have completed the normal installation of the pipeline in which the commands 
are placed into your virtualenv's bin directory.

It is fine to complete the installation(``python setup.py install``) prior
to setting up your databases, but just know that ``verifydatabases`` will report
errors until after you have completed setting up your databases.

Directory Setup
===============

You need to first decide where you will be placing all of your databases.
You will need to choose a place that has at least 300 Gigabytes of space.

A common scenario is to create the databases directory on a shared drive
that has a lot of free space and then symlink that shared path to your
home directory as follows:

.. code-block:: bash

    ln -s /path/to/shared/databases ~/databases

Putting the databases on a shared drive will allow others to use your 
databases as well, however, it can slow down your analysis as the
databases will need to be read across the network every time.

The instructions here have you creating the structure
under your home directory inside of a directory called
databases. If you choose a different location you will need to
modify all the instructions below replacing ``~/databases`` with
the path you choose.

You will also need to change the ``pathdiscov/files/config.yaml``
file to point to that location as well. It by default also points to
~/databases.

Create databases directory structure
------------------------------------

.. code-block:: bash
    
    mkdir -p ~/databases/{humandna,humanrna,ncbi,diamond}
    mkdir -p ~/databases/ncbi/blast/{nt,nr}
    mkdir -p ~/databases/ncbi/taxonomy

Blast
=====

Blast databases correspond to :doc:`stages/iterative_blast_phylo`'s 
``blast_db_list`` and ``blast_pro_db``

In general you just need to unpack the nt/nr databases from ncbi(or wherever) 
into ~/databases/ncbi/blast/nt,nr

There is a shell script included that you can use to do this for you.
This may take a long time depending on your network connection.

.. code-block:: bash

    pathdiscov/scripts/get_blast_dbs.sh ~/databases/ncbi/blast nt nr taxdb

Taxonomy
========

Taxonomy databases correspond to :doc:`stages/iterative_blast_phylo`'s 
``taxonomy_nodes`` and ``taxonomy_names``

You need to download and extract the taxonomy databases as well so the pipeline
can extract taxonomy names for each of the blast results

.. code-block:: bash

    pushd ~/databases/ncbi/taxonomy
    wget http://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz -O - | tar xzvf -
    popd

Diamond
=======

Diamond corresponds to :doc:`stages/iterative_blast_phylo`'s 
``blast_db_list``

Download and index protein database for diamond blastx.

You should make sure to view the `makedb`_ command, specifically the part about
how the ``-b`` option controls how much memory diamond will take when it runs
using the database you build. In the example below we use ``-b 2`` which should
consume about 12GB of RAM so if you don't have that much RAM on your computer where
you will run diamond, or if you have way more RAM you may want to change that
number as the higher it is, the faster diamond will run but consume more memory.

You can use the ``free -hg`` command which will show you how many GB you have(listed
under the Total column)

.. code-block:: bash
      
    pushd ~/databases/diamond
    wget ftp://ftp.ncbi.nih.gov/blast/db/FASTA/nr.gz
    gunzip nr.gz
    diamond makedb -d diamondnr -v --log --in nr -b 2
    popd

Alternatively you can generate the diamond database from an already downloaded
blast nr database

.. code-block:: bash

    pushd ~/databases/diamond
    blastdbcmd -db ~/databases/ncbi/blast/nr/nr -entry all > blastnr.fasta
    diamond makedb -d diamondnr --log --in blastnr.fasta -b 2
    rm blastnr.fasta

Host Genome Setup
=================

The host genome setup corresponds to the :doc:`stages/host_map`'s
``mapper_db_list``

General steps to build host genome
----------------------------------

#. Download
#. Unpack download
#. build index

Links to different genome sites to download from
------------------------------------------------

* UCSC
    http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/
* Ensemble
    http://www.ensembl.org/info/data/ftp/index.html
* NCBI
    ftp://ftp.ncbi.nih.gov/genomes/

Building the Genome Indexes
---------------------------

The instructions below default to downloading and building the Human Genome
DNA and RNA databases.

If you want to build different host genomes you can download the fasta file from
one of the sources listed above and index them using the steps below
(replacing the hg38 fasta file with the path to the fasta file you download).

Ensure you are in the pathdiscov git cloned directory then proceed.

DNA
^^^

#. Download and unpack

    .. code-block:: bash

        pushd ~/databases/humandna
        wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromFa.tar.gz
        tar -xzvf hg38.chromFa.tar.gz

#. Concatenate all host fasta [Optional]

    If you have multiple hosts, you may download the fasta files of all 
    hosts to same folder ('chroms/') and concatenate as show below.
    You may also modify the names accordingly, example instead of hg38, you may 
    name 'allHost.fa'

    .. code-block:: bash

        cat chroms/\*.fa > hg38_all.fa

#. Index the downloaded fasta

    * Bowtie

        .. code-block:: bash

            bowtie2-build hg38_all.fa hg38

    * Snap

        .. code-block:: bash

            snap index hg38_all.fa hg38 -s 20

#. Clean up download

    .. code-block:: bash

        rm chroms/\*_random.fa
        rm chroms/\*alt.fa
        rm -rf chroms
        rm hg38.chromFa.tar.gz
        rm hg38_all.fa

#. Setup config.yaml to utilize indexed database

    .. code-block:: bash

        popd
        sed -i 's%host_dna:.*%host_dna: humandna/hg38%' pathdiscov/files/config.yaml

RNA
^^^

Download human rna from the same URL, the version of the genome might be different.

#. Download and unpack

    .. code-block:: bash
       
        pushd ~/databases/humanrna
        wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/mrna.fa.gz
        gunzip mrna.fa.gz

#. Index the downloaded fasta

    * Bowtie

        .. code-block:: bash

            bowtie2-build mrna.fa hg38_mrna

    * Snap

        .. code-block:: bash

            snap index mrna.fa hg38_mrna -s 20

#. Cleanu up download

    .. code-block:: bash

        rm mrna.fa

#. Setup config.yaml to utilize indexed database

    .. code-block:: bash

        popd
        sed -i 's%host_rna:.*%host_rna: humanrna/hg38_mrna%' pathdiscov/files/config.yaml

Verify Databases
================

Note: This command is only available after you install. Unfortunately at this point you cannot use verifydatabases until after you have finished the entire installation.

You will probably want to ensure that the pipeline can find all of your databases. There is now a handy script that you can use to do this prior to installing.

:doc:`scripts/verifydatabases`

.. _makedb: https://github.com/bbuchfink/diamond#makedb-options
