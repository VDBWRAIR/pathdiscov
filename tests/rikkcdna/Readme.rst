Rikkcdna
========

These need to be built in order for tests to run

Small rikketsia database to test with::

    rikkcdna.fasta
    rikkrna.fasta

rikkrna.fasta is just a translated version of rikkcdna.fasta

Building Blast DB's
-------------------

.. code-block:: bash

    cd tests/rikkcdna
    makeblastdb -dbtype nucl -title rikkcdna -in rikkcdna.fasta -out rikkcdna
    makeblastdb -dbtype nucl -title rikkrna -in rikkrna.fasta -out rikkrna

Building Diamond DB
-------------------

.. code-block:: bash

    cd tests/rikkcdna
    diamond makedb -p 12 -d rikknr -v --log --in rikkrna.fasta -b 0.5

Building Host DB
----------------

Just use top 2 reads from both dna and rna as host

.. code-block:: bash

    head -4 rikkcdna.fasta > rikkhostdna.fasta
    head -4 rikkrna.fasta > rikkhostrna.fasta

Build Bowtie
^^^^^^^^^^^^

.. code-block:: bash

    bowtie2-build rikkhostdna.fasta rikkdna
    bowtie2-build rikkhostrna.fasta rikkrna

Build Snap
^^^^^^^^^^

.. code-block:: bash

    snap index rikkhostdna.fasta rikkdna -s 20 -O1000
    snap index rikkhostrna.fasta rikkrna -s 20 -O1000
