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
