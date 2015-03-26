Rikkcdna
========

These need to be built in order for tests to run

Small rikketsia database to test with::

    rikkcdna.fasta

Building Blast DB
-----------------

.. code-block:: bash

    cd tests/rikkcdna
    makeblastdb -dbtype nucl -title rikkcdna -in rikkcdna.fasta

Building Diamond DB
-------------------

.. code-block:: bash

    cd tests/rikkcdna
    diamond makedb -p 12 -d rikknr -v --log --in rikkcdna.fasta -b 0.5
