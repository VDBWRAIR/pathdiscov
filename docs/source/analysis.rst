====================
Running the Pipeline
====================

.. _activate:

Activating the Pipeline
=======================

You must always ensure that the virtualenv that you installed into during the
:doc:`install` is activated prior to running any portion of the pipeline.

You only have to do this if you have closed your terminal since the last time your
activated.

You need to specify the full path to the activate script such as

.. code-block:: bash

    $> . /path/to/usamriidPathDiscov/usamriidPathDiscov/bin/activate

When you ran the :doc:`install` the full path to activate should have been printed
during step #2

You can read more about what virtualenv is `here <https://virtualenv.pypa.io/en/latest/>`_

To get help
===========

.. code-block:: bash

    usamriidPathDiscov_cli -h 

If your fastq file has a `.fq` extension, make sure to rename to `.fastq` extension.
The name of the fastq file doesn't matter.

General use
===========

.. code-block:: bash

    usamriidPathDiscov_cli -R1 F.fastq -R2 R.fastq --outdir testoutDir 

Don't forget to give the full path to your forward and reverse files if the reads 
are not in your current analysis directory that you will be running the pipeline in

To run with default **param.txt** file and default host database (human)

.. code-block:: bash

    usamriidPathDiscov_cli -R1 $(pwd)/testData/F.fastq -R2 $(pwd)/testData/R.fastq --outdir testoutDir

**Note**: The pipeline, by default, runs in the following order::

    step1 host_map quality_filter ray2_assembly iterative_blast_phylo orf_filter

Checking error logs
===================

If it fails then an error is reported that generally suggest where it failed by
checking the key files created at each stage. Most likely, the error occurs on the 
suggested stage or the stage before it. You will likely have to check the log files
to get an idea what went wrong and go from there.

To check the log for example under host_map

.. code-block:: bash

    cat testoutDir/results/host_map_1/logs/*.e

Creating custom param.txt
=========================

You need to first initialize a project directory with the input reads and the template
param.txt file.

.. code-block:: bash

    usamriidPathDiscov_cli -R1 $(pwd)/testData/F.fastq -R2 $(pwd)/testData/R.fastq --param --outdir testoutDir

Then, open :ref:`testoutDir/input/param.txt <paramtxt>` and manually edit the databases and 
paramaters you would like to change.

Execute the following line to use the :ref:`paramtxt` you have edited to complete the analysis

.. code-block:: bash

    usamriidPathDiscov_cli -R1 $(pwd)/testData/F.fastq -R2 $(pwd)/testData/R.fastq --noparam --outdir testoutDir
