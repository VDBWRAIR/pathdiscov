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

    $> . /path/to/pathdiscov/pathdiscov/bin/activate

When you ran the :doc:`install` the full path to activate should have been printed
during step #2

You can read more about what virtualenv is `here <https://virtualenv.pypa.io/en/latest/>`_

To get help
===========

.. code-block:: bash

    pathdiscov_cli -h 

If your fastq file has a `.fq` extension, make sure to rename to `.fastq` extension.
The name of the fastq file doesn't matter.

Usage Examples
==============

You may use ``.fastq``, ``.fastq.gz``, ``.sff`` and ``.sff.gz`` files for ``-R1`` and ``-R2``

Typical Usage
-------------

To run with default **param.txt** file created during ``python setup.py``

.. code-block:: bash

    pathdiscov_cli -R1 F.fastq -R2 R.fastq --outdir testoutDir 

Single Read file
----------------

   .. code-block:: bash
           
       pathdiscov_cli -R1 testData/454Reads.sff --outdir testoutDir


Creating custom param.txt
-------------------------

You need to first initialize a project directory with the input reads and the template
param.txt file.

.. code-block:: bash

    pathdiscov_cli -R1 testData/F.fastq -R2 testData/R.fastq --param --outdir testoutDir

Then, open :ref:`testoutDir/input/param.txt <paramtxt>` and manually edit the databases and 
paramaters you would like to change.

Execute the following line to use the :ref:`paramtxt` you have edited to complete the analysis

.. code-block:: bash

    pathdiscov_cli -R1 testData/F.fastq -R2 testData/R.fastq --noparam --outdir testoutDir

Using Sun Grid Engine
---------------------
    
If your cluster support SGE, use ``--use-sge`` to activate sge_iterative_blast_phylo instead of iterative_blast_phylo

.. code-block:: bash

     pathdiscov_cli -R1 testData/F.fastq.gz -R2 testData/R.fastq.gz --outdir testoutDir --use-sge

**Note**: The pipeline, by default, runs in the following order::

    step1 host_map quality_filter ray2_assembly iterative_blast_phylo orf_filter

Specifying number of unassembled reads to blast
-----------------------------------------------

During the Ray assembly stage, all reads that passed quality filter and host_map
are Denovo Assembled using Ray. Then the same reads are mapped back to the
contigs so that the reads that mapped and those that did not map can be 
separated. This is where the unassembled reads are pulled from.

The iterative_blast_phylo_2 step blasts these unassembled reads and generates
reports to tell you what your unassembled reads are composed of. In some cases
you can end up with too many reads to blast(which may take an indefinate amount
of time).

The option ``--blast-unassembled`` will allow you to limit the amount of unassembled
reads that will be sent into the interative_blast_phylo_2 step.

By default, the pipeline uses 1000(configured in config.yaml) but you can
override this as follows:

Use top 100 reads from unassembled reads files.

.. code-block:: bash

    pathdiscov_cli --R1 testData/F.fastq --R2 testData/R.fastq --outdir testoutDir --blast-unassembled 100

Checking error logs
===================

If it fails then an error is reported that generally suggest where it failed by
checking the key files created at each stage. Most likely, the error occurs on the 
suggested stage or the stage before it. You will likely have to check the log files
to get an idea what went wrong and go from there.

To check the log for example under host_map

.. code-block:: bash

    cat testoutDir/results/host_map_1/logs/*.e

Running the Test Data Set
=========================

The pipeline comes with some test data that is essentially a very small set of reads
for a rickettsia dataset.

There is also a pre-built rickettsia nt database to make running the test dataset
much faster(2 minutes vs 45 minutes with the full nt database)

Create our custom param.txt file that utilizes rikkcdna nt database
-------------------------------------------------------------------

#. Copy the rikkcdna nt database into your databases

    .. code-block:: bash

        $> mkdir -p ~/databases/ncbi/blast/rikkcdna
        $> cp tests/rikkcdna/rikkcdna.{nhr,nin,nsq} ~/databases/ncbi/blast/rikkcdna/

#. Generate default param.txt so we can modify

    .. code-block:: bash

        $> pathdiscov_cli --outdir testrun --R1 testData/F.fastq --R2 testData/R.fastq --param
        Starting project: testrun

        ________________________________________
        Tasks which will be run:


        Task enters queue = "mkdir('testrun',   'testrun/input',   'testrun/results',   'testrun/logs')   before pathdiscov.main.createPram " 
        Completed Task = "mkdir('testrun',   'testrun/input',   'testrun/results',   'testrun/logs')   before pathdiscov.main.createPram " 
        Task enters queue = 'pathdiscov.main.createPram' 
        Completed Task = 'pathdiscov.main.createPram' 
        Time to complete the task .....0:00:00

    **Note**: We added ``--param`` to the end to have the pipeline generate only
    the param.txt and setup the project directory structure

    At this point there is a param.txt inside of testrun/input that contains all
    the defaults pulled from the config.yaml you edited during the install.
    We can now edit that to point the iterative_blast_phylo stages to point to
    the rikkcdna nt database instead of the entire nt database

#. Customize param.txt for this run

    Here we will use set to replace all occurances of ``nt/nt`` with 
    ``rikkcdna/rikkcdna`` in our param.txt

    .. code-block:: bash

        $> sed -i 's%nt/nt%rikkcdna/rikkcdna%g' testrun/input/param.txt

Run the data with the pipeline
------------------------------

Now that the param.txt is edited we can simply run the pipeline using the 
``--noparam`` argument to tell the pipeline to use the param.txt that is already
there(otherwise it would just move the existing testrun directory to testrun.bk and
create a new testrun project and run with the default param.txt.

.. code-block:: bash

    $> pathdiscov_cli --outdir testrun --R1 testData/F.fastq --R2 testData/R.fastq --noparam
    Starting project: testrun

    ________________________________________
    Tasks which will be run:


    Task enters queue = "mkdir('testrun/results/quality_analysis')   before pathdiscov.main.fastQC " 
    Completed Task = "mkdir('testrun/results/quality_analysis')   before pathdiscov.main.fastQC " 
    Task enters queue = 'pathdiscov.main.fastQC' 
    fastqc /home/AMED/tyghe.vallard/Projects/pathdiscov/testData/F.fastq -o testrun/results/quality_analysis 2>&1 | tee -a testrun/results/quality_analysis/analysis_quality.log
    Completed Task = 'pathdiscov.main.fastQC' 
    Task enters queue = 'pathdiscov.main.priStage' 
    run_standard.pl --sample testrun --paramfile testrun/input/param.txt --outputdir testrun/results --R1 /home/AMED/tyghe.vallard/Projects/pathdiscov/testData/F.fastq --blast_unassembled 1000 2>&1 | tee -a testrun/results/analysis.log
    [cmd] pathogen.pl --sample testrun --command step1 host_map quality_filter ray2_assembly iterative_blast_phylo --paramfile testrun/input/param.txt --outputdir testrun/results --R1 /home/AMED/tyghe.vallard/Projects/pathdiscov/testData/F.fastq  --SGE 0
    [module] step1
    [cmd] step1 --sample testrun --paramfile testrun/input/param.txt --outputdir /home/AMED/tyghe.vallard/Projects/pathdiscov/testrun/results/step1 --logs /home/AMED/tyghe.vallard/Projects/pathdiscov/testrun/results/step1/logs --timestamp 20150723-08.56 --R1 /home/AMED/tyghe.vallard/Projects/pathdiscov/testData/F.fastq --R2 none
    [module] host_map
    [cmd] host_map.pl --sample testrun --paramfile testrun/input/param.txt --outputdir /home/AMED/tyghe.vallard/Projects/pathdiscov/testrun/results/host_map_1 --logs /home/AMED/tyghe.vallard/Projects/pathdiscov/testrun/results/host_map_1/logs --timestamp 20150723-08.56 --R1 /home/AMED/tyghe.vallard/Projects/pathdiscov/testrun/results/step1/step1.R1 --R2 none --fastafile no --wellpaired 1 --run_iteration 1
    [module] quality_filter
    [cmd] quality_filter.pl --sample testrun --paramfile testrun/input/param.txt --outputdir /home/AMED/tyghe.vallard/Projects/pathdiscov/testrun/results/quality_filter --logs /home/AMED/tyghe.vallard/Projects/pathdiscov/testrun/results/quality_filter/logs --timestamp 20150723-08.56 --R1 /home/AMED/tyghe.vallard/Projects/pathdiscov/testrun/results/host_map_1/host_map_1.R1 --R2 none
    [module] ray2_assembly
    [cmd] ray2_assembly.pl --sample testrun --paramfile testrun/input/param.txt --outputdir /home/AMED/tyghe.vallard/Projects/pathdiscov/testrun/results/ray2_assembly_1 --logs /home/AMED/tyghe.vallard/Projects/pathdiscov/testrun/results/ray2_assembly_1/logs --timestamp 20150723-08.56 --R1 /home/AMED/tyghe.vallard/Projects/pathdiscov/testrun/results/quality_filter/quality_filter.R1 --R2 none --fastafile no --run_iteration 1
    [module] iterative_blast_phylo
    [cmd] iterative_blast_phylo.pl --sample testrun --paramfile testrun/input/param.txt --outputdir /home/AMED/tyghe.vallard/Projects/pathdiscov/testrun/results/iterative_blast_phylo_1 --logs /home/AMED/tyghe.vallard/Projects/pathdiscov/testrun/results/iterative_blast_phylo_1/logs --timestamp 20150723-08.56 --R1 /home/AMED/tyghe.vallard/Projects/pathdiscov/testrun/results/ray2_assembly_1/ray2_assembly_1.fasta --fastafile yes --run_iteration 1 --contig 1
    [cmd] cat testrun/results/ray2_assembly_1/1.R1.unmap.fastq | head -4000 > testrun/results/ray2_assembly_1/head.1.R1.unmap.fastq
    [cmd] pathogen.pl --sample testrun --command iterative_blast_phylo_2 --paramfile testrun/input/param.txt --outputdir testrun/results --R1 testrun/results/ray2_assembly_1/head.1.R1.unmap.fastq --R2 none
    [module] iterative_blast_phylo
    [cmd] iterative_blast_phylo.pl --sample testrun --paramfile testrun/input/param.txt --outputdir /home/AMED/tyghe.vallard/Projects/pathdiscov/testrun/results/iterative_blast_phylo_2 --logs /home/AMED/tyghe.vallard/Projects/pathdiscov/testrun/results/iterative_blast_phylo_2/logs --timestamp 20150723-08.58 --R1 /home/AMED/tyghe.vallard/Projects/pathdiscov/testrun/results/ray2_assembly_1/head.1.R1.unmap.fastq --R2 /home/AMED/tyghe.vallard/Projects/pathdiscov/none --fastafile no --run_iteration 2 --contig 0
    [cmd] readcount.pl --sample testrun --outputdir testrun/results/output --projdir testrun/results --dirlist "step1,quality_filter,host_map_1,ray2_assembly_1,iterative_blast_phylo_1,iterative_blast_phylo_2" --trackread
    [cmd] process_counts.pl --sample testrun --outputdir testrun/results/output > testrun/results/output/stats.txt
    [cmd] augment_report.sh testrun/results testrun
    [cmd] pathogen.pl --checkerror --outputdir testrun/results
    Completed Task = 'pathdiscov.main.priStage' 
    Task enters queue = 'pathdiscov.main.symlink' 
    Completed Task = 'pathdiscov.main.symlink' 
    Time to complete the task .....0:02:20

**Note**: The actual command that runs the pipeline and generates the param.txt 
are almost identical we just change the ``--param`` and ``--noparam``

So here the pipeline ran and you can see the status as it ran and what it was doing
as the pipeline lists every command(``[cmd]``) that is run.
