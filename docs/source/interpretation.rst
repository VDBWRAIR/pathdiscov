=======================
Interpreting the Output
=======================

While running the pipeline is quite easy once it is installed and all the databases
are setup, interpreting the output can be confusing if you are not sure what to look
for.

Here we will go through an example project that was run and all of the relevant output
files and what they mean.

Pipeline Stages
===============

First you should get aquainted with the stages the pipeline runs.
There are 5 stages available for the pipeline to run in the order they are run:

* :ref:`quality_analysis <quality-analysis-output-interpretation>`
* :ref:`step1 <step1-output-interpretation>`
* :ref:`host_map <host-map-output-interpretation>`
* :ref:`quality_filter <quality-filter-output-interpretation>`
* :ref:`ray2_assembly <ray2-assembly-output-interpretation>`
* :ref:`iterative_blast_phylo_1 <iterative-blast-phylo-output-interpretation>`
* :ref:`iterative_blast_phylo_2 <iterative-blast-phylo-output-interpretation>`

That is, the pipeline will always produce these 6 directories under the ``results``
directory.

General Interpretation
======================

Here we will briefly list each stage and files that are most important. For a more
detailed explenation of the stage and any relevant files generated, see the links
above for each pipeline.

* ray2_assembly_1

    This stage is where all reads that did not get filtered due to host_map or
    quality_filter will be denovo assembled into contigs.

    Relevant files:

        * 1.R1.unmap.fastq
        * 1.R2.unmap.fastq
        * out.cap.fa
        * out.ray.fa

    Relevant directories:

        * reads_by_contig/

* iterative_blast_phylo_1

    This stage is where your ray assembled contigs will get blasted

    Relevant files:

        * reports/contig.*smallreport*.txt

* iterative_blast_phylo_2

    This is where the unassembled reads that passed host-map and quality filtering,
    but did not map to any contigs in the ray assembly

    Relevant files:

        * reports/{R1,R2}.*smallreport*.txt

Stage Output Read File
----------------------

When the pipeline runs, it will run these pipeline stages in the order they are 
listed above. For the most part, each stage will contain at least one symbolic
link that contains the name of the stage. Usually, this is the output file that
will be used by the next stage in the pipeline as input.

.. _count-files-interpretation:

Count Files
-----------

Each stage should contain a ``.count`` file for each of the input reads. These
count files contain the number of input reads that were present as well as how
many reads there were after each item listed.

Here are the results from when the pipeline was run during development for
comparison. If you run the testData your results should be very similar if not 
identical(mostly depending on your database setup)

.. code-block:: bash

    $ for stage in step1 host_map_1 quality_filter ray2_assembly_1 iterative_blast_phylo_1 iterative_blast_phylo_2; do grep -H '.' testoutDir/results/$stage/\*.count; done
    testoutDir/results/step1/R1.count:rawfile   250
    testoutDir/results/step1/R2.count:rawfile   250
    testoutDir/results/host_map_1/R1.count:input    250
    testoutDir/results/host_map_1/R1.count:bowtie2_genome_local 193
    testoutDir/results/host_map_1/R1.count:bowtie2_transcript_local 193
    testoutDir/results/host_map_1/R2.count:input    250
    testoutDir/results/host_map_1/R2.count:bowtie2_genome_local 197
    testoutDir/results/host_map_1/R2.count:bowtie2_transcript_local 197
    testoutDir/results/quality_filter/R1.count:input    193
    testoutDir/results/quality_filter/R1.count:cut_adapt    183
    testoutDir/results/quality_filter/R1.count:prinseq  158
    testoutDir/results/quality_filter/R2.count:input    197
    testoutDir/results/quality_filter/R2.count:cut_adapt    184
    testoutDir/results/quality_filter/R2.count:prinseq  158
    testoutDir/results/ray2_assembly_1/assembly.count:ray_contigs   87
    testoutDir/results/ray2_assembly_1/assembly.count:cap_contigs   87
    testoutDir/results/ray2_assembly_1/R1.count:input   158
    testoutDir/results/ray2_assembly_1/R1.count:unassembled_reads   66
    testoutDir/results/ray2_assembly_1/R2.count:input   158
    testoutDir/results/ray2_assembly_1/R2.count:unassembled_reads   69
    testoutDir/results/iterative_blast_phylo_1/contig.count:input   87
    testoutDir/results/iterative_blast_phylo_1/contig.count:megablast   5
    testoutDir/results/iterative_blast_phylo_1/contig.count:dc-megablast    4
    testoutDir/results/iterative_blast_phylo_2/R1.count:input   66
    testoutDir/results/iterative_blast_phylo_2/R1.count:megablast   4
    testoutDir/results/iterative_blast_phylo_2/R1.count:dc-megablast    3
    testoutDir/results/iterative_blast_phylo_2/R2.count:input   69
    testoutDir/results/iterative_blast_phylo_2/R2.count:megablast   7
    testoutDir/results/iterative_blast_phylo_2/R2.count:dc-megablast    4

How to read this output:

* step1

    * 250 input reads in your F.fastq and R.fastq

* host_map

    * 250+250 reads entered the stage
    * 193+197 reads exited the stage
    
* quality_filter

    * 193+197 reads entered the stage
    * 158+158 reads exited the stage

* ray2_assembly

    * 158+158 reads entered the stage
    * 87 ray contigs were built
    * 0 contigs were extended by cap3 as the count was the same as ray
    * 66+69 reads were left over that were not assembled into contigs

* iterative_blast_phylo_1 -- contig blast

    * 87 contigs entered the stage
    * 4 contigs were left over that did not blast to anything

* iterative_blast_phylo_2 -- unassembled read blast

    * 66+69 reads entered the stage
    * 3+4 reads were left over that did not blast to anything

You will notice that the last count in each stage should match the input line of
the count file in the next stage in the line.
