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

Count Files
-----------

Each stage should contain a ``.count`` file for each of the input reads. These
count files contain the number of input reads that were present as well as how
many reads there were after each item listed.

Example iterative_blast_phylo_1 ``contig.count`` file::
    
    input   304
    megablast   176
    dc-megablast    117

You can see that the input read file(``1.contig.fasta``) contained 304 reads. This
number should match the number of contigs that were generated in the ray2_assembly_1
stage.

We can see that it does here::

    ray_contigs 309
    cap_contigs 304

If we want to do another check, we can look at the iterative_blast_phylo_2 stage
where we would expect to see the input counts for R1 and R2 match the 
unassembled_reads count in the ray2_assembly_1/R1.count and ray2_assembly_1/R2.count

ray2_assembly_1/R1.count::

    input   46474
    unassembled_reads   27779

ray2_assembly_1/R1.count::

    input   46474
    unassembled_reads   29091

iterative_blast_phylo_2/R1.count::

    input   1000
    megablast   95
    dc-megablast    88

iterative_blast_phylo_2/R1.count::

    input   1000
    megablast   115
    dc-megablast    98

You can see here that the input counts do not match. Input was only 1000, for
R1 and R2. Why? Because the pipeline only uses the first 1000 unassembled reads
from each of R1 and R2 otherwise the unassembled blast would take too long to 
complete.
