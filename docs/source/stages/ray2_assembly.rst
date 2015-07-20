============
ray_assembly
============

#. Break R1 and R2 (if exists) into paired and unpaired fastq files
#. Assemble paired and unpaired files together using Ray2
#. Runs the outputted fasta file from Ray2(``out.ray.cap``) through Cap3 to try 
   and assemble the contigs better
#. Build bowtie index on outputted cap3 concatenated singlet + contig
#. Map original paired and unpaired reads to contig index 
#. Use mapping to find unmapped reads

Configuration Options
=====================

command ray2_assembly

* kmer
    Example: ``25``
* ninst
    Number of mpiexec instances to use.

    You can specify NUMINST to have it replaced with NODE_NUM from :ref:`config-yaml-base`
    
    Example: ``NUMINST``
* cap
    Boolean to specify to run CAP3 after Ray completes
    
    Choices: ``1`` indicates True, ``0`` indicates False
* cap_options
    Options to pass to CAP3
* map2contigs
    Boolean to specify to map reads back to assembly

    Choices: ``1`` to map reads, ``0`` to skip
* bowtie2_options
    Options for bowtie when map2contigs is ``1``

    Example: ``--local``
* parse_contigs
    If set to ``yes`` then runs parse_contigs on bowtie_mapping/out.sam and creates
    reads_by_contig directory with fastq files for each contig that contain reads
    that mapped to each contig

    Choices: ``yes`` or ``1``, ``no`` or ``0``

.. _ray2-assembly-output-interpretation:

Output Interpretation
=====================

This stage pulls in the reads that were unmapped during the host_map stage and
splits them into Paired and unpaired read files. Ray2 requires that you supply
paired and unpaired read files seapartely, which is why you will notice in the
ray2_assembly directory you will likely have 4 fastq files:

* R1.paired.fastq
* R1.single.fastq
* R2.paired.fastq
* R2.single.fastq

out.ray.fa
----------

Once the files are split, Ray2 is run and directed to place all output inside of the
``ray2_assembly/results`` directory and all Ray output is placed inside of the
``ray2_assembly/logs_assembly`` directory. Ray produces quite a few metric files
that really are not that much of use in the results directory. We are only concered
with the Contigs.fasta file which contains the fasta contigs built by Ray.
This file is formatted such that each sequence may span multiple lines, so the
pipeline reformats it and places each contig sequence onto its own line which is
where the ``ray2_assembly/out.ray.fa`` file comes from.

out.cap.fa
----------

If your param.txt has ``cap 1`` (Which it is by default) then the pipeline will
run ``out.ray.fa`` through the cap3 program to try and join the Ray contigs together
to form larger contigs. cap3 produces quite a few files all prefixed with
``out.ray.fa.cap``. The important files from cap3 are ``.singlets`` and ``.contigs``
``out.ray.fa.cap.singlets`` contains all the ``out.ray.fa`` contigs that cap3 could
not form into larger contigs, where ``out.ray.fa.cap.contigs`` contains all the
larger contigs cap3 was able to form. Thus, ``out.cap.fa`` is simply the ``.singlets``
file concattenated with the ``.contigs`` file.

ray2_assembly_1.fasta
---------------------

This is simply a symbolic link that points to either ``out.ray.fa`` or ``out.cap.fa``

1.R1.unmap.fastq and 1.R2.unmap.fastq
-------------------------------------

If your param.txt has ``map2contigs 1`` (Which it is by default) then the pipeline
will use bowtie2 to build an index from ``ray2_assembly_1.fasta``. This is where the 
``bowtie_index`` directory comes from. Once the index is build it maps all input 
reads(unmapped reads from host_map same input as Ray2) onto this index. This produces
the ``bowtie_mapping`` directory. The end result is a bam file containing the 
assembly. This file is then used to extract all reads that mapped to the index 
which produces ``1.R1.unmap.fastq`` and ``1.R2.unmap.fastq``.

reads_by_contig
---------------

If your paramt.xt has both ``map2contigs 1`` and ``parse_contigs yes`` then the
pipeline will utilize the ``bowtie_mapping/out.bam`` assembly to parse out each
contig's reads and place them into ``reads_by_contig``. Each file in this directory
is named after the contig inside of the ``ray2_assembly_1.fasta``

All output files
================

* 1.R1.unmap.fastq, 1.R2.unmap.fastq
    Unmapped reads found from bowtie mapping
* cap3.out
    cap3 output stats
* bowtie2_index
    Directory that contains bowtie index built from cap3 contigs
* bowtie2_mapping
    Directory that contains all bowtie mapping results
* R1.paired.fastq, R2.paired.fastq
    1-to-1 paired reads
* R1.single.fastq, R2.single.fastq
    Unpaired reads
* ray2_assembly_1.fasta
    Symlink to out.cap.fa
* out.ray.fa
    Ray contigs
* out.cap.fa
    Cap contigs + cap singlets with ids replaced with numbers
* out.ray.fa.cap.concat
    Cap contigs + cap singlets
* out.ray.fa.cap.contigs
    Consensus sequences
* out.ray.fa.cap.singlets
    Reads not used in cap3 assembly
* out.ray.fa.cap.ace
    Ace assembly file
* out.ray.fa.cap.contigs.qual
    Contig quality file
* assembly.count
    Number of contigs+singlets generated by CAP3
* contig.id
    Mapping file of new contig number to original Ray2/CAP3 id
* contig_len.txt
    Length of each contig
* contig_numreads.txt
    How many reads support each contig
* R1.count, R2.count
    Count of resulting R1 and R2 unmapped reads
* results
    Directory containing Ray2 assembly
* reads_by_contig
    Directory containing fastq files for each contig that contain reads that
    mapped to each contig 

Undocumented Ouput
------------------

* out.ray.fa.cap.contigs.links
* head.1.R1.unmap.fastq
* head.1.R2.unmap.fastq
* logs
* logs_assembly
* out.ray.fa.cap.info
