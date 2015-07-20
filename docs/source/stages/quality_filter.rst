==============
quality_filter
==============

Run cutadapt and possibly prinseq on input reads

#. Run cutadapt on R1 and R2 if exists
#. Count reads after each step
#. Run prinseq if specified
#. Count reads

Configuration Options
=====================

command quality_filter

* cutadapt_options_R1
    Options to pass to cutadapt for R1 read file

    Should not contain the input/output options

    Example: ``-a GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -g GCCGGAGCTCTGCAGATATC -a GATATCTGCAGAGCTCCGGC -m 50 --match-read-wildcards``
* cutadapt_options_R2
    Options to pass to cutadapt for R1 read file

    Should not contain the input/output options

    Example: ``-a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -g GCCGGAGCTCTGCAGATATC -a GATATCTGCAGAGCTCCGGC -m 50 --match-read-wildcards``
* prinseq_options
    Options to pass to prinseq

    Should not include ``-out_good``, ``-out_bad`` or ``fastq``

    Example: ``-min_len 50 -derep 14 -lc_method dust -lc_threshold 3 -trim_ns_left 1 -trim_ns_right 1 -trim_qual_right 15``

.. _quality-filter-output-interpretation:

Output Interpretation
=====================

This stage runs input files from step1 through cutadapt and/or prinseq depending on
what is set in your param.txt

By default, quality_filter will run your original input read files first through
cutadapt which produces the ``R1.cut.fastq`` and ``R2.cut.fastq`` files. These two
files contain all reads that have been trimmed due to low quality or found primers..

After cutadapt runs, prinseq is run(if ``prinseq_options`` exists in param.txt) on 
``R1.cut.fastq`` and ``R2.cut.fastq`` files. This produces ``R1.prinseq.fastq`` and
``R2.prinseq.fastq`` which contain the finalized trimmed and filtered read files
that will go on to the next stages.

quality_filter.R1 and quality_filter.R2
---------------------------------------

symlink to either the cutadapt or prinseq output files dending on how param.txt is
configured.

All output files
================

* quality_filter.R1, quality_filter.R2
    Symlink to resulting trimmed reads file
* R1.count, R2.count
    Counts for every stage
* R2.discard, R1.discard
    Read ids that were discarded
* R1.cut.fastq, R2.cut.fastq
    Resulting read files
* R1.prinseq.fastq, R2.prinseq.fastq
  Reads passed prinseq filter
* R1.prinseq.bad.fastq, R2.prinseq.bad.fastq
  Reads that failed pinseq filter
