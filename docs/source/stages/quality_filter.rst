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

Output
======

* quality_filter.R1, quality_filter.R2
    Symlink to resulting trimmed reads file
* R1.count, R2.count
    Counts for every stage
* R2.discard, R1.discard
    Read ids that were discarded
* R1.cut.fastq, R2.cut.fastq
    Resulting read files
