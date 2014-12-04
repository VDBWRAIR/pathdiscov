==========
orf_filter
==========

#. Run getorf on each input reads
    #. Run getorf
    #. Generate unique hash/set of the first part of read names from getorf output
    #. Extract only reads that match the unique hash/set from input read file

Configuration Options
=====================

command orf_filter

* getorf_options
    Options to pass to getorf except ``-sequence`` and ``-outseq``

    Example: ``-minsize 60 -find 0``

Output
======

* out.cap.fa, 1.R1.unmap.fastq, 1.R2.unmap.fastq 
    Symlink to original input reads
* R1.orfout.fa, R2.orfout.fa
    Output from getorf
* orf_filter.R1, orf_filter.R2
    Resulting reads that have passed filter
