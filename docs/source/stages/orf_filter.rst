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

Interpreting Output
===================

orf_filter is intended to extract only the reads/contigs in the supplied input
file that contain orfs.

orf_filter.R1, orf_filter.R2, orf_filter.contig
-----------------------------------------------

Resulting reads that only contain orfs

R1.orfout.fa, R2.orfout.fa and contig.orfout.fa
-----------------------------------------------

This is the outputted fasta sequence from the getorf program
