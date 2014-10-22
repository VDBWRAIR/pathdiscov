#!/bin/bash

# check quality score bounds for a fastq file ( phred 33 or phred 64? )

cat $1 | head -1000 | perl -ne 'BEGIN{my $nr = 1;}{ if ($nr%4==0) { chomp ($_); my $j; my $inputline=$_; for ($j= 0; $j < length($inputline); $j++) {print ord(substr($inputline,$j,1))."\n"} } $nr++;}' | sort -u -n | awk '{a=$1; if (NR==1) {print "qual bounds"; print $0;}}END{print a}'

