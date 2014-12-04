=====
step1
=====

Handles converting sff -> fastq and renaming read names to easy integers

#. Convert sff -> fastq if seq platform is set to 454
#. Unzip gzip files if necessary
#. Change all fastq id's to sequenctial integers

Output
======

* R1.fastq, R2.fastq
    Same as input reads except with ID's replaced
* R1.count, R2.count
    Count of input reads for each input file
* R1.id, R2.id
    Map of old id's to numeric ids
* step1.R1, step1.R2
    Symlink that points to output files(R1.fastq, R2.fastq)
