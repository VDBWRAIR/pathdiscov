=====
step1
=====

Handles converting sff -> fastq and renaming read names to easy integers

#. Convert sff -> fastq if SEQUENCE_PLATFORM is set to 454
    Option is configured through `pathdiscov/files/config.yaml.base <../../../../pathdiscov/files/config.yaml.base>`_ during install
#. Unzip gzip files if necessary
#. Change all fastq id's to sequenctial integers
#. Count each mate

Configuration Options
=====================

command step1

* seq_platform
    The sequencing platform used. SEQPLATFORM is replaced with SEQUENCE_PLATFORM from :ref:`config-yaml-base`

    Choices: ``illumina``, ``454``

    Default: ``SEQPLATFORM``

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
