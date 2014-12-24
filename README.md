usamriidPathDiscov
==================

About usamriidPathDiscov_cli
---------

``usamriidPathDiscov_cli`` is a command line computational pipeline for pathogen discovery. The application contains all the necessary tools to do the analyis in one go. It is easy to setup, nothing to edit. Install and use it!!

Clone the repository
--------------------

```
git clone $(eval echo https://$(read -p "Gitub username: " gu; echo $gu)@github.com/VDBWRAIR/usamriidPathDiscov.git)
```
```
cd usamriidPathDiscov
```

Installation
------------

See [Installation](docs/source/install.rst)

Configuration
-------------

See [Configuration](docs/source/configuration.rst)

Using  usamriidPathDiscov
------------------------

To get help::

```
usamriidPathDiscov_cli   -h 
```

If your fastq file has a `.fq` extension, make sure to rename to `.fastq` extension. The name of the fastq file doesn't matters.

General  use::

```
usamriidPathDiscov_cli -R1 F.fastq  -R2 R.fastq  --outdir  testoutDir 
```

Don't forget to give the full path for your forward and reverse files if the reads are not in your current analysis directory that you will be running the pipeline in::

####  To run with  default *`param.txt`* file and default host database (human)::

```
usamriidPathDiscov_cli -R1 $(pwd)/testData/F.fastq  -R2 $(pwd)/testData/R.fastq  --outdir  testoutDir
```
Note: The pipeline runs in the following order: `step1,host_map,quality_filter,ray2_assembly,iterative_blast_phylo,orf_filter`.If it fails the application  generally suggest where it fails  by
checking the key files created at each stage. Most likely, the error occurs on the suggested stage or the stage before it,  hence you may check the log file to get a clue and fix it.

To check the log for example under `host_map`

```
cat testoutDir/results/host_map_1/logs/*.e
```

#### To create  example *`param.txt`* file, edit the parameters and run::

```
usamriidPathDiscov_cli -R1 $(pwd)/testData/F.fastq  -R2 $(pwd)/testData/R.fastq --param  --outdir  testoutDir
```
Then, open `$(pwd)/testoutDir/input/param.txt` and manually edit the databases and paramaters you would like to change.

Execute the following line to use the `param.txt` you have edited to complete the analysis::

```
usamriidPathDiscov_cli -R1 $(pwd)/testData/F.fastq  -R2 $(pwd)/testData/R.fastq --noparam  --outdir  testoutDir
```

Supported Python Versions
=========================

This package  supports the following versions out of the box:

* CPython 2.6, 2.7, 3.3
* PyPy 1.9


Licenses
========

The code that makes up this Python project is licensed under the GPL license. Feel free to use it in your free software/open-source or proprietary projects.

Issues
======

Please report any bugs or requests that you have using the GitHub issue tracker!

Development
===========

Authors
=======

* Micheal Wiley
* Jason Ladner
* Dereje Jima
* Tyghe Vallard
