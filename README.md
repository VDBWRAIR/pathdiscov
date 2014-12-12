usamriidPathDiscov
==================

About usamriidPathDiscov_cli
---------

``usamriidPathDiscov_cli`` is a command line computational pipeline for  pathogen discovery.The application contain all the necessary tools to do the analyis in one go. It is easy to setup, nothing to edit. Install and use it!!

Clone the repository
--------------------

```
git clone $(eval echo https://$(read -p "Gitub username: " gu; echo $gu)@github.com/VDBWRAIR/usamriidPathDiscov.git)
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

Don't forget to give the full path for your forward and reverse lanes if the reads isn't in your analysis directory::

####  To run with  default *`param.txt`* file and default host database (human)::

```
usamriidPathDiscov_cli -R1 $(pwd)/testData/F.fastq  -R2 $(pwd)/testData/R.fastq  --outdir  testoutDir
```

#### To create  example *`param.txt`* file, edit the parameters and run::

```
usamriidPathDiscov_cli -R1 $(pwd)/testData/F.fastq  -R2 $(pwd)/testData/R.fastq --param  --outdir  testoutDir
```
Then, open `$(pwd)/testoutDir/input/param.txt` and manually edit  the databases and paramaters you like.

Followed by executing the following line to use the `param.txt` you have edited to complete the analysis::

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

* Mickeal Wiley
* Jason
* Dereje Jima
* Tyghe Vallard
