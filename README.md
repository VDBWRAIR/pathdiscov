usamriidPathDiscov
==================

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

To use::

```
usamriidPathDiscov_cli -R1 F.fastq  -R2 R.fastq  --outdir  testoutDir 
```

Don't forget to give the full path for your forward and reverse lanes if the reads isn't in your analysis directory::

```
usamriidPathDiscov_cli -R1 $(pwd)/testData/F.fastq  -R2 $(pwd)/testData/R.fastq  --outdir  testoutDir
```

Using usamriidPathDiscov_cli
---------

``usamriidPathDiscov_cli`` is a command line computational pipeline for pathogen discovery.The application contain all the necessary tools to do the analyis in one go. It is easy to setup, nothing to edit. Install and use it!!



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
