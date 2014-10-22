=========================
 Python Project usamriidPathDiscov
=========================


Project Setup
=============

Instructions
------------

# Install System Packages

  ```
  yum install openmpi openmpi-devel git python-devel
  yum groupinstall Development tools
  ```
  
# To install

  ```
  git clone  https://github.com/VDBWRAIR/usamriidPathDiscov.git
  cd  usamriidPathDiscov
  wget --no-check-certificate https://pypi.python.org/packages/source/v/virtualenv/virtualenv-1.11.6.tar.gz -O- | tar xzf -
  python virtualenv-1.11.6/virtualenv.py usamriidPathDiscov
  source usamriidPathDiscov/bin/activate
  pip install paver
  python setup.py install
  deactivate
  source ~/.bashrc
  ```

Using  usamriidPathDiscov
------------------------

To get help::
```
   usamriidPathDiscov_cli   -h 
```
Make sure you have indexed database under your  home directory,

example::
```
   
  ~/databases
  ```

If you extracted the databse to a  different location, make a symbolic link at your home directory::

```
   ln -s  path_to/databases    ~/databases
   ```

That is all it needs, the databases are forced to be at your directory
to make the setting easier.'

If your fastq file has a `.fq` extension, make sure to rename to `.fastq` extension. The name of the fastq file doesn't matters.



To use::

```

   usamriidPathDiscov_cli -R1 F.fastq  -R2 R.fastq  --outdir  testoutDir 

```

Don't forget to give the full path for your forward and reverse lanes if the reads isn't in your analysis directory::
```

    usamriidPathDiscov_cli -R1 ~/testData/F.fastq  -R2 ~/testData/R.fastq  --outdir  testoutDir

```

Using  usamriidPathDiscov_cli
---------

``usamriidPathDiscov_cli`` is a command line computational pipeline for pathogen discovery.The application contain all the necessary tools to do the analyis in one go. It is easy to setup, nothing to edit. Install and use it!!



Supported Python Versions
=========================

This package  supports the following versions out of the box:

* CPython 2.6, 2.7, 3.3
* PyPy 1.9


Licenses
========

The code that makes up this Python project is licensed under the GPL license (fake). Feel free to use it in your free software/open-source or proprietary projects.

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
* Tyghe
