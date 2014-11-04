=========================
 Python Project usamriidPathDiscov
=========================


Project Setup
=============

Instructions
------------

# Install System Packages

  ```
  yum install openmpi openmpi-devel git python-devel zlib-devel ncurses-devel
  ```

  ```
  yum groupinstall Development tools
  ```
  
# To install

  1. Download Emboss

   At this time you have to manually download the EMBOSS package manually using your browser.
   
   Download the EMBOSS-6.6.0.tar.gz from EMBOSS ftp site or from the github repo into usamriidPathDiscov/downloads/ using one of the links below:
   - ftp://emboss.open-bio.org/pub/EMBOSS/EMBOSS-6.6.0.tar.gz
   - https://github.com/VDBWRAIR/usamriidPathDiscov/releases/download/v4.0.3/EMBOSS-6.6.0.tar.gz

  2. Run installation instructions(you should be able to copy paste this entire section)

    ```
    git clone https://$(read -p "Gitub username: " gu; echo $gu)@github.com/VDBWRAIR/usamriidPathDiscov.git
    cd usamriidPathDiscov
    cp ~/Downloads/EMBOSS-6.6.0.tar.gz usamriidPathDiscov/download/
    wget --no-check-certificate https://pypi.python.org/packages/source/v/virtualenv/virtualenv-1.11.6.tar.gz -O- | tar xzf -
    python virtualenv-1.11.6/virtualenv.py usamriidPathDiscov
    source usamriidPathDiscov/bin/activate
    pip install paver
    python setup.py install
    deactivate
    source ~/.bashrc
    ```
    
  3. Setup databases under your home directory

    1. Setup databases directory
    
      ```
      mkdir -p ~/databases/{humandna,humanrna,ncbi}
      mkdir -p ~/databases/ncbi/blast/{nt,nr}
      ```
      
    2. You need to have the human dna/rna genome indexed by bowtie
    
      You can download pre-built indexes from Illumina's iGenomes page(http://support.illumina.com/sequencing/sequencing_software/igenome.html)
      You need to put the human dna under ~/databases/humandna and the human rna under ~/databases/humanrna
      
    3. You need both the dna and rna ncbi databases setup under ~/databases/ncbi/blast
    
       ```
       get_blast_dbs.sh ~/databases/ncbi/blast nt nr taxdb
       ```

  4. Quick Verify of all components

    ```
    # These should now all be in your path so should work
    apps=( bwa samtools bowtie2 Ray Ray2 cutadapt getorf run_standard_stable4.pl )
    for p in ${apps[@]}; do $p --help 2>&1 | grep -qiE '[main]|usage|useage|qualifiers' && echo "$p runs" || echo "$p broken?"; done
    ```

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
* Tyghe Vallard
