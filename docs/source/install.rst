============
Installation
============

.. _install-system-packages:

System Packages
===============

In order to installed developer tools you will need root priveledges; that is, somebody who can use
su or sudo

CentOS
------

.. code-block:: bash

    #> yum install openmpi openmpi-devel git python-devel zlib-devel ncurses-devel freetype-devel libpng-devel wget
    #> yum groupinstall Development Tools
    
Ubuntu
------

.. code-block:: bash

    #> apt-get install openmpi-bin libopenmpi-dev python-dev git zlib1g-dev build-essential libncurses5	libncurses5-dev libpng12-dev libfreetype6-dev


Installation
============

#. Clone the repository

    .. code-block:: bash

        git clone $(eval echo https://$(read -p "Gitub username: " gu; echo $gu)@github.com/VDBWRAIR/usamriidPathDiscov.git)
        cd usamriidPathDiscov

#. Setup a :ref:`virtualenv <activate>` to install into

    .. code-block:: bash

        wget --no-check-certificate https://pypi.python.org/packages/source/v/virtualenv/virtualenv-1.11.6.tar.gz -O- | tar xzf -
        python virtualenv-1.11.6/virtualenv.py usamriidPathDiscov
        source usamriidPathDiscov/bin/activate
        pip install paver
        cat <<EOF



        **************************************************************
        Remember this path as you will need it every time you activate
        the pipeline:


        $(pwd)/usamriidPathDiscov/bin/activate


        **************************************************************


        EOF
        echo

#. Databases setup

    You must refer to built documentation to set up these databases. These databases must be built before you can verify below.

    See :doc:`databases` or `<databases.rst>`_ if you have not built the docs

#. Setup usamriidPathDiscov/files/config.yaml.base

    #. Copy config.yaml.base to config.yaml

        .. code-block:: bash

            cp usamriidPathDiscov/files/config.yaml{.base,}

    #. Edit config.yaml to suite your setup
    
        .. code-block:: bash

            vim usamriidPathDiscov/files/config.yaml

        Example edits:

        .. code-block:: bash

            SEQUENCE_PLATFORM: illumina #choices are: illumina,454
            NODE_NUM: 10 # number of blast partition depending on the number of CPU on your computer. If you have 12 CPU on on your workstation, '10' works, if you have more CPU increase this number

#. Install the pipeline into the virtualenv

    .. code-block:: bash

        python setup.py install

#. Build and view the complete documentation

    This will open a new firefox window that will display the built documentation
    that you can continue on where you left off here

    .. code-block:: bash

        cd docs
        make clean && make html
        firefox build/html/install.html &
        cd ..

#. Quick verify of a few things

    * See if required executables are available

        .. code-block:: bash

            # These should now all be in your path so should work
            apps=( bwa samtools bowtie2 blastx blastn Ray Ray2 cutadapt getorf run_standard_stable4.pl fastqc prinseq-lite.pl )
            for p in ${apps[@]}; do $p --help 2>&1 | grep -qiE '[main]|usage|useage|qualifiers' && echo "$p ok" || echo "$p broken?"; done

    * See if your databases are available as specified in config

        .. code-block:: bash

            verifydatabases usamriidPathDiscov/files/config.yaml

#. Optional: Run a Paired-end  dataset

    Anytime you run the pipeline you need to activate the pipeline first. If the pipeline is activated you will see ```
    (usamriidPathDiscov)``` in front of your prompt.
    
    If it is not activated:
    
    .. code-block:: bash
    
        source ~/usamriidPathDiscov/usamriidPathDiscov/bin/activate

    If your blast database is quite large (like the default nt database) this could take up to 2 hours...
    It is recommended that you trim down your nt databases to just the things that you are interested in.
    
    You may input `*.fastq` or `*.fastq.gz` file

    .. code-block:: bash

         usamriidPathDiscov_cli -R1 $(pwd)/testData/F.fastq.gz -R2     $(pwd)/testData/R.fastq.gz --outdir testoutDir
 
OR

    .. code-block:: bash

        usamriidPathDiscov_cli -R1 $(pwd)/testData/F.fastq -R2 $(pwd)/testData/R.fastq --outdir testoutDir

#. Option 2: Run a Single-end dataset

   .. code-block:: bash
           
       samriidPathDiscov_cli -R1 $(pwd)/testData/F.fastq.gz  --outdir      testoutDir

OR

    .. code-block:: bash
        
       samriidPathDiscov_cli -R1 $(pwd)/testData/F.fastq  --outdir testoutDir

#. Option 3: Using Sun Grid Engine
    
   If your cluster support SGE, pass `--sge 1`, default is  `--sge 0`
    
   .. code-block:: bash

         usamriidPathDiscov_cli -R1 $(pwd)/testData/F.fastq.gz -R2  $(pwd)/testData/R.fastq.gz --outdir testoutDir --sge 1

Offline Installation
====================

There may be some instances where you need to install onto an offline workstation. You can achieve this by the following method

#. Clone the usamriidPathDiscov project from github
#. Download all of the required software prior to installation and place in usamriidPathDiscov/download
    * `htslib <https://github.com/samtools/htslib>`_
    * `samtools <https://github.com/samtools/samtools>`_
    * `bwa <https://github.com/lh3/bwa>`_
    * `fastqc <http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.2.zip>`_
#. Download all of the required python packages

    .. code-block:: bash

        mkdir -p usamriidPathDiscov/download/python_packages; pip install --no-use-wheel -d usamriidPathDiscov/download/python_packages -r requirements-dev.txt 
        pip install --no-use-wheel -d usamriidPathDiscov/download/python_packages virtualenv paver

#. Once downloaded make sure all of the files are extracted if needed and the following directories/files exist
    * usamriidPathDiscov/download/htslib
    * usamriidPathDiscov/download/samtools
    * usamriidPathDiscov/download/bwa
    * usamriidPathDiscov/download/fastqc_v0.11.2.zip
#. Now you can copy the git cloned usamriidPathDiscov directory to your offline workstation to kick off the install

    .. code-block:: bash

        cd usamriidPathDiscov

#. Install virtualenv and python packages into that virtualenv

    .. code-block:: bash

        tar xzf usamriidPathDiscov/download/python_packages/virtualenv*
        python virtualenv*/virtualenv.py usamriidPathDiscov
        . usamriidPathDiscov/bin/activate
        pip install --no-index --find-links=usamriidPathDiscov/download/python_packages six argparse numpy paver
        pip install --no-index --find-links=usamriidPathDiscov/download/python_packages -r requirements-dev.txt
#. Now you can start the normal installation process from the Databases setup step

