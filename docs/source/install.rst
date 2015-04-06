============
Installation
============

.. _install-system-packages:

System Packages
===============

In order to installed developer tools you will need root privileges; that is, somebody who can use
su or sudo

CentOS
------

.. code-block:: bash

    #> yum install openmpi openmpi-devel git python-devel zlib-devel ncurses-devel freetype-devel libpng-devel wget java-1.6.0 dejavu*
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
        
    .. code-block:: bash
    
        cd usamriidPathDiscov

#. Setup a :ref:`virtualenv <activate>` to install into and build documentation

    #. Install virtualenv python environment

        .. code-block:: bash

            wget --no-check-certificate https://pypi.python.org/packages/source/v/virtualenv/virtualenv-1.11.6.tar.gz -O- | tar xzf -
            python virtualenv-1.11.6/virtualenv.py usamriidPathDiscov

    #. Activate the virtualenv to install everything into

        .. code-block:: bash

            source usamriidPathDiscov/bin/activate
            pip install paver

    #. If you want to view/install the built html documentation

        .. code-block:: bash

            paver doc_html
            firefox docs/build/html/install.html#id1

    #. If you want to view/install the man page documentation

        .. code-block:: bash

            paver doc_man
            mkdir -p usamriidPathDiscov/man/man1
            cp docs/build/man/* usamriidPathDiscov/man/man1
            man usamriidPathDiscov


#. Setup usamriidPathDiscov/files/config.yaml.base

        .. code-block:: bash

            mkdir -p ~/tmp  # or change the location  of  `tmp` dir in `usamriidPathDiscov/files/config.yaml.base` required for `diamond`
            cp usamriidPathDiscov/files/config.yaml{.base,}


#. Edit config.yaml to suite your setup

    .. code-block:: bash

        vim usamriidPathDiscov/files/config.yaml

    Example edits:

    .. code-block:: bash

        SEQUENCE_PLATFORM: illumina #choices are: illumina,454

#. Databases setup

    You must refer to built documentation to set up these databases. These databases must be built before you can verify below.

    See :doc:`databases` or `<databases.rst>`_ if you have not built the docs


#. Install the pipeline into the virtualenv

    .. code-block:: bash

        python setup.py install

#. Quick verify of a few things

    * See if required executables are available

        .. code-block:: bash

            # These should now all be in your path so should work
            apps=( bwa samtools bowtie2 blastx blastn Ray Ray2 cutadapt getorf run_standard_stable4.pl fastqc prinseq-lite.pl diamond snap usamriidPathDiscov_cli)
            for p in ${apps[@]}; do $p --help 2>&1 | grep -qiE '[main]|usage|useage|qualifiers' && echo "$p ok" || echo "$p broken?"; done

    * See if your databases are available as specified in config

        .. code-block:: bash

            verifydatabases usamriidPathDiscov/files/config.yaml

#. Optional: Run a Paired-end dataset

    Anytime you run the pipeline you need to activate the pipeline first. If the pipeline is activated you will see 
    ```(usamriidPathDiscov)``` in front of your prompt.
    
    If it is not activated:
    
    .. code-block:: bash
    
        source ~/usamriidPathDiscov/usamriidPathDiscov/bin/activate

    You may change the number of CPU based on the resource in your
    system.

    .. code-block:: bash

        usamriidPathDiscov_cli -R1 testData/F.fastq -R2 testData/R.fastq --outdir testoutDir --cpuNum 12

    If your blast database is quite large (like the default nt database) this could take up to 2 hours...
    It is recommended that you trim down your nt databases to just the things that you are interested in.

Offline Installation
====================

There may be some instances where you need to install onto an offline workstation. You can achieve this by the following method

#. Clone the usamriidPathDiscov project from github
#. Download all of the required software prior to installation and place in usamriidPathDiscov/download
    * `htslib <https://github.com/samtools/htslib>`_
    * `samtools <https://github.com/samtools/samtools>`_
    * `bwa <https://github.com/lh3/bwa>`_
    * `fastqc <http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.2.zip>`_
    * `snap <https://github.com/amplab/snap.git>`_
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

