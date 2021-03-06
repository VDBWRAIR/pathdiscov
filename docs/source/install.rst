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

        git clone $(eval echo https://$(read -p "Gitub username: " gu; echo $gu)@github.com/VDBWRAIR/pathdiscov.git)
        
    .. code-block:: bash
    
        cd pathdiscov

#. Setup a :ref:`virtualenv <activate>` to install into and build documentation

    #. Install virtualenv python environment

        .. code-block:: bash

            wget https://raw.githubusercontent.com/necrolyte2/bootstrap_vi/master/bootstrap_vi.py -O- | python - pathdiscov --prompt="(pathdiscov)"

    #. Activate the virtualenv to install everything into

        .. code-block:: bash

            source pathdiscov/bin/activate
            pip install paver

    #. If you want to view/install the built html documentation

        .. code-block:: bash

            paver doc_html
            firefox docs/build/html/install.html#id1

    #. If you want to view/install the man page documentation

        .. code-block:: bash

            paver doc_man
            mkdir -p pathdiscov/man/man1
            cp docs/build/man/* pathdiscov/man/man1
            man pathdiscov


#. Setup pathdiscov/files/config.yaml.base

        .. code-block:: bash

            mkdir -p ~/tmp  # or change the location  of  `tmp` dir in `pathdiscov/files/config.yaml.base` required for `diamond`
            cp pathdiscov/files/config.yaml{.base,}


#. Edit config.yaml to suite your setup

    .. code-block:: bash

        vim pathdiscov/files/config.yaml

    Example edits:

    .. code-block:: bash

        SEQUENCE_PLATFORM: illumina #choices are: illumina,454


#. Install the pipeline into the virtualenv

    .. code-block:: bash

        python setup.py install

#. Databases setup

    You must refer to built documentation to set up these databases. These databases must be built before you can verify below.

    See :doc:`databases` or `<databases.rst>`_ if you have not built the docs

#. Quick verify of a few things

    * See if required executables are available

        .. code-block:: bash

            # These should now all be in your path so should work
            apps=( bwa samtools bowtie2 blastx blastn Ray Ray2 cutadapt getorf run_standard.pl fastqc prinseq-lite.pl diamond snap pathdiscov_cli)
            for p in ${apps[@]}; do $p --help 2>&1 | grep -qiE '\[main\]|usage|useage|qualifiers|DESCRIPTION|Syntax' && echo "$p ok" || echo "$p broken?"; done

    * See if your databases are available as specified in config

        .. code-block:: bash

            verifydatabases pathdiscov/files/config.yaml

#. Optional: Run a Paired-end dataset

    Anytime you run the pipeline you need to activate the pipeline first. If the pipeline is activated you will see 
    ```(pathdiscov)``` in front of your prompt.
    
    If it is not activated:
    
    .. code-block:: bash
    
        source ~/pathdiscov/pathdiscov/bin/activate

    .. code-block:: bash

        pathdiscov_cli --R1 testData/F.fastq --R2 testData/R.fastq --outdir testoutDir

    You can check this project against when it was run during development by heading
    over to :ref:`Inspect Stage Counts<count-files-interpretation>`

Offline Installation
====================

There may be some instances where you need to install onto an offline workstation. You can achieve this by the following method

#. Clone the pathdiscov project from github
#. Download all of the required software prior to installation and place in pathdiscov/download
    * `htslib <https://github.com/samtools/htslib>`_
    * `samtools <https://github.com/samtools/samtools>`_
    * `bwa <https://github.com/lh3/bwa>`_
    * `fastqc <http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.2.zip>`_
    * `snap <https://github.com/amplab/snap.git>`_
#. Download all of the required python packages

    .. code-block:: bash

        mkdir -p pathdiscov/download/python_packages; pip install --no-use-wheel -d pathdiscov/download/python_packages -r requirements-dev.txt 
        pip install --no-use-wheel -d pathdiscov/download/python_packages virtualenv paver

#. Once downloaded make sure all of the files are extracted if needed and the following directories/files exist
    * pathdiscov/download/htslib
    * pathdiscov/download/samtools
    * pathdiscov/download/bwa
    * pathdiscov/download/fastqc_v0.11.2.zip
#. Now you can copy the git cloned pathdiscov directory to your offline workstation to kick off the install

    .. code-block:: bash

        cd pathdiscov

#. Install virtualenv and python packages into that virtualenv

    .. code-block:: bash

        tar xzf pathdiscov/download/python_packages/virtualenv*
        python virtualenv*/virtualenv.py pathdiscov
        . pathdiscov/bin/activate
        pip install --no-index --find-links=pathdiscov/download/python_packages six argparse numpy paver
        pip install --no-index --find-links=pathdiscov/download/python_packages -r requirements-dev.txt
#. Now you can start the normal installation process from the Databases setup step

