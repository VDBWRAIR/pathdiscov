============
Installation
============

.. _install-system-packages:

System Packages
===============

Install developmental tools like compilers and make and such
You need an administrator to run the following commands. That is, somebody who can use
su or sudo

CentOS
------

.. code-block:: bash

    #> yum install openmpi openmpi-devel git python-devel zlib-devel ncurses-devel freetype-devel libpng-devel
    #> yum groupinstall Development Tools
    
Ubuntu
------

.. code-block:: bash

    #> apt-get install openmpi-bin libopenmpi-dev python-dev git zlib1g-dev build-essential libncurses5	libncurses5-dev libpng12-dev libfreetype6-dev


Installation
============

#. Setup a virtualenv to install into

    .. code-block:: bash

        wget --no-check-certificate https://pypi.python.org/packages/source/v/virtualenv/virtualenv-1.11.6.tar.gz -O- | tar xzf -
        python virtualenv-1.11.6/virtualenv.py usamriidPathDiscov
        source usamriidPathDiscov/bin/activate
        pip install paver

#. Build and view the complete documentation

    This will open a new firefox window that will display the built documentation
    that you can continue on where you left off here

    .. code-block:: bash

        cd docs
        make clean && make html
        firefox build/html/install.html
        cd ..

#. Setup `usamriidPathDiscov/files/config.yaml.base <../../../usamriidPathDiscov/files/config.yaml.base>`_

    #. Copy config.yaml.base to config.yaml

        .. code-block:: bash

            cp usamriidPathDiscov/files/config.yaml{.base,}

    #. Edit config.yaml to suite your setup

        Example:

        .. code-block:: bash

            SEQUENCE_PLATFORM: illumina  #choices are: illumina,454
            NODE_NUM: 10  # number of blast partition depending on the number of CPU on your computer. If you have 12 CPU on on your workstation, '10' works, if you have more CPU increase this number

#. Install the pipeline into the virtualenv

    .. code-block:: bash

        python setup.py install
#. Blast/Bowtie databases setup

    See :doc:`databases`

#. Quick verify of necessary executables

    .. code-block:: bash

        # These should now all be in your path so should work
        apps=( bwa samtools bowtie2 blastx blastn Ray Ray2 cutadapt getorf run_standard_stable4.pl fastqc )
        for p in ${apps[@]}; do $p --help 2>&1 | grep -qiE '[main]|usage|useage|qualifiers' && echo "$p ok" || echo "$p broken?"; done

#. Optional: Run a sample dataset

    If your blast database is quite large(like the default nt database) this could take up to 2 hours...
    It is recommended that you trim down your nt databases to just the things that you are interested in

    .. code-block:: bash

        usamriidPathDiscov_cli -R1 $(pwd)/testData/F.fastq  -R2 $(pwd)/testData/R.fastq  --outdir  testoutDir
