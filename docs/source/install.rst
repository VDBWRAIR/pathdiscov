Installation
============

System Packages
---------------

Install developmental tools like compilers and make and such
You need an administrator to run the following commands. That is, somebody who can use
su or sudo

.. code-block:: bash

    #> yum install openmpi openmpi-devel git python-devel zlib-devel ncurses-devel
    #> yum groupinstall Development Tools

Installation
------------

1. Clone the repository

    .. code-block:: bash

        git clone $(eval echo https://$(read -p "Gitub username: " gu; echo $gu)@github.com/VDBWRAIR/usamriidPathDiscov.git)
        cd usamriidPathDiscov

2. Setup a virtualenv to install into

    .. code-block:: bash

        wget --no-check-certificate https://pypi.python.org/packages/source/v/virtualenv/virtualenv-1.11.6.tar.gz -O- | tar xzf -
        python virtualenv-1.11.6/virtualenv.py usamriidPathDiscov
        source usamriidPathDiscov/bin/activate
        pip install paver


3. Edit "usamriidPathDiscov/files/config.yaml.base" if necessary. For example edit the following lines...

    .. code-block:: bash

        SEQUENCE_PLATFORM: illumina  #choices are: illumina,454
        NODE_NUM: 10  # number of blast partition depending on the number of CPU on your computer. If you have 12 CPU on on your workstation, '10' works, if you have more CPU increase this number

4. Install the pipeline into the virtualenv

    .. code-block:: bash

        python setup.py install
5. Blast/Bowtie databases setup

    These databases will be dependent on your specific situation, but by default we will use the blast nt/nr databases as well as the human genome to
    do host mapping

    1. Setup directories

        .. code-block:: bash
        
            mkdir -p ~/databases/{humandna,humanrna,ncbi}
            mkdir -p ~/databases/ncbi/blast/{nt,nr}

    2. You also need to get bowtie indexed host genomes

        You should be able to fetch these from Illumina, Ensamble or NCBI. Illumina seems to have the easiest ones as they are pre-built where ncbi and ensamble require you to build them manually with bowtie
        Illumina's iGenomes page can be found at http://support.illumina.com/sequencing/sequencing_software/igenome.html

        Whichever ones you decide on, you need to put them under ~/databases/humandna and ~/databases/humanrna
        We will add more documentation on how to do other hosts later, but in general you can check out the configuration.rst file for
        more information on how to configure the pipeline to use different indexes.

    3. You need to then setup the blast databases

        In general you just need to unpack the nt/nr databases from ncbi(or wherever) into ~/databases/ncbi/blast/nt,nr,taxdb
        There is a shell script you should be able to use to do this for you as well.
        This may take longer time depending on your connection. 

        .. code-block:: bash

            usamriidPathDiscov/scripts/get_blast_dbs.sh ~/databases/ncbi/blast nt nr taxdb

6. Quick verify of necessary executables

    .. code-block:: bash

        # These should now all be in your path so should work
        apps=( bwa samtools bowtie2 Ray Ray2 cutadapt getorf run_standard_stable4.pl fastqc )
        for p in ${apps[@]}; do $p --help 2>&1 | grep -qiE '[main]|usage|useage|qualifiers' && echo "$p runs" || echo "$p broken?"; done

7. Optional: Run a sample dataset

    If your blast database is quite large(like the default nt database) this could take up to 2 hours...
    It is recommended that you trim down your nt databases to just the things that you are interested in

    .. code-block:: bash

        usamriidPathDiscov_cli -R1 $(pwd)/testData/F.fastq  -R2 $(pwd)/testData/R.fastq  --outdir  testoutDir
