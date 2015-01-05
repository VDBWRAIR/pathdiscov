====
Help
====

Eventually you will run across some errors. No application/software is without bugs. Here we will compile all of the most common errors and what to look for to find out what is going on

.. _faq:

Frequently Asked Questions
==========================

#. How many CPU's does my computer have?

    This will print how many CPU Cores your computer has

    .. code-block:: bash

        lscpu | awk -F':' 'BEGIN {cpu=1} /(Core|Socket)/ {gsub(/ /,"",$0); cpu *= $2;} END {print cpu}'

#. I'm not sure if the pipeline completed sucessfully. How do I check the log files?

    Each stage has its own log directory where you can look for errors and messages. An easy way to just look at all log output is
    
    .. code-block:: bash

        cat results/*/logs/*

    If you are looking for something specific, you can even search for it. Say if you want to find any logs that have the word error in them

    .. code-block:: bash

        grep -i 'error' results/*logs/*
