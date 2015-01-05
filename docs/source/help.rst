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
