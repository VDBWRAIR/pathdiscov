========
make_pie
========

Creates a graphic for each supplied project based on all of the .top.blast.phylo files found in both the
iterative_blast_phylo output directories.

The graphic dipicts sample compisition in the form of 4 pie graphics depicting:

* Vector composition
* Host composition
* Pathogen composition
* Sample composition

Example
=======

This example was built using the testoutDir created via the example in the :doc:`install`

.. code-block:: bash

    make_pie testoutDir

This will create a directory called host_vector_pathogen inside of the current directory. Inside of that directory
you will find the file testoutDir.png that looks like this:

.. image:: _static/example_make_pie.png
    :width: 100%

Vector Composition
==================

Based on the class field and only includes Insecta

Host Composition
================

Based on the class field and only includes Mammalia

Pathogen composition
====================

Based on the superkingdom field and only includes Bacteria and Viruses

API
===

usamriidPathDiscov.make_pie
---------------------------

.. automodule:: usamriidPathDiscov.make_pie
    :members:
    :undoc-members:
