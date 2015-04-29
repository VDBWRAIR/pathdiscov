========
make_pie
========

Creates a graphic for each supplied project based on all of the .top.blast.phylo 
files found in both the iterative_blast_phylo output directories.
Likely these counts will not reflect all your input reads as the host_map stage 
typically runs prior to these and will remove reads that mapped to your host

The graphic depicts shows sample composition broken up in 4 ways:

* Vector composition
* Host composition
* Pathogen composition
* Sample composition

Vector, Host and Pathogen are broken up into both a plot and a pie graphic
because some times there are a lot of results that overwhelm the pie graphic.
The plot shows all results found while the pie graphic shows you quickly the top
hits(where counts are greater than 5) and combines all results that are less 
than 5 into a single 'other' wedge.

Example
=======

This example was built using the testoutDir created via the example in the :doc:`../install`

.. code-block:: bash

    make_pie testoutDir

This will create a directory called host_vector_pathogen inside of the current 
directory. Inside of that directory you will find the file testoutDir.png that 
looks like this:

.. image:: ../_static/example_make_pie.png
    :width: 100%

Vector Composition
==================

Based on the class field and only includes Insecta by default.
It is configurable by using the ``--vectorclasses`` argument.

Host Composition
================

Based on the class field and only includes Mammalia by default.
It is configurable by using the ``--hostclasses`` argument.

Pathogen Composition
====================

Based on the superkingdom field and only includes Bacteria and Viruses by
default.
It is configurable by using the ``--pathogenclasses`` argument.

Sample Composition
==================

Overview that shows Host, Vector and Pathogen composition as a whole for the 
sample.

Files Used from analysis
========================

* results/iterative_blast_phylo*/\*.top.blast.phylo

How phylo files are used
========================

* Grab all results/iterative_blast_phylo_*/\*.top.blast.phylo files
* Each File is opened and read line by line
* Line is split by tab and the very top row is used as the column names.
  The following are the columns that are used
  * count
  * superkingdom
  * class
  * species
* A mapping of vectors, hosts and pathogens is created based on each found 
  species
* If class matches Mammalia
  * increment count for host
  * increment mapping count for hosts[species]
* If class matches Insecta
  * increment count for vectors
  * increment mapping count for vectors[species]
* If superkingdom matches Viruses or Bacteria
  * increment count for pathogens
  * increment mapping count for pathogens[species]
* You now have counts for the following
  * total hosts and total hosts broken down by species
  * total vectors, total vectors broken down by species
  * total pathogens, total pathogens broken down by species
  * hosts + vectors + pathogens = total overall

*Note*: Sometimes there is no species set and there is a dash. If this happens
then columns to the left will be searched until a non dash column 
is found and that will be used for the species name. This means you may end up
with a superkingom for species name or similar.

Each graphic is created by looking through the created mappings and using
the keys(species) for the pie slice labels and the counts associated to those 
species for the values.
