============
make_summary
============

Creates a summary report for all given projects by compiling together information from various output files from each project.
The summary generated is in tab separated(tsv) format and is by default printed to the terminal.
To save the output, the user needs to use the common Unix redirect operator `>` such as

.. code-block:: bash

    make_summary testoutDir > summary.tsv

Then you can easily open the report using your spreadsheet program and using only the Tab as the delimiter when asked.

Reading the Report
==================

The report is split up into two sections side by side and the column widths may need to be adusted in your spreadsheet program to view it properly.

The column on the left contains the Contig report while the column on the right contains the unasemmbled reads report

Contig Report
-------------

The contig report contains the following columns

* Sample Name
    Name of the sample(project directory name)
* Num Reads
    Total number of reads
* Non-Host Num reads
    Reads not mapped to host genome index
* Num Ctg
    Number of contigs generated
* Num blast0 Ctg
    Not sure
* N50
    N50 of all contig lengths
* Ctg#
    The name of the contig in the project
* Ctg bp
    The length of the contig
* numReads
    Number of reads supporting contig
* Accession, Family, Genus, description
    Blast information for contig

Unassembled Report
------------------

The unasembled report contains information on the reads that did not map
The rows are compiled from multiple files

* Num unassem
    Number of reads that were not mapped
* Num blast0 Unassem
    Not sure
* num reads
    Number of reads matching blast row
* Accession, Family, Virus Genus, descrip
    Blast information for reads
