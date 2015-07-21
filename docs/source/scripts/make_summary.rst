============
make_summary
============

Creates a summary report for all given projects by compiling together information 
from various output files from each project. The summary generated is in tab 
separated(tsv) format and is by default printed to the terminal.
To save the output, the user needs to use the common Unix redirect 
operator `>` such as

.. code-block:: bash

    make_summary testoutDir > summary.tsv

Then you can easily open the report using your spreadsheet program and using only 
the Tab as the delimiter when asked.

Reading the Report
==================

Make sure to open the report using only the tab character as the field separator
or you will not be able to view the report correctly.

The report is split up into two sections side by side and the column widths may 
need to be adusted in your spreadsheet program to view it properly.

The column on the left contains the Contig report while the column on the right 
contains the unasemmbled reads report

Contig Report
-------------

The contig report contains the following columns

* Sample Name
    Name of the sample(project directory name)
* Num Reads
    Total number of input reads
    Sum of dc-megablast for R1 and R2 inside of iterative_blast_phylo_1/\*.count
* Non-Host Num reads
    Total of reads left over after host mapping and quality filter
    Sum of input lines from ray2_assembly_1/\*.count
* Num Ctg
    Number of contigs generated
    last line in ray2_assembly_1/assembly.count(probably cap_contigs)
* Num Unblasted Ctg
    Number of contigs that had 0 blast results
    The last value in iterative_blast_phylo_1/contig.counts
* N50
    N50 of all contig lengths
* Assembly Length
    The length of all contigs concatted together
* Ctg#
    The name of the contig in the project
* Ctg bp
    The length of the contig
    Gathered from ray2_assembly_1/contig_len.txt which is generated from 
    augment_report.sh
* numReads
    Number of reads supporting contig
    Gathered from ray2_assembly_1/contig_numreads.txt which comes
   from the ray2_assembly_1/bowtie_mapping/out.bam by counting all mapped reads
* Accession, Family, Genus, description
    Blast information for contig

Unassembled Report
------------------

The unasembled report contains information on the reads that did not map
The rows are compiled from multiple files

* Num input unassem
    Number of reads selected from unassembled reads from ray2_assembly_1/bowtie_mapping
    This number may be a lot smaller than you expect due to the 
    ``--blast_unassembled`` argument given to run_standard.pl(Default 1000)
* Num Unblasted Unassem
    Number of unassembled reads that had 0 blast results
    Essentially the last value in iterative_blast_phylo_2/\*.counts
* Sum of reads that group to the same ``--group-by``
    iterative_blast_phylo_2/reports/{R1,R2}.\*.top.smallreport.txt are both used
    and all results are grouped on the ``--group-by`` supplied and all unique
    values are summed.
* Accession, Family, Virus Genus, descrip
    Blast information for reads
