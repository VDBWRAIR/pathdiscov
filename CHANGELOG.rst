Changelog
=========

Version 4.3.0
-------------

* Added get_blast_reads script which allows user to retrieve all reads for a specific
  blast col=>val pair
* Fixed a minor bug where par_block_blast.pl in iterative_blast_phylo would spawn
  an extra blast process in some cases
* Fixed issue with make_pie where debug lines were printed
* Removed/Consolidated redundant scripts
* pathdiscov/scripts/* now copied into virtualenv's bin and are referenced
  without need for path_scripts
* All stages can now be run independently with minimal effort. That is, you can run
  host_map.pl, iterative_blast_phylo.pl, ray2_assembly.pl, quality_filter.pl
  orf_filter.pl, or any script that came with the pipeline inside of the
  pathdiscov/scripts directory directly.
  This makes debugging the pipeline much easier as you can simply ensure the 
  virtualenv is activated and copy/paste any command from the log files into your
  terminal and it should run.
* Count files for every stage now contain input line to verify that each stage
  is receiving the same amount of reads from the prior stage.
* Fixed bug related to orf_filter incorrectly being run inside of 
  iterative_blast_phylo_1
* iterative_blast_phylo_1/contig.count now contains orf_filer count
* linecount now automatically detects file format for fasta/fastq and will fallback
  to counting all lines in file
* Fixed a bug with make_summary where ``--group-by`` was not being utilized
* pathdiscov_cli now logs any [cmd] or [module] line from analysis.log as the
  pipeline runs to give more feedback of what stage the pipeline is on as it runs.
* pathdiscov_cli's output is now more concise and does not print extra statements
  that are not useful.
* Running the Pipeline docs now contain an actual full example to help show
  how to use the pipeline
* Fixes a bug where blastx was never being run

Version 4.2.3
-------------

* Fixed a bug with make_summary looking for missing order column
* Fixed documentation for database generation

Version 4.2.2
-------------

* Fixed a bug with verifydatabases introduced when config.yaml had host_dna and
  host_rna implemented
* Fixed host_map documentation to reflect above bug fix
* ray2_assembly_1 now produces reads_by_contig that contains fastq files for each
  contig containing the reads that mapped to it

Version 4.2.1
-------------

* Fixed the version of snap to 0.15
* Fixed some errors with prinseq perl scripts not using correct perl path
* Fixed documentation in reference to some commands not existing during
  database setup steps

Version 4.2
-----------

* Added diamond as option for blast_db in param.txt(iterative_blast_phylo)
* Added snap as option for host_map stage
* Added error message if system is not 64 bit
* Added instructions on how to create diamond indexed database
* Changed tests/rikkcdna such that developer builds databases instead of shipping
  prebuilt diamond and blast dbs
* Added pathdiscov_cli to verify script to make sure everything installed
* Fixed typo in installation.rst
* Updated documentation for iterative_blast_phylo on how to use diamond
* Fixed issue where fastqc was not logging stderr to analysis_quaility.log
* Renamed run_standard_stable4.pl -> run_standard.pl
* Removed duplicated run_standard_*.pl scripts
* Pipeline no longer creates absolute path'd symlinks which makes the output
  directory portable
* Broke up the databases documentation a bit for each of the indexing steps
  such that the user can more clearly see how each database is created
* Added references in databases documentation to the corresponding 
  stage's documentation to link the two together.
* Added more documentation about the databases directory structure
* Added ability to select snap as host_map aligner
* Added documentation inside of host_map stage about snap as well as a more
  thorough description of each output file
* Added -amos flag for Ray2 so that AMOS.afg file is created inside of
  results/ray2_assembly_X/results/
* Added --blast-unassembled option to main script
* Added documentation for --blast-unassembled
* Main script now echos run_standard.pl command used into the analysis.log
* run_standard.pl now echos more about unassembled read numbers that are used
  into analysis.log
* orf_filter can now work on files that are not specifically formatted with
  @1 or >1. That is, default miseq fastq files work now
* added count files for orf_filter
* orf_filter now accepts --contig argument to work on single contig file and 
  will produce contig.* files instead of R1.*/R2.*
* no longer runs orf_filter as separate stage. Now runs only if blastx or
  diamond are in the blast_task_list of iterative_blast_phylo
* fixed a bug where ``paver doc_man`` didn't work
* Fixed a bug with make_summary where contig may have wrong length/numreads
* Fixed a bug where having trailing slash on a project path when using
  make_summary would not produce sample names in the summary
* Fixed a bug when there were more unassembled reads than contigs that the
  unassembled reads would have an extra column
* Changed config.yaml such that human_dna and human_rna are now host_dna
  and host_rna
* Fixed error in install documentation where the verify executables step
  would always return 'ok'
* make_pie now will parse pretty much any tab delimited file to get results
* make_pie now places total count next to percentages

Version 4.1
-----------

* Fixed issue with Database setup where sed command did not correctly replace
  path to bowtie index in config.yaml
* Fixed issue with bowtie2-build not being in path during database setup
* Fixed issue with iterative_blast_phylo handling fasta config files
* Implemented functional tests on testData
* Added java as dependency for pipeline as fastqc requires it
* Fixed issue where absolute path to none was being passed around
* Merged in a cleaner param.txt given to us from usamriid
* Fixed an issue where iterative_blast_phylo_2 was trying to use out.cap.fa from
  ray2_assembly step instead of unmapped reads
* Fixed issue where error logs had -num_descriptions ignored in them
* Added arguments to make_pie that allow you to specify host, vector and pathogen
  so the graphic can be customized
* Fixed an issue where step1 was not handling bad --R2 input(such as absolute path
  to none)
* Fixed an issue with make_pie where when the column in the blast report it was
  parsing was empty, it would return an empty string as the label. It now 
  looks for the next column to the left for the label until it finds a non empty
  field.
* Fixed a bug where if --outdir was given an absolute path many bad side-effects
  happened.
* Added SGE/PBS support for iterative_blast_phylo
* Added --cpu option to pathdiscov_cli
* Fixed a bug where --R1 and --R2 input files were copied into the input folder
  of the analysis creating duplicate files. Now input will only have input files
  if they are unpacked(.gz) or are sff files and converted to .fastq
* Added gzip compressed input file support(.fastq.gz and .sff.gz)
* Added single read support. Previously both --R1 and --R2 needed to be supplied
* Added support to use --R1/--R2 as well as -R1/-R2 to pathdiscov_cli
* Various improvements to the documentation
* Documentation now has instructions for offline installation
