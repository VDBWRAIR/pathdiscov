Changelog
=========

Version 4.2
-----------

* Added diamond as option for blast_db in param.txt(iterative_blast_phylo)
* Added snap as option for host_map stage
* Added error message if system is not 64 bit
* Added instructions on how to create diamond indexed database
* Changed tests/rikkcdna such that developer builds databases instead of shipping
  prebuilt diamond and blast dbs
* Added usamriidPathDiscov_cli to verify script to make sure everything installed
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
* Added --cpu option to usamriidPathDiscov_cli
* Fixed a bug where --R1 and --R2 input files were copied into the input folder
  of the analysis creating duplicate files. Now input will only have input files
  if they are unpacked(.gz) or are sff files and converted to .fastq
* Added gzip compressed input file support(.fastq.gz and .sff.gz)
* Added single read support. Previously both --R1 and --R2 needed to be supplied
* Added support to use --R1/--R2 as well as -R1/-R2 to usamriidPathDiscov_cli
* Various improvements to the documentation
* Documentation now has instructions for offline installation
