#!/usr/bin/env perl

# ---------------------------------------------------------------
#																#
# 	pathogen.pl													#
#																#
# 	Oliver Elliott oe2118@columbia.edu							#
# 	Jason Ladner jason.t.ladner.ctr@us.army.mil					#
# 	Francesco Abate francesco.abate@gmail.com					#
#																#
#	date:			20121203									#
#	version:		1											#
#																#
#	notes			first modular version						#
#																#
#	This wrapper calls a series of commands and 				#
#	their attendant parameter files (or a single,				#
#	master parameter file), and runs them in the 				#
#	order you choose. It mangages parameter files 				#
#	and log directories and gives you a readout (deltat)		#
#	in seconds of how long your module took to run. 			#
#	The scripts for each module should be stored in				# 
#	a unique folder in the the main scripts directory 			#
#	When a module is called, it will look for scripts in		#
#	its own folder. This wrapper can be used as a template		#
#	for large pipelines.										#
#																#
# ---------------------------------------------------------------

# example:
# pathogen.pl --command initial_mapping --paramfile sample_param_files/pathogen.param --outputdir tmp --R1 input_R1.fastq --R2 input_R2.fastq --sample mysample

use strict;
use Getopt::Long;
use Pod::Usage;
use File::Copy;
use File::Basename;
use Cwd 'abs_path';

use FindBin qw($RealBin);
use lib "$RealBin/Local_Module";

# local modules:
use Verbose_Sys;
# use Parse_ParameterFile;

# TO DO:
# make module hookup automatic
# put stuff in functions where possible
# make error check and count automatic options w flags

# declare global variables:
my ($outputdir);													# output dir
my (@command);														# array of commands
my (@paramfile);													# array of param files
my ($sample);														# sample name
my ($example);														# bool for param file (if 1, print example param file)
my ($filekey);														# bool for file key (if 1, print all possible files created plus description)
my ($help);															# bool for help (if 1, print help page then exit)
my ($man);															# bool for man page (if 1, print man page then exit)
my ($path_scripts, $path_output, $path_logs);						# path variables (scripts, output dir, logs)
my ($r1, $r2);														# paths to r1 and r2, if present
my ($abs_r1, $abs_r2);												# absolute paths to r1 r2
my ($out_r1, $out_r2);												# these variables are used to hold the output from the previous stage == input for next stage
my ($main_start_time, $start_date, $end_date, $main_end_time);		# benchmark variables
my ($pfile);														# parameter file
my ($isfasta);														# a bool for fasta (1 if fasta, 0 if fastq)
my ($is_fasta);														# same as isfasta but a string ("yes" "no") rather than a bool (redundant - bad style!)
my ($sge);															# bool for sun grid engine (if 1,  use SGE qsub parallelization; if 0, no SGE; default is 0)
my ($wellpaired);													# bool for well-pairedness (1 if R1 and R2 are well-paired --- i.e., they have the same read IDs; 0 if not --- i.e., the files are asymmetrical)
my ($checkerror);													# check error
my ($verbose);														# turn on verbose
my ($boolrun);														# bool for running a command (the idea is, this defaults to false but if the output of a command does not already exist, it becomes true and the command is run)
my ($boolcleanup);													# bool for cleaning up (i.e., removing unwanted files after the run has finished)
my $boolcontigname = 0;													# bool for changing the name of the iterative blast module to "contig" rather than "R1"

# set defaults
$example = 0;
$filekey = 0;
$help = 0;
$man = 0;
$sample = "sample";
$abs_r1="none";
$abs_r2="none";
$out_r1 = "none";
$out_r2 = "none";
$is_fasta = "no";
$wellpaired=1;
$sge = 0;

GetOptions ('outputdir=s' => \$outputdir,
            'command=s{1,}' => \@command,			# 1 or more
            'paramfile=s{1,}' => \@paramfile,		# 1 or more
            'R1=s' => \$r1,							# mate 1 or single end reads file
            'R2=s' => \$r2,							# mate 2 file
            'sample=s' => \$sample,    				# same name            
            'fastafile=s' => \$is_fasta,			# is fasta (default assumes fastq) (allow values: "yes" "no")
            'fasta' => \$isfasta,            		# is fasta (default assumes fastq) (bad style - this is a boolean that does the same thing as the above fastafile option but uses a bool instead of a string)            
            'SGE' => \$sge,							# use qsub (default: no)
            'example' => \$example,              	# boolean for example parameter file
            'key' => \$filekey,  	            	# boolean for key to output files
            'checkerror' => \$checkerror,			# check error after running the program
			'verbose' => \$verbose,					# check error after running the program but print all error logs            
			'cleanup' => \$boolcleanup,				# remove intermediate files
            'help' => \$help, 						# show help
            'man' => \$man);						# show man page
            
pod2usage(1) if ($help);
pod2usage(-verbose => 2) if ($man);

# -------------------- main --------------------

# set paths

# scripts paths 
$path_scripts=$RealBin;

# get abs path to output directory:
$path_output=abs_path($outputdir);

# logs directory
$path_logs = $path_output."/logs";

# get abs path to r1 r2:
$abs_r1=abs_path($r1) if (defined($r1));
$abs_r2=abs_path($r2) if (defined($r2));

# set redundant fasta bool (bad style!)
if ($isfasta)
{
	$is_fasta="yes";
}

# print various files, then exit
if ($example)
{
    system("cat $path_scripts/lib/python*/site-packages/usamriidPathDiscov/files/sample.param");
	exit;
}

if ($filekey)
{
	system("cat $path_scripts/files/filekey.txt");
	exit;
}

if ($checkerror)
{
	if ($verbose) { system("$path_scripts/scripts/checkerror.sh $path_output 1"); } else { system("$path_scripts/scripts/checkerror.sh $path_output 0"); }
	exit;
}

if ($boolcleanup)
{
	system("$path_scripts/scripts/cleanup.sh $path_output");
	exit;
}

# check for user errors 
&check_error;		

# print title 
&echo_title;	

# get main starting time in seconds
$main_start_time = time();

# echo start date
$start_date = get_date();
print "[START] $start_date\n";

# create project directory, logs directory
print "[echo] create project logs directory\n";

my $cmd = "mkdir -p $path_logs";
print_system($cmd);

print "\n";		

# if there's only one param file, use it
if (scalar(@paramfile)==1)
{
	$pfile = $paramfile[0];
	
	die "[error] parameter file $pfile not found.\n\n" if (!( -e $pfile ));

	# save a time stamped copy of the parameter file in the main logs directory
	my $logcopy=sprintf($path_logs."/".$sample.".".$start_date."_".basename($pfile));
	copy($pfile, $logcopy);	
} 

# run the modules sequentially
for (my $i = 0; $i < scalar(@command); $i++) 
{	
	# if more than one param file, use the appropriate one in the list
	if (scalar(@paramfile) > 1)
	{
		$pfile = $paramfile[$i];

		# save a time stamped copy of the parameter file in the main logs directory
		my $logcopy=sprintf($path_logs."/".$start_date."_".basename($pfile));
		copy($pfile, $logcopy) if ( -e $pfile );			 		
	} 

	if ( -e $pfile )	# make sure the parameter file exists
	{						
		if ( $command[$i] eq "step1" ) # step1 must be run first ! 
		{
			print "[module] $command[$i] \n";
			$out_r2 = &check_exist($out_r2, "$path_output/$command[$i]/$command[$i].R2");
			# set boolrun to false, but if the output of a command does not already exist, it becomes true and the command is run. 
			# do output R1 second because this is the one whose existence we want to detemine boolrun (for non-paired data, R2 will always not exist)
			$boolrun = 0;			
			$out_r1 = &check_exist($out_r1, "$path_output/$command[$i]/$command[$i].R1");			
			print "\n" if (!($boolrun));
			
			# if output doesnt exist, run module
			if ($boolrun)
			{
				# make module's logs directory
				my $cmd = "mkdir -p $path_output/$command[$i]/logs";
				print_system($cmd);
				
				my $cmd = "$path_scripts/step1/step1.pl --sample $sample --paramfile $pfile --outputdir $path_output/step1 --logs $path_output/step1/logs --timestamp $start_date --R1 $abs_r1 --R2 $abs_r2";
				verbose_system($cmd);
				print "\n";
						
				$out_r1 = "$path_output/$command[$i]/$command[$i].R1" if ( -e "$path_output/$command[$i]/$command[$i].R1" );	
				$out_r2 = "$path_output/$command[$i]/$command[$i].R2" if ( -e "$path_output/$command[$i]/$command[$i].R2" );
			}
		}
		elsif ( $command[$i] eq "quality_filter" )
		{
			print "[module] $command[$i] \n";
			$out_r2 = &check_exist($out_r2, "$path_output/$command[$i]/$command[$i].R2");
			$boolrun = 0;
			$out_r1 = &check_exist($out_r1, "$path_output/$command[$i]/$command[$i].R1");
			print "\n" if (!($boolrun));
			
			# if output doesnt exist, run module
			if ($boolrun)
			{			
				# if more than one module, output from prev stage becomes input for next stage
				if (scalar(@command)>1 && $i>0) {$abs_r1=$out_r1; $abs_r2=$out_r2;}
				
				my $cmd = "mkdir -p $path_output/$command[$i]/logs";
				print_system($cmd);
					
				my $cmd = "$path_scripts/quality_filter/quality_filter.pl --sample $sample --paramfile $pfile --outputdir $path_output/quality_filter --logs $path_output/quality_filter/logs --timestamp $start_date --R1 $abs_r1 --R2 $abs_r2";
				verbose_system($cmd);
				print "\n";

				$out_r1 = "$path_output/$command[$i]/$command[$i].R1" if ( -e "$path_output/$command[$i]/$command[$i].R1" );	
				$out_r2 = "$path_output/$command[$i]/$command[$i].R2" if ( -e "$path_output/$command[$i]/$command[$i].R2" );								
			}
			
			# the output of this stage, for the next stage, is not well paired so set this:
			$wellpaired = 0;	
		}		
		elsif ( $command[$i] eq "host_map" || $command[$i] =~ m/^host_map_(\d){1}$/ )
		{
			my ($num); # run_iteration
			if ($command[$i] =~ m/^host_map_(\d){1}$/ )
			{
				$num=$1;
			}
			else
			{
				$num=1;				
			}
			
			print "[module] host_map\n";
			print "[iteration] $num\n";		
			$out_r2 = &check_exist($out_r2, "$path_output/host_map_$num/host_map_$num.R2");
			$boolrun = 0;
			$out_r1 = &check_exist($out_r1, "$path_output/host_map_$num/host_map_$num.R1");
			print "\n" if (!($boolrun));								
		
			# if output doesnt exist, run module
			if ($boolrun)
			{
				# if more than one module, output from prev stage becomes input for next stage
				if (scalar(@command)>1 && $i>0) {$abs_r1=$out_r1; $abs_r2=$out_r2;}
				
				# make module's logs directory			
				my $cmd = "mkdir -p $path_output/host_map_$num/logs";
				print_system($cmd);
						
				my $cmd = "$path_scripts/host_map/host_map.pl --sample $sample --paramfile $pfile --outputdir $path_output/host_map_$num --logs $path_output/host_map_$num/logs --timestamp $start_date --R1 $abs_r1 --R2 $abs_r2 --fastafile $is_fasta --wellpaired $wellpaired --run_iteration $num";
				verbose_system($cmd);
				print "\n";
				
				$out_r1 = "$path_output/host_map_$num/host_map_$num.R1" if ( -e "$path_output/host_map_$num/host_map_$num.R1" );	
				$out_r2 = "$path_output/host_map_$num/host_map_$num.R2" if ( -e "$path_output/host_map_$num/host_map_$num.R2" );		
			}
			
			# the output of this stage, for the next stage, is not well paired so set this:
			$wellpaired = 0;			
		}
		elsif ( $command[$i] eq "orf_filter" )
		{
			print "[module] $command[$i] \n";
			$out_r2 = &check_exist($out_r2, "$path_output/$command[$i]/$command[$i].R2");
			$boolrun = 0;
			$out_r1 = &check_exist($out_r1, "$path_output/$command[$i]/$command[$i].R1");
			print "\n" if (!($boolrun));
			
			# if output doesnt exist, run module
			if ($boolrun)
			{				
				# if more than one module, output from prev stage becomes input for next stage
				if (scalar(@command)>1 && $i>0) {$abs_r1=$out_r1; $abs_r2=$out_r2;}
				
				# make module's logs directory
				my $cmd = "mkdir -p $path_output/$command[$i]/logs";
				print_system($cmd);
						
				my $cmd = "$path_scripts/orf_filter/orf_filter.pl --sample $sample --paramfile $pfile --outputdir $path_output/$command[$i] --logs $path_output/$command[$i]/logs --timestamp $start_date --R1 $abs_r1 --R2 $abs_r2 --fastafile $is_fasta";
				verbose_system($cmd);
				print "\n";

				$out_r1 = "$path_output/$command[$i]/$command[$i].R1" if ( -e "$path_output/$command[$i]/$command[$i].R1" );	
				$out_r2 = "$path_output/$command[$i]/$command[$i].R2" if ( -e "$path_output/$command[$i]/$command[$i].R2" );												
			}
			
			# the output of this stage, for the next stage, is a fasta file so set this:
			$is_fasta = "yes";	
		}
		elsif ( $command[$i] eq "nohost_blast" )
		{
			print "[module] $command[$i] \n";
			$out_r2 = &check_exist($out_r2, "$path_output/$command[$i]/$command[$i].R2");
			$boolrun = 0;
			$out_r1 = &check_exist($out_r1, "$path_output/$command[$i]/$command[$i].R1");
			print "\n" if (!($boolrun));
			
			# if output doesnt exist, run module
			if ($boolrun)
			{				
				# if more than one module, output from prev stage becomes input for next stage
				if (scalar(@command)>1 && $i>0) {$abs_r1=$out_r1; $abs_r2=$out_r2;}
				
				# make module's logs directory
				my $cmd = "mkdir -p $path_output/$command[$i]/logs";
				print_system($cmd);
						
				my $cmd = "$path_scripts/nohost_blast/nohost_blast.pl --sample $sample --paramfile $pfile --outputdir $path_output/nohost_blast --logs $path_output/nohost_blast/logs --timestamp $start_date --R1 $abs_r1 --R2 $abs_r2 --fastafile $is_fasta";
				verbose_system($cmd);
				print "\n";

				$out_r1 = "$path_output/$command[$i]/$command[$i].R1" if ( -e "$path_output/$command[$i]/$command[$i].R1" );	
				$out_r2 = "$path_output/$command[$i]/$command[$i].R2" if ( -e "$path_output/$command[$i]/$command[$i].R2" );												
			}
			
			# the output of this stage, for the next stage, is a fasta file so set this:
			$is_fasta = "yes";	
		}
		elsif ( $command[$i] eq "iterative_blast_phylo" || $command[$i] =~ m/^iterative_blast_phylo_(\d){1}$/ )
		{
			my ($num); # run_iteration
			if ($command[$i] =~ m/^iterative_blast_phylo_(\d){1}$/ )
			{
				$num=$1;
			}
			else
			{
				$num=1;				
			}
						
			print "[module] iterative_blast_phylo\n";
			print "[iteration] $num\n";
			$out_r2 = &check_exist($out_r2, "$path_output/iterative_blast_phylo_$num/iterative_blast_phylo_$num.R2");
			$boolrun = 0;
			$out_r1 = &check_exist($out_r1, "$path_output/iterative_blast_phylo_$num/iterative_blast_phylo_$num.R1");
			print "\n" if (!($boolrun));
			
			# if output doesnt exist, run module
			if ($boolrun)			
			{
				# if more than one module, output from prev stage becomes input for next stage
				if (scalar(@command)>1 && $i>0) {$abs_r1=$out_r1; $abs_r2=$out_r2;}
				
				# make module's logs directory			
				my $cmd = "mkdir -p $path_output/iterative_blast_phylo_$num/logs";
				print_system($cmd);
				
				my ($command_prefix);
				if ($sge)
				{
					$command_prefix="sge_iterative_blast_phylo";
				}
				else
				{
					$command_prefix="iterative_blast_phylo";
				}
						
				my $cmd = "$path_scripts/$command_prefix/$command_prefix.pl --sample $sample --paramfile $pfile --outputdir $path_output/iterative_blast_phylo_$num --logs $path_output/iterative_blast_phylo_$num/logs --timestamp $start_date --R1 $abs_r1 --R2 $abs_r2 --fastafile $is_fasta --run_iteration $num --contig $boolcontigname";
				verbose_system($cmd);
				print "\n";
				
				$out_r1 = "$path_output/iterative_blast_phylo_$num/iterative_blast_phylo_$num.R1" if ( -e "$path_output/iterative_blast_phylo_$num/iterative_blast_phylo_$num.R1" );	
				$out_r2 = "$path_output/iterative_blast_phylo_$num/iterative_blast_phylo_$num.R2" if ( -e "$path_output/iterative_blast_phylo_$num/iterative_blast_phylo_$num.R2" );
			}
			
			# the output of this stage, for the next stage, is a fasta file so set this:
			$is_fasta = "yes";				
					
		}		
		elsif ( $command[$i] eq "ray2_assembly" || $command[$i] =~ m/^ray2_assembly_(\d){1}$/ )		
		{
			my ($num); # run_iteration
			if ($command[$i] =~ m/^ray2_assembly_(\d){1}$/ )
			{
				$num=$1;
			}
			else
			{
				$num=1;				
			}
						
			print "[module] ray2_assembly\n";			
			print "[iteration] $num\n";
			$boolrun = 0;			
			$out_r1 = &check_exist($out_r1, "$path_output/ray2_assembly_$num/ray2_assembly_$num.fasta");
			print "\n" if (!($boolrun));
			
			# if output already exists, amend out_r2 to be "none", else run module as usual
			if (!($boolrun))
			{
				$out_r2 = "none";
			}
			else
			{
				# if more than one module, output from prev stage becomes input for next stage
				if (scalar(@command)>1 && $i>0) {$abs_r1=$out_r1; $abs_r2=$out_r2;}
				
				# make module's logs directory			
				my $cmd = "mkdir -p $path_output/ray2_assembly_$num/logs";
				print_system($cmd);
					
				my $cmd = "$path_scripts/ray2_assembly/ray2_assembly.pl --sample $sample --paramfile $pfile --outputdir $path_output/ray2_assembly_$num --logs $path_output/ray2_assembly_$num/logs --timestamp $start_date --R1 $abs_r1 --R2 $abs_r2 --fastafile $is_fasta --run_iteration $num";
				verbose_system($cmd);
				print "\n";
				
				# assembly produces a single output file:
				$out_r1 = "$path_output/ray2_assembly_$num/ray2_assembly_$num.fasta" if ( -e "$path_output/ray2_assembly_$num/ray2_assembly_$num.fasta" );			
				$out_r2 = "none";
			}
			
			# the output of this stage, for the next stage, is a fasta file so set this:
			$is_fasta = "yes";				
			# use name contig instead of R1 in iterative blast
			$boolcontigname = 1;
		}								
		# ... ADD NEW MODULES HERE ! ...
		else
		{
			print "[error] command \"".$command[$i]."\" not found. Continuing to run ...\n\n";
		}
	}
	else
	{
		print "[error] parameter file \"".$pfile."\" not found. Continuing to run ...\n\n";
	}
}

$end_date = get_date();
print "[END] $end_date";

# get main ending time in seconds
$main_end_time = time();
print "\n[DELTAT] ",$main_end_time-$main_start_time,"\n";

# ------------------------------------------------------------------------------------------------
# SUBROUTINES

# ---------------------- print title  ----------------------
sub echo_title
{
	print "\n";
	print "#########################################################################################","\n";
	print "############################## PATHOGEN DISCOVERY PIPELINE ##############################","\n";	
	print "#########################################################################################","\n";
	print "\n";	
}

# -------------------- check for error  --------------------
# check for various errors
sub check_error 
{
	if (!(defined($outputdir)))
	{
		print "\n[error] no output directory specified\n\n"; 
		pod2usage(1);		
	}
	if (!( -e $outputdir ))
	{
		print "\n[error] output directory does not exist\n\n"; 
		pod2usage(1);		
	}	
	if (!(defined(@command)))
	{
		print "\n[error] no command option specified\n\n"; 
		pod2usage(1);		
	}
	if (!(defined(@paramfile)))
	{
		print "\n[error] no param file specified\n\n"; 
		pod2usage(1);			 		
	}
	if ( scalar(@paramfile) > 1 && scalar(@paramfile) != scalar(@command) )
	{
		print "\n[error] You must specify EITHER a single parameter file OR a parameter file for every command\n\n"; 
		pod2usage(1);			 		
	}
	if ( $path_scripts eq $path_output )
	{
		print "\n[error] You cannot run this project in the scripts directory. Choose a different output directory\n\n"; 
		pod2usage(1);			 		
	}
	if (!(eval "require Parse_ParameterFile"))
	{
		print "\n[error] The required perl modules cannot be found. Make sure you have added the script\'s local modules to your \$PERL5LIB path and exported it.\n\nPlease run:\nPERL5LIB=\$PERL5LIB:$path_scripts/Local_Module; export PERL5LIB\n\n";
		exit; 
	}
		
}

# -------------------- check for file existence  --------------------
# check if file exists (if so, return it; if not, return orignal file)
sub check_exist 
{
	# When the subroutine is called any parameters are passed as a list in the special @_ (array form of $_)
	# @_ exists magically - anything passed in

	my $orig = shift;	# original file
	my $file = shift;  	# file whose existence is checked
	
	if ( -e $file )
	{
		print "[file exists] $file \n";
		return $file;
	}
	else
	{	
		# output of command doesnt not already exist, so set global var boolrun to true to run the command
		$boolrun = 1;
		return $orig;
	}		
}

# -------------------- get current date -------------------			
# get current date
sub get_date 
{
	(my $sec, my $min, my $hour, my $mday, my $mon, my $year, my $wday, my $yday, my $isdst)=localtime(time);

	# return string in the form, e.g., 20121204-14.36	
	return sprintf("%4d%02d%02d-%02d.%02d", $year+1900,$mon+1,$mday,$hour,$min);
}

# ------------------------ Perl pod -----------------------

# this produces: the help and man page

__END__

=head1 NAME

pathogen.pl - run modules in the pathogen discovery pipeline 

=head1 SYNOPSIS

pathogen.pl [B<--sample> sample] [B<--outputdir> output_directory] [B<--command> command(s)] [B<--paramfile> parameter_file(s)] [options]

=head1 ARGUMENTS

=over 8

=item B<--sample>

Sample name.

=item B<--outputdir>

The output directory of the project. All of the output will be written in this directory. If you give a directory that does not exist, the program will create it.

=item B<--command>

The command(s) (or modules) you want to run.

=item B<--paramfile>

The parameter file for the command. If you are running multiple commands, you can either specify a parameter file for each command (in the same order as the commands) or a single parameter file for all of them.

=back

=head1 OPTIONS

=over 8

=item B<--R1>

.fastq or .sff input file. If paired reads, this file represents the R1 reads. Can be zipped or not. Required only for the "initial_mapping" command.

=item B<--R2>

.fastq or .sff input file of the R2 reads. Can be zipped or not. Required only for the "initial_mapping" command with paired reads.

=item B<--fasta>

The program should know what your input files are, but use this in the special case you want to override it and are using fasta input files.

=item B<--SGE>

turn on SGE parallelization (i.e., qsub jobs) for certain modules.

=item B<--example>

Prints example parameter files for each command.

=item B<--key>

Prints a list of the output files for each command.

=item B<--help>

Prints the help message and exits.

=item B<--man>

Prints the full manual page (including examples) and exits.

=item B<--checkerror>

For a run that's already complete: check the logs files. Must specify outputdir for this option.

=item B<--verbose>

turn on verbose mode

=item B<--cleanup>

Remove intermediate files after a run has finished. Must specify outputdir for this option.

=back

=head1 DESCRIPTION

pathogen.pl runs modules in the pathogen discovery pipeline. It also automates organizing directories, parameter files, and logs files for the modules.

The following modules are available:

B<step1> maps the fastq IDs into simple numerical IDs and processes the .sff file(s) if 454 input files. step1 must be run first. Hence, it is titled step1. 

B<quality_filter> does quality filtering. Specifically, it runs up to two iterations of cutadapt, if specified, and one of prinseq, if specified. 

B<host_map> perform a chain of alignments with a choice of bwa or bowtie2. At each stage, collect the sequences that didn't map, and use these to map to the next reference.
 
B<nohost_blast> for filtering a sample without a known reference. Blasts a small subset of your input file to NCBI NT, creates a new blast db based on the taxid of what it hits, and blasts your whole input file to this new db.

B<iterative_blast> choose an arbitrary number of blast references and algorithmic flavors. Blast to the first reference. Collect the sequences that didn't blast. Use these to blast to the next reference, and so on.

B<iterative_blast_phylo> the same as "iterative_blast" but only for use with the NCBI NT and NR blast databases. This module also annotates the blast output with taxid and makes counts of superkingdom, class, order, genus, species etc.

B<ray2_assembly> do an assembly with ray2.

NOTE: The commands "host_map", "iterative_blast", "iterative_blast_phylo", "ray2_assembly" can be run up to 10 times. If you wish to run a second iteration of "iterative_blast", for example, use the command "iterative_blast_2" and make a corresponding entry in your parameter file.

=head1 EXAMPLE

run "host_map" with paired end .fastq files:

B<pathogen.pl --sample sample --command step1 host_map --paramfile pathogen.param --outputdir /my/output/directory --R1 input_R1.fastq --R2 input_R2.fastq>

run "host_map" with single end .sff files:

B<pathogen.pl --sample sample --command step1 host_map --paramfile pathogen.param --outputdir /my/output/directory --R1 input.sff>

run "host_map" followed by "iterative_blast" followed by "ray2_assembly" (if you've already run "step1" in the same outputdir, you need not give inputs):

B<pathogen.pl --sample sample --command step1 host_map iterative_blast ray2_assembly --paramfile pathogen.param --outputdir /my/output/directory>

run "quality_filter" followed by "host_map":

B<pathogen.pl --sample sample --command step1 quality_filter host_map --paramfile pathogen.param --outputdir /my/output/directory>

run "host_map" followed by "quality_filter" followed by a second iteration of "host_map":

B<pathogen.pl --sample sample --command step1 host_map quality_filter host_map_2 --paramfile pathogen.param --outputdir /my/output/directory>

run "host_map" followed by "iterative_blast" using a separate parameter file for each:

B<pathogen.pl --sample sample --command step1 host_map iterative_blast --paramfile step1.param host_map.param iterative_blast.param --outputdir /my/output/directory>

run "iterative_blast" as a stand alone module:

B<pathogen.pl --sample sample --command iterative_blast --paramfile pathogen.param --outputdir /my/output/directory --R1 input_R1.fastq --R2 input_R2.fastq>

run "nohost_blast":

B<pathogen.pl --sample sample --command step1 nohost_blast --paramfile nohost_blast.param --outputdir /my/output/directory --R1 input_R1.fastq --R2 input_R2.fastq>

=head1 MORE INFOMRATION

This pipeline aspires to be modular. "step1", however, must always be run first. If you want to run a number of modules, you can string them together in any order:

pathogen.pl --sample sample --command step1 host_map iterative_blast ray2_assembly --paramfile pathogen.param --outputdir /my/output/directory --R1 input_R1.fastq --R2 input_R2.fastq

What's happening here? "host_map" maps your reads to a host (or a number of hosts), and what doesn't map becomes the input for "iterative_blast". What doesn't blast, in turn, becomes the input for "ray2_assembly".

You can also run the modules one at a time as follows:

pathogen.pl --sample sample --command step1 host_map --paramfile pathogen.param --outputdir /my/output/directory --R1 input_R1.fastq --R2 input_R2.fastq

pathogen.pl --sample sample --command host_map iterative_blast --paramfile pathogen.param --outputdir /my/output/directory

pathogen.pl --sample sample --command iterative_blast ray2_assembly --paramfile pathogen.param --outputdir /my/output/directory

Note in the second case and third cases, the previous module was also included in the command chain. Why? The program will not re-run a module after its output directory has already been created, but it needs to know what the input is for the current module so, even if you're only running "iterative_blast", using "--command host_map iterative_blast" tells the program that the input for iterative_blast comes from the previous stage, "host_map". If you want to re-run a module over again, delete or rename its output directory.

Here's another example: let's add another blast stage onto the end of what we have above. This would be run as:

pathogen.pl --sample sample --command step1 host_map iterative_blast ray2_assembly iterative_blast_2 --paramfile pathogen.param --outputdir /my/output/directory --R1 input_R1.fastq --R2 input_R2.fastq

Note that we've run "iterative_blast_2". You can add "_n" (where n is a number between 1 and 9) to commands like "iterative_blast" and "host_map" to run them an nth time.

This scripts serves as a wrapper for many programs, which themselves have numerous parameters. Because it would be too cumbersome to pass all of these things from the command line, it uses a parameter file to hold all of these options and paths. You should modify this to suit your purposes. Run B<pathogen.pl --example> to get a template parameter file.

=head1 SAMPLE OUTPUT

Here's sample output for the case where you map, quality filter, then map again:

$ pathogen.pl --sample mysample --command step1 host_map quality_filter host_map_2 --paramfile param.txt --outputdir /my/path --R1 /my/data/1.fastq --R2 /my/data/2.fastq

[module] step1

[cmd] mkdir -p /my/path/step1/logs

[cmd] /my/scripts/step1/step1.pl --sample mysample --paramfile param.txt --outputdir /my/path/step1 --logs /my/path/step1/logs --timestamp 20130112-15.59 --R1 /my/data/1.fastq --R2 /my/data/2.fastq

[deltat] 100


[module] host_map

[iteration] 1

[cmd] mkdir -p /my/path/host_map_1/logs

[cmd] /my/scripts/host_map/host_map.pl --sample mysample --paramfile param.txt --outputdir /my/path/host_map_1 --logs /my/path/host_map_1/logs --timestamp 20130112-15.59 --R1 /my/path/step1/step1.R1 --R2 /my/path/step1/step1.R2 --fastafile no --wellpaired 1 --run_iteration 1

[deltat] 200


[module] quality_filter

[cmd] mkdir -p /my/path/quality_filter/logs

[cmd] /my/scripts/quality_filter/quality_filter.pl --sample mysample --paramfile param.txt --outputdir /my/path/quality_filter --logs /my/path/quality_filter/logs --timestamp 20130112-15.59 --R1 /my/path/host_map_1/host_map_1.R1 --R2 /my/path/host_map_1/host_map_1.R2

[deltat] 50


[module] host_map

[iteration] 2

[cmd] mkdir -p /my/path/host_map_2/logs

[cmd] /my/scripts/host_map/host_map.pl --sample mysample --paramfile param.txt --outputdir /my/path/host_map_2 --logs /my/path/host_map_2/logs --timestamp 20130112-15.59 --R1 /my/path/quality_filter/quality_filter.R1 --R2 /my/path/quality_filter/quality_filter.R2 --fastafile no --wellpaired 0 --run_iteration 2

[deltat] 50

=head1 DIAGNOSTICS

The output of each module, saved in its own logs directory, employs the following tags for easy interpretation:

[module] the name of a command

[echo] a message or comment

[cmd] a line of code the program is running

[error] an error

[deltat] the time a something took to run in seconds

[file exists] a message to note a file has been found and will not be overwritten

=head1 AUTHOR

Written by 

Oliver Elliott (oe2118@columbia.edu), Jason Ladner (jason.t.ladner.ctr@us.army.mil), Francesco Abate (francesco.abate@gmail.com)

=head1 VERSION

Version 1.0; January 3, 2013

=cut
