#!/usr/bin/perl

use strict;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use Cwd 'abs_path';

use FindBin qw($RealBin);
use lib "$RealBin/../Local_Module";
# local modules:
use Verbose_Sys;
use Parse_ParameterFile;

# TO DO:
# put R1 R2 in loop instead of repeat & make parallel
# links use rel not abs path
# improve syntax of string commands
# add fasta option
# re write file functions in C++

# example
# step1.pl --sample $sample --paramfile $pfile --outputdir $path_output/step1 --logs $path_output/step1/logs --R1 input_R1.fastq --R2 input_R2.fastq --timestamp 23672

# declare vars
my ($outputdir);		# output dir
my ($logs);			# logs path
my ($pfile);			# parameter file 
my ($r1);			# sample R1
my ($r2);			# sample R2
my ($sample);			# sample name
my ($help);			# help bool
my ($man);			# man bool
my ($iszip);			# boolean for zipped
my $timestamp="0";		# time stamp (default is 0)
my $path_scripts=$RealBin;	# get abs path to directory in which script resides (where to look for sister scripts)
my @mates=("R1","R2");		# strings for mates
my $command="step1";		# set command

GetOptions (	'outputdir=s' => \$outputdir,	# outputdir
		'logs=s' => \$logs,		# logs
		'paramfile=s' => \$pfile,	# parameter
		'R1=s' => \$r1,			# R1
		'R2=s' => \$r2,			# R2
		'sample=s' => \$sample,         # sample
		'timestamp=s' => \$timestamp,   # time stamp
		'help' => \$help, 		# help
		'man' => \$man);		# man
            
pod2usage(1) if ($help);
pod2usage(-verbose => 2) if ($man);            

die "[error] required input parameters not found" if (!( defined($outputdir) && defined($logs) && defined($pfile) && defined($r1) ));
            
# get hash ref
my $href = &parse_param($pfile);
# get hash of parameters
my %hoh=%$href;

# --------------------------------------
# param 

# allowed parameters in hash after parsing: (modify for your module)
# command
# seq_platform

# --------------------------------------

# get absolute path to outputdir
$outputdir=abs_path($outputdir);

# open output and error logs 
open my $olog, '>', $logs."/".$sample.".".$timestamp."-out.o";
open my $elog, '>', $logs."/".$sample.".".$timestamp."-out.e";

# redirect standard output and error into logs
open STDOUT, '>&', $olog;
open STDERR, '>&', $elog; 	

print "[START]\n";
print "[hash] ";
print Dumper \ %hoh;

chdir($outputdir) or die "[error] cannot chdir to $outputdir";

# --------------------------------------
# files	produced

# "R1.fastq";		# after turning IDs into numbers
# "R1.id";		# map old IDs to numerical IDs
# "R2.fastq";
# "R2.id";
# "R1.unzip";		# variable for unzipped file if nec
# "R1.unzip.fastq";	# another variable used if seq platform == 454
# "R2.unzip";		
# "R2.unzip.fastq";	
# "R1.count";		# a file to track read counts
# "R2.count";

# --------------------------------------
# main

# unzip convert sff
# the idea here is that $hoh{$command}{$mate} will hold whatever the input is, after the specified level of processing

$hoh{$command}{"R1"}=$r1;
$hoh{$command}{"R2"}=$r2 if (defined($r2) && $r2 ne "none");

foreach my $mate (@mates) 
{
	# check if defined, not "none", and non-zero
	if ( defined($hoh{$command}{$mate}) && $hoh{$command}{$mate} ne "none" && -s $hoh{$command}{$mate} )
	{
		# boolean for zipped
		$iszip=0;

		# if zipped, unzip
		if ($hoh{$command}{$mate} =~ m/\.gz$/)
		{
			$iszip=1;
			
			my $cmd = "cat $hoh{$command}{$mate} | gunzip > $mate.unzip";
			verbose_system($cmd);

			# input file now unzipped	
			$hoh{$command}{$mate} = "$mate.unzip";
		}
		
		# check seq platform: if 454, change sff to fastq
		if ($hoh{$command}{"seq_platform"} eq "454")
		{
			print "[echo] convert sff to fastq $mate\n";
			# args: input output
			my $cmd = "$path_scripts/sff2fastq.sh $hoh{$command}{$mate} $mate.sff.fastq";
			verbose_system($cmd);
			
			$hoh{$command}{$mate} = "$mate.sff.fastq";	
		}
		
		# change IDs to numbers, truncate 3rd line of fastq to "+", replace "." in seq with "N"
		print "[echo] change fastq IDs to numeric IDs $mate\n";							
		# args: function, inputfile, outputfile, key				
		my $cmd = "$path_scripts/perlscripts_wrapper.pl change_fastq_id $hoh{$command}{$mate} $mate.fastq $mate.id";
		verbose_system($cmd);
		
		# rm upzipped copy
		system("rm $mate.unzip") if ($iszip);							
		
		# make link with prefix as command
		system("ln -sf $mate.fastq $command.$mate");	
		
		# count lines in file
		# args: input file, filtering_program_name, output file, isfastq, concat
		my $cmd = "$path_scripts/linecount.sh $mate.fastq rawfile $mate.count 1 0";
		print "[cmd] ",$cmd,"\n";
		system($cmd);
	}	  
}
			
# --------------------------------------

print "[END]\n";
# print "[hash] ";
# print Dumper \ %hoh;

close(STDOUT);
close(STDERR);  
close($olog);
close($elog);

# ------------------------ Perl pod -----------------------

# this produces: the help and man page

__END__

=head1 NAME

step1.pl - converst fastq ids to numbers (and convert sff if 454)

=head1 SYNOPSIS

step1.pl [B<--sample> sample] [B<--outputdir> output_directory] [B<--logs> log_directory] [B<--paramfile> parameter_file(s)]  [B<--R1> R1_input] [B<--timestamp> a time stamp] [options] 

=head1 ARGUMENTS

=over 8

=item B<--sample>

Sample name.

=item B<--outputdir>

The output directory. All of the output will be written in this directory.

=item B<--logs>

The logs directory. stdout and sterr will be written in this folder.

=item B<--paramfile>

The parameter file for the command.

=item B<--R1>

.fastq or .sff input file. If paired reads, this file represents the R1 reads. Can be zipped or not.

=item B<--timestamp>

A time stamp. This ensures your log files will be unique.

=back

=head1 OPTIONS

=over 8

=item B<--R2>

.fastq or .sff input file of the R2 reads. Can be zipped or not. Required only for paired reads.

=item B<--example>

Prints example parameter file.

=item B<--help>

Prints the help message and exits.

=item B<--man>

Prints the full manual page (including examples) and exits.

=back

=head1 DESCRIPTION

step1.pl first takes your input file (.sff or .fastq) (which can be gzipped) and coverts the ugly identifiers into easy-to-manipulate numbers. 

=head1 EXAMPLE

run "step1" with paired end .fastq files:

B<step1.pl --sample sample --paramfile step1.param --outputdir /my/output/directory --logs /my/logs/directory --R1 input_R1.fastq --R2 input_R2.fastq --timestamp 121207-11.30>

=cut
