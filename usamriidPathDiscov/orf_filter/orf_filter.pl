#!/usr/bin/env perl 
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

# example
# orf_filter.pl --sample $sample --paramfile $pfile --outputdir $path_output/orf_filter --logs $path_output/orf_filter/logs --R1 input_R1.fasta --R2 input_R2.fasta --timestamp 23672

# declare vars
my ($outputdir);				# output dir
my ($logs);					# logs path
my ($pfile);					# parameter file 
my ($r1);					# sample R1
my ($r2);					# sample R2
my ($sample);					# sample name
my $timestamp="0";				# time stamp (default is 0)
my ($output_r1);				# output R1
my ($output_r2);				# output R2
my $path_scripts=$RealBin;			# get abs path to directory in which script resides (where to look for sister scripts)
my @mates=("R1","R2");				# strings for mates
my $command="orf_filter";			# set command
my $fastafile=1;				# is fasta


GetOptions (	'outputdir=s' => \$outputdir,	# outputdir
		'logs=s' => \$logs,		# logs
		'paramfile=s' => \$pfile,	# parameter
		'fastafile=i' => \$fastafile,	# is fasta 
		'R1=s' => \$r1,			# R1
		'R2=s' => \$r2,			# R2
		'sample=s' => \$sample,		# sample
		'timestamp=s' => \$timestamp);	# time stamp
            
die "[error] required input parameters not found" if (!( defined($outputdir) && defined($logs) && defined($pfile) && defined($r1) ));
die "[error] input must be fasta" if (!( $fastafile ));
            
# get hash ref
my $href = &parse_param($pfile);
# get hash of parameters
my %hoh=%$href;

# --------------------------------------
# param 

# allowed parameters in hash after parsing: (modify for your module)
# command
# cutadapt_options
# cutadapt_options2
# prinseq_options

# --------------------------------------

# get absolute path to outputdir
$outputdir=abs_path($outputdir);

# open output and error logs 
open my $olog, '>', $logs."/".$sample.".".$timestamp."-out.o";
open my $elog, '>', $logs."/".$sample.".".$timestamp."-out.e";

# redirect standard output and error into logs
open STDOUT, '>&', $olog;
open STDERR, '>&', $elog;

if ( $r1 ne "none" && defined($r1) )
{
	$hoh{$command}{"R1"}=abs_path($r1);
}

if ( $r2 ne "none" && defined($r2) )
{
	$hoh{$command}{"R2"}=abs_path($r2);
} 	

print "[START]\n";
print "[hash] ";
print Dumper \ %hoh;

chdir($outputdir) or die "[error] cannot chdir to $outputdir";

# --------------------------------------
# main

foreach my $mate (@mates) 
{
	# check if defined, and non-zero
	if ( defined($hoh{$command}{$mate}) && -s $hoh{$command}{$mate} )
	{
		# system("ln -sf $hoh{$command}{$mate} $command.$mate");		
		system("ln -sf $hoh{$command}{$mate}");		

		if ($hoh{$command}{"getorf_options"})
		{
			print "[echo] get orf $mate \n";
			# E.g., getorf -sequence 3.contig.noblast.fasta -outseq testtest -minsize 60 -find 0 
			my $cmd = "getorf -sequence $hoh{$command}{$mate} -outseq $mate.orfout.fa $hoh{$command}{\"getorf_options\"}";	
			verbose_system($cmd);

			# my $cmd = "cat $mate.orfout.fa | $RealBin/../scripts/fastajoinlines > $mate.orfout.join.fa";
			# verbose_system($cmd);
			
			# my $href = fastaid_firstword_hash("$mate.orfout.fa");
			# print Dumper $href;

			my %h_fasta = map {/>(\w+)_(\w+)\s(.*)/; $1 => 1} split(/\n/, `cat $mate.orfout.fa`);
			# print Dumper \ %h_fasta;

			print "[echo] filter $hoh{$command}{$mate} by orfs in $mate.orfout.join.fa\n";
			get_subset_by_fastaid($hoh{$command}{$mate}, "$command.$mate", \%h_fasta);
		}
	} # defined
	elsif ($mate eq "R1")
	{
		print "input not defined\n";
	}
} # mate

print "[END]\n";
# print "[hash] ";
# print Dumper \ %hoh;

close(STDOUT);
close(STDERR);  
close($olog);
close($elog);

# return hash of the IDs of a fasta file to "1"
sub fastaid_firstword_hash  
{
	# arg 1 - file name

	my $infile = shift;
	
	my %h = ();				# hash	

	if ( -s $infile )		# if file nonzero
	{	
		open(my $fh, '<', $infile);
	
		while (<$fh>)
		{
			# ID is every line starting with ">"
			if ( $_ =~ m/^>/ )
			{
				chomp $_;
				
				# header looks like >c3_1 [3 - 89] 
				if ($_ =~ m/>(\S+)(\s+)(.*)/)
				{
					# don't print leading ">"
					my $key = $1;
							
					# hash key is fastq id
					$h{$key} = 1;
				}	
			}
		}
			
		close($fh);
	}

	return \%h;
}

# get subset of file
sub get_subset_by_fastaid
{
	my $infile = shift;
	my $outfile = shift;
	my $hash_fasta_ref = shift;
	
	my %h = %$hash_fasta_ref;	# hash	

	# print Dumper \ %h;

	if ( -s $infile )		# if file nonzero
	{	
		open(my $fh, '<', $infile);
		open(my $fh2, '>', $outfile);
	
		while (<$fh>)
		{
			# ID is every line starting with ">"
			if ( $_ =~ m/^>/ )
			{
				chomp $_;
				
				if ($_ =~ m/>(\S+)(.*)/)
				{
					# don't print leading ">"
					my $key = $1;
					# print ($key,"\n");			
					if ($h{$key})
					{
						print $fh2 $_,"\n";
						$_ = <$fh>;
						print $fh2 $_;
					}
				}		
			}
		}
			
		close($fh);
		close($fh2);
	}

	return \%h;
}
