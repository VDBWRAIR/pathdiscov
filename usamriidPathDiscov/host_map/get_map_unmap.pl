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
use Sam_Functions;

# take a sam file piped in from std:in and print IDs into various files depending on mate pair and if mapped

# example
# samtools view out.bam | get_map_unmap.pl --R1map R1.map.id --R2map R2.map.id --R1unmap R1.unmap.id --R2unmap R2.unmap.id --singletonmap singleton.map.id --singletonunmap singleton.unmap.id --paired --unpaired

# to do:
# put code into functions

# declare vars
my ($output1);					# R1 paired map
my ($output2);					# R2 paired unmap
my ($output3);					# R1 paired map
my ($output4);					# R2 paired unmap
my ($output5);					# nonpaired map
my ($output6);					# nonpaired unmap
my $boolpair=0;					# bool for paired reads 
my $boolunpair=0;				# bool for unpaired reads
my $boolsid=0;					# bool for printing subject id

GetOptions ('R1map=s' => \$output1,					# R1 paired map
            'R1unmap=s' => \$output2,				# R1 paired unmap
            'R2map=s' => \$output3,					# R2 paired map
            'R2unmap=s' => \$output4,				# R2 paired unmap
            'singletonmap=s' => \$output5,			# singleton map
            'singletonunmap=s' => \$output6,		# singleton unmap
            'printsid' => \$boolsid,				# bool for printing subject ID
            'paired' => \$boolpair,					# bool for paired reads 
            'unpaired' => \$boolunpair);			# bool for unpaired reads
            										# if combining paired and unpaired reads, use both these flags
# --------------------------------------
# main

open my $fh1, '>', $output1 if ($boolpair);		# R1 paired map
open my $fh2, '>', $output2 if ($boolpair);		# R1 paired unmap
open my $fh3, '>', $output3 if ($boolpair);		# R2 paired map
open my $fh4, '>', $output4 if ($boolpair);		# R2 paired unmap
open my $fh5, '>', $output5 if ($boolunpair);	# singleton map
open my $fh6, '>', $output6 if ($boolunpair);	# singleton unmap

# loop thro what's piped in
while (<STDIN>)  
{
	# not header
	if ($_ !~ m/^@/) 
	{		
		chomp $_;		
		my $wholeline = $_;
		my @line = split;

		my $qid=$line[0];		# query id
		my $flag=$line[1];		# flag
		my $sid=$line[2];		# subject id		
	
		# R1
		if (is_firstread($flag) == 1 && $boolpair)
		{
			if (is_unmapped($flag) == 1)
			{
				print $fh2 $qid,"\n";
				# print "id: $qid flag: $flag r1 unmap \n";
			}
			else
			{
				# print "id: $qid flag: $flag r1 map \n";
				if ($boolsid)
				{
					print $fh1 $qid,"\t",$sid,"\n";
				}
				else
				{
					print $fh1 $qid,"\n";
				}
			}
		}
		# R2
		elsif (is_secondread($flag) == 1 && $boolpair)
		{
			if (is_unmapped($flag) == 1)
			{
				print $fh4 $qid,"\n";
				# print "id: $qid flag: $flag r2 unmap \n";
			}
			else
			{			
				# print "id: $qid flag: $flag r2 map \n";
				if ($boolsid)
				{
					print $fh3 $qid,"\t",$sid,"\n";
				}
				else
				{
					print $fh3 $qid,"\n";
				}				
			}
		}
		# singletons	
		elsif (is_firstread($flag) == 0 && is_secondread($flag) == 0 && $boolunpair)
		{
			if (is_unmapped($flag) == 1)
			{
				print $fh6 $qid,"\n";
				# print "id: $qid flag: $flag singleton unmap \n";
			}
			else
			{
				# print "id: $qid flag: $flag singleton map \n";
				if ($boolsid)
				{
					print $fh5 $qid,"\t",$sid,"\n";
				}
				else
				{
					print $fh5 $qid,"\n";
				}			
			}
		}
	#	else
	#	{
	#		print STDERR "[error] cannot interpret flag \n";
	#		exit;
	#	}
	}
}

# close file handles
close $fh1 if ($boolpair);
close $fh2 if ($boolpair);
close $fh3 if ($boolpair);
close $fh4 if ($boolpair);
close $fh5 if ($boolunpair);
close $fh6 if ($boolunpair);
