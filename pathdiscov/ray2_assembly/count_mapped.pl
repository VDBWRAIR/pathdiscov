#!/usr/bin/perl

use strict;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use Cwd 'abs_path';

use FindBin qw($RealBin);
use lib "$RealBin/../lib/Local_Module";
# local modules:
use Verbose_Sys;
use Sam_Functions;

# take a sam file piped in from std:in and count reads mapping to given elt (already in header?)

# example
# samtools view out.bam | count_mapped.pl

# hash
my %h = ();

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

		# if not unmapped
		if (is_unmapped($flag) == 0)
		{
			# if hash entry not exist, make entry
			if (!$h{$sid})
			{
				$h{$sid}=1;
			}
			# else increment
			else
			{
				$h{$sid}++;
			}
		}
	}
}

foreach my $key ( keys %h )
{
	print $key,"\t",$h{$key},"\n";
}
