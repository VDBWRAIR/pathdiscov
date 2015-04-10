#!/usr/bin/perl

# this stores the query IDs of a blast file in a hash table, and uses it to extract whatever IDs of a fasta dont match them (i.e., if this fasta file was the input to the blast, it gets the reads that dont map)

# $ARGV[0]	BLAST output
# $ARGV[1]	fasta input that was used for BLAST

use Data::Dumper;

open my $infile1, '<', $ARGV[0];	# blast output
open my $infile2, '<', $ARGV[1];	# fasta input

# open the second file and create a hash based on IDs

my %myh = ();

# create a hash of query IDs	
while (<$infile1>)
{
	split;							# split $_ on default ws
	$mykey=$_[0];					# get first elt (i.e., query ID)
	$myh{$mykey}=1;					# hash to 1
}

# print Dumper \ %myh;

# get reads from fasta which dont match the query IDs 
# IMPORTANT: this assumes that the fasta sequences are on a single line!
while (<$infile2>) 
{

	chomp($_);						# take off \n	
	
	if ($_ =~ m/>(\S+)(\s+)(\S+)/)	# if id contains ws, take the first part
	{
		# this ID will be the key to your hash
		$id=$1;
	}
	else
 	{
		# else take whole thing
	 	$id=substr($_,1);			# take off leading ">"
	}	
	
	if ($myh{$id})					# if match, do nothing
	{
		$_ = <$infile2>;			# cycle through one more line of the file (the sequence)
	}
	else							# if not match, print id & sequence
	{
		print ">",$id,"\n";
		$_ = <$infile2>;
		print $_;
	}
}

close $infile1;
close $infile2;
