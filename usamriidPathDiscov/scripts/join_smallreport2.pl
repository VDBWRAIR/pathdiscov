#!/usr/bin/perl

use strict;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use Cwd 'abs_path';		# akin to bash readlink

# add info about the contigs in mate file

# example
# join_smallreport.pl --outputdir tmp --prefix sample --R1report report1.txt --R1qualdiscard R1.qual.discard.txt --R1hostdiscard R1.host.discard.txt --R2report report1.txt --R2qualdiscard R2.qual.discard.txt --R2hostdiscard R2.host.discard.txt

# -----------------------------------

# args 
# $ARGV[0] 
# $ARGV[1] 
# ...

# learning perl ! 

# \d [0-9] Any digit
# \D [^0-9] Any character not a digit
# \w [0-9a-zA-z_] Any "word character"
# \W [^0-9a-zA-z_] Any character not a word character
# \s [ \t\n\r\f] whitespace (space, tab, newline, carriage return, form feed)
# \S [^ \t\n\r\f] Any non-whitespace character

# *      Match 0 or more times
# +      Match 1 or more times
# ?      Match 1 or 0 times
# {n}    Match exactly n times
# {n,}   Match at least n times
# {n,m}  Match at least n but not more than m times

# -----------------------------------

# declare vars
my ($outputdir);	# output dir
my ($prefix);		# sample name
my ($report1);		# small report R1
my ($report2);		# small report R2
my ($qdiscard1);	# qual discard R1
my ($qdiscard2);	# qual discard R2
my ($hdiscard1);	# host discard R1
my ($hdiscard2);	# host discard R2

GetOptions ('outputdir=s' => \$outputdir,				# output directory
			'prefix=s' => \$prefix,						# prefix
			'R1report=s' => \$report1,					# R1 report file
			'R2report=s' => \$report2,					# R2 report file
			'R1qualdiscard=s' => \$qdiscard1,			# qual discard R1
			'R2qualdiscard=s' => \$qdiscard2,			# qual discard R2
			'R1hostdiscard=s' => \$hdiscard1,			# host discard R1
			'R2hostdiscard=s' => \$hdiscard2);			# host discard R2			
                       
# get abs paths

$outputdir=abs_path($outputdir);
$report1=abs_path($report1);
$report2=abs_path($report2);
$qdiscard1=abs_path($qdiscard1);
$qdiscard2=abs_path($qdiscard2);
$hdiscard1=abs_path($hdiscard1);
$hdiscard2=abs_path($hdiscard2);

chdir($outputdir) or die "[error] cannot chdir to ".$outputdir;

# --------------------------------------

my (%h1, %h2, %h3, %h4, %h5, %h6);		# hashes for the input files

# get hash reference 
my $href = &hash_file1($report1);
# deref to get hash
%h1=%$href;

my $href = &hash_file1($report2);
# deref to get hash
%h2=%$href;

my $href = &hash_file2($qdiscard1);
%h3=%$href;

my $href = &hash_file2($qdiscard2);
%h4=%$href;

my $href = &hash_file2($hdiscard1);
%h5=%$href;

my $href = &hash_file2($hdiscard2);
%h6=%$href;

# this hash will be the union of keys
my %key_union = ();

# get the union of the sets (i.e., all elts represented uniquely):
foreach my $key ( keys %h1 )
{
	$key_union{$key} = 1;
}

foreach my $key ( keys %h2 )
{
	$key_union{$key} = 1;
}

# print Dumper \ %h1;
# print Dumper \ %h2;
# print Dumper \ %key_union;

my $header="qseqid"."\t"."pair2sametax"."\t"."R1.discard"."\t"."R1.blast_alg"."\t"."R1.qseqid"."\t"."R1.sseqid"."\t"."R1.pident"."\t"."R1.length"."\t"."R1.mismatch"."\t"."R1.gapopen"."\t"."R1.qstart"."\t"."R1.qend"."\t"."R1.sstart"."\t"."R1.send"."\t"."R1.evalue"."\t"."R1.bitscore"."\t"."R1.qlen"."\t"."R1.slen"."\t"."R1.taxid"."\t"."R1.superkingdom"."\t"."R1.order"."\t"."R1.family"."\t"."R1.genus"."\t"."R1.descrip"."\t"."R2.discard"."\t"."R2.blast_alg"."\t"."R2.qseqid"."\t"."R2.sseqid"."\t"."R2.pident"."\t"."R2.length"."\t"."R2.mismatch"."\t"."R2.gapopen"."\t"."R2.qstart"."\t"."R2.qend"."\t"."R2.sstart"."\t"."R2.send"."\t"."R2.evalue"."\t"."R2.bitscore"."\t"."R2.qlen"."\t"."R2.slen"."\t"."R2.taxid"."\t"."R2.superkingdom"."\t"."R2.order"."\t"."R2.family"."\t"."R2.genus"."\t"."R2.descrip"."\n";

open(my $outfile, ">", "$outputdir/$prefix.joinR1R2.smallreport.txt");

print $outfile $header;

my $spacer="\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t";	

my ($outline);		# line of output
my ($sametax);		# bool true if R1 R2 map to same taxid (default: false)
my ($taxid1);		# a holder for taxid

foreach my $key ( keys %key_union )
{
	$outline = "";
	$sametax = 0;
	$taxid1 = 0;

	# R1
	# in report
	if ($h1{$key})
	{
		# dicard column (put 0 for not discarded)
		$outline = $outline."0"."\t".$h1{$key}."\t";

		# get taxid
		my @line = split("\t",$h1{$key}); 
		$taxid1 = $line[15];
	}
	# in qual discards
	elsif ($h3{$key})
	{
		$outline = $outline."qual"."\t".$spacer;	
	}
	# in host discards
	elsif ($h5{$key})
	{
		$outline = $outline."host"."\t".$spacer;	
	}
	# not found	
	else 
	{
		$outline = $outline."notfound"."\t".$spacer;	
	}

	# R2
	# in report	
	if ($h2{$key})
	{
		$outline = $outline."0"."\t".$h2{$key}."\n";
		
		# get taxid and see if it matches R1
		my @line = split("\t",$h2{$key}); 
		$sametax = 1 if ($line[15] == $taxid1);
				
	}
	# in qual discards
	elsif ($h4{$key})
	{
		$outline = $outline."qual"."\t".$spacer."\n";			
	}
	# in host discards
	elsif ($h6{$key})
	{
		$outline = $outline."host"."\t".$spacer."\n";			
	}
	# not found	
	else 
	{
		$outline = $outline."notfound"."\t".$spacer."\n";			
	}
	
	# add qid in front
	$outline = $key."\t".$sametax."\t".$outline;	
	
	# print
	print $outfile $outline;
}

close($outfile);

# -------------------------------------------------------------

# Parse smallreport and return hash of 2nd col => whole line
sub hash_file1  
{
	
	# skip first line and then hash whole line to the 2nd column as keys (which is query IDs)
	
	# When the subroutine is called any parameters are passed as a list in the special @_ (array form of $_)
	# @_ exists magically - anything passed in

	my $infile = shift;  
	open(my $fh, '<', $infile);

	my %h = ();						# hash of hashes - to store {command}{parameters}
	
	my $skipline = 1; 				# skip first line (header)
		
	# parse configuration file
	while (<$fh>)
	{
		if (!$skipline)
		{
			chomp $_;
				
			# line is the whole line
			my $line = $_;
			
			# split on ws
			my @row = split;
			
			# hash key - 2nd column (recall: array index starts from 0)
			$h{$row[1]} = $line;
		}	
		$skipline = 0;
	}

	close($fh);

	return \%h;
}

# Parse discard file and return hash of 1st col => "1"
sub hash_file2  
{
	
	# When the subroutine is called any parameters are passed as a list in the special @_ (array form of $_)
	# @_ exists magically - anything passed in

	my $infile = shift;  
	open(my $fh, '<', $infile);

	my %h = ();						# hash of hashes - to store {command}{parameters}
			
	# parse configuration file
	while (<$fh>)
	{
		chomp $_;
			
		# line is the whole line
		my $line = $_;
		
		# hash key - 2nd column (recall: array index starts from 0)
		$h{$line} = 1;
	}

	close($fh);

	return \%h;
}
