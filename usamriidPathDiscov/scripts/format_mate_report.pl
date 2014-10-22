#!/usr/bin/perl

use strict;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;

# add info about the contigs in mate file

# example
# format_mate_report.pl --outputdir tmp --projdir . --prefix sample --report1 report1.txt --report2 report2.txt --command_list host_map,quality_filter

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
my ($outputdir, $blastdir, $blist, $prefix, $bigreport);

GetOptions ('outputdir=s' => \$outputdir,		# output directory
			'report1=s' => \$report1,			# report 1 file
			'report2=s' => \$report2,			# report 2 file
			'bigreport' => \$bigreport,			# a boolean to control whether or not to make the full report (it's large - usu unnec - the top report should suffice)
			'prefix=s' => \$prefix,				# output files name prefix
            'command_list=s' => \$clist,		# list of module dirs in which to check for filtered reads 
			'projdir=s' => \$projdir);			# output directory
            
die "[error] required input parameters not found" if (!( defined($outputdir) && defined($blastdir) && defined($blist) ));
           
# get abs path to directory in which script resides:
# my $mytmp=`readlink -m $0`; 					# perl's readlink function is lousy, so use bash
# my $path_scripts=`dirname $mytmp`; 
# chomp $path_scripts;

# get absolute paths
my $mytmp=`readlink -m $outputdir`;
chomp $mytmp;
$outputdir = $mytmp;

my $mytmp=`readlink -m $projdir`;
chomp $mytmp;
$projdir = $mytmp;

my $mytmp=`readlink -m $report1`;
chomp $mytmp;
$report1 = $mytmp;

my $mytmp=`readlink -m $report2`;
chomp $mytmp;
$report2 = $mytmp;

my @command_list = split(',',$clist);

chdir($outputdir) or die "[error] cannot chdir to ".$outputdir;

# --------------------------------------

open(my $outfile, ">", "$outputdir/$prefix.R1.R2.report.txt");

# map taxid to phylo stuff
open(my $infile1, "<", "$report1"); 				
my %h1=();
my $skipline=1;	# skip the first header line
while (<$infile1>)
{
	print $outphylo $blast_task_list[$i],"\t",$_ if ($skipline==0);					
	$skipline=0;
		
	# $infile3 - e.g., 
	# taxid	count	superkingdom	kingdom	class	order	family	genus	species
	# 565995	49.00	Viruses	-	-	Mononegavirales	Filoviridae	Ebolavirus	Bundibugyo ebolavirus	
	# 11036	3.00	Viruses	-	-	-	Togaviridae	Alphavirus	Venezuelan equine encephalitis virus

	my @ln=split; 
	# map: taxid --> superkingdom order family genus 
	$h3{$ln[0]}=$ln[2]."\t".$ln[5]."\t".$ln[6]."\t".$ln[7];			
}
close($infile1);	

close($outfile)