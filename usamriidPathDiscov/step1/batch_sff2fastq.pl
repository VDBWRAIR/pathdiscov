#!/usr/bin/perl

# batch convert sff to fastq. You can use in the following ways:
# $ ./batch_sff2fastq.pl 1.sff 2.sff
# $ ./batch_sff2fastq.pl *.sff
# $ ./batch_sff2fastq.pl --list list.txt

# dependencies:
# fastq_convert.pl
# sff2fastq.sh
# (these must be in the directory from which this script is run)

use Getopt::Long;

GetOptions ('list=s' => \$list);	# if list option use a test file list of .sff files

# get path to scripts
$tmp=`readlink -m $0`; 
$path_scripts=`dirname $tmp`; 
chomp $path_scripts;

$aref = \@ARGV;				# aref will hold the reference to either the array of arguments OR the array of of files given in a list (if it's provided) 

if (defined($list))
{
	open my $infile1, '<', $list;	# inputfile ( a list of files)

	@filelist = ();			# initialize empty array

	# loop through list - store in @filelist
	while (<$infile1>)
	{
		chomp;
		push(@filelist, $_);
	}

	close $infile1;

	$aref = \@filelist;		# reset array reference pointer

}

# loop through either argv OR list (if defined)
foreach (@{$aref})
{
	if ( $_ =~ m/\.sff$/)
	{
		$myfile=$_;
		$myfastq=$myfile;
		$myfastq =~ s/\.sff/\.fastq/;
		# print "converting $myfile to $myfastq \n";
		my $cmd = "$path_scripts/sff2fastq.sh $myfile $myfastq";
		print "[cmd] ",$cmd,"\n";
		system($cmd);
	}
	else
	{
		print "$_ is not an .sff file \n";
	}
}