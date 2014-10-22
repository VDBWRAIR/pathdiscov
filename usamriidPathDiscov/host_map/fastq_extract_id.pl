#!/usr/bin/perl

# this extracts the 4 line chunks from a fastq file, given in argument 1 file, for IDs give in argument 2 file
# for example, your second file might look like this:
# 2
# 4
# 7

open my $infile1, '<', $ARGV[0];	# fastq file
open my $infile2, '<', $ARGV[1];	# file of IDs

# open the second file and create a hash based on IDs

my %myh = ();
	
while (<$infile2>)
{
	chomp($_);
	$tmpvar=$_;

	# if there's the pattern word space word, take only the first word
	if ($tmpvar =~ m/(\S+)(\s+)(\S+)/) 
	{
		# this ID will be the key to your hash
		$mykey="@".$1;
	}
	else
	{
		# else take whole thing
		$mykey="@".$tmpvar;
	}

#	print ("key: ",$mykey,"\n");
		
	# map the ID to 1
	$myh{$mykey}=1;
}

# loop through the first file and, if ID matches hash ID, print 4 lines as well as chunk from the second file

while (<$infile1>) 
{

	chomp($_);
		
	my $line=$_;

	if ($line =~ m/(\S+)(\s+)(\S+)/)
	{
		# this ID will be the key to your hash
		$mykey=$1;
	}
	else
 	{
		# else take whole thing
		$mykey=$line;
	}

	if ($myh{$mykey})
	{
		print $line,"\n";
		$_ = <$infile1>;
		print $_;
		$_ = <$infile1>;
		print $_;
		$_ = <$infile1>;
		print $_;
	}
}

close $infile1;
close $infile2;
