#!/usr/bin/perl

# given a list of taxids you pass as an argument, this will output all the NCBI GIs for those taxid s

# $ARGV[0] list of taxids

# to run, e.g.,
# blastdbcmd -db /data/db/nt -entry all -outfmt '%g %T' > gi2taxid.txt
# cat gi2taxid.txt | get_gi_from_taxid.pl my_taxid_list.txt

use Data::Dumper;

open my $infile1, '<', $ARGV[0];

%myh = ();

# hash your taxids	
while (<$infile1>)
{
	chomp($_);
		
	# this ID will be the key to your hash
	my $mykey=$_;
	
	# map the ID to a four line chunk
	$myh{$mykey}=1;
}

# print Dumper \ %myh;

while (<STDIN>)
{
	my @myrow=split;
	
	# if taxid exists in the hash, print GI number
	print $myrow[0],"\n" if $myh{$myrow[1]};
}