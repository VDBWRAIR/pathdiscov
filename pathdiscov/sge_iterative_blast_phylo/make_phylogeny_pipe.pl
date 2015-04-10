#!/usr/bin/perl

# tally phylogenies from an "augmented" BLAST output file (i.e. --- one with the NCBI taxids added on )

# $ARGV[0]	input file
# $ARGV[1]	ncbi taxonomy nodes.dmp
# $ARGV[2]	ncbi taxonomy names.dmp

# useage -
# cat taxid_file | make_phylogeny_pipe.pl /data/taxonomy/nodes.dmp /data/taxonomy/names.dmp --verbose

# RECALL:
# nodes.dmp
# ---------
# 
# This file represents taxonomy nodes. The description for each node includes 
# the following fields:
# 
# 	tax_id								-- node id in GenBank taxonomy database
#  	parent tax_id						-- parent node id in GenBank taxonomy database
#  	rank								-- rank of this node (superkingdom, kingdom, ...) 
#  	embl code							-- locus-name prefix; not unique
#  	division id							-- see division.dmp file
#  	inherited div flag  (1 or 0)		-- 1 if node inherits division from parent
#  	genetic code id						-- see gencode.dmp file
#  	inherited GC  flag  (1 or 0)		-- 1 if node inherits genetic code from parent
#  	mitochondrial genetic code id		-- see gencode.dmp file
#  	inherited MGC flag  (1 or 0)		-- 1 if node inherits mitochondrial gencode from parent
#  	GenBank hidden flag (1 or 0)      	-- 1 if name is suppressed in GenBank entry lineage
#  	hidden subtree root flag (1 or 0) 	-- 1 if this subtree has no sequence data yet
#  	comments							-- free-text comments and citations
# 
# names.dmp
# ---------
# Taxonomy names file has these fields:
# 
# 	tax_id					-- the id of node associated with this name
# 	name_txt				-- name itself
# 	unique name				-- the unique variant of this name if name not unique
# 	name class				-- (synonym, common name, ...)
# 

use Data::Dumper;
use Getopt::Long;

GetOptions ('verbose' => \$verbose);	# if verbose print using @rank_header else print using @short_header	

open my $infile2, '<', $ARGV[0];		# nodes.dmp
open my $infile3, '<', $ARGV[1];		# names.dmp

@rank_header = ("superkingdom", "kingdom", "phylum", "superclass",
          		"infraclass", "class", "subclass", "superorder",
          		"order", "suborder", "infraorder", "parvorder",
          		"superfamily", "family", "subfamily","tribe",
          		"subtribe", "genus", "subgenus", "species group",
          		"species");
          		
@short_header = ("superkingdom", "kingdom", "class", "order", "family", "genus", "species");

# loop through nodes.dmp
# make a map of tax_id to parent tax_id & tax_id to rank
while (<$infile2>)
{
	chomp;
	$_ =~ s/\t//g;						# eliminate tabs
	my @a = split (/\|/);				# split		
	$parent[$a[0]] = $a[1];				# map tax_id to parent
	$rank[$a[0]] = $a[2];				# map tax_id to rank
}

# loop through names.dmp
# make a map of tax_id to name
while (<$infile3>)
{
	if ($_ =~ m/scientific name/)		# get only scientific name
	{
		chomp;
		$_ =~ s/\t//g;					# eliminate tabs
		my @a = split (/\|/);			# split
		$sciname[$a[0]] = $a[1];		# map tax_id to scientific name
	}	
}

# print header
if ($verbose)
{
	print join("\t",@rank_header),"\n";
}
else
{
	print join("\t",@short_header),"\n";
}

%hoh = ();								# hash of hashes - to store a taxid along with everything in the rank header
										# e.g., hoh{taxid}{superkingdom} = Viruses OR hoh{taxid}{family} = Hominidae 

# loop through input file (taxids)
# for each taxid, make a map of rank --> name, if taxid not yet in hash
while (<STDIN>)
{
	chomp;
	$taxid = $_; 						# get taxid
	my $a = $taxid;	 					# set a = taxid
	# print "tax: ",$taxid,"\n";
	
	# $rank = "species";
	# print "$a[0] \n";

	# accumulate rank information
	if ( !(defined $hoh{$a}) )
	{
		# print "in loop \n";

		$i = 1;
		while ( $i < 32 && $taxid != 1 ) # there's only about 30 possible links max; if taxid == 1, it means we've hit the root and we're at the end
		{
	
			# print "iteration: $i\ttaxid: $taxid\tname: $sciname[$taxid]\tparentid: $parent[$taxid]\trank: $rank[$taxid]\n";
	
			if ( $rank[$taxid] ne "no rank" )
			{
				# print "in loop \n";
				$hoh{$a}{$rank[$taxid]} = $sciname[$taxid];			
			}
			
			$taxid=$parent[$taxid];		# set taxid to tax id of parent
		    $i++;
	    }
	}
	
	# print 
	if ($verbose)
	{
		foreach (@rank_header) 
		{
			if (defined $hoh{$a}{$_})
			{
				print $hoh{$a}{$_},"\t";
			}
			else
			{
				print "-","\t";
			}
		}
		print "\n";
	}
	else
	{
		foreach (@short_header) 
		{
			if (defined $hoh{$a}{$_})
			{
				print $hoh{$a}{$_},"\t";
			}
			else
			{
				print "-","\t";
			}
		}
		print "\n";
	}	
	
}

# print Dumper \ %hoh;

close $infile1;
close $infile2;