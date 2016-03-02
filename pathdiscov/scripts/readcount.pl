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
use Easy_Hash;

# add info about the contigs in mate file
# does 3 things:
# (1) produces a file of the complete read id mapping
# (2) produces a count report
# (3) makes a graph of this

# example
# readcount.pl --sample sample --outputdir . --projdir ../results --dirlist "step1,quality_filter,host_map_1,ray2_assembly_1,iterative_blast_phylo_1,iterative_blast_phylo_2" --trackread 
# example 2 
# $d/scripts/readcount.pl --sample sample --outputdir . --projdir ../../results/VLP-4 --dirlist "step1,quality_filter,host_map_1,ray2_assembly_1,iterative_blast_phylo_1,iterative_blast_phylo_2" --trackread 

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
my ($outputdir);			# output dir
my ($projdir);				# output dir
my ($dlist);				# comma-delimited list of proj dirs
my ($sample);				# sample name
my ($boolgraph);			# bool make graph
my ($booltrackread);		# bool make full read tracking report
my ($header);				# header for output file

GetOptions ('outputdir=s' => \$outputdir,		# output directory
			'projdir=s' => \$projdir,			# project directory
			'dirlist=s' => \$dlist,				# comma-delimited list of proj dirs
			'sample=s' => \$sample,				# sample
			'graph' => \$boolgraph,				# bool make graph
			'trackread' => \$booltrackread);	# bool make full read tracking report
                       
# get abs paths
$outputdir=abs_path($outputdir);
$projdir=abs_path($projdir);

chdir($outputdir) or die "[error] cannot chdir to ".$outputdir;

# --------------------------------------
# main

my $outfile_R1 = "$outputdir/$sample.R1.count.txt";				# count file for R1
my $outfile_R2 = "$outputdir/$sample.R2.count.txt";				# count file for R2
my $outfile_contig = "$outputdir/$sample.contig.count.txt";		# count file for contigs

my @dir_list = split(',',$dlist);								# list of dirs in the prj folder

my (%h1, %h2);													# hashes R1 and R2 reads step1
my (%h3, %h4);													# hashes R1 and R2 reads qual
my (%h5, %h6);													# hashes R1 and R2 reads host_map
my (%h7, %h8);													# hashes R1 and R2 unassembled reads from assembly
my (%r2c_R1, %r2c_R2, %r2c_single);								# reads 2 contig (from the assembly step)
my (%c2t);														# contig 2 taxid - iterative_blast_phylo_1 step
my (%r2t_R1, %r2t_R2);											# read 2 taxid - iterative_blast_phylo_2 step


# loop over each dir
foreach my $d (@dir_list) 
{						
	if ( $d eq "step1" ) # step1 must be run first ! 
	{
		if ( -e "$projdir/$d/R1.count" )
		{
			system("cat $projdir/$d/R1.count > $outfile_R1");
		}

		if ( -e "$projdir/$d/R2.count" )
		{
			system("cat $projdir/$d/R2.count > $outfile_R2");
		}
		
		if ($booltrackread)
		{
			# get hash reference (args: file, column number)
			my $href = &column_hash("$projdir/$d/R1.id", 1);
			# deref to get hash
			%h1=%$href;
			
			if ( -e "$projdir/$d/R2.id" )
			{
				my $href = &column_hash("$projdir/$d/R2.id", 1);
				%h2=%$href;
			}								
		}		
	}
	elsif ( $d eq "quality_filter" )
	{
		if ( -e "$projdir/$d/R1.count" )
		{
			system("cat $projdir/$d/R1.count >> $outfile_R1");
		}

		if ( -e "$projdir/$d/R2.count" )
		{
			system("cat $projdir/$d/R2.count >> $outfile_R2");
		}
		
		if ($booltrackread)
		{
			# get hash reference (args: file, column number)
			my $href = &column_hash("$projdir/$d/R1.discard", 1);
			# deref to get hash
			%h3=%$href;
			
			if ( -e "$projdir/$d/R2.discard" )
			{
				my $href = &column_hash("$projdir/$d/R2.discard", 1);
				%h4=%$href;
			}								
		}		
	}		
	elsif ( $d eq "host_map" || $d =~ m/^host_map_(\d){1}$/ )
	{
		if ( -e "$projdir/$d/R1.count" )
		{
			system("cat $projdir/$d/R1.count | grep -v input >> $outfile_R1");
		}

		if ( -e "$projdir/$d/R2.count" )
		{
			system("cat $projdir/$d/R2.count | grep -v input >> $outfile_R2");
		}
		
		if ($booltrackread)
		{
			# get hash reference (args: file, column number)
			my $href = &column_hash("$projdir/$d/R1.discard", 1);
			# deref to get hash
			%h5=%$href;
			
			if ( -e "$projdir/$d/R2.discard" )
			{
				my $href = &column_hash("$projdir/$d/R2.discard", 1);
				%h6=%$href;
			}								
		}		
	}
	elsif ( $d eq "nohost_blast" )
	{
		;
	}
	elsif ( $d eq "ray2_assembly_1")		
	{
		if ( -e "$projdir/$d/R1.count" )
		{
			system("cat $projdir/$d/R1.count | grep -v input >> $outfile_R1");
		}

		if ( -e "$projdir/$d/R2.count" )
		{
			system("cat $projdir/$d/R2.count | grep -v input >> $outfile_R2");
		}

		if ( -e "$projdir/$d/R2.count" )
		{
			system("cat $projdir/$d/assembly.count > $outfile_contig");
		}
		
		if ($booltrackread)
		{
			# get hash reference (args: file, column number)
			my $href = &fastqid_hash("$projdir/$d/1.R1.unmap.fastq", 1);
			# deref to get hash
			%h7=%$href;
			
			my $href = &col1_to_col2_hash("$projdir/$d/bowtie2_mapping/R1.map.id");
			%r2c_R1=%$href;

			my $href = &col1_to_col2_hash("$projdir/$d/bowtie2_mapping/singleton.map.id");
			%r2c_single=%$href;

			if ( -e "$projdir/$d/1.R2.unmap.fastq" )
			{
				my $href = &fastqid_hash("$projdir/$d/1.R2.unmap.fastq", 1);
				%h8=%$href;

				my $href = &col1_to_col2_hash("$projdir/$d/bowtie2_mapping/R2.map.id");
				%r2c_R2=%$href;
			}								
		}		

	}
	elsif ( $d eq "iterative_blast_phylo_1")
	{
		# ASSUME iterative_blast_phylo_1 IS A BLAST OF THE CONTIGS WHILE iterative_blast_phylo_2 IS A BLAST OF THE UNASSEM READS
		
		if ( -e "$projdir/$d/R1.count" )
		{
			system("cat $projdir/$d/R1.count | grep -v input >> $outfile_contig");

			# get to contig to taxid mapping
			system("cat $projdir/$d/reports/R1.*.smallreport.txt | cut -f2,16 | sed '1d' > tmp_contig2taxid.txt");
			my $href = &col1_to_col2_hash("tmp_contig2taxid.txt");
			%c2t=%$href;
		}
	}
	elsif ( $d eq "iterative_blast_phylo_2")
	{
		# ASSUME iterative_blast_phylo_1 IS A BLAST OF THE CONTIGS WHILE iterative_blast_phylo_2 IS A BLAST OF THE UNASSEM READS
		
		if ( -e "$projdir/$d/R1.count" )
		{
			system("cat $projdir/$d/R1.count | grep -v input >> $outfile_R1");

			# get to contig to taxid mapping
			system("cat $projdir/$d/reports/R1.*.smallreport.txt | cut -f2,16 | sed '1d' > tmp_R1_2_taxid.txt");
			my $href = &col1_to_col2_hash("tmp_R1_2_taxid.txt");
			%r2t_R1=%$href;			
					
		}

		if ( -e "$projdir/$d/R2.count" )
		{
			system("cat $projdir/$d/R2.count | grep -v input >> $outfile_R2");
			
			# get to contig to taxid mapping
			system("cat $projdir/$d/reports/R2.*.smallreport.txt | cut -f2,16 | sed '1d' > tmp_R2_2_taxid.txt");
			my $href = &col1_to_col2_hash("tmp_R2_2_taxid.txt");
			%r2t_R2=%$href;			
		}
	}	
}

system("rm -f tmp_R1_2_taxid.txt tmp_R2_2_taxid.txt tmp_contig2taxid.txt");

open(my $fh1, ">", "$outputdir/reads.txt");

$header="id"."\t"."R1";
$header="id"."\t"."R1"."\t"."R2" if ( -e "$projdir/step1/R2.id" );

print $fh1 $header."\n";

foreach my $key ( keys %h1 )
{
	my $outline = $key;
		
	if ($h3{$key})
	{
		$outline=$outline."\t"."qual";
	}
	elsif ($h5{$key})
	{
		$outline=$outline."\t"."host";
	}
	elsif ($r2c_R1{$key})
	{
		$outline=$outline."\t"."assem:".$r2c_R1{$key};

		# if contig blasts, print taxid, else print 0
		if ($c2t{$r2c_R1{$key}})
		{
			$outline=$outline.",blast:".$c2t{$r2c_R1{$key}};
		}
		else
		{
			$outline=$outline.",blast:0";
		}
	}
	elsif ($r2c_single{$key})
	{
		$outline=$outline."\t"."assem(single):".$r2c_single{$key};

		# if contig blasts, print taxid, else print 0
		if ($c2t{$r2c_single{$key}})
		{
			$outline=$outline.",blast:".$c2t{$r2c_single{$key}};
		}
		else
		{
			$outline=$outline.",blast:0";
		}
	}	
	elsif ($h7{$key})
	{
		$outline=$outline."\t"."unassem";

		# if read blasts, print taxid, else print 0
		if ($r2t_R1{$key})
		{
			$outline=$outline.",blast:".$r2t_R1{$key};
		}
		else
		{
			$outline=$outline.",blast:0";
		}		
	}
	else
	{
		$outline=$outline."\t"."0";
	}
	
	# if R2
	if ($h2{$key})
	{
		if ($h4{$key})
		{
			$outline=$outline."\t"."qual";
		}
		elsif ($h6{$key})
		{
			$outline=$outline."\t"."host";
		}
		elsif ($r2c_R2{$key})
		{
			$outline=$outline."\t"."assem:".$r2c_R2{$key};
			
			# if contig blasts, print taxid, else print 0
			if ($c2t{$r2c_R2{$key}})
			{
				$outline=$outline.",blast:".$c2t{$r2c_R2{$key}};
			}
			else
			{
				$outline=$outline.",blast:0";
			}			
		}
		elsif ($r2c_single{$key})
		{
			$outline=$outline."\t"."assem(single):".$r2c_single{$key};			
			# if contig blasts, print taxid, else print 0
			if ($c2t{$r2c_single{$key}})
			{
				$outline=$outline.",blast:".$c2t{$r2c_single{$key}};
			}
			else
			{
				$outline=$outline.",blast:0";
			}			
		}		
		elsif ($h8{$key})
		{
			$outline=$outline."\t"."unassem";

			# if read blasts, print taxid, else print 0
			if ($r2t_R2{$key})
			{
				$outline=$outline.",blast:".$r2t_R2{$key};
			}
			else
			{
				$outline=$outline.",blast:0";
			}			
		}		
		else
		{
			$outline=$outline."\t"."0";
		}		
	}

	# print
	print $fh1 $outline."\n";
}

close($fh1);
