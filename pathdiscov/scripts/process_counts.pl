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
use Easy_Hash;

# process the count file further

# example
# process_counts.pl --sample sample --outputdir .

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
my ($sample);				# sample name

GetOptions ('outputdir=s' => \$outputdir,		# output directory
			'sample=s' => \$sample);			# sample
                  
system("mkdir -p $outputdir");

# get abs paths
$outputdir=abs_path($outputdir);

chdir($outputdir) or die "[error] cannot chdir to ".$outputdir;

# --------------------------------------
# main

my (%h1, %h2, %h3, %h4, %h5);	# hashes
			
my $href = col1_to_col2_hash("$outputdir/$sample.R1.count.txt");
%h1=%$href;

my $href = col1_to_col2_hash("$outputdir/$sample.contig.count.txt");
%h3=%$href;

# unblasted contigs
my $href = fastaid_hash("$outputdir/../iterative_blast_phylo_1/iterative_blast_phylo_1.R1");
%h4=%$href;

my $href = col1_to_col2_hash("$outputdir/../ray2_assembly_1/contig_numreads.txt");
%h5=%$href;

my $unblasted_read_sum = 0;		# sum of reads in contigs that dont blast
my $blasted_read_sum = 0;		# sum of reads in contigs that blast	

# get number of reads in blasted contigs and unblasted contigs
foreach my $key ( keys %h5 )
{		
	if ($h4{$key})
	{
		$unblasted_read_sum = $unblasted_read_sum + $h5{$key};
	}
	else
	{
		$blasted_read_sum = $blasted_read_sum + $h5{$key};
	}
}

# sum over second column (args: file, column)
my $tot_assemb = column_sum("$outputdir/../ray2_assembly_1/contig_numreads.txt", 2);

my $num_contigs = $h3{"ray_assembly_1"};

my $num_contigs_unblast = $h3{"blastx"};

my $percent_unblast = int(100*$num_contigs_unblast/$num_contigs);

# note: the number of unassembled reads here may not nec be the total (it will only be the total if the user chose to blast ALL unassembled reads)
my $r1_unassemb_reads=fastaid_count("$outputdir/../iterative_blast_phylo_2/1.contig.fasta"); #Dereje: changed R1 to contig
my $r1_unassemb_unblast_reads=fastaid_count("$outputdir/../iterative_blast_phylo_2/iterative_blast_phylo_2.contig");#Dereje: changed R1 to contig
my $percent_r1_unassemb_unblast_reads = int(100*$r1_unassemb_unblast_reads/$r1_unassemb_reads);

# if mate exists
if ( -e "$outputdir/../step1/R2.id" )
{
	my $href = col1_to_col2_hash("$outputdir/$sample.R2.count.txt");
	%h2=%$href;

	my $r2_unassemb_reads=fastaid_count("$outputdir/../iterative_blast_phylo_2/2.contig.fasta");# changed 1.R2.fasta to 2.contig.fasta
	my $r2_unassemb_unblast_reads=fastaid_count("$outputdir/../iterative_blast_phylo_2/iterative_blast_phylo_2.contig");
	my $percent_r2_unassemb_unblast_reads = int(100*$r2_unassemb_unblast_reads/$r2_unassemb_reads);

	my $percent_tot_unassemb_unblast_reads = int(100*($r1_unassemb_unblast_reads + $r2_unassemb_unblast_reads)/($r1_unassemb_reads + $r2_unassemb_reads));

	my $tot= $h1{"rawfile"} + $h2{"rawfile"};
	my $tot_unassemb= $h1{"ray_1"} + $h2{"ray_1"};
	
	print("R1: $h1{\"rawfile\"} reads, R2: $h2{\"rawfile\"} reads, R1+R2: $tot reads for raw input\n");
	print("input to assembly: ",$tot_assemb + $tot_unassemb,"\n");
	print("R1: $h1{\"ray_1\"} unassembled reads, R2: $h2{\"ray_1\"} unassembled reads, R1+R2: $tot_unassemb unassembled reads\n");
	print("R1+R2: $tot_assemb assembled reads\n");
	my $percent_assemb = int(100*$tot_assemb/($tot_assemb + $tot_unassemb));
	
	print("percent of reads assembled: $percent_assemb % \n");
	print("number of contigs: $num_contigs \n");
	print("number of contigs: unblasted $num_contigs_unblast \n");
	print("percent of contigs: unblasted $percent_unblast %\n");

	print("R1+R2: reads comprising unblasted contigs: $unblasted_read_sum \n");
	print("R1+R2: reads comprising blasted contigs: $blasted_read_sum \n");

	print("R1: number of unassembled reads that were blasted: ",$r1_unassemb_reads,"\n");
	print("R2: number of unassembled reads that were blasted: ",$r2_unassemb_reads,"\n");
	print("tot: number of unassembled reads that were blasted: ",$r1_unassemb_reads + $r2_unassemb_reads,"\n");
		
	print("R1: number and percent of unassembled reads that do not blast: ",$r1_unassemb_unblast_reads,", $percent_r1_unassemb_unblast_reads % \n");
	print("R2: number and percent of unassembled reads that do not blast: ",$r2_unassemb_unblast_reads,", $percent_r2_unassemb_unblast_reads % \n");
	print("tot: number and percent of unassembled reads that do not blast: ",$r1_unassemb_unblast_reads + $r2_unassemb_unblast_reads,", $percent_tot_unassemb_unblast_reads % \n");

#	my $percent_x = int(100*($unblasted_read_sum + $tot_unassemb)/($tot));
#	print("reads comprising unblasted contigs plus total unassembled reads divided by tot input reads $percent_x %\n");

	my $percent_y = int(100*($unblasted_read_sum + $r1_unassemb_unblast_reads + $r2_unassemb_unblast_reads)/($tot));
	print("number and percent of reads comprising unblasted contigs plus total unblasted unassembled reads divided by tot input reads: ",$unblasted_read_sum + $r1_unassemb_unblast_reads + $r2_unassemb_unblast_reads,", $percent_y %\n");
}
else
{
	my $percent_tot_unassemb_unblast_reads = int(100*($r1_unassemb_unblast_reads)/($r1_unassemb_reads));

	my $tot= $h1{"rawfile"};
	my $tot_unassemb= $h1{"ray_1"};
	
	print("tot: $tot reads for raw input\n");
	print("input to assembly: ",$tot_assemb + $tot_unassemb,"\n");
	print("tot: $tot_unassemb unassembled reads\n");
	print("tot: $tot_assemb assembled reads\n");
	my $percent_assemb = int(100*$tot_assemb/($tot_assemb + $tot_unassemb));
	
	print("percent of reads assembled: $percent_assemb % \n");
	print("number of contigs: $num_contigs \n");
	print("number of contigs: unblasted $num_contigs_unblast \n");
	print("percent of contigs: unblasted $percent_unblast %\n");

	print("tot: reads comprising unblasted contigs: $unblasted_read_sum \n");
	print("tot: reads comprising blasted contigs: $blasted_read_sum \n");

	print("tot: number of unassembled reads that were blasted: ",$r1_unassemb_reads,"\n");
		
	print("tot: number and percent of unassembled reads that do not blast: ",$r1_unassemb_unblast_reads,", $percent_r1_unassemb_unblast_reads % \n");

	my $percent_y = int(100*($unblasted_read_sum + $r1_unassemb_unblast_reads)/($tot));
	print("number and percent of reads comprising unblasted contigs plus total unblasted unassembled reads divided by tot input reads: ",$unblasted_read_sum + $r1_unassemb_unblast_reads,", $percent_y %\n");
}

# print Dumper \ %h4;

# open(my $fh1, ">", "$outputdir/numbers.txt");
# close($fh1);
