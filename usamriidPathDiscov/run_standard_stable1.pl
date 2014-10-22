#!/usr/bin/perl

# example:
# run_standard.pl --sample mysample --outputdir ../results --R1 R1.fastq --R2 R2.fastq --blast_unassembled 1000

use strict;
use Getopt::Long;
use Cwd 'abs_path';

use FindBin qw($RealBin);
use lib "$RealBin/Local_Module";

# local modules:
use Verbose_Sys;

# declare global variables:
my ($outputdir);													# output dir
my ($sample);														# sample name
my ($numreads);														# number of unassembled reads to blast
my ($path_scripts);													# path variables (scripts)
my ($r1, $r2);														# paths to r1 and r2, if present

# default:
$r2="none";
$numreads=0;

GetOptions ('outputdir=s' => \$outputdir,
			'blast_unassembled=i' => \$numreads,
            'R1=s' => \$r1,
            'R2=s' => \$r2,
            'sample=s' => \$sample);
            
# -------------------- main --------------------

# scripts paths 
$path_scripts=$RealBin;

print("-|-------pathogen pipeline-------|-\n");
print_system("$path_scripts/pathogen.pl --sample $sample --command step1 quality_filter host_map ray2_assembly iterative_blast_phylo --paramfile param.txt --outputdir $outputdir --R1 $r1 --R2 $r2");

# fastq files come in 4-line chunks
$numreads=$numreads*4;

# default is use all unassembled reads
if ($numreads == 0)
{
	$r1="$outputdir/ray2_assembly_1/1.R1.unmap.fastq";
	$r2="$outputdir/ray2_assembly_1/1.R2.unmap.fastq" if ($r2 ne "none");
}
# else if num reads specified
else
{
	# R1
	system("cat $outputdir/ray2_assembly_1/1.R1.unmap.fastq | head -$numreads > $outputdir/ray2_assembly_1/head.1.R1.unmap.fastq");
	$r1="$outputdir/ray2_assembly_1/head.1.R1.unmap.fastq";

	# R2
	system("cat $outputdir/ray2_assembly_1/1.R2.unmap.fastq | head -$numreads > $outputdir/ray2_assembly_1/head.1.R2.unmap.fastq") if ($r2 ne "none");
	$r2="$outputdir/ray2_assembly_1/head.1.R2.unmap.fastq" if ($r2 ne "none");
}

print("\n-|-------unassembled reads-------|-\n");
print_system("$path_scripts/pathogen.pl --sample $sample --command iterative_blast_phylo_2 --paramfile param.txt --outputdir $outputdir --R1 $r1 --R2 $r2");

print("\n-|-------read counts-------|-\n");
system("mkdir -p $outputdir/output");
print_system("$path_scripts/scripts/readcount.pl --sample $sample --outputdir $outputdir/output --projdir $outputdir --dirlist \"step1,quality_filter,host_map_1,ray2_assembly_1,iterative_blast_phylo_1,iterative_blast_phylo_2\" --trackread");

print_system("$path_scripts/scripts/process_counts.pl --sample $sample --outputdir $outputdir/output > $outputdir/output/stats.txt");

if ($r2 ne "none")
{
	print_system("$path_scripts/scripts/join_smallreport.pl --outputdir $outputdir/iterative_blast_phylo_2/reports --prefix $sample --R1report $outputdir/iterative_blast_phylo_2/reports/R1.$sample.top.smallreport.txt --R2report $outputdir/iterative_blast_phylo_2/reports/R2.$sample.top.smallreport.txt --R1qualdiscard $outputdir/quality_filter/R1.discard --R1hostdiscard $outputdir/host_map_1/R1.discard --R2qualdiscard $outputdir/quality_filter/R2.discard --R2hostdiscard $outputdir/host_map_1/R2.discard");

}

print("\n-|-------checkerror-------|-\n");
print_system("$path_scripts/pathogen.pl --checkerror --outputdir $outputdir");

print("\n-|-------cleanup-------|-\n");
print_system("$path_scripts/pathogen.pl --cleanup --outputdir $outputdir");
