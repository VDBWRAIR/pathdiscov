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
use Parse_ParameterFile;

# declare vars
my ($command, $j, $i, $k, $output, $pfile, $r, $run_iteration, $path_scripts, $qrls, $boolphylo);



GetOptions ('iteration_j=i' => \$j,					# j
			'iteration_i=i' => \$i,					# i
			'iteration_k=i' => \$k,					# k
            'paramfile=s' => \$pfile,				# parameter
            'output=s' => \$output,					# output            
            'mate=s' => \$r,						# mate ("R1" or "R2")
            'phylo' => \$boolphylo,            		# boolean which controls if you want to do phylo stuff             
            'run_iteration=i' => \$run_iteration,  	# global # of times you run the script                       
            'path_scripts=s' => \$path_scripts,		# path to scripts
            'qrls=s' => \$qrls);					# qrls id              
            
if (!(eval "require Parse_ParameterFile"))
{
	print "\n[error] The required perl modules cannot be found.\n";
	exit;
}

if (!( -s $j.".R1.fasta" )) # if not nonzero
{
	print "\n[echo] file ".$j.".R1.fasta empty.\n";
	exit;	            
}

# get hash ref
my $href = &parse_param($pfile);
# get hash of parameters
my %hoh=%$href;

# set command
if ($run_iteration>1)
{	
	$command="iterative_blast_phylo_".$run_iteration;	
}
else
{
	$command="iterative_blast_phylo";
} 

# get arrays
my @blast_db_list = split(',',$hoh{$command}{"blast_db_list"});
my @blast_task_list = split(',',$hoh{$command}{"blast_task_list"});
my @blast_options_list = split(',',$hoh{$command}{"blast_options_list"});

# print $command,"\n";
# print $hoh{$command}{"taxonomy_nodes"},"\n";

print "[START]\n";
print "[echo] get phylogeny counts\n";
# args: outputdir, input_file(form: query_id gi_number ...), outputfile_annotate, outputfile taxid2queryid, outputfile, nodes.dmp, names.dmp, ntdb
my $cmd = "$path_scripts/phylogeny_wrapper.sh tmp_".$r."_$j $j.$r.blast $j.$r.blast.ann $j.$r.blast.t2q $j.$r.blast.phylo $hoh{$command}{\"taxonomy_nodes\"} $hoh{$command}{\"taxonomy_names\"} $blast_db_list[$i]";
if ($boolphylo)
{
	print "[cmd] ",$cmd,"\n";
	system($cmd);	
}

my %tmph=();
open(my $infile, "<", "$j.$r.blast"); 
open(my $outfile, ">", "$j.$r.top.blast");
while (<$infile>)
{	
	my @ln=split; 
	if (!($tmph{$ln[0]})) {print $outfile $_;} # if no preexisting entry in the hash, print it 
	$tmph{$ln[0]}=1;			
}				 		 
close($infile); 
close($outfile); 

# args: outputdir, input_file(form: query_id gi_number ...), outputfile_annotate, outputfile taxid2queryid, outputfile, nodes.dmp, names.dmp, ntdb
my $cmd = "$path_scripts/phylogeny_wrapper.sh tmp_".$r."_$j $j.$r.top.blast $j.$r.top.blast.ann $j.$r.top.blast.t2q $j.$r.top.blast.phylo $hoh{$command}{\"taxonomy_nodes\"} $hoh{$command}{\"taxonomy_names\"} $blast_db_list[$i]";
if ($boolphylo)
{
	print "[cmd] ",$cmd,"\n";
	system($cmd);	
}
				
# get reads that didnt blast
# args: blast output, fasta input
my $cmd = "$path_scripts/get_unblast_reads.pl $j.$r.blast $j.$r.fasta > $j.$r.noblast.fasta";
print "[cmd] ",$cmd,"\n";
system($cmd);

# args: input file, filtering_program_name, output file, 2->fasta, concat
my $cmd = "linecount $j.$r.noblast.fasta $blast_task_list[$i] $r.count 2 1";
print "[cmd] ",$cmd,"\n";
system($cmd);			

system("ln -sf $j.$r.noblast.fasta $k.$r.fasta");

# get absolute path to output file
my $mytmp=`readlink -m $j.$r.noblast.fasta`;
chomp $mytmp;

my $cmd = "ln -sf $mytmp $output";
print "[cmd] ",$cmd,"\n";
system($cmd);

# get counts per superclass:
# args: outputdir, R1 or R2
my $cmd = "$path_scripts/prepare_plot_phylo_percents.sh . $r";
if ($boolphylo)
{
	print "[cmd] ",$cmd,"\n";
	system($cmd);	
}

# release hold on next job if jid > 0
if ( $qrls ne "0" )
{
	my $qcmd="qsub -V -N qrls.$i -e ./logs -o ./logs -l mem=1G,time=1:: -S /bin/sh -cwd $path_scripts/job_release.sh $qrls";	
	print "[qcmd] ",$qcmd,"\n";
	my $qsub_message=`$qcmd`;	
	print($qsub_message); 
}

print "[END]\n";


