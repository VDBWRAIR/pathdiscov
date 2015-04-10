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
use Parse_ParameterFile;

# example
# sge_iterative_blast_phylo.pl --sample $sample --paramfile $pfile --outputdir $path_output/iterative_blast_phylo --logs $path_output/iterative_blast_phylo/logs --R1 input_R1.fastq --R2 input_R2.fastq --timestamp 23672

# this is the same as iterative blast but it also makes a "phylogeny count." Unlike iterative_blast, this only works with the blast db s nt or nr
# this one is optimized for the sge job scheduling system

# TO DO:
# make $qtime and $qmem params
# make report generation automatic
# make sure boolphylo is working - TEST W NON NT DB

# declare vars				

my ($outputdir);			# output dir
my ($logs);					# logs path
my ($pfile);				# parameter file 
my ($r1);					# sample R1
my ($r2);					# sample R2
# my ($abs_r1);				# abs path to R1
# my ($abs_r2);				# abs path to R2
my ($sample);				# sample name
my ($isfasta);				# is fasta bool
my ($fastafile);			# is fasta bool
my ($run_iteration);		# run iteration
my ($help);					# help bool
my ($man);					# man bool
my ($timestamp);			# time stamp
my ($command);				# sample name 
my ($boolphylo);			# a boolean to determine whether or not to do phylo stuff
my ($output_r1);			# output file R1
my ($output_r2);			# output file R2
my $path_scripts=$RealBin;	# scripts
my @mates=("R1","R2");		# strings for mates
my $jid_list="";			# comma seperated list of job IDs (used to hold on)
my $qtime=4;				# the time for qsub-ing the jobs first and second
my $qmem=4;					# the mem for qsub-ing the jobs first and second

# default
$timestamp="0";
$run_iteration=1;

GetOptions ('outputdir=s' => \$outputdir,			# outputdir
			'logs=s' => \$logs,						# logs
            'paramfile=s' => \$pfile,				# parameter
            'R1=s' => \$r1,							# R1
            'R2=s' => \$r2,							# R2
            'sample=s' => \$sample,            		# sample            
            'timestamp=s' => \$timestamp,      		# time stamp
            'fastafile=s' => \$fastafile,			# is fasta (default assumes fastq)
            'fasta' => \$isfasta,            		# is fasta (default assumes fastq)
            'nophylo' => \$boolphylo,            	# boolean which controls if you want to do phylo stuff            
            'run_iteration=i' => \$run_iteration,  	# global # of times you run the script             
            'help' => \$help, 						# help
            'man' => \$man);						# man
            
pod2usage(1) if ($help);
pod2usage(-verbose => 2) if ($man);            

die "[error] required input parameters not found" if (!( defined($outputdir) && defined($logs) && defined($pfile) ));

$boolphylo=!($boolphylo);	# default - boolean to do phylogeny stuff is true
           
# get hash ref
my $href = &parse_param($pfile);
# get hash of parameters
my %hoh=%$href;

# --------------------------------------
# param 

# allowed parameters in hash after parsing: (modify for your module)
# blast_db_list
# blast_task_list
# blast_options_list
# ninst_list
# taxonomy_nodes
# taxonomy_names
# qtime_list
# qmem_list

# --------------------------------------

# get absolute path to outputdir
$outputdir=abs_path($outputdir);

# get absolute path to param file
$pfile=abs_path($pfile);

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
my @ninst_list = split(',',$hoh{$command}{"ninst_list"});
my @qtime_list = split(',',$hoh{$command}{"qtime_list"});
my @qmem_list = split(',',$hoh{$command}{"qmem_list"});

# die "[error] blast param lists different sizes" if (!( scalar(@blast_db_list)==scalar(@blast_task_list) && scalar(@blast_db_list)==scalar(@blast_options_list) && scalar(@blast_db_list)==scalar(@chunk_list) && scalar(@blast_db_list)==scalar(@ninst_list) ));
die "[error] blast param lists different sizes" if (!( scalar(@blast_db_list)==scalar(@blast_task_list) && scalar(@blast_db_list)==scalar(@ninst_list) ));

# open output and error logs 
open my $olog, '>', $logs."/".$sample.".".$timestamp."-out.o";
open my $elog, '>', $logs."/".$sample.".".$timestamp."-out.e";

# redirect standard output and error into logs
open STDOUT, '>&', $olog;
open STDERR, '>&', $elog;

if ( $r1 ne "none" && defined($r1) )
{
	$hoh{$command}{"R1"}=abs_path($r1);
}

if ( $r2 ne "none" && defined($r2) )
{
	$hoh{$command}{"R2"}=abs_path($r2);
}

die "[error] input file not found" if (!( -e $hoh{$command}{"R1"} ));

print "[START]\n";
print "[hash] ";
print Dumper \ %hoh;

chdir($outputdir) or die "[error] cannot chdir to $outputdir";

# make special logs dir for blast
my $cmd = "mkdir -p logs_blast";
print_system($cmd);

# --------------------------------------
# final output files - These will be the names of links that point to whatever the last file is				

$output_r1 = $outputdir."/iterative_blast_phylo_".$run_iteration.".R1";	# this will point to the last output stage after the specified level of processing
$output_r2 = $outputdir."/iterative_blast_phylo_".$run_iteration.".R2";	# if paired data, this will point to the last output stage after the specified level of processing		

# --------------------------------------
# main

# loop over each mate
foreach my $mate (@mates) 
{
	# check if defined, and non-zero
	if ( defined($hoh{$command}{$mate}) && -s $hoh{$command}{$mate} )
	{

		# if input fasta file, link to it. O.w., convert fastq to fasta
		if ($isfasta || $fastafile eq "yes")
		{
			system("ln -sf $hoh{$command}{$mate} 1.".$mate.".fasta");
		}
		else
		{
			print "[echo] convert input fastq to fasta $mate \n";
			my $cmd = "cat $hoh{$command}{$mate} | $path_scripts/fastq2fasta.awk > 1.".$mate.".fasta";
			print_system($cmd);
		}

		print "[echo] run a series of blasts $mate \n";	
		
		# count lines in file
		# args: input file, filtering_program_name, output file, 2->fasta, concat
		my $cmd = "$path_scripts/linecount.sh 1.".$mate.".fasta input $mate.count 2 0";
		print_system($cmd);
		
		# here's the larger picture: submit order: iteration n ... 3,2,1; submit order: second first; pass jid2 to first; first releases second when finished; when second itself finishes, it releases the next iteration of first
		# part of this could perhaps be implemented more simply with -hold_jid
		
		my $jid_next = 0;
		
		# for (my $i = 0; $i < scalar(@blast_db_list); $i++) 
		for (my $i = scalar(@blast_db_list) - 1; $i>=0; $i--)	# go backwards, so you can grab jid s
		{
			print "[echo] $blast_task_list[$i] \n";	
			
			my $start_time = time();
			
			# use 1-based, not 0-based, counting
			my $j = $i + 1;
			my $k = $i + 2;
	
			# e.g., if 1.R1.fasta nonzero
			# this line wont work here because files dont exist yet so cant check their size!
#			if ( -s $j.".".$mate.".fasta" )
#			{
				# make tmp dirs (e.g., mkdir tmp_R1_1)
				my $cmd = "mkdir -p tmp_".$mate."_$j";
				print_system($cmd);			
					
				# submit second step first, and immediately hold it. second must release next iteration of first
				my $cmd = "$path_scripts/sge_iterative_blast_phylo_second.pl --iteration_j $j --iteration_i $i --iteration_k $k --paramfile $pfile --mate $mate --path_scripts $path_scripts --run_iteration $run_iteration --output $outputdir/iterative_blast_phylo_".$run_iteration.".$mate --qrls $jid_next --phylo $boolphylo";
				
				my $qcmd="qsub -V -N R1.secnd$j -e ./logs -o ./logs -l mem=".$qmem."G,time=".$qtime.":: -S /usr/bin/perl -cwd $cmd";
				
				my $jid2=qsub_system($qcmd);
				print "second jid: ",$jid2,"\n";
				`qhold $jid2`;			
			
				# submit first step, and release second step when done
			
				# blast in chunks	
				my $cmd = "$path_scripts/sge_block_blast.pl --inputfasta $outputdir/$j.$mate.fasta --outputdir tmp_".$mate."_$j --outfile $outputdir/$j.$mate.blast --outheader $outputdir/blast.header --scripts $path_scripts --ninst $ninst_list[$i] --db $blast_db_list[$i] --blast_type $blast_task_list[$i] --task $blast_task_list[$i] --blast_options \'".$blast_options_list[$i]."\' --qtime $qtime_list[$i] --qmem $qmem_list[$i] --qrls $jid2 --iteration_j $j --mate $mate";		
				
				my $qcmd="qsub -V -N R1.first$j -e ./logs -o ./logs -l mem=".$qmem."G,time=".$qtime.":: -S /usr/bin/perl -cwd $cmd";	
				
				my $jid1=qsub_system($qcmd);		 
				print "first jid: ",$jid1,"\n";
				# hold the job unless it s the 1st iteration 
				`qhold $jid1` if ($i>0);
			
				$jid_next = $jid1;
				
				# update jid_list to include jobs submitted
				$jid_list=$jid_list.",".$jid1.",".$jid2;
				
#			} # if fasta nonzero
		} # iterative blast			
	} # defined, and non-zero
} # loop over each mate

# make reports
my $cmd = "mkdir -p reports";
print_system($cmd);

# format_iterative_blast_phylo.pl --outputdir tmp --prefix sample --blastdir iterative_blast_phylo_1 --blast_list megablast,dc-megablast,blastx
my $cmd = "$path_scripts/format_iterative_blast_phylo.pl --outputdir reports --prefix $sample --blastdir . --blast_list $hoh{$command}{\"blast_task_list\"}";
my $qcmd="qsub -hold_jid $jid_list -V -N reports -e ./logs -o ./logs -l mem=".$qmem."G,time=".$qtime.":: -S /usr/bin/perl -cwd $cmd";
print_system($qcmd);

# --------------------------------------

print "[END]\n";
# print "[hash] ";
# print Dumper \ %hoh;

close(STDOUT);
close(STDERR);  
close($olog);
close($elog);

# ------------------------ Perl pod -----------------------

# this produces: the help and man page

__END__

=head1 NAME

iterative_blast_phylo.pl - run a number, specified by the user, of chained blasts to nt or nr. At each stage, only the sequences that did not blast pass to the next stage.

=head1 SYNOPSIS

iterative_blast_phylo.pl [B<--sample> sample] [B<--outputdir> output_directory] [B<--logs> log_directory] [B<--paramfile> parameter_file(s)]  [B<--R1> R1_input] [B<--timestamp> a time stamp] [options] 

=head1 ARGUMENTS

=over 8

=item B<--sample>

Sample name.

=item B<--outputdir>

The output directory. All of the output will be written in this directory.

=item B<--logs>

The logs directory. stdout and sterr will be written in this folder.

=item B<--paramfile>

The parameter file for the command.

=item B<--R1>

.fastq or .sff input file. If paired reads, this file represents the R1 reads. Can be zipped or not.

=item B<--timestamp>

A time stamp. This ensures your log files will be unique.

=back

=head1 OPTIONS

=over 8

=item B<--R2>

file of the R2 sequences. Required only for paired reads.

=item B<--fasta>

Input is fasta file (default assumes fastq).

=item B<--help>

Prints the help message and exits.

=item B<--man>

Prints the full manual page (including examples) and exits.

=back

=head1 DESCRIPTION

iterative_blast_phylo.pl takes a series of blast db's and runs a series of blasts on them. For example, you might run megablast, followed by dc-megablast, followed by blastx. At each stage, the input for the blast is what didn't blast in the previous stage. Because this module assumes the blast db's are either NCBI nt or NCBI nr, it is also able to annotate the output with taxid etc.

=head1 EXAMPLE

run "iterative_blast_phylo" with paired end .fastq files:

B<iterative_blast_phylo.pl --sample sample --paramfile iterative_blast_phylo.param --outputdir /my/output/directory --logs /my/logs/directory --R1 input_R1.fastq --R2 input_R2.fastq --timestamp 121207-11.30>

=cut

