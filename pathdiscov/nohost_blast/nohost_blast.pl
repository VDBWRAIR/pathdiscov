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
# nohost_blast.pl --sample $sample --paramfile $pfile --outputdir $path_output/nohost_blast --logs $path_output/nohost_blast/logs --R1 input_R1.fastq --R2 input_R2.fastq --timestamp 23672

# TO DO:
# put r1 r2 in a loop and error check ! 

# declare vars
my ($outputdir);			# output dir
my ($logs);					# logs path
my ($pfile);				# parameter file 
my ($r1);					# sample R1
my ($r2);					# sample R2
my ($sample);				# sample name
my ($isfasta);				# is fasta bool
my ($fastafile);			# is fasta bool
my ($help);					# help bool
my ($man);					# man bool
my ($timestamp);			# time stamp
my $path_scripts=$RealBin;	# get abs path to directory in which script resides
my @mates=("R1","R2");		# strings for mates
my $command="nohost_blast";	# set command

# default
$timestamp="0";

GetOptions ('outputdir=s' => \$outputdir,		# outputdir
			'logs=s' => \$logs,					# logs
            'paramfile=s' => \$pfile,			# parameter
            'R1=s' => \$r1,						# R1
            'R2=s' => \$r2,						# R2
            'sample=s' => \$sample,            	# sample            
            'timestamp=s' => \$timestamp,      	# time stamp
            'fastafile=s' => \$fastafile,		# is fasta (default assumes fastq)            
            'fasta' => \$isfasta,            	# is fasta (default assumes fastq)
            'help' => \$help, 					# help
            'man' => \$man);					# man
            
pod2usage(1) if ($help);
pod2usage(-verbose => 2) if ($man);            

die "[error] required input parameters not found" if (!( defined($outputdir) && defined($logs) && defined($pfile) ));
            
# get hash ref
my $href = &parse_param($pfile);
# get hash of parameters
my %hoh=%$href;

# --------------------------------------
# param 

# allowed parameters in hash after parsing: (modify for your module)
# command
# ncbi_nt_db
# gi2taxid
# num_subset_seq
# blast_type
# blast_task
# blast_options
# ninst

# --------------------------------------

# get absolute path to outputdir
$outputdir=abs_path($outputdir);

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

# --------------------------------------
# main

foreach my $mate (@mates) 
{
	# check if defined, and non-zero
	if ( defined($hoh{$command}{$mate}) && -s $hoh{$command}{$mate} )
	{
		# if input fasta file, link to it. O.w., convert fastq to fasta
		if ($isfasta || $fastafile eq "yes")
		{
			system("ln -sf $hoh{$command}{$mate} $mate.fasta");
		}
		else
		{
			print "[echo] convert input fastq to fasta $mate \n";
			my $cmd = "cat $hoh{$command}{$mate} | $path_scripts/fastq2fasta.awk > $mate.fasta";
			print_system($cmd);
		}
	}

	# make the subset blast db (only once - for the first mate) 	
	if ( $mate eq "R1" )
	{
		print "[echo] blast a subset of your input file to the ncbi nt db and use the taxids to make an indexed subset db\n";
		# args: outputdir, input fasta file, output blast_db name, ntdb, gi2taxid
		# output BLAST prefix: outputdir/subset_db
		my $cmd = "$path_scripts/subset_db_blast.sh . $mate.fasta subset_db $hoh{$command}{\"ncbi_nt_db\"} $hoh{$command}{\"gi2taxid\"} $hoh{$command}{\"num_subset_seq\"}";
		verbose_system($cmd);
	}
	
	print "[echo] count\n";
	# args: input file, filtering_program_name, output file, 2->fasta, concat
	my $cmd = "$path_scripts/linecount.sh $mate.fasta input $mate.count 2 0";
	print_system($cmd);
		
	print "[echo] blast input against subset db $mate \n";
	# dont try a regular blast - it will take ages! 
	# my $cmd = $path_scripts."/blast_wrapper.pl --type blastn --query $mate.fasta --db subset_db --task megablast --out 1.blast --options \'-evalue 1e-4 -word_size 28\'";
	my $cmd = "mkdir -p tmp_$mate";
	print_system($cmd);
	# note: $hoh{$command}{blast_options} can be a string with spaces in it, so you have to be careful when passing it as an argument
	my $cmd = "$path_scripts/par_block_blast.pl --outputdir tmp_$mate --inputfasta $outputdir/$mate.fasta --db $outputdir/subset_db --blast_type $hoh{$command}{\"blast_type\"} --task $hoh{$command}{\"blast_task\"} --ninst $hoh{$command}{\"ninst\"} --outfile $outputdir/$mate.blast --outheader $outputdir/blast.header --blast_options \"$hoh{$command}{\"blast_options\"}\"";		
	verbose_system($cmd);
	
	print "[echo] get what didnt blast $mate \n";
	my $cmd = "$path_scripts/get_unblast_reads.pl $mate.blast $mate.fasta > $mate.noblast.fasta";
	print_system($cmd);
		
	system("ln -sf $mate.noblast.fasta nohost_blast.$mate");
	
	print "[echo] count\n";
	# args: input file, filtering_program_name, output file, 2->fasta, concat
	my $cmd = "$path_scripts/linecount.sh $mate.noblast.fasta ".$hoh{$command}{"blast_task"}."_subset $mate.count 2 1";
	print_system($cmd);	
}
			
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

nohost_blast.pl - blast a subset of your input file to the ncbi nt db and use the taxids to make an indexed subset db

=head1 SYNOPSIS

nohost_blast.pl [B<--sample> sample] [B<--outputdir> output_directory] [B<--logs> log_directory] [B<--paramfile> parameter_file(s)]  [B<--R1> R1_input] [B<--timestamp> a time stamp] [options] 

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

nohost_blast.pl first blasts a subset of your input fastq file (the initial 50 sequences) to the ncbi nt db. The sequences that blast return ncbi gi (GenInfo Identifier) numbers.
The script finds the taxids associated with these gi numbers and then, in turn, finds all the gi numbers associated with these taxids (a much larger set of gi numbers than the
original one). Next, it extracts all the sequences from the nt db with these gi numbers and makes an indexed subset db, which is much smaller than nt. Finally, the script blasts 
your full input fastq file against this subset db.

=head1 EXAMPLE

run "nohost_blast" with paired end .fastq files:

B<nohost_blast.pl --sample sample --paramfile nohost_blast.param --outputdir /my/output/directory --logs /my/logs/directory --R1 input_R1.fastq --R2 input_R2.fastq --timestamp 121207-11.30>

=cut
