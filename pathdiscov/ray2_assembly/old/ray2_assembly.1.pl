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
# ray2_assembly.pl --sample $sample --paramfile $pfile --outputdir $path_output/ray2_assembly --logs $path_output/ray2_assembly/logs --R1 input_R1.fastq --R2 input_R2.fastq --timestamp 23672

# to do
# input_R1 -> R1
# count reads mapping to contigs
# √ put Ray logs in seperate place
# add discards
# count case where contig file empty

my ($outputdir);			# output dir
my ($logs);					# logs path
my ($pfile);				# parameter file 
my ($r1);					# sample R1
my ($r2);					# sample R2
my ($abs_r1);				# abs path to R1
my ($abs_r2);				# abs path to R2
my ($sample);				# sample name
my ($isfasta);				# is fasta bool
my ($fastafile);			# is fasta bool
my ($run_iteration);		# run iteration
my ($help);					# help bool
my ($man);					# man bool
my ($timestamp);			# time stamp
my ($command);				# sample name 
my ($output);				# output of assembly
my ($output_R1_unmap);		# reads which dont map to assembly (unassembled), R1
my ($output_R2_unmap);		# reads which dont map to assembly (unassembled), R2
my $path_scripts=$RealBin;	# path to scripts

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
            'run_iteration=i' => \$run_iteration,  	# global # of times you run the script           
            'help' => \$help, 						# help
            'man' => \$man);						# man
            
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
# kmer
# map2contigs
# bowtie2_options

# --------------------------------------

# get absolute path to outputdir
$outputdir=abs_path($outputdir);

# set command
if ($run_iteration>1)
{	
	$command="ray2_assembly_".$run_iteration;	
}
else
{
	$command="ray2_assembly";
}

# open output and error logs 
open my $olog, '>', $logs."/".$sample.".".$timestamp."-out.o";
open my $elog, '>', $logs."/".$sample.".".$timestamp."-out.e";

# redirect standard output and error into logs
open STDOUT, '>&', $olog;
open STDERR, '>&', $elog;

if ( $r1 ne "none" && defined($r1) )
{
	$abs_r1 = abs_path($r1);
	$hoh{$command}{"R1"}=$abs_r1;
}

if ( $r2 ne "none" && defined($r2) )
{
	$abs_r2 = abs_path($r2);
	$hoh{$command}{"R2"}=$abs_r2;		
}

die "[error] input file not found" if (!( -e $hoh{$command}{"R1"} ));

print "[START]\n";
print "[hash] ";
print Dumper \ %hoh;

chdir($outputdir) or die "[error] cannot chdir to $outputdir";

# make special logs dir for blast
my $cmd = "mkdir -p logs_assembly";
print "[cmd] ",$cmd,"\n";		
system($cmd);

# --------------------------------------
# output files				

$output = $outputdir."/ray2_assembly_".$run_iteration.".fasta";			
$output_R1_unmap = $outputdir."/".$run_iteration.".R1.unmap.fastq";	
$output_R2_unmap = $outputdir."/".$run_iteration.".R2.unmap.fastq";			

my $filecount = $outputdir."/assembly.count";		# a file to track read counts

# add counts ! 

# --------------------------------------
# main

# do assembly
# if fasta
if ($isfasta || $fastafile eq "yes")
{
	# if R2 exists
	if ($hoh{$command}{"R2"})
	{
		print "[echo] find mate pairs and singletons\n";
		my $cmd = "$path_scripts/perlscripts_wrapper.pl get_common_uneven_fastas $hoh{$command}{\"R1\"} $hoh{$command}{\"R2\"} R1.single.fasta R1.paired.fasta R2.single.fasta R2.paired.fasta";
		verbose_system($cmd);		
		
		print "[echo] ray2 assembly - paired end fasta\n";
		my $cmd = "Ray2 -o results -k $hoh{$command}{\"kmer\"} -s R1.single.fasta -s R2.single.fasta -p R1.paired.fasta R2.paired.fasta > logs_assembly/assembly.o 2> logs_assembly/assembly.e";
		# And you can have as many -s and -p flags as you want, so you can have multiple paired end datasets; plus, Ray2 uses suffix to determine fastq or fasta
		verbose_system($cmd);		
	}
	else
	{
		print "[echo] ray2 assembly - single end fasta\n";

		my $cmd = "ln -sf $hoh{$command}{\"R1\"} R1.single.fasta";
		system($cmd);		
		
		my $cmd = "Ray2 -o results -k $hoh{$command}{\"kmer\"} -s R1.single.fasta > logs_assembly/assembly.o 2> logs_assembly/assembly.e";
		print "[cmd] ",$cmd,"\n";
		verbose_system($cmd);
	}
	
	print "[echo] map2contigs option turned off for fasta formatted files\n" if ( $hoh{$command}{"map2contigs"} eq "yes" || $hoh{$command}{"map2contigs"} eq "1" );	
	
}
# if fastq
else
{
	# if R2 exists
	if ($hoh{$command}{"R2"})
	{
		print "[echo] find mate pairs and singletons\n";
		my $cmd = "$path_scripts/perlscripts_wrapper.pl get_common_uneven_files $hoh{$command}{\"R1\"} $hoh{$command}{\"R2\"} R1.single.fastq R1.paired.fastq R2.single.fastq R2.paired.fastq";
		verbose_system($cmd);	

		print "[echo] ray2 assembly - paired end fastq\n";
		
		# do assembly according to which files exist 
		if ( -s "R1.single.fastq" && -s "R2.single.fastq" )
		{
			my $cmd = "Ray2 -o results -k $hoh{$command}{\"kmer\"} -s R1.single.fastq -s R2.single.fastq -p R1.paired.fastq R2.paired.fastq > logs_assembly/assembly.o 2> logs_assembly/assembly.e";
			verbose_system($cmd);	
		}
		elsif ( -s "R1.single.fastq" )
		{
			my $cmd = "Ray2 -o results -k $hoh{$command}{\"kmer\"} -s R1.single.fastq -p R1.paired.fastq R2.paired.fastq > logs_assembly/assembly.o 2> logs_assembly/assembly.e";
			verbose_system($cmd);	
		}
		elsif ( -s "R2.single.fastq" )
		{
			my $cmd = "Ray2 -o results -k $hoh{$command}{\"kmer\"} -s R2.single.fastq -p R1.paired.fastq R2.paired.fastq > logs_assembly/assembly.o 2> logs_assembly/assembly.e";
			verbose_system($cmd);			
		}
		else
		{
			my $cmd = "Ray2 -o results -k $hoh{$command}{\"kmer\"} -p R1.paired.fastq R2.paired.fastq > logs_assembly/assembly.o 2> logs_assembly/assembly.e";
			verbose_system($cmd);		
		}
	}
	else
	{
		print "[echo] ray2 assembly - single end fastq\n";
		
		my $cmd = "ln -sf $hoh{$command}{\"R1\"} R1.single.fastq";
		system($cmd);		
		
		my $cmd = "Ray2 -o results -k $hoh{$command}{\"kmer\"} -s R1.single.fastq > logs_assembly/assembly.o 2> logs_assembly/assembly.e";
		verbose_system($cmd);		
	}
}

# put contigs seqs on a single line and replace ID strings w numbers
# my $cmd = "cat $outputdir/results/Contigs.fasta | $path_scripts/../scripts/fastajoinlines > $output";
# args: input output
my $cmd = "$path_scripts/joinlines.sh $outputdir/results/Contigs.fasta $output";
print "[cmd] ",$cmd,"\n";
system($cmd);

# args: input file, filtering_program_name, output file, 2->fasta, concat
my $cmd = "$path_scripts/linecount.sh $output ray_assembly_".$run_iteration." $filecount 2 0";
print "[cmd] ",$cmd,"\n";
system($cmd);

# if map to contigs option:
if ( $hoh{$command}{"map2contigs"} eq "yes" || $hoh{$command}{"map2contigs"} eq "1" )
{
	# if fasta, exit
	if ($isfasta || $fastafile eq "yes")
	{
		exit;		
	}
	
	my $cmd = "mkdir -p bowtie2_index";
	print "[cmd] ",$cmd,"\n";
	system($cmd);	

	# index contigs
	$cmd = "bowtie2-build $output bowtie2_index/contigs";
	verbose_system($cmd);

	my $cmd = "mkdir -p bowtie2_mapping";
	print "[cmd] ",$cmd,"\n";
	system($cmd);
	
	# if R2 exists
	if ($hoh{$command}{"R2"})
	{	
		# map according to which mates and singletons are non-zero	
		if ( -s "R1.single.fastq" && -s "R2.single.fastq" )
		{
			my $cmd = "bowtie2 -q -x bowtie2_index/contigs -1 R1.paired.fastq -2 R2.paired.fastq -U R1.single.fastq,R2.single.fastq -S bowtie2_mapping/out.sam $hoh{$command}{\"bowtie2_options\"}";
			verbose_system($cmd);		
		}
		elsif ( -s "R1.single.fastq" )
		{
			my $cmd = "bowtie2 -q -x bowtie2_index/contigs -1 R1.paired.fastq -2 R2.paired.fastq -U R1.single.fastq -S bowtie2_mapping/out.sam $hoh{$command}{\"bowtie2_options\"}";
			verbose_system($cmd);		
		}
		elsif ( -s "R2.single.fastq" )
		{
			my $cmd = "bowtie2 -q -x bowtie2_index/contigs -1 R1.paired.fastq -2 R2.paired.fastq -U R2.single.fastq -S bowtie2_mapping/out.sam $hoh{$command}{\"bowtie2_options\"}";
			verbose_system($cmd);		
		}
		else
		{
			my $cmd ="bowtie2 -q -x bowtie2_index/contigs -1 R1.paired.fastq -2 R2.paired.fastq -S bowtie2_mapping/out.sam $hoh{$command}{\"bowtie2_options\"}";
			verbose_system($cmd);		
		}
	}
	else
	{
		my $cmd = "bowtie2 -q -x bowtie2_index/contigs -U R1.single.fastq -S bowtie2_mapping/out.sam $hoh{$command}{\"bowtie2_options\"}";
		verbose_system($cmd);		
	}

	# get number of reads mapping to each contig, as well as number of split reads
	$cmd = "cat bowtie2_mapping/out.sam | $path_scripts/count_reads2contigs.awk | sed \'s|\*|unmapped|\' | sort -k2,2nr > contig_numreads.txt";
	verbose_system($cmd);	
	
	# get unmapped reads
	# args: outputdir scripts r1 r2
	$cmd = "$path_scripts/get_unmap.sh bowtie2_mapping $path_scripts $hoh{$command}{\"R1\"} $hoh{$command}{\"R2\"}";
	verbose_system($cmd);	
	
	# create link
	my $cmd = "ln -sf bowtie2_mapping/R1.unmap.fastq $output_R1_unmap";
	system($cmd);

	# count lines
	# args: input file, filtering_program_name, output file, 2->fasta, concat
	my $cmd = "$path_scripts/linecount.sh $output_R1_unmap ray_".$run_iteration.".R1.unassembled $filecount 1 1";
	print "[cmd] ",$cmd,"\n";
	system($cmd);	

	# if R2 exists
	if ( -s "$outputdir/bowtie2_mapping/R2.unmap.fastq" )
	{
		# create link
		my $cmd = "ln -sf bowtie2_mapping/R2.unmap.fastq $output_R2_unmap";
		system($cmd);
	
		# count lines
		# args: input file, filtering_program_name, output file, 2->fasta, concat
		my $cmd = "$path_scripts/linecount.sh $output_R2_unmap ray_".$run_iteration.".R2.unassembled $filecount 1 1";
		print "[cmd] ",$cmd,"\n";
		system($cmd);		
	}

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

ray2_assembly.pl - perform a Ray2 assembly.

=head1 SYNOPSIS

ray2_assembly.pl [B<--sample> sample] [B<--outputdir> output_directory] [B<--logs> log_directory] [B<--paramfile> parameter_file(s)]  [B<--R1> R1_input] [B<--timestamp> a time stamp] [options] 

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

ray2_assembly.pl serves as a wrapper for ray2.

=head1 EXAMPLE

run "ray2_assembly" with paired end .fastq files:

B<ray2_assembly.pl --sample sample --paramfile ray2_assembly.param --outputdir /my/output/directory --logs /my/logs/directory --R1 input_R1.fastq --R2 input_R2.fastq --timestamp 121207-11.30>

=cut
