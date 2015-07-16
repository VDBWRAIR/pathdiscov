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
# quality_filter.pl --sample $sample --paramfile $pfile --outputdir $path_output/quality_filter --logs $path_output/quality_filter/logs --R1 input_R1.fastq --R2 input_R2.fastq --timestamp 23672

# TO DO:
# put R1 R2 in loop instead of repeat & * * * make parallel * * * 
# links use rel not abs path
# improve syntax of string commands
# add $mate.discard

# declare vars
my ($outputdir);				# output dir
my ($logs);					# logs path
my ($pfile);					# parameter file 
my ($r1);					# sample R1
my ($r2);					# sample R2
my ($sample);					# sample name
my ($help);					# help bool
my ($man);					# man bool
my $timestamp="0";				# time stamp (default is 0)
my ($output_r1);				# output R1
my ($output_r2);				# output R2
my $path_scripts=$RealBin;			# get abs path to directory in which script resides (where to look for sister scripts)
my @mates=("R1","R2");				# strings for mates
my $command="quality_filter";			# set command


GetOptions (	'outputdir=s' => \$outputdir,	# outputdir
		'logs=s' => \$logs,		# logs
		'paramfile=s' => \$pfile,	# parameter
		'R1=s' => \$r1,			# R1
		'R2=s' => \$r2,			# R2
		'sample=s' => \$sample,		# sample
		'timestamp=s' => \$timestamp,	# time stamp
		'help' => \$help,		# help
		'man' => \$man);		# man
            
pod2usage(1) if ($help);
pod2usage(-verbose => 2) if ($man);            

die "[error] required input parameters not found" if (!( defined($outputdir) && defined($logs) && defined($pfile) && defined($r1) ));
            
# get hash ref
my $href = &parse_param($pfile);
# get hash of parameters
my %hoh=%$href;

# --------------------------------------
# param 

# allowed parameters in hash after parsing: (modify for your module)
# command
# cutadapt_options
# cutadapt_options2
# prinseq_options

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
    my $cmd = "linecount $hoh{$command}{$mate} input $mate.count 1 1";
    print "[cmd] ",$cmd,"\n";
    system($cmd);
}

if ( $r2 ne "none" && defined($r2) )
{
	$hoh{$command}{"R2"}=abs_path($r2);
    my $cmd = "linecount $hoh{$command}{$mate} input $mate.count 1 1";
    print "[cmd] ",$cmd,"\n";
    system($cmd);
} 	

print "[START]\n";
print "[hash] ";
print Dumper \ %hoh;

chdir($outputdir) or die "[error] cannot chdir to $outputdir";

# --------------------------------------
# files				

# my $file1out = $outputdir."/R1.cut.fastq";		# after cutting adapters (if requested)
# my $file1out1 = $outputdir."/R1.cut2.fastq";		# after cutting adapters second time (if requested)
# my $file1out2 = $outputdir."/R1.prinseq.fastq";	# after prinseq (if requested)
			
# my $file2out = $outputdir."/R2.cut.fastq";
# my $file2out1 = $outputdir."/R2.cut2.fastq";		
# my $file2out2 = $outputdir."/R2.prinseq.fastq";

# my $file1count = $outputdir."/R1.count";		# a file to track read counts
# my $file2count = $outputdir."/R2.count";

# --------------------------------------
# main

foreach my $mate (@mates) 
{
	# check if defined, and non-zero
	if ( defined($hoh{$command}{$mate}) && -s $hoh{$command}{$mate} )
	{
		# not the best style, but copy the input entry in the hash for use at the end
		$hoh{$command}{$mate."input"} = $hoh{$command}{$mate};

		system("ln -sf $hoh{$command}{$mate} $command.$mate");		

		# ---------------------------------------------------------
		# if cut adapter, do that now
		if ( -e "$mate.cut.fastq" )
		{	
			# if file already exists, put it in the hash
			print "[file exists] $mate.cut.fastq \n";			
			$hoh{$command}{$mate} = "$mate.cut.fastq";	
		}
		else
		{			
			# if ($hoh{$command}{"cutadapt_options"})
			# {
			# 	# run cut adapt
			# 	print "***cut_adapter***\n";
			# 	print "[echo] filter reads w adapters $mate iteration:1\n";
			# 	# E.g., cutadapt R1.fastq -g GATGCATGCAAGACC -o out2.fastq -m 20 --match-read-wildcards
			# 	my $cmd = "cutadapt $hoh{$command}{$mate} -o $mate.cut.fastq $hoh{$command}{\"cutadapt_options\"}";	
			# 	verbose_system($cmd);
				
			# 	# link current output to file
			# 	$hoh{$command}{$mate} = "$mate.cut.fastq";
			# 	system("ln -sf $mate.cut.fastq $command.$mate");		
		
			# 	# count lines
			# 	my $cmd = "linecount $hoh{$command}{$mate} cut_adapt $mate.count 1 1";
			# 	print "[cmd] ",$cmd,"\n";
			# 	system($cmd);		
			# }
			if ($hoh{$command}{"cutadapt_options_R1"} and $mate eq "R1")
			{
				# run cut adapt
				print "***cut_adapter***\n";
				print "[echo] filter reads w adapters $mate iteration:1\n";
				# E.g., cutadapt R1.fastq -g GATGCATGCAAGACC -o out2.fastq -m 20 --match-read-wildcards
				my $cmd = "cutadapt $hoh{$command}{$mate} -o $mate.cut.fastq $hoh{$command}{\"cutadapt_options_R1\"}";	
				verbose_system($cmd);
				
				# link current output to file
				$hoh{$command}{$mate} = "$mate.cut.fastq";
				system("ln -sf $mate.cut.fastq $command.$mate");		
		
				# count lines
				my $cmd = "linecount $hoh{$command}{$mate} cut_adapt $mate.count 1 1";
				print "[cmd] ",$cmd,"\n";
				system($cmd);		
			}
			if ($hoh{$command}{"cutadapt_options_R2"} and $mate eq "R2")
			{
				# run cut adapt
				print "***cut_adapter***\n";
				print "[echo] filter reads w adapters $mate iteration:1\n";
				# E.g., cutadapt R1.fastq -g GATGCATGCAAGACC -o out2.fastq -m 20 --match-read-wildcards
				my $cmd = "cutadapt $hoh{$command}{$mate} -o $mate.cut.fastq $hoh{$command}{\"cutadapt_options_R2\"}";	
				verbose_system($cmd);
				
				# link current output to file
				$hoh{$command}{$mate} = "$mate.cut.fastq";
				system("ln -sf $mate.cut.fastq $command.$mate");		
		
				# count lines
				my $cmd = "linecount $hoh{$command}{$mate} cut_adapt $mate.count 1 1";
				print "[cmd] ",$cmd,"\n";
				system($cmd);		
			}
		}
		
		# # if cut adapter a second time, do that now
		# if ( -e "$mate.cut2.fastq" )
		# {	
		# 	# if file already exists, put it in the hash
		# 	print "[file exists] $mate.cut2.fastq \n";			
		# 	$hoh{$command}{$mate} = "$mate.cut2.fastq";	
		# }
		# else
		# {			
		# 	if ($hoh{$command}{"cutadapt_options2"})
		# 	{
		# 		# run cut adapt
		# 		print "[echo] filter reads w adapters $mate iteration:2\n";
		# 		my $cmd = "cutadapt $hoh{$command}{$mate} -o $mate.cut2.fastq $hoh{$command}{\"cutadapt_options2\"}";	
		# 		verbose_system($cmd);
				
		# 		# link current output to file
		# 		$hoh{$command}{$mate} = "$mate.cut2.fastq";
		# 		system("ln -sf $mate.cut2.fastq $command.$mate");		
		
		# 		# count lines
		# 		my $cmd = "linecount $hoh{$command}{$mate} cut_adapt2 $mate.count 1 1";
		# 		print "[cmd] ",$cmd,"\n";
		# 		system($cmd);		
		# 	}
		# }		
		
		# ---------------------------------------------------------
		# if prinseq, do that now
		if ( -e "$mate.prinseq.fastq" )
		{	
			# if file already exists, put it in the hash
			print "[file exists] $mate.prinseq.fastq \n";
			$hoh{$command}{$mate} = "$mate.prinseq.fastq";
		}
		else
		{			
			if ($hoh{$command}{"prinseq_options"})
			{
				print "***prinseq***\n";
				print "[echo] quality filter $mate\n";
				# E.g., prinseq-lite.pl -log -verbose -fastq 1.cut.fastq -min_len 50 -ns_max_p 10 -derep 12345 -out_good 1.prinseq -out_bad 1.bad
				my $cmd = "prinseq-lite.pl -fastq $hoh{$command}{$mate} -out_good $mate.prinseq -out_bad $mate.prinseq.bad $hoh{$command}{\"prinseq_options\"}";				
				verbose_system($cmd);
		
				# link current output to file
				$hoh{$command}{$mate} = "$mate.prinseq.fastq";
				system("ln -sf $mate.prinseq.fastq $command.$mate");
				
				my $cmd = "linecount $hoh{$command}{$mate} prinseq $mate.count 1 1";
				print "[cmd] ",$cmd,"\n";
				system($cmd);
			}
		}
		
		# --------------------------------------
		
		# get discard IDs - i.e., the reads filtered in this step
		my $cmd = "fastaq_tools_diff.exe --fastq $hoh{$command}{$mate.\"input\"} --fastq $hoh{$command}{$mate} > $mate.discard";
		verbose_system($cmd);		

	} # defined
	else
	{
		print "input not defined\n";
	}
} # mate

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

quality_filter.pl - does adapter and quality filtering

=head1 SYNOPSIS

quality_filter.pl [B<--sample> sample] [B<--outputdir> output_directory] [B<--logs> log_directory] [B<--paramfile> parameter_file(s)]  [B<--R1> R1_input] [B<--timestamp> a time stamp] [options] 

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

.fastq or .sff input file of the R2 reads. Can be zipped or not. Required only for paired reads.

=item B<--example>

Prints example parameter file.

=item B<--help>

Prints the help message and exits.

=item B<--man>

Prints the full manual page (including examples) and exits.

=back

=head1 DESCRIPTION

quality_filter.pl runs up to 2 iterations of cutadapt, if specified, and prinseq, if specified.

=head1 EXAMPLE

run "quality_filter" with paired end .fastq files:

B<quality_filter.pl --sample sample --paramfile quality_filter.param --outputdir /my/output/directory --logs /my/logs/directory --R1 input_R1.fastq --R2 input_R2.fastq --timestamp 121207-11.30>

=cut
