#!/usr/bin/env perl 
use strict;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use Cwd 'abs_path';
use File::Spec;

use FindBin qw($RealBin);
use lib "$RealBin/../lib/Local_Module";
# local modules:
use Verbose_Sys;
use Parse_ParameterFile;

# example
# orf_filter.pl --sample $sample --paramfile $pfile --outputdir $path_output/orf_filter --logs $path_output/orf_filter/logs --R1 input_R1.fasta --R2 input_R2.fasta --timestamp 23672

# declare vars
my ($outputdir);				# output dir
my ($logs);					# logs path
my ($pfile);					# parameter file 
my ($r1);					# sample R1
my ($r2);					# sample R2
my ($sample);					# sample name
my $timestamp="0";				# time stamp (default is 0)
my ($output_r1);				# output R1
my ($output_r2);				# output R2
my $contig = 0;			# contig bool
my $path_scripts=$RealBin;			# get abs path to directory in which script resides (where to look for sister scripts)
my @mates=("R1","R2");				# strings for mates
my $command="orf_filter";			# set command


GetOptions (	'outputdir=s' => \$outputdir,	# outputdir
		'logs=s' => \$logs,		# logs
		'paramfile=s' => \$pfile,	# parameter
		'R1=s' => \$r1,			# R1
		'R2=s' => \$r2,			# R2
		'sample=s' => \$sample,		# sample
        'contig=i' => \$contig, 			# contig bool
		'timestamp=s' => \$timestamp);	# time stamp
            
die "[error] required input parameters not found" if (!( defined($outputdir) && defined($logs) && defined($pfile) && defined($r1) ));

@mates=("contig") if ($contig);  # Reset mates to the name contig
            
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
	$hoh{$command}{$mates[0]} = abs_path($r1);
}

if ( $r2 ne "none" && defined($r2) )
{
	$hoh{$command}{$mates[1]} = abs_path($r2);
} 	

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
        my $mate_rel = File::Spec->abs2rel($hoh{$command}{$mate});
		verbose_system("ln -sf $mate_rel");		

		if ($hoh{$command}{"getorf_options"})
		{
            print "[echo] Count input\n";
            verbose_system("linecount $hoh{$command}{$mate} input $mate.count fastq 0");
			print "[echo] get orf $mate \n";
			# E.g., getorf -sequence 3.contig.noblast.fasta -outseq testtest -minsize 60 -find 0 
			my $cmd = "getorf -sequence $hoh{$command}{$mate} -outseq $mate.orfout.fa $hoh{$command}{\"getorf_options\"}";	
			verbose_system($cmd);

			# my $cmd = "cat $mate.orfout.fa | fastajoinlines > $mate.orfout.join.fa";
			# verbose_system($cmd);
			
			# my $href = fastaid_firstword_hash("$mate.orfout.fa");
			# print Dumper $href;

            # orf fasta file identifiers look like this
            # >6_1 [109 - 192]
            # So make a hash of all the unique first digits
			#my %h_fasta = map {/>(\w+)_(\w+)\s(.*)/; $1 => 1} split(/\n/, `cat $mate.orfout.fa`);
            # Should be >(anything)_\d\s(.*
			my %h_fasta = map {/>(.*?)_(\d+)\s(.*)/; $1 => 1} split(/\n/, `cat $mate.orfout.fa`);
			#print Dumper \ %h_fasta;

			print "[echo] filter $hoh{$command}{$mate} by orfs into $command.$mate\n";
			get_subset_by_fastaid($hoh{$command}{$mate}, "$command.$mate", \%h_fasta);
            verbose_system("linecount $command.$mate orf_filter $mate.count fastq 1");
		}
        else {
            print "Doing nothing because there are no getorf_options set\n";
        }
	} # defined
	elsif ($mate eq "R1")
	{
        if(! -s $hoh{$command}{$mate} ) {
            print("$hoh{$command}{$mate} is empty\n");
            # just create an empty output file for this step then
            print("Creating empty $command.$mate\n");
            open(my $fh, '>', $command.".".$mate);
            close($fh);
        } else {
            print "$mate not defined\n";
        }
	}
} # mate

print "[END]\n";
# print "[hash] ";
# print Dumper \ %hoh;

close(STDOUT);
close(STDERR);  
close($olog);
close($elog);

# return hash of the IDs of a fasta file to "1"
sub fastaid_firstword_hash  
{
	# arg 1 - file name

	my $infile = shift;
	
	my %h = ();				# hash	

	if ( -s $infile )		# if file nonzero
	{	
		open(my $fh, '<', $infile);
	
		while (<$fh>)
		{
			# ID is every line starting with ">"
			if ( $_ =~ m/^>/ )
			{
				chomp $_;
				
				# header looks like >c3_1 [3 - 89] 
				if ($_ =~ m/>(\S+)(\s+)(.*)/)
				{
					# don't print leading ">"
					my $key = $1;
							
					# hash key is fastq id
					$h{$key} = 1;
				}	
			}
		}
			
		close($fh);
	}

	return \%h;
}

# get subset of file
sub get_subset_by_fastaid
{
	my $infile = shift;
	my $outfile = shift;
	my $hash_fasta_ref = shift;
	
	my %h = %$hash_fasta_ref;	# hash	

    # Default to fasta but detect later
    my $format = "fasta";

	# print Dumper \ %h;

	if ( -s $infile )		# if file nonzero
	{	
        # Detect fasta or fastq file
		open(my $fh, '<', $infile);
        $_ = <$fh>;
        if($_ =~ /^@/) {
            $format = "fastq";
            print("Detected $infile as fastq\n");
        } else {
            print("Detected $infile as fasta\n");
        }
        close($fh);

		open(my $fh, '<', $infile);
		open(my $fh2, '>', $outfile);
	
		while (<$fh>)
		{
            chomp $_;
            
            # parse out id line
            $_ =~ m/[>@](\S+)(.*)/;
            #print("$format ID: $_\n");
            # don't print leading ">"
            my $key = $1;
            # print ($key,"\n");			
            if ($h{$key})
            {
                #print("Keeping $key\n");
                # Print identifier line
                print $fh2 $_,"\n";
                # Print sequence line
                $_ = <$fh>;
                print $fh2 $_;
                if($format == "fastq") {
                    # Input file is fastq so read 4 lines total
                    $_ = <$fh>;
                    print $fh2 $_;
                    $_ = <$fh>;
                    print $fh2 $_;
                }
			} else {
                # Not keeping this sequence
                # need to burn off the sequence
                $_ = <$fh>;
                if($format == "fastq") {
                    # Fastq so burn off + line and qual line
                    $_ = <$fh>;
                    $_ = <$fh>;
                }
            }
		}
			
		close($fh);
		close($fh2);
	}

	return \%h;
}
