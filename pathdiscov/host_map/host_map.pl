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

# example
# host_map.pl --sample $sample --paramfile $pfile --outputdir $path_output/host_map --logs $path_output/host_map/logs --R1 input_R1.fastq --R2 input_R2.fastq --timestamp 23672

# TO DO:
# deal w fasta input

# declare vars
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
my ($wellpaired);			# wellpaired - 1 if paired data "well-paired" (i.e., mate files have exact same IDs), 0 if not; default: well-paired. The idea is, if we know data is well paired it's faster b/c no need to check
my $path_scripts=$RealBin;	# get abs path to directory in which script resides
my @mates=("R1","R2");		# strings for mates

# default
$timestamp="0";
$run_iteration=1;
$wellpaired=1;			# wellpaired - 1 if paired data "well-paired" (i.e., mate files have exact same IDs), 0 if not; default: well-paired

GetOptions ('outputdir=s' => \$outputdir,			# outputdir
			'logs=s' => \$logs,						# logs
            'paramfile=s' => \$pfile,				# parameter
            'R1=s' => \$r1,							# R1
            'R2=s' => \$r2,							# R2
            'sample=s' => \$sample,            		# sample            
            'timestamp=s' => \$timestamp,      		# time stamp
            'fastafile=s' => \$fastafile,			# is fasta (default assumes fastq)            
            'fasta' => \$isfasta,            		# is fasta (default assumes fastq). Note: $isfasta & $fastafile do the same thing. The former is for use from the command line, while the later is useful because it's not a bool)
            'run_iteration=i' => \$run_iteration,  	# global # of times you run the script
            'wellpaired=i' => \$wellpaired,  		# 1 if R1 and R2 are well-paired (i.e., they have the same read IDs), 0 if not (i.e., the files are asymmetrical)                          
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
# mapper_program_list		(e.g., "bwa,bowtie2,bowtie2")
# mapper_db_list			
# mapper_name_list			(e.g., bwa_genome,bowtie2_genome,bowtie2_transcriptome)
# mapper_options_list		(e.g., ",--local,--local")

# --------------------------------------

# get absolute path to outputdir
$outputdir=abs_path($outputdir);

# set command
if ($run_iteration>1)
{	
	$command="host_map_".$run_iteration;	
}
else
{
	$command="host_map";
}

# get arrays
my @mapper_program_list = split(',',$hoh{$command}{"mapper_program_list"});
my @mapper_db_list = split(',',$hoh{$command}{"mapper_db_list"});
my @mapper_name_list = split(',',$hoh{$command}{"mapper_name_list"});
my @mapper_options_list = split(',',$hoh{$command}{"mapper_options_list"});

# die "[error] param lists different sizes" if (!( scalar(@mapper_db_list)==scalar(@mapper_program_list) && scalar(@mapper_db_list)==scalar(@mapper_name_list) && scalar(@mapper_db_list)==scalar(@mapper_options_list) ));
die "[error] param lists different sizes" if (!( scalar(@mapper_db_list)==scalar(@mapper_program_list) ));


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
	$hoh{$command}{"R1input"} = $abs_r1;	# not the best style, but copy the input entry in the hash for use at the end
}

if ( $r2 ne "none" && defined($r2) )
{
	$abs_r2 = abs_path($r2);
	$hoh{$command}{"R2"}=$abs_r2;
	$hoh{$command}{"R2input"} = $abs_r2;		
} 

die "[error] input file not found" if (!( -e $hoh{$command}{"R1"} ));

print "[START]\n";
print "[hash] ";
print Dumper \ %hoh;

chdir($outputdir) or die "[error] cannot chdir to $outputdir";

# --------------------------------------
# main

# if input fasta file, link to it. O.w., convert fastq to fasta
if ($isfasta || $fastafile eq "yes")
{
	print "[error] currently, this module cannot handle fasta format\n";
	exit;		
}

print "[echo] run a series of alignments\n";

foreach my $mate (@mates) 
{
	# check if defined, and non-zero
	if ( defined($hoh{$command}{$mate}) && -s $hoh{$command}{$mate} )
	{
		# count lines in file
		# args: input file, filtering_program_name, output file, 1->fastq, concat
		my $cmd = "linecount $hoh{$command}{$mate} input $mate.count 1 0";
		print_system($cmd);	
	}
}

for (my $i = 0; $i < scalar(@mapper_db_list); $i++) 
{
	print "[echo] $mapper_name_list[$i] \n";	
	
	my $j = $i + 1;
	my $k = $i + 2;

	print "[iteration] $j\n";	
	
	# make tmp dirs for each iteration
	my $cmd = "mkdir -p map_$j";
	print_system($cmd);	
	
    # Pull out aligner
    my $aligner = $mapper_program_list[$i];
    if($aligner eq "bowtie2" || $aligner eq "snap") {
        print "[echo] Aligner to use: $aligner\n";
    } else {
        print "[ERROR] Aligner specified $aligner is not a valid mapper program\n";
    }
	
	# - - - - - - - 
	# do the mapping 
	
	# paired data
	if ( defined($hoh{$command}{"R2"}) && -s $hoh{$command}{"R2"} )
	{
        print "[echo] Paired data detected\n";
        print "[echo] wellpaired is set as: $wellpaired\n";
		# well paired - the idea is, if we know well paired in advance, get speed up b/c dont have to check
		if ($wellpaired)
		{
            print "[echo] Doing wellpaired mapping\n";
			# map
            my $cmd = "";
            if($aligner eq "bowtie2") {
                $cmd = "bowtie2 -q -x $mapper_db_list[$i] -1 $hoh{$command}{\"R1\"} -2 $hoh{$command}{\"R2\"} -S map_$j/out.sam $mapper_options_list[$i]";
            } elsif($aligner eq "snap") {
                $cmd = "snap paired $mapper_db_list[$i] $hoh{$command}{\"R1\"} $hoh{$command}{\"R2\"} -o map_$j/out.sam $mapper_options_list[$i]";
            }
			verbose_system($cmd);
			
			# get unmapped read IDs
			my $cmd = "cat map_$j/out.sam | get_map_unmap.pl --R1map map_$j/R1.map.id --R2map map_$j/R2.map.id --R1unmap map_$j/R1.unmap.id --R2unmap map_$j/R2.unmap.id --paired --printsid";
			verbose_system($cmd);
			
			# get unmap reads
			my $cmd="fastq_extract_id.pl $hoh{$command}{\"R1\"} map_$j/R1.unmap.id > map_$j/R1.unmap.fastq";
			print_system($cmd);

			# get unmap reads
			my $cmd="fastq_extract_id.pl $hoh{$command}{\"R2\"} map_$j/R2.unmap.id > map_$j/R2.unmap.fastq";
			print_system($cmd);	
    
            # The next time we map, likely that the input files are no longer wellpaired
            $wellpaired=0;
		}
		# paired but not well paired
		else
		{
            print "[echo] Doing non wellpaired mapping\n";
			# get well paired reads and non-well paired reads
			my $cmd="perlscripts_wrapper.pl get_common_uneven_files $hoh{$command}{\"R1\"} $hoh{$command}{\"R2\"} map_$j/R1.single.fastq map_$j/R1.paired.fastq map_$j/R2.single.fastq map_$j/R2.paired.fastq";
			verbose_system($cmd);

            if($aligner eq "bowtie2") {
                # map
                # check what's zero and what's non zero and run bowtie accordingly
                if ( -s "map_$j/R1.single.fastq" && -s "map_$j/R2.single.fastq" )
                {
                    $cmd="bowtie2 -q -x $mapper_db_list[$i] -1 map_$j/R1.paired.fastq -2 map_$j/R2.paired.fastq -U map_$j/R1.single.fastq,map_$j/R2.single.fastq -S map_$j/out.sam $mapper_options_list[$i]";
                }
                elsif ( -s "map_$j/R1.single.fastq" )
                {
                    $cmd="bowtie2 -q -x $mapper_db_list[$i] -1 map_$j/R1.paired.fastq -2 map_$j/R2.paired.fastq -U map_$j/R1.single.fastq -S map_$j/out.sam $mapper_options_list[$i]";
                }
                elsif ( -s "map_$j/R2.single.fastq" )
                {
                    $cmd="bowtie2 -q -x $mapper_db_list[$i] -1 map_$j/R1.paired.fastq -2 map_$j/R2.paired.fastq -U map_$j/R2.single.fastq -S map_$j/out.sam $mapper_options_list[$i]";
                }
                else
                {
                    $cmd="bowtie2 -q -x $mapper_db_list[$i] -1 map_$j/R1.paired.fastq -2 map_$j/R2.paired.fastq -S map_$j/out.sam $mapper_options_list[$i]";	
                }
            } elsif($aligner eq "snap") {
                # Just map paired data and let snap handle non-wellpaired data
                $cmd = "snap paired $mapper_db_list[$i] $hoh{$command}{\"R1\"} $hoh{$command}{\"R2\"} -o map_$j/out.sam $mapper_options_list[$i]";
            }
			verbose_system($cmd);
			
			# get unmapped read IDs
			my $cmd = "cat map_$j/out.sam | get_map_unmap.pl --R1map map_$j/R1.map.id --R2map map_$j/R2.map.id --R1unmap map_$j/R1.unmap.id --R2unmap map_$j/R2.unmap.id --singletonmap map_$j/singleton.map.id --singletonunmap map_$j/singleton.unmap.id --paired --unpaired --printsid";
			verbose_system($cmd);
			
			# get unmap reads (R1, paired)
			my $cmd="fastq_extract_id.pl $hoh{$command}{\"R1\"} map_$j/R1.unmap.id > map_$j/R1.unmap.fastq";
			print_system($cmd);

			# get unmap reads (R2, paired)
			my $cmd="fastq_extract_id.pl $hoh{$command}{\"R2\"} map_$j/R2.unmap.id > map_$j/R2.unmap.fastq";
			print_system($cmd);

			# get unmap reads (R1, singleton)
			my $cmd="fastq_extract_id.pl $hoh{$command}{\"R1\"} map_$j/singleton.unmap.id >> map_$j/R1.unmap.fastq";
			print_system($cmd) if ( -s "map_$j/singleton.unmap.id" );

			# get unmap reads (R2, singleton)
			my $cmd="fastq_extract_id.pl $hoh{$command}{\"R2\"} map_$j/singleton.unmap.id >> map_$j/R2.unmap.fastq";
			print_system($cmd) if ( -s "map_$j/singleton.unmap.id" );	
		}
	}
	# non paired data
	else
	{
        print "[echo] Doing single read mapping\n";
		# map
        my $cmd = "";
        if($aligner eq "bowtie2") {
            $cmd = "bowtie2 -q -x $mapper_db_list[$i] -U $hoh{$command}{\"R1\"} -S map_$j/out.sam $mapper_options_list[$i]";
        } elsif($aligner eq "snap") {
            $cmd = "snap single $mapper_db_list[$i] $hoh{$command}{\"R1\"} $hoh{$command}{\"R2\"} -o map_$j/out.sam $mapper_options_list[$i]";
        }
		verbose_system($cmd);
		
		# get unmapped read IDs
		my $cmd = "cat map_$j/out.sam | get_map_unmap.pl --singletonmap map_$j/singleton.map.id --singletonunmap map_$j/singleton.unmap.id --unpaired --printsid";
		verbose_system($cmd);
		
		# get unmap reads (R1, singleton)
		my $cmd="fastq_extract_id.pl $hoh{$command}{\"R1\"} map_$j/singleton.unmap.id > map_$j/R1.unmap.fastq";
		print_system($cmd);		
	}
	
	my $cmd="samtools view -bS map_$j/out.sam > map_$j/out.bam";
	verbose_system($cmd);		

	my $cmd="rm map_$j/out.sam";
	system($cmd);

	print "[stats] paired stats\n";
	my $cmd="samtools flagstat map_$j/out.bam";	
	print_system($cmd);
	# - - - - - - - 	
	
	# make links and do read counting 
	foreach my $mate (@mates) 
	{
		# check if defined, and non-zero
		if ( defined($hoh{$command}{$mate}) && -s $hoh{$command}{$mate} )
		{
			# update mate to new output
			$hoh{$command}{$mate} = "$outputdir/map_$j/$mate.unmap.fastq.tmp";										
		
			my $cmd = "linecount $hoh{$command}{$mate} $mapper_name_list[$i] $mate.count 1 1";
			print "[cmd] ",$cmd,"\n";
			system($cmd);
			
			# (take this out of loop for efficiency - only needs to be done on the final iteration)
			# get discard IDs - i.e., the reads filtered in this step
			my $cmd = "fastaq_tools_diff.exe --fastq $hoh{$command}{$mate.\"input\"} --fastq $hoh{$command}{$mate} | sort > $mate.discard";
			verbose_system($cmd);
		}
	}

    my $cmd = "drop_mapped --saved $hoh{$command}{"R1"} --dropped "R2".discard --out $mapout" 

    verbose_system($cmd);

    my $cmd = "drop_mapped --saved $hoh{$command}{"R2"} --dropped "R1".discard --out $mapout" 

    verbose_system($cmd);

  foreach my $mate (@mates) 
  {
	system("ln -sf map_$j/$mate.unmap.fastq $outputdir/host_map_$run_iteration.$mate");
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

host_map.pl - run a number, specified by the user, of chained alignments to a host. At each stage, only the sequences that did not map to a host pass to the next stage.

=head1 SYNOPSIS

host_map.pl [B<--sample> sample] [B<--outputdir> output_directory] [B<--logs> log_directory] [B<--paramfile> parameter_file(s)]  [B<--R1> R1_input] [B<--timestamp> a time stamp] [options] 

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

host_map.pl takes a series of indexed db's and runs a series of alignments on them. For example, you might run bowtie2 with db=genome, followed bowtie2 with db=transcriptome. At each stage, the input is what didn't align in the previous stage.

=head1 EXAMPLE

run "host_map" with paired end .fastq files:

B<host_map.pl --sample sample --paramfile host_map.param --outputdir /my/output/directory --logs /my/logs/directory --R1 input_R1.fastq --R2 input_R2.fastq --timestamp 121207-11.30>

=cut
