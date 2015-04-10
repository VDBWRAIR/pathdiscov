#!/usr/bin/perl

use strict;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use POSIX;
use Cwd 'abs_path';

use FindBin qw($RealBin);
use lib "$RealBin/../Local_Module";
# local modules:
use Verbose_Sys;
use Parse_ParameterFile;

# example
# sge_block_blast.pl --inputfasta x --outputdir x --outfile x --outheader x --scripts x --ninst 5 --db x --blast_type x --task x --blast_options x --qtime x --qmem x --qrls x

# declare vars
my ($outputdir, $inputfasta, $outfile, $outheader, $path_scripts, $j, $r, $ninst, $db, $blast_type, $task, $blast_options, $qtime, $qmem, $qrls);

GetOptions ('outputdir=s' => \$outputdir,					# outputdir
			'inputfasta=s' => \$inputfasta,					# input fasta file
			'outfile=s' => \$outfile,						# output file
			'outheader=s' => \$outheader,					# output header
			'scripts=s' => \$path_scripts,					# path to scripts
			'iteration_j=i' => \$j,							# j
            'mate=s' => \$r,								# mate ("R1" or "R2")												
            'ninst=i' => \$ninst,							# # instsance (or qsub) 
            'db=s' => \$db,									# blast db
            'blast_type=s' => \$blast_type,					# blast type
            'task=s' => \$task,            					# blast task
            'blast_options=s' => \$blast_options,			# blast options
            'qtime=s' => \$qtime,							# qsub time
            'qmem=s' => \$qmem,								# qsub mem
            'qrls=s' => \$qrls);							# qrls id            

# it's better to pass to scripts as a variable b/c o.w., you risk an error "/opt/gridengine/default/spool/b15b/job_scripts/blast_wrapper.pl not found" if run as a job
if (!defined($path_scripts))
{
	# get abs path to directory in which script resides:
	# perl's readlink function is lousy, so use bash
	my $mytmp=`readlink -m $0`;
	$path_scripts=`dirname $mytmp`; 
	chomp $path_scripts;
}

# get absolute path to outputdir
my $mytmp=`readlink -m $outputdir`;
chomp $mytmp;
$outputdir = $mytmp;


die "[error] input file not found" if (!( -e $inputfasta ));

if (!( -s $inputfasta )) # if not nonzero
{
	print "\n[echo] file $inputfasta empty.\n";
	exit;	            
}

chdir($outputdir) or die "[error] cannot chdir to ".$outputdir;

# ----------------

my $prefix = "tmpsplit";
my $minchunk = 20;

my $len=`cat $inputfasta | wc -l`;
chomp($len);

die "[error] number of lines in fasta file must be even (use fastajoinlines if nec)" if ( $len % 2 == 1 );
die "[error] attempting to run too many instances" if ( $ninst > 1000 );

# assume even input --- get size of split files + remainder
my $chunk = floor($len/$ninst);
my $remainder = $len % $ninst;

print "wc_input_length = $len\n";
print "instances = $ninst\n";
print "wc_split_length = $chunk\n";
print "wc_remainder = $remainder\n";

if ( $len <= $minchunk )		# if file size less than or equal to min chunk  
{
	print "\nlength less than min chunk - run a single instance\n\n";
	$chunk = $len;
	$remainder = 0;
	$ninst = 1;		
	
	print "wc_input_length = $len\n";
	print "instances = $ninst\n";
	print "wc_split_length = $chunk\n";
	print "wc_remainder = $remainder\n";		
}
elsif ( $chunk < $minchunk )	# if chunk too small  
{
	print "\nchuck too small, decrement #instance\n\n";
	
	# loop through split files
	while ( $len/$ninst < $minchunk ) 
	{		
	    $ninst--;
	}

	$chunk = floor($len/$ninst);
	$remainder = $len % $ninst;		
	
	print "wc_input_length = $len\n";
	print "instances = $ninst\n";
	print "wc_split_length = $chunk\n";
	print "wc_remainder = $remainder\n";		
}

# if chunk odd, unacceptable because fasta file - must keep even. Hence, increment chunk size to make it even. 
if ( $chunk % 2 == 1 )
{
	print "\nchuck odd, increment\n\n";
	$chunk++; 
	$ninst = floor($len/$chunk);
	$remainder = $len % $chunk;	
	
	print "wc_input_length = $len\n";
	print "instances = $ninst\n";
	print "wc_split_length = $chunk\n";
	print "wc_remainder = $remainder\n";	
		
}

# check that assumptions still satisfied
die "[error] chunk must be even" if ( $chunk % 2 == 1 );
die "[error] chunk remainder must be even" if ( $remainder % 2 == 1 );

`split -a 3 -d -l $chunk $inputfasta $prefix`;

# if remainder is greater than 0, add another instance for the subfile that contains the remainder 
if ( $remainder != 0 )
{
	print "\nremainder nonzero, increment ninst\n\n";
	$ninst++;

	print "wc_input_length = $len\n";
	print "instances = $ninst\n";
	print "wc_split_length = $chunk\n";
	print "wc_remainder = $remainder\n";
	
}

my $i = 0;

print "\n";

my $jid_list="";

# loop through split files
while ( $i < $ninst ) 
{
	
	# pad i to match files
	my ($suffix);	
	if (length($i) == 1 )
	{
		$suffix="00".$i;
	}
	elsif (length($i) == 2 )
	{
		$suffix="0".$i;
	}	
	elsif (length($i) == 3 )
	{
		$suffix=$i;		
	}	
	
# 	print $prefix.$suffix,"\n";
	
	my $cmd = "$path_scripts/blast_wrapper.pl --type ${blast_type} --query $prefix$suffix --db ${db} --task ${task} --out blastout${suffix} --options \"${blast_options}\"";

	my $qcmd="qsub -V -N blast".$i.".".$r.".".$j." -e ../logs_blast -o ../logs_blast -l mem=".$qmem."G,time=".$qtime.":: -S /usr/bin/perl -cwd $cmd";	
	my $tmpjid=qsub_system($qcmd);
	$jid_list=$jid_list.$tmpjid.",";	
	
    $i++;
}

print "jids: ",$jid_list,"\n";

my $qcmd="qsub -V -hold_jid $jid_list -N blastfin.".$r.".".$j." -e ../logs_blast -o ../logs_blast -l mem=1G,time=1:: -S /bin/sh -cwd $path_scripts/job_release_concat_blast.sh $qrls $outfile";	
print "[qcmd] ",$qcmd,"\n";
my $qsub_message=`$qcmd`;	
print($qsub_message); 

