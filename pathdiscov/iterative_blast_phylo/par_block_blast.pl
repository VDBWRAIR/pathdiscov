#!/usr/bin/perl

use strict;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use POSIX;

# example
# par_block_blast.pl --inputfasta x --outputdir x --outfile x --outheader x --ninst 5 --db x --blast_type x --task x --blast_options x

# to do: delete logsdir or use it

# declare vars
my ($outputdir, $logsdir, $inputfasta, $outfile, $outheader, $ninst, $db, $blast_type, $task, $blast_options);

GetOptions ('outputdir=s' => \$outputdir,					# outputdir
			'logs=s' => \$logsdir,							# logs dir
			'inputfasta=s' => \$inputfasta,					# input fasta file
			'outfile=s' => \$outfile,						# output file
			'outheader=s' => \$outheader,					# output header			
            'ninst=i' => \$ninst,							# # instsance (or qsub) 
            'db=s' => \$db,									# blast db
            'blast_type=s' => \$blast_type,					# blast type
            'task=s' => \$task,            					# blast task
            'blast_options=s' => \$blast_options);			# blast options

# get abs path to directory in which script resides:
my $mytmp=`readlink -m $0`; 					# perl's readlink function is lousy, so use bash
my $path_scripts=`dirname $mytmp`; 
chomp $path_scripts;

# get absolute path to outputdir
my $mytmp=`readlink -m $outputdir`;
chomp $mytmp;
$outputdir = $mytmp;


die "[error] input file not found" if (!( -e $inputfasta ));

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

my @children = ();	# pids of children
my ($pid);

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

	$pid = fork();
	 
	if ( $pid ) # pid nonzero
	{
		# parent
		print "[echo] parent process: submit child PID $pid \n";
		push @children, $pid;
	}
	elsif ( $pid == 0 ) 
	{
		# child
		# from the web: "real action happens here"
#		print "[echo] child process: child PID $pid \n";
		my ($cmd);
		$cmd = "$path_scripts/blast_wrapper.pl --type ${blast_type} --query $prefix$suffix --db ${db} --task ${task} --out blastout${suffix} --options \"${blast_options}\"";
		print "[cmd] ",$cmd,"\n";
		system($cmd);
		exit 0;
	}
	else 
	{
	     die "Couldn't fork: $!\n";
	}
			
    $i++;
}

foreach my $child (@children) 
{
#	print "PID ",$pid,": wait on ".$child."\n";
	print "[echo] wait on ".$child."\n";		
	waitpid($child,0);
	system("cat blastout* > $outfile");
}

print "[echo] all processes done\n";	
