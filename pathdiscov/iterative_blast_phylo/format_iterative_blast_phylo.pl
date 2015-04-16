#!/usr/bin/perl

use strict;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;

# stitch the report files in the iterative_blast dir into one unified piece

# example
# format_iterative_blast_phylo.pl --outputdir tmp --prefix sample --blastdir iterative_blast_phylo_1 --blast_list megablast,dc-megablast,blastx

# -----------------------------------

# args 
# $ARGV[0] 
# $ARGV[1] 
# ...

# learning perl ! 

# \d [0-9] Any digit
# \D [^0-9] Any character not a digit
# \w [0-9a-zA-z_] Any "word character"
# \W [^0-9a-zA-z_] Any character not a word character
# \s [ \t\n\r\f] whitespace (space, tab, newline, carriage return, form feed)
# \S [^ \t\n\r\f] Any non-whitespace character

# *      Match 0 or more times
# +      Match 1 or more times
# ?      Match 1 or 0 times
# {n}    Match exactly n times
# {n,}   Match at least n times
# {n,m}  Match at least n but not more than m times

# -----------------------------------

# declare vars
my ($outputdir);		# output dir
my ($blastdir);			# dir where blast files reside
my ($blist);			# blast task list (e.g., megablast, blastx)
my ($prefix);			# prefix name
my ($bigreport);		# a bool that's "1" if you want to make the big report (default: "0")
my ($header);			# string for report header 
my ($header_phylo);		# string for phylo header
my @mates=("contig","R1","R2"); # strings for mates 

GetOptions ('outputdir=s' => \$outputdir,		# output directory
	    'blastdir=s' => \$blastdir,			# blast directory
	    'bigreport' => \$bigreport,			# a boolean to control whether or not to make the full report (it's large - usu unnec - the top report should suffice)
	    'prefix=s' => \$prefix,			# output files name prefix
            'blast_list=s' => \$blist);			# blast algorithm list
            
die "[error] required input parameters not found" if (!( defined($outputdir) && defined($blastdir) && defined($blist) ));
           
# get abs path to directory in which script resides:
# my $mytmp=`readlink -m $0`; 					# perl's readlink function is lousy, so use bash
# my $path_scripts=`dirname $mytmp`; 
# chomp $path_scripts;

# get absolute paths
my $mytmp=`readlink -m $outputdir`;
chomp $mytmp;
$outputdir = $mytmp;

my $mytmp=`readlink -m $blastdir`;
chomp $mytmp;
$blastdir = $mytmp;

my @blast_task_list = split(',',$blist);

chdir($outputdir) or die "[error] cannot chdir to ".$outputdir;

# --------------------------------------

foreach my $mate (@mates) 
{	
	if ( -e "$blastdir/1.$mate.fasta")
	{
		
		print "***",$mate,"***","\n";
	
		open(my $outfile, ">", "$outputdir/$mate.$prefix.report.txt") if ($bigreport);	# report
		open(my $outfile_top, ">", "$outputdir/$mate.$prefix.top.report.txt");			# top hit report
		open(my $outphylo, ">", "$outputdir/$mate.$prefix.phylo.txt");					# phylo report
		open(my $outphylo_top, ">", "$outputdir/$mate.$prefix.top.phylo.txt");			# phylo top hit report 
		
		# again, not great style, because header is dependent on blast format
		$header="qseq\tblast_alg\tqseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\ttaxid\tsuperkingdom\torder\tfamily\tgenus\tdescrip\n";
		
		$header_phylo="blast_alg\ttaxid\tcount\tsuperkingdom\tkingdom\tclass\torder\tfamily\tgenus\tspecies\tdescrip\n";	
		
		print $outfile $header if ($bigreport);
		print $outfile_top $header;	
		print $outphylo $header_phylo;
		print $outphylo_top $header_phylo;	
			
		open(my $infile2, "<", "$blastdir/1.$mate.fasta"); 	
		
		# map the qid to the fasta seq
		my %h2=();
		while (<$infile2>)
		{
			# dont chomp - the reason is that for IDs like ">3" the regex will chomp for you --- ">3\n" matches m/>(\S*)(\s)(.*)/
			
			# this ID will be the key to your hash
			my $mykey=$_;
	
			# format is ">contig-0 5828 nucleotides" or ">3"
			$mykey =~ m/>(\S*)(\s)(.*)/;
			$mykey = $1;
		
			$_ = <$infile2>;
			chomp($_);
		
			$h2{$mykey}=$_;
		}	
		
		# loop thro blast iterations
		for (my $i = 0; $i < scalar(@blast_task_list); $i++) 
		{
			my $j = $i + 1;			
			print "[echo] $j.$mate.blast \n";
	
			if ( -s "$blastdir/$j.$mate.blast")
			{
				print "[echo] $blast_task_list[$i] \n";	
			
				# map taxid to phylo stuff
				open(my $infile3, "<", "$blastdir/$j.$mate.blast.phylo"); 				
				my %h3=();
				my $skipline=1;	# skip the first header line
				while (<$infile3>)
				{
					print $outphylo $blast_task_list[$i],"\t",$_ if ($skipline==0);					
					$skipline=0;
						
					# $infile3 - e.g., 
					# taxid	count	superkingdom	kingdom	class	order	family	genus	species
					# 565995	49.00	Viruses	-	-	Mononegavirales	Filoviridae	Ebolavirus	Bundibugyo ebolavirus	
					# 11036	3.00	Viruses	-	-	-	Togaviridae	Alphavirus	Venezuelan equine encephalitis virus
		
					my @ln=split; 
					# map: taxid --> superkingdom order family genus 
					$h3{$ln[0]}=$ln[2]."\t".$ln[5]."\t".$ln[6]."\t".$ln[7];			
				}
				close($infile3);			
	
				# now write to top file
				my $skipline=1;	# skip the first header line
				open(my $infile3, "<", "$blastdir/$j.$mate.top.blast.phylo"); 				
				while (<$infile3>)
				{
					print $outphylo_top $blast_task_list[$i],"\t",$_ if ($skipline==0);
					$skipline=0;										
				}
				close($infile3);
				
				# print Dumper \ %h3;
		
				if ($bigreport)		
				{	
					print "$j.$mate.blast\n";	
			
					open(my $infile1, "<", "$blastdir/$j.$mate.blast"); 
					open(my $infile4, "<", "$blastdir/$j.$mate.blast.ann"); 			
					
					if ( -s "$blastdir/$j.$mate.blast" ) # if $infile1 is nonzero
					{
						while (<$infile1>)
						{	
							chomp($_);
							my $f1ln=$_;			# line from $j.$mate.blast
							
							$_ = <$infile4>;
							chomp($_);			
							my $f2ln=$_;			# line from $j.$mate.blast.ann (USE THE FACT $j.$mate.blast and $j.$mate.blast.ann MUST HAVE SAME # OF LINES)
							
							# $infile4 - e.g., 
							# contig-0	10844	Bacteriophage S13 circular DNA, complete genome
							# contig-1000000	11036	Venezuelan equine encephalitis virus strain V3526, complete genome
							
							$f2ln =~ m/(\S*)(\s)(\S*)(\s*)(.*)/;
							my $taxid = $3;
							my $descrip = $5;			
							$descrip =~ s/ /_/g;
							
							my @ln=split(/\s/, $f1ln);			
				
							my ($field1, $field2);
				
							# look for qid to get full sequence
							if ($h2{$ln[0]})
							{
								$field1=$h2{$ln[0]};
							}
							else
							{
								$field1="-";	
							}
				
							# look for taxid
							if ($h3{$taxid})
							{
								$field2=$h3{$taxid};
							}
							else
							{
								$field2="-"."\t"."-"."\t"."-"."\t"."-";	
							}
				
							print $outfile $field1,"\t",$blast_task_list[$i],"\t",$f1ln,"\t",$taxid,"\t",$field2,"\t",$descrip,"\n";		
						} # loop thro blast output file
					} # nonzero
		
					close($infile1);
					close($infile4);
				}
				
				# could do this with one line of awk but ... do the same for top.blast
							
				print "$j.$mate.top.blast\n";
		
				open(my $infile1_top, "<", "$blastdir/$j.$mate.top.blast"); 
				open(my $infile4_top, "<", "$blastdir/$j.$mate.top.blast.ann"); 
	
				if ( -s "$blastdir/$j.$mate.top.blast" ) # if $infile1 is nonzero
				{
					while (<$infile1_top>)
					{	
						chomp($_);
						my $f1ln=$_;			# line from $j.$mate.blast
						
						$_ = <$infile4_top>;
						chomp($_);			
						my $f2ln=$_;			# line from $j.$mate.blast.ann (USE THE FACT $j.$mate.blast and $j.$mate.blast.ann MUST HAVE SAME # OF LINES)
						
						# $infile4 - e.g., 
						# contig-0	10844	Bacteriophage S13 circular DNA, complete genome
						# contig-1000000	11036	Venezuelan equine encephalitis virus strain V3526, complete genome
						
						$f2ln =~ m/(\S*)(\s)(\S*)(\s*)(.*)/;
						my $taxid = $3;
						my $descrip = $5;			
						$descrip =~ s/ /_/g;
						
						my @ln=split(/\s/, $f1ln);			
			
						my ($field1, $field2);
			
						# look for qid to get full sequence
						if ($h2{$ln[0]})
						{
							$field1=$h2{$ln[0]};
						}
						else
						{
							$field1="-";	
						}
			
						# look for taxid
						if ($h3{$taxid})
						{
							$field2=$h3{$taxid};
						}
						else
						{
							$field2="-"."\t"."-"."\t"."-"."\t"."-";	
						}
			
						print $outfile_top $field1,"\t",$blast_task_list[$i],"\t",$f1ln,"\t",$taxid,"\t",$field2,"\t",$descrip,"\n";		
					} # loop thro blast output file
				} # nonzero				 		 
		
				close($infile1_top);	
				close($infile4_top);	
		
			} # blast file nonzero
		} # blast iteration
	
		close($infile2);
			
		close($outfile) if ($bigreport); 	
		close($outfile_top);
		close($outphylo); 	
		close($outphylo_top);

		# this is hack-ish: cut off space-consuming seq cols to make smaller report
		system("cat $outputdir/$mate.$prefix.report.txt | cut -f2-16,18- > $outputdir/$mate.$prefix.smallreport.txt") if ($bigreport);		# report
		system("cat $outputdir/$mate.$prefix.top.report.txt | cut -f2-16,18- > $outputdir/$mate.$prefix.top.smallreport.txt");			# top hit report		
				
	} # mate exists
} # mate
