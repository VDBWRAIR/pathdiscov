#!/usr/bin/perl

# $ARGV[0] input 
# $ARGV[1] output

# takes input like this:

# contig-0	10844	Bacteriophage S13 circular DNA, complete genome
# contig-0	10844	Bacteriophage S13 circular DNA, complete genome
# contig-0	10847	Coliphage phi-X174, complete genome
# contig-0	10847	Enterobacteria phage phiX174 strain beta4, complete genome
# contig-0	10847	Enterobacteria phage phiX174 isolate XC+mbD10im2, complete genome
# contig-1000000	11036	Venezuelan equine encephalitis virus strain TC-83, complete genome
# contig-1000000	364092	VEEV replicon vector YFV-C3opt, complete sequence
# contig-1000000	364089	VEEV replicon vector YFV-C2opt, complete sequence
# contig-1000000	364086	VEEV replicon vector YFV-C2, complete sequence
# contig-1000000	364087	VEEV replicon vector YFV-C1, complete sequence
# contig-1000000	364086	VEEV replicon vector YFV-C2, complete sequence
# contig-1000000	11036	Venezuelan equine encephalitis virus, complete genome
# contig-1000000	364091	VEEV replicon vector YFV-C2opt-NS2mut, complete sequence
# contig-1000000	364090	VEEV replicon vector YFV-C3opt-NS2mut, complete sequence
# contig-2000000	565995	Bundibugyo ebolavirus, complete genome
# contig-3000000	565995	Bundibugyo ebolavirus, complete genome
# contig-3000000	186541	Cote d'Ivoire ebolavirus, complete genome
# contig-4000000	565995	Bundibugyo ebolavirus, complete genome
# contig-5000000	565995	Bundibugyo ebolavirus, complete genome
# contig-6000000	565995	Bundibugyo ebolavirus, complete genome
# contig-7000000	565995	Bundibugyo ebolavirus, complete genome
# contig-8000000	565995	Bundibugyo ebolavirus, complete genome
# contig-8000000	565995	Bundibugyo ebolavirus, complete genome
# contig-9000000	565995	Bundibugyo ebolavirus, complete genome
# contig-10000000	565995	Bundibugyo ebolavirus, complete genome
# contig-11000000	565995	Bundibugyo ebolavirus, complete genome
# contig-12000000	6279	Brugia malayi rRNA promoter binding protein partial mRNA
# contig-12000000	6279	Brugia malayi unspecific monooxygenase  partial mRNA
# contig-12000000	6279	Brugia malayi alpha-L1 nicotinic acetyl choline receptor, putative partial mRNA
# contig-12000000	6279	Brugia malayi alpha-L1 nicotinic acetyl choline receptor, putative partial mRNA
# contig-12000000	6279	Brugia malayi alpha-L1 nicotinic acetyl choline receptor, putative partial mRNA

# and produces output like this:

# 565995	contig-6000000,contig-4000000,contig-8000000,contig-3000000,contig-7000000,contig-2000000,contig-11000000,contig-10000000,contig-5000000,contig-9000000
# 10847	contig-0
# 11036	contig-1000000
# 364090	contig-1000000
# 364092	contig-1000000
# 186541	contig-3000000
# 6279	contig-12000000
# 364089	contig-1000000
# 364086	contig-1000000
# 10844	contig-0
# 364087	contig-1000000
# 364091	contig-1000000

# update: counts added

my %h=();

open(my $infile, "<", $ARGV[0]); 
open(my $outfile, ">", $ARGV[1]);

while (<$infile>)
{	
	my @a=split; 
	$h{$a[1]}{$a[0]}=1;			
}				 		 

foreach my $key1 ( keys %h ) 
{ 
	print $outfile $key1,"\t"; 

	my $counter=1; 
	
	foreach my $key2 ( keys %{$h{$key1}} ) 
	{ 
		if ($counter==1) {print $outfile $key2; $counter++;} else {print  $outfile ",",$key2; $counter++;} 
	} 
	print $outfile "\t$counter\n";
}

close($infile); 
close($outfile); 
