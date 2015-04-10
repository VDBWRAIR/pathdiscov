#!/usr/bin/perl

# system("echo [Directory] `pwd`");
# system("echo [Machine] `uname -n`");
# system("echo [Start] `date`");
# print "[args] ".join(" ", @ARGV)."\n";

if ($ARGV[0] eq "change_fastq_id")
{
	&change_fastq_id ( \$ARGV[1], \$ARGV[2], \$ARGV[3] );	
}
elsif ($ARGV[0] eq "get_common_uneven_files")
{
	# input_R1 input_R2 R1.single.fastq R1.paired.fastq R2.single.fastq R2.paired.fastq
	&get_common_uneven_files ( \$ARGV[1], \$ARGV[2], \$ARGV[3], \$ARGV[4], \$ARGV[5], \$ARGV[6] );
}
elsif ($ARGV[0] eq "get_common_uneven_fastas")
{
	# input_R1 input_R2 R1.single.fasta R1.paired.fasta R2.single.fasta R2.paired.fasta
	&get_common_uneven_fastas ( \$ARGV[1], \$ARGV[2], \$ARGV[3], \$ARGV[4], \$ARGV[5], \$ARGV[6] );
}
elsif ($ARGV[0] eq "interleave_uneven_files")
{
	&interleave_uneven_files ( \$ARGV[1], \$ARGV[2], \$ARGV[3], \$ARGV[4], \$ARGV[5] );	
}
elsif ($ARGV[0] eq "interleave_uneven_files2")
{
	&interleave_uneven_files2 ( \$ARGV[1], \$ARGV[2], \$ARGV[3], \$ARGV[4], \$ARGV[5], \$ARGV[6], \$ARGV[7] );	
}
elsif ($ARGV[0] eq "interleave_even_files")
{	
	&interleave_even_files ( \$ARGV[1], \$ARGV[2], \$ARGV[3] );	
}
	
# system("echo [End] `date`");
	
# -------------------- change fastq IDs --------------------			
# change fastq ids to numbers 
# added later: make sure the third line is only a single char ("+")
# added later: replace any "." s in the sequence with "N" s
sub change_fastq_id 
{

	(my $ref1, my $ref2, my $ref3) = @_;

	my $input = ($$ref1);
	my $output = ($$ref2);
	my $output2 = ($$ref3);

#	print $input,"\n";
	open my $infile1, '<', $input;
	open my $outfile1, '>', $output;
	open my $outfile2, '>', $output2;
	
	my $counter=1;	# count rows

	while (<$infile1>)
	{	
		if ($counter%4==1) # ID line
		{
			print $outfile1 ("@".($counter+3)/4,"\n");
			# don't print leading "@"
			print $outfile2 (($counter+3)/4,"\t",substr($_,1));
		}
		elsif ($counter%4==2) # seq line
		{
			$_ =~ s/\./N/g; # if dots, replace with "N"s
			print $outfile1 $_;
		}		 
		elsif ($counter%4==3) # third line ("+" or "+ID")
		{
			print $outfile1 substr($_,0,1),"\n";
		} 		
		else # qual 
		{
			print $outfile1 $_;
		} 	
		$counter++;
	}
	
	close $infile1;
	close $outfile1;
	close $outfile2;	
}

# -------------------- find common IDs --------------------			
# find common IDs among uneven files in chunks of 4 (nec if they've been cleaned with fastx, thus removing some reads in an asymmetrical manner). produces 4 files
sub get_common_uneven_files
{

	(my $ref1, my $ref2, my $ref3, my $ref4, my $ref5, my $ref6) = @_;
	
	my $input1 = ($$ref1);
	my $input2 = ($$ref2);
	my $output1 = ($$ref3);		# 1.clean_single
	my $output2 = ($$ref4);		# 1.clean_paired
	my $output3 = ($$ref5);		# 2.clean_single
	my $output4 = ($$ref6);		# 2.clean_paired

	open my $infile1, '<', $input1;
	open my $infile2, '<', $input2;
	
	open my $outfile1, '>', $output1;	# 1.clean_single
	open my $outfile2, '>', $output2;	# 1.clean_paired
	open my $outfile3, '>', $output3;	# 2.clean_single
	open my $outfile4, '>', $output4;	# 2.clean_paired
			
	# open the second file and create a hash based on IDs

	my %myh = ();
	
	while (<$infile2>)
	{
		chomp($_);
		
		# this ID will be the key to your hash
		my $mykey=$_;
	
		$chunk=$_."\n";
		$_ = <$infile2>;
		$chunk=$chunk.$_;
		$_ = <$infile2>;
	    $chunk=$chunk.$_;
	    $_ = <$infile2>;
	    $chunk=$chunk.$_;
	
		# map the ID to a four line chunk
		$myh{$mykey}=$chunk;
	}
	
	# loop through the first file and, if ID matches hash ID, print 4 lines as well as chunk from the second file
	while (<$infile1>) 
	{
		chomp($_);
		
		my $mykey=$_;
		
		if ($myh{$mykey})
		{
			# print 4 lines of 1st file to 1.clean_paired	
			print $outfile2 $_,"\n";
			$_ = <$infile1>;	
			print $outfile2 $_; 
	        $_ = <$infile1>;
	        print $outfile2 $_;
	        $_ = <$infile1>;
	        print $outfile2 $_;
	        
			# print 4 lines of 2nd file to 2.clean_paired	
			print $outfile4 $myh{$mykey};
			# then delete the val from the hash table
			delete($myh{$mykey});        
		} 
		else
		{
			# print 4 lines	to 1.clean_single	
			print $outfile1 $_,"\n";
			$_ = <$infile1>;	
			print $outfile1 $_; 
	        $_ = <$infile1>;
	        print $outfile1 $_;
	        $_ = <$infile1>;
	        print $outfile1 $_;
		}
	}
	
	# print the remaining vals in 1st file to 2.clean_single
	my @myvals = values %myh;
	print $outfile3 join("", @myvals);
		
	close $infile1;
	close $infile2;
	close $outfile1;
	close $outfile2;
	close $outfile3;
	close $outfile4;
}

# -------------------- find common IDs --------------------			
# same as above but for fastas rather than fastqs. produces 4 files
sub get_common_uneven_fastas
{

	(my $ref1, my $ref2, my $ref3, my $ref4, my $ref5, my $ref6) = @_;
	
	my $input1 = ($$ref1);
	my $input2 = ($$ref2);
	my $output1 = ($$ref3);		# 1.clean_single
	my $output2 = ($$ref4);		# 1.clean_paired
	my $output3 = ($$ref5);		# 2.clean_single
	my $output4 = ($$ref6);		# 2.clean_paired

	open my $infile1, '<', $input1;
	open my $infile2, '<', $input2;
	
	open my $outfile1, '>', $output1;	# 1.clean_single
	open my $outfile2, '>', $output2;	# 1.clean_paired
	open my $outfile3, '>', $output3;	# 2.clean_single
	open my $outfile4, '>', $output4;	# 2.clean_paired
			
	# open the second file and create a hash based on IDs

	my %myh = ();
	
	while (<$infile2>)
	{
		chomp($_);
		
		# this ID will be the key to your hash
		my $mykey=$_;
	
		$chunk=$_."\n";
		$_ = <$infile2>;
		$chunk=$chunk.$_;
	
		# map the ID to a 2 line chunk
		$myh{$mykey}=$chunk;
	}
	
	# loop through the first file and, if ID matches hash ID, print 4 lines as well as chunk from the second file
	while (<$infile1>) 
	{
		chomp($_);
		
		my $mykey=$_;
		
		if ($myh{$mykey})
		{
			# print 4 lines of 1st file to 1.clean_paired	
			print $outfile2 $_,"\n";
			$_ = <$infile1>;	
			print $outfile2 $_; 
	        
			# print 2 lines of 2nd file to 2.clean_paired	
			print $outfile4 $myh{$mykey};
			# then delete the val from the hash table
			delete($myh{$mykey});        
		} 
		else
		{
			# print 4 lines	to 1.clean_single	
			print $outfile1 $_,"\n";
			$_ = <$infile1>;	
			print $outfile1 $_;
		}
	}
	
	# print the remaining vals in 1st file to 2.clean_single
	my @myvals = values %myh;
	print $outfile3 join("", @myvals);
		
	close $infile1;
	close $infile2;
	close $outfile1;
	close $outfile2;
	close $outfile3;
	close $outfile4;
}

# -------------------- interleave uneven files --------------------			
# interleave uneven files in chunks of 4 (nec if they've been cleaned with fastx, thus removing some reads in an asymmetrical manner). produces 3 files
# Velvet requires interleaving and it works for ABySS, too. 
# format fastq files for VELVET: For paired-end reads, the assumption is that each read is next to its mate read. In other words if the reads are indexed from 0, then reads 0 & 1 are paired, 2 & 3, 4 $ 5, etc.# format fastq files for ABySS: paired reads must have same IDs (interleaving ok)
sub interleave_uneven_files
{

	(my $ref1, my $ref2, my $ref3, my $ref4, my $ref5) = @_;
	
	my $input1 = ($$ref1);
	my $input2 = ($$ref2);
	my $output1 = ($$ref3);		# 1.clean_single
	my $output2 = ($$ref4);		# 2.clean_single
	my $output3 = ($$ref5);		# merge

	open my $infile1, '<', $input1;
	open my $infile2, '<', $input2;
	
	open my $outfile1, '>', $output1;	# 1.clean_single
	open my $outfile2, '>', $output2;	# 2.clean_single
	open my $outfile3, '>', $output3;	# merge
			
	# open the second file and create a hash based on IDs

	my %myh = ();
	
	while (<$infile2>)
	{
		chomp($_);
		
		# this ID will be the key to your hash
		my $mykey=$_;
	
		$chunk=$_."\n";
		$_ = <$infile2>;
		$chunk=$chunk.$_;
		$_ = <$infile2>;
	    $chunk=$chunk.$_;
	    $_ = <$infile2>;
	    $chunk=$chunk.$_;
	
		# map the ID to a four line chunk
		$myh{$mykey}=$chunk;
	}
	
	# loop through the first file and, if ID matches hash ID, print 4 lines as well as chunk from the second file
	while (<$infile1>) 
	{
		chomp($_);
		
		my $mykey=$_;
		
		if ($myh{$mykey})
		{
			# print 4 lines of 1st file to merge	
			print $outfile3 $_,"\n";
			$_ = <$infile1>;	
			print $outfile3 $_; 
	        $_ = <$infile1>;
	        print $outfile3 $_;
	        $_ = <$infile1>;
	        print $outfile3 $_;
	        
			# print 4 lines of 2nd file to merge	
			print $outfile3 $myh{$mykey};
			# then delete the val from the hash table
			delete($myh{$mykey});        
		} 
		else
		{
			# print 4 lines	to 1.clean_single	
			print $outfile1 $_,"\n";
			$_ = <$infile1>;	
			print $outfile1 $_; 
	        $_ = <$infile1>;
	        print $outfile1 $_;
	        $_ = <$infile1>;
	        print $outfile1 $_;
		}
	}
	
	# print the remaining vals in 1st file to 2.clean_single
	my @myvals = values %myh;
	print $outfile2 join("", @myvals);
		
	close $infile1;
	close $infile2;
	close $outfile1;
	close $outfile2;
	close $outfile3;

}

# -------------------- interleave uneven files2 --------------------			
# version 2 produces an interleaved, as well as a non-interleaved, output
sub interleave_uneven_files2
{
	
	(my $ref1, my $ref2, my $ref3, my $ref4, my $ref5, my $ref6, my $ref7) = @_;
	
	my $input1 = ($$ref1);
	my $input2 = ($$ref2);
	my $output1 = ($$ref3);		# 1.clean_single
	my $output2 = ($$ref4);		# 2.clean_single
	my $output3 = ($$ref5);		# merge
	my $output4 = ($$ref6);		# 1.clean_paired
	my $output5 = ($$ref7);		# 2.clean_paired
		
	open my $infile1, '<', $input1;
	open my $infile2, '<', $input2;
	
	open my $outfile1, '>', $output1;	# 1.clean_single
	open my $outfile2, '>', $output2;	# 2.clean_single
	open my $outfile3, '>', $output3;	# merge
	open my $outfile4, '>', $output4;	# 1.clean_paired
	open my $outfile5, '>', $output5;	# 2.clean_paired	
	
			
	# open the second file and create a hash based on IDs

	my %myh = ();
	
	while (<$infile2>)
	{
		chomp($_);
		
		# this ID will be the key to your hash
		my $mykey=$_;
	
		$chunk=$_."\n";
		$_ = <$infile2>;
		$chunk=$chunk.$_;
		$_ = <$infile2>;
	    $chunk=$chunk.$_;
	    $_ = <$infile2>;
	    $chunk=$chunk.$_;
	
		# map the ID to a four line chunk
		$myh{$mykey}=$chunk;
	}
	
	# loop through the first file and, if ID matches hash ID, print 4 lines as well as chunk from the second file
	while (<$infile1>) 
	{
		chomp($_);
		
		my $mykey=$_;
		
		if ($myh{$mykey})
		{
			# print 4 lines of 1st file to merge	
			print $outfile3 $_,"\n";	
			print $outfile4 $_,"\n";			
			$_ = <$infile1>;	
			print $outfile3 $_;
			print $outfile4 $_; 			
	        $_ = <$infile1>;
	        print $outfile3 $_;
			print $outfile4 $_; 
	        $_ = <$infile1>;
	        print $outfile3 $_;
			print $outfile4 $_; 
				        
			# print 4 lines of 2nd file to merge	
			print $outfile3 $myh{$mykey};
			print $outfile5 $myh{$mykey};			
			# then delete the val from the hash table
			delete($myh{$mykey});        
		} 
		else
		{
			# print 4 lines	to 1.clean_single	
			print $outfile1 $_,"\n";
			$_ = <$infile1>;	
			print $outfile1 $_; 
	        $_ = <$infile1>;
	        print $outfile1 $_;
	        $_ = <$infile1>;
	        print $outfile1 $_;
		}
	}
	
	# print the remaining vals in 1st file to 2.clean_single
	my @myvals = values %myh;
	print $outfile2 join("", @myvals);
		
	close $infile1;
	close $infile2;
	close $outfile1;
	close $outfile2;
	close $outfile3;
	close $outfile4;
	close $outfile5;	
}

# -------------------- interleave even files --------------------			
# merge two paired fastqs
sub interleave_even_files
{

	(my $ref1, my $ref2, my $ref3) = @_;
	
	my $input1 = ($$ref1);
	my $input2 = ($$ref2);
	my $output1 = ($$ref3);

	open my $infile1, '<', $input1;
	open my $infile2, '<', $input2;
	
	open my $outfile1, '>', $output1;				


	while(<$infile1>) 
	{
		# only want first field - quake adds stupid words to the fastq id
		chomp($_); 
		my @myrow = split(/\s+/, $_);	
		print $outfile1 $myrow[0],"\n";
		$_ = <$infile1>;
		print $outfile1 $_; 
		$_ = <$infile1>;
		print $outfile1 $_; 
		$_ = <$infile1>;
		print $outfile1 $_; 
	
		$_ = <$infile2>;
		chomp($_); 
		my @myrow2 = split(/\s+/, $_);	
		print $outfile1 $myrow2[0],"\n";
		$_ = <$infile2>;
		print $outfile1 $_;
		$_ = <$infile2>;
		print $outfile1 $_;
		$_ = <$infile2>;
		print $outfile1 $_;
	}
		
	close $infile1;
	close $infile2;
	close $outfile1;
}
