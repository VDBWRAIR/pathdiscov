use strict;

# pass in samtools flag: return true if the read is the first read in a pair, else false
sub is_firstread  
{
	# When the subroutine is called any parameters are passed as a list in the special @_ (array form of $_)
	# @_ exists magically - anything passed in

	# first argument is the sam file flag
	my $flag = shift; 
	
	# turn flag into binary 
	my $binarynum=sprintf ("%b", $flag);	
	
	# reverse it
	my $rbinarynum=reverse($binarynum);

	# padding, so string will be const length
	my $padding="";

	for (my $count = 1; $count <= 11-length($rbinarynum); $count++) 
	{
 		$padding=$padding."0";
 	}
 	
 	# pad
	$rbinarynum=$rbinarynum.$padding;
	
	return substr($rbinarynum, 6, 1);

	# bit fields: 
		
#	print substr($rbinarynum, 0, 1),"\t","the read is paired in sequencing","\n";
#	print substr($rbinarynum, 1, 1),"\t","the read is mapped in a proper pair","\n";
#	print substr($rbinarynum, 2, 1),"\t","the query sequence itself is unmapped","\n";
#	print substr($rbinarynum, 3, 1),"\t","the mate is unmapped","\n";
#	print substr($rbinarynum, 4, 1),"\t","strand of the query (1 for reverse)","\n";
#	print substr($rbinarynum, 5, 1),"\t","strand of the mate","\n";
#	print substr($rbinarynum, 6, 1),"\t","the read is the first read in a pair","\n";
#	print substr($rbinarynum, 7, 1),"\t","the read is the second read in a pair","\n";
#	print substr($rbinarynum, 8, 1),"\t","the alignment is not primary","\n";
#	print substr($rbinarynum, 9, 1),"\t","QC failure","\n";
#	print substr($rbinarynum, 10, 1),"\t","optical or PCR duplicate","\n";

}

# pass in samtools flag: return true if the read is the second read in a pair, else false
sub is_secondread  
{
	# When the subroutine is called any parameters are passed as a list in the special @_ (array form of $_)
	# @_ exists magically - anything passed in

	# first argument is the sam file flag
	my $flag = shift; 
	
	# turn flag into binary 
	my $binarynum=sprintf ("%b", $flag);	
	
	# reverse it
	my $rbinarynum=reverse($binarynum);

	# padding, so string will be const length
	my $padding="";

	for (my $count = 1; $count <= 11-length($rbinarynum); $count++) 
	{
 		$padding=$padding."0";
 	}
 	
 	# pad
	$rbinarynum=$rbinarynum.$padding;
	
	return substr($rbinarynum, 7, 1);
}

# pass in samtools flag: return true if the read is unmapped, else false
sub is_unmapped  
{
	# When the subroutine is called any parameters are passed as a list in the special @_ (array form of $_)
	# @_ exists magically - anything passed in

	# first argument is the sam file flag
	my $flag = shift; 
	
	# turn flag into binary 
	my $binarynum=sprintf ("%b", $flag);	
	
	# reverse it
	my $rbinarynum=reverse($binarynum);

	# padding, so string will be const length
	my $padding="";

	for (my $count = 1; $count <= 11-length($rbinarynum); $count++) 
	{
 		$padding=$padding."0";
 	}
 	
 	# pad
	$rbinarynum=$rbinarynum.$padding;
	
	return substr($rbinarynum, 2, 1);
}

1;
