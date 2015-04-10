use strict;

# return hash of each line, as a key, hashed to one (line => "1")
sub line_hash  
{	
	# When the subroutine is called any parameters are passed as a list in the special @_ (array form of $_)
	# @_ exists magically - anything passed in

	# arg 1 - file name
	my $infile = shift; 	
	
	my %h = ();				# hash	

	if ( -s $infile )		# if file nonzero
	{	 
		open(my $fh, '<', $infile);
				
		# loop thro file
		while (<$fh>)
		{
			chomp $_;
				
			# line is the whole line
			my $line = $_;
			
			# hash line to "1"
			$h{$line} = 1;
		}
	
		close($fh);
	}

	return \%h;
}

# return hash of a specified column hashed to one (column => "1")
sub column_hash  
{		
	# When the subroutine is called any parameters are passed as a list in the special @_ (array form of $_)
	# @_ exists magically - anything passed in

	# arg 1 - file name
	# arg 2 - column number

	my $infile = shift; 
	
	my %h = ();				# hash	

	if ( -s $infile )		# if file nonzero
	{	
		open(my $fh, '<', $infile);
	
		my $colnum = shift; 	#  column number
		$colnum--;				# shift from 1-based to 0-based counting (recall: in perl, array index starts from 0)
		
		my $skipline = 0; 		# skip first line (header) ... this is false now, so the first line wont be skipped
			
		# loop thro file
		while (<$fh>)
		{
			if (!$skipline)
			{
				chomp $_;
					
				# line is the whole line
				my $line = $_;
				
				# split on ws
				my @row = split;
				
				# hash key is specified column (recall: array index starts from 0)
				# $h{$row[$colnum]} = $line;
				$h{$row[$colnum]} = 1;
			}	
			$skipline = 0;
		}
	
		close($fh);
	}

	return \%h;
}

# return hash of column1 to column2
sub col1_to_col2_hash  
{		
	# When the subroutine is called any parameters are passed as a list in the special @_ (array form of $_)
	# @_ exists magically - anything passed in

	# arg 1 - file name

	my $infile = shift; 
	
	my %h = ();				# hash	

	if ( -s $infile )		# if file nonzero
	{	
		open(my $fh, '<', $infile);
		
		my $skipline = 0; 		# skip first line (header) ... this is false now, so the first line wont be skipped
			
		# loop thro file
		while (<$fh>)
		{
			if (!$skipline)
			{
				chomp $_;
					
				# line is the whole line
				my $line = $_;
				
				# split on ws
				my @row = split;
				
				$h{$row[0]} = $row[1];
			}	
			$skipline = 0;
		}
		
		close($fh);
	}

	return \%h;
}

# return hash of the IDs of a fastq file to "1")
sub fastqid_hash  
{
	# When the subroutine is called any parameters are passed as a list in the special @_ (array form of $_)
	# @_ exists magically - anything passed in

	# arg 1 - file name

	my $infile = shift;
	
	my %h = ();				# hash	

	if ( -s $infile )		# if file nonzero
	{	
		open(my $fh, '<', $infile);
	
		my $nr = 1;		# row number
				
		# loop thro file
		while (<$fh>)
		{
			# ID is every fourth line
			if ( $nr%4 == 1 )
			{
				chomp $_;
					
				# don't print leading "@"
				my $key = substr($_,1);
						
				# hash key is fastq id
				$h{$key} = 1;
			}
			
			$nr++;
		}
	
		close($fh);
	}

	return \%h;
}

# return hash of the IDs of a fasta file to "1")
sub fastaid_hash  
{
	# When the subroutine is called any parameters are passed as a list in the special @_ (array form of $_)
	# @_ exists magically - anything passed in

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
					
				# don't print leading ">"
				my $key = substr($_,1);
						
				# hash key is fastq id
				$h{$key} = 1;
			}
		}
			
		close($fh);
	}

	return \%h;
}

# sum up over a column, specified as input
sub column_sum  
{		
	# When the subroutine is called any parameters are passed as a list in the special @_ (array form of $_)
	# @_ exists magically - anything passed in

	# arg 1 - file name
	# arg 2 - column number

	my $infile = shift; 
	
	my %h = ();				# hash	

	my $sum = 0;

	if ( -s $infile )		# if file nonzero
	{	
		open(my $fh, '<', $infile);
	
		my $colnum = shift; 	#  column number
		$colnum--;				# shift from 1-based to 0-based counting (recall: in perl, array index starts from 0)
		
		my $skipline = 0; 		# skip first line (header) ... this is false now, so the first line wont be skipped
			
		# loop thro file
		while (<$fh>)
		{
			if (!$skipline)
			{
				chomp $_;
					
				# line is the whole line
				my $line = $_;
				
				# split on ws
				my @row = split;
				
				$sum=$sum+$row[$colnum];
			}	
			$skipline = 0;
		}
	
		close($fh);
	}

	return $sum;
}

# count the IDs in a fasta file
sub fastaid_count  
{
	# When the subroutine is called any parameters are passed as a list in the special @_ (array form of $_)
	# @_ exists magically - anything passed in

	# arg 1 - file name

	my $infile = shift;
	
	my $sum = 0;

	if ( -s $infile )		# if file nonzero
	{	
		open(my $fh, '<', $infile);
	
		while (<$fh>)
		{
			# ID is every line starting with ">"
			if ( $_ =~ m/^>/ )
			{
				$sum++;
			}
		}
			
		close($fh);
	}

	return $sum;
}

1;

