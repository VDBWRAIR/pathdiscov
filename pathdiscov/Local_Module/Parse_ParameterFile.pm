use strict;

# Parse parameter file and return params and settings in a hash
sub parse_param  
{
	# When the subroutine is called any parameters are passed as a list in the special @_ (array form of $_)
	# @_ exists magically - anything passed in

	my $infile = shift;  
	open(my $fh, '<', $infile);

	my %hoh = ();					# hash of hashes - to store {command}{parameters}
	
	my $current_command = "";
	
	# parse configuration file
	while (<$fh>)
	{
		# ignore comments and empty lines
		if ($_ !~ m/^#/ && $_ !~ m/^(\s*)$/) 
		{
			chomp $_;
					
			# split on hash to get uncommented section
			my @row0 = split(/#/, $_);
			 
			# line is the uncommented portion of the row
			my $line=$row0[0];

			my $param="";			
			my $setting="no";
	
			# nonspace space anything
			if ($line =~ m/(\S+)(\s+)(.*)/)
			{
				$param=$1;	
				$setting=$3;
				$setting=~s/\s*$//;		# removing trailing whitespace
			}
			else
			{
				# if can't parse, skip by setting setting to "_"
				$setting="-";
			}
											
			# -------------------- for samples --------------------
			# if setting is a "no", dont make a hash key for it
			if (!($setting eq "-" || $setting eq "no" || $setting eq "0"))
			{
				if ($param eq "command") 
				{			
					# name of the current command 
					$current_command = $setting;								
				}
				else
				{
					$hoh{$current_command}{$param} = $setting;	
				}					
			} # not ignore setting
		} # not comment
	}

	close($fh);

	return \%hoh;
}

1;

