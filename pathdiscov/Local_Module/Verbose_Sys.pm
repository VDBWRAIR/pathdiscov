use strict;
use File::Basename;

# run a system commmand "verbosely", echo-ing the command as well as giving the time
sub verbose_system  
{
	# When the subroutine is called any parameters are passed as a list in the special @_ (array form of $_)
	# @_ exists magically - anything passed in

	my $cmd = shift;  
    my @parts = split(' ', $cmd);
    my $exe = basename($parts[0]);

	my $start_time = time();		
    print_system($cmd);
 	my $end_time = time();
 	print "[$exe deltat] ",$end_time-$start_time,"\n";
}

# run a system commmand but print the command first
sub print_system  
{
	my $cmd = shift;  

	print "[cmd] ",$cmd,"\n";
	system($cmd);
}

# run a system commmand to qsub and return job id
# ASSUME qsub message is of the form, e.g., "Your job 3137692 ("my_job") has been submitted"
sub qsub_system  
{
	my $qcmd = shift;
	 
	print "[qcmd] ",$qcmd,"\n";
	
	my $qsub_message=`$qcmd`;
	
	print($qsub_message);	
	chomp($qsub_message); 
	
	my @h=split(/ /,$qsub_message);
	
	# return job ID
	return $h[2];
}

1;
