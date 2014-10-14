#!/usr/bin/perl

# perform a weighted count of second column. The first column acts as the weighting factor. So, for example, if there is a streak of 6 1's in the first column, 
# every number in the second column gets a weight of 1/6th. At the end, all of the numbers in the second column are counted (every unique # appears once together 
# with its count). 

# In practice this is used to provided a count of taxid's as weighted by the queryids

# ASSUME: col1 is qid, col2 is taxid

# e.g., 

# $ cat junk.txt
# 1	10
# 1	10
# 1	10
# 1	20
# 1	20
# 1	10
# 3	5
# 3	5
# 4	8
# 4	9
# 4	9
# 4	10
# 4	4
# 6	40
# 60	405

# useage: 
# cat tmp.annotate.blast | weighted_count.pl > tmp.count	

use Getopt::Long;

GetOptions ('nosum' => \$nosum);		# dont sum taxid's from different qid's	


$counter=1;								# counter for steak of query IDs


%total_sum = (); 						# this hash is the total sum for every taxid


# loop through input file (taxids)
while (<STDIN>)
{
	@a=split; 							# split line on whitespace
	if ($.==1)							# if NR==1 
	{
		$qid=$a[0];						# get qid 
		
		# b will be a hash table that's kept per streak of qid's; after, it's cleared
		
		$b{$a[1]}=1;					# get first taxid
	} 
#	elsif ($a[0]==$qid) 				# if steak of qid's
	elsif ($a[0] eq $qid) 				# if steak of qid's - use strings instead of numbers - that way, can accommodate non numeral IDs
	{
		$qid=$a[0]; 					# set qid to first column (not really nec since already equal)
		
		if (!defined($b{$a[1]})) 		# if entry not defined in hash, add one
		{
			$b{$a[1]}=1;			
		} 
		else 					
		{
			$b{$a[1]}++;				# if defined, increment
		} 
		
		$counter++;						# increment qid streak counter
	} 
	else 								# if change in qid
	{
		# print "***",$qid,"\t",$counter,"\n"; 

		foreach my $key (sort(keys %b)) 
		{
#			print $key,"\t",$b{$key},"\n";
			print $key,"\t",$b{$key}/$counter,"\n" if $nosum;

			if (!defined($total_sum{$key})) 					# if entry not defined in hash, add one
			{
				$total_sum{$key}=$b{$key}/$counter;			
			} 
			else 					
			{
				$total_sum{$key}+=$b{$key}/$counter;			# if defined, add on to it
			}
		}

		$qid=$a[0];						# set new qid 
		
		$counter=1; 

		%b = (); 						# clear hash
		$b{$a[1]}=1;					# add first entry in hash
	}
}

# do the last one:

# print "***",$qid,"\t",$counter,"\n"; 

foreach my $key (sort(keys %b)) 
{
#		print $key,"\t",$b{$key},"\n";
		print $key,"\t",$b{$key}/$counter,"\n" if $nosum;

		if (!defined($total_sum{$key})) 					# if entry not defined in hash, add one
		{
			$total_sum{$key}=$b{$key}/$counter;			
		} 
		else 					
		{
			$total_sum{$key}+=$b{$key}/$counter;			# if defined, add on to it
		}
}

# now print results - the total sum

if (!$nosum)
{
	foreach my $key (sort(keys %total_sum)) 
	{
#		print $key,"\t",$total_sum{$key},"\n";
		printf("%d\t%.2f\n", $key, $total_sum{$key});
	}
}
