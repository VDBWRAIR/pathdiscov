#!/usr/bin/env perl


use strict;
use warnings;
use Getopt::Long;



my $field_file1=1;
my $field_file2=1;
my $lower;
my $delim;

GetOptions(
    'one=i' => \$field_file1,
    'two=i' => \$field_file2,
    'delim=s' => \$delim,
    'lower' => \$lower,
);

if(defined $delim){
    # my $compiled= eval{qr/$delim/};
    # die "$0: bad pattern ($delim):\n$@" unless $compiled;
}
else {
    $delim="\t";
}

# tableX.pl --one 3 --two 2 <FILE1> <FILE2>
# tableX.pl <FILE1> <FILE2>   # same as tableconcatlines
# tableX.pl --delim ',' <FILE1> <FILE2>   # same as tableconcatline.csv

open (my $fh, '<', $ARGV[1]);
my %h;
while (<$fh>) {
    chomp;
    my @F=split /$delim/;
    my $idx= $field_file2 - 1;
    my @excluded=exclude_from_array(\@F, $idx);
    if($lower){
	$h{lc($F[$idx])}=join($delim, @excluded);
    }
    else {
	$h{$F[$idx]}=join($delim, @excluded);
    }
}
close $fh;

open (my $fh2, '<', $ARGV[0]);
while (<$fh2>) {
    chomp;
    my @F=split /$delim/;
    my $idx= $field_file1 - 1;
    my $is_defined;
    if($lower){
	$is_defined= 1 if(defined $h{lc($F[$idx])});
    }
    else{
	$is_defined= 1 if(defined $h{$F[$idx]});
    }
    
    if($is_defined) {
	if($lower){
	    print $_ . $delim . $h{lc($F[$idx])} ;
	}
	else {
	    print $_ . $delim . $h{$F[$idx]} ;
	}
    }
    else {
	print $_ ;
    }
    print "\n";

}

close $fh2;

# exclude_from_array(\@arr, 2);
# returns the array with the item in the second index removed;
sub exclude_from_array {
    my $arr_ref = shift;
    my $idx = shift;
    my @arr;
    my $counter=0;
    for(my $i=0; $i < scalar(@{$arr_ref}) ; $i++) {
	if($i == $idx) {
	    next;
	}
	else {
	    $arr[$counter] = $arr_ref->[$i];
	}
	$counter++;
    }

    return @arr;
}


