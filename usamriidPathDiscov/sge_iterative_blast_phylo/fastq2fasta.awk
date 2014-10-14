#!/bin/awk -f

{ 
	if (NR%4==1) 
	{
		print ">"substr($0,2,length);
	}

	if (NR%4==2) 
	{
		print;
	}
}

