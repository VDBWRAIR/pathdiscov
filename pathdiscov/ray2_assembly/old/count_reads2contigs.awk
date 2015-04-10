#!/bin/awk -f

# this counts of the number of pairs + singletons that map to a particular contig - i.e., a pair counts as 1 and a singleton counts as 1
# also counts the number of split reads that map to a contig (obviously, applicable for paired reads only)

# input: a sam file piped in from std:in

# example: cat bowtie2_mapping/out.sam | $path_scripts/count_reads2contigs.awk | sed \'s|\*|unmapped|\' | sort -k2,2nr > contig_numreads.txt

{
	if ($0 !~ /^@/) # not header
	{

		# count the split reads mapping to contig
		if ($1==prev1 && $3!=prev3 && NR>1) # id same but mapping destination different
		{
		        # increment contig
		        if (!y[$3]) {y[$3]=1} else {y[$3]++}
		        # increment prev contig
		        if (!y[prev3]) {y[prev3]=1} else {y[prev3]++}
		}

		# count pairs and singletons that map to contig
		if ($1==prev1 && $3==prev3 && NR>1) # id same and mapping destination same
		{
			# do nothing
			; 
		} 
		else 
		{
			prev1=$1; prev3=$3; 
			# if there exists no entry in the hash, create one; else increment count
			if (!x[$3]) {x[$3]=1} else {x[$3]++}
		}

		# print $0; print $1; print $3; print prev1; print prev3; printf "\n";

	}
}
END{ for (i in x){ printf i"\t"x[i]"\t"; if (y[i]) {print y[i]} else {print "0"} } }
