#!/usr/bin/perl

# a wrapper for BLAST: choose either blastn or blastx, choose algorithms, choose options, choose to make a report

# note: defaults:
# -num_descriptions=10
# -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"

# e.g., blast_wrapper.pl --type blastn --query input.fasta --db /my/db --task megablast --out output.blast --options '-evalue 1e-4 -word_size 28'

use Getopt::Long;


# default

$outfmt="\"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen\"";

GetOptions ('query=s' => \$query,		# inputfile
            'db=s' => \$db,				# db
            'task=s' => \$task,			# (options are: megablast dc-megablast blastn)
            'type=s' => \$type,			# (options are: blastn blastx) 
            'out=s' => \$out,			# outputfile
            'options=s' => \$options,	# (options are: any BLAST options)                                    
            'outfmt=s' => \$outfmt);	# (options are: BLAST outfmt options) 

# e.g., 
# blastn -query tmpsplit${suffix} -db $db -task $task -out blastout${suffix} -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" -num_threads=1 -num_descriptions=10 -evalue 1e-4 -word_size 28

if ($type eq "blastx" || !(defined($task)))
{
        $task_option="";
}
else
{
        $task_option="-task $task";
}

my $cmd = "$type -query $query -db $db $task_option -out $out -outfmt $outfmt -num_descriptions=10 $options";				
print "[cmd] ",$cmd,"\n";
system($cmd);
