#!/usr/bin/perl

# a wrapper for BLAST: choose either blastn or blastx, choose algorithms, choose options, choose to make a report

# note: defaults:
# -num_descriptions=10
# -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"

# e.g., blast_wrapper.pl --type blastn --query input.fasta --db /my/db --task megablast --out output.blast --options '-evalue 1e-4 -word_size 28'

use Getopt::Long;
# use FindBin;
# use lib "$FindBin::Bin/../Local_Module";
use FindBin qw($RealBin);
use lib "$RealBin/../Local_Module";
use Verbose_Sys;

# default

$outfmt="\"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qseq\"";

GetOptions ('query=s' => \$query,		# inputfile
            'db=s' => \$db,				# db
            'task=s' => \$task,			# (options are: megablast dc-megablast blastn)
            'type=s' => \$type,			# (options are: blastn blastx) 
            'out=s' => \$out,			# outputfile
            'options=s' => \$options,	# (options are: any BLAST options)                                    
            'outfmt=s' => \$outfmt);	# (options are: BLAST outfmt options) 

# e.g., 
# blastn -query tmpsplit${suffix} -db $db -task $task -out blastout${suffix} -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" -num_threads=1 -num_descriptions=10 -evalue 1e-4 -word_size 28

# this is awkwardly coded. re-do ? 

if ($type eq "blastx" || !(defined($task)))
{
        $task_option="";
}
else
{
        $task_option="-task $task";
}

# if task defined, force type to be blastn
if ($task eq "megablast" || $task eq "dc-megablast" || $task eq "blastn")
{
        $type="blastn";
}

print "[start]\n";

my $cmd = "$type -query $query -db $db $task_option -out $out -outfmt $outfmt -num_descriptions=10 $options";				
verbose_system($cmd);

print "[end]\n";
