from os.path import dirname,basename,abspath,exists
import os
from functools import partial
import argparse
import shlex

try:
    from collections import OrderedDict
except ImportError:
    from ordereddict import OrderedDict

import sh

# Staticly set options for blast
MAX_TARGET_SEQS = 10
BLAST_FORMAT = "\"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore\""

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--outputdir',
        help='path to put blast results'
    )
    parser.add_argument(
        '--logs',
        help='path to place log files'
    )
    parser.add_argument(
        '--inputfasta',
        help='fasta sequence file to blast as -query'
    )
    parser.add_argument(
        '--outfile',
        help='blast output file to write that will contain all results'
    )
    parser.add_argument(
        '--outheader',
        help='Deprecated'
    )
    parser.add_argument(
        '--ninst',
        type=int,
        help='Number of total cpus to use'
    )
    parser.add_argument(
        '--db',
        help='Blast db path'
    )
    parser.add_argument(
        '--blast_type',
        choices=('diamond','blastn','blastx'),
        help='Executable to call'
    )
    parser.add_argument(
        '--task',
        choices=('megablast','dc-megablast','blastn','blastx'),
        help='-task to use for blast/diamond'
    )
    parser.add_argument(
        '--blast_options',
        help='Options to pass on to blast'
    )
    return parser.parse_args()

def parallel_blast(inputfile, outfile, ninst, db, blasttype, task, blastoptions):
    '''
    Runs blast commands in parallel on a given fasta file

    :param file inputfile: Input fasta file handle
    :param file outfile: Output file handle
    :param int ninst: number of cpus to use if not in PBS or SGE job
    :param str db: Database path to blast against
    :param str blasttype: Blast exe to use
    :param str task: Blast task to run with -task option for blasttype
    :param str blastoptions: other options to pass to blast
    '''
    parallel = partial(
        sh.Command('parallel'), u=True, pipe=True, block='100k', recstart='>',
    )
    sshlogins = ' '.join(generate_sshlogins(ninst))
    parallelcmd = sshlogins
    parallelcmd += " \"{blastcmd}\""
    
    #my $cmd = "$type -query $query -db $db $task_option -out $out -outfmt $outfmt -max_target_seqs 10 $options";
    blastcmd = "$(which {blasttype}) -task {task} -db {db} " \
        "-max_target_seqs {max_target_seqs} -outfmt {blastfmt} {blastoptions} " \
        "-query - "
    blastcmd = blastcmd.format(
        blasttype=blasttype, task=task, db=db, max_target_seqs=MAX_TARGET_SEQS,
        blastfmt=BLAST_FORMAT, blastoptions=blastoptions
    )
    cmd = parallelcmd.format(
        blastcmd=blastcmd,
    )
    parallel(cmd, _in=inputfile, _out=outfile)

    '''
    cat input.fa | /usr/bin/time parallel -u --sshloginfile ${PBS_NODEFILE} --pipe --block 100k --recstart '>' -P ${PBS_NUM_PPN} "$(which blastn) -max_target_seqs 10 -db /media/VD_Research/databases/ncbi/blast/nt/nt -evalue 0.01 -outfmt  \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore\" -query -" > results.blast
    '''

def get_hostfile():
    '''
    Return the path to the hostfile|machinefile depending on SGE or PBS
    If local, returns ''
    '''
    return os.environ.get('PBS_NODEFILE', os.environ.get('PE_HOSTFILE', ''))

def parse_hostfile(hostfile_fh):
    '''
    Parse job queue system's hostfile for either PBS style or SGE style

    These are available with $PBS_NODEFILE and $PE_HOSTFILE during jobs that are run
    Always assumes first column is hostname
    If second column exists, then it is cpu count
    Any other columns are ignored

    :param file hostfile_fh: file like object that contains host file information
    '''
    hosts = OrderedDict()
    for line in hostfile_fh:
        x = line.split()
        hostname = x[0]
        if hostname not in hosts:
            hosts[hostname] = 0
        cpus = 1 # Default is to increment 1 for ever found(machinefile)
        if len(x) > 1: #Non PBS machinefile
            try:
                cpus = int(x[1])
            except ValueError as e:
                raise ValueError(
                    "Invalid second column value found in hostfile: {0}".format(x[1])
                )
        hosts[hostname] += cpus
    return hosts.items()

def generate_sshlogins(ninst=None):
    '''
    Generate a list of --sshlogins compatible with GNU Parallel such that they
    match PBS_NODEFILE, PE_HOSTFILE or just the local host

    :param int ninst: Number of cpus to use if not in SGE or PBS job
    :returns: ['--sshlogin cpu/host1', ..., '--sshlogin cpu/hostN']
    '''
    path = get_hostfile()
    if not path:
        if ninst is None:
            ninst = ''
        else:
            ninst = '{0}/'.format(ninst)
        return ['--sshlogin {0}:'.format(ninst)]
    else:
        with open(path) as fh:
            hosts = parse_hostfile(fh)
            return ['--sshlogin {1}/{0}'.format(*x) for x in hosts]

def main():
    args = parse_args()
    args.outputdir = abspath(args.outputdir)
    assert exists(args.outputdir), '[error] {0} does not exist'.format(args.outputdir)
    assert exists(args.inputfasta), '[error] {0} does not exist'.format(args.inputfasta)
    with open(args.inputfasta) as infile:
        with open(args.outfile, 'w') as outfile:
            parallel_blast(
                infile, outfile, args.ninst, args.db, args.blast_type, args.task,
                args.blast_options 
            )

if __name__ == '__main__':
    main()
