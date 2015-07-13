#!/usr/bin/env python

'''
You supply a pathdiscov project directory and a value to search(optionally also
can specify a different blast column) and the following occurs:

Contig reads
------------

projdir/results/iterative_blast_phylo_1/reports/*smallreport* is opened and all
rows where the matching supplied ``--blastcol`` matches the ``blastvalue`` are
used. Since each row contains contig name for result, then the resulting
fastq files inside of projdir/results/ray2_assembly_1/reads_by_contig can be
concatenated together for matching contig names and placed in requested output
fastq file.

Unassembled Reads
-----------------

Iterate over both R1 and R2 smallreport files inside of 
projdir/results/iterative_blast_phylo_2/reports and all
rows where the matching supplied ``--blastcol`` matches the ``blastvalue`` are
used.
The resulting matched columns will then be used to match reads inside of the input
R1 and R2 fastq files. Reads are extracted and concatenated into the requested
output file.


The end result is a single fastq file that contains all input reads that compose
all matching results for the supplied ``blastvalue``
'''

import sh
import sys
import itertools
from Bio import SeqIO
from os.path import *
import os
import argparse
from glob import glob

from pathdiscov import make_summary

def get_fastq_records_for_unassembled(projpath, blastcol, blastval, writefqidslist):
    '''
    Return a list of SeqRecord generators for all records that match from each
    input fastq file for unassembled reads
    Iterates over projpath/results/iterative_blast_phylo_2/reports/*.*smallreport.txt
    (which should be 1 or 2 files) and then uses resulting matching read names.
    These read names are then used to produce generators that only return the
    records from projpath/step1/R1.fastq and/or projpath/step1/R2.fastq that match
    
    Essentially only returning seqrecs for unassembled reads that match requested
    blastcol/blastval
    '''
    _blastfiles = join(
        projpath, 'results/iterative_blast_phylo_2/reports/*smallreport.txt'
    )
    _inputfiles = join(
        projpath, 'results/step1/*.fastq'
    )
    reportpaths = glob(_blastfiles)
    if len(reportpaths) == 0:
        raise ValueError('Missing {0} files'.format(_blastfiles))
    inputfiles = glob(_inputfiles)
    if len(inputfiles) == 0:
        raise ValueError('Missing {0} files'.format(_inputfiles))

    fqrecs = []
    with open('fastqids.lst', 'w') as fh:
        # Iterate over the input and report together
        for inputfile, smallreport in zip(inputfiles, reportpaths):
            fqids = get_qseqids_from_blastreport(smallreport, blastcol, blastval)
            fqg = filter_fastq_by_ids(inputfile, fqids)
            fqrecs.append(fqg)
            for _id in fqids:
                fh.write(_id + '\n')
    return fqrecs

def filter_fastq_by_ids(seqrecs, idlist):
    '''
    Return a generator of seqrecs that only returns seqrec ids from idlist
    '''
    for seqrec in SeqIO.parse(seqrecs, 'fastq'):
        if seqrec.id in idlist:
            yield seqrec

def get_fastq_records_for_contigs(projpath, blastcol, blastval, writecontigfilelist):
    '''
    Return list of SeqIO.parse for all contig fastq files

    for matching blast records
    Looks in projpath/results/iterative_blast_phylo_1/reports/contig.*smallreport.txt
    for blast results
    Looks in projpath/results/ray2_assembly_1/reads_by_contig
    for matching contig fastq files in form of <contigname>.group.fq
    '''
    _blastfile = join(
        projpath, 'results/iterative_blast_phylo_1/reports/contig.*smallreport.txt'
    )
    blastfile = glob(_blastfile)
    contigspath = join(projpath, 'results/ray2_assembly_1/reads_by_contig')
    # Make sure blast file exists
    if len(blastfile) != 1:
        raise ValueError('{0} is missing from project {1}'.format(
            _blastfile, projpath
        ))
    blastfile = blastfile[0]
    # Only want unique names so convert to set
    contigids = set(get_qseqids_from_blastreport(blastfile, blastcol, blastval))
    contig_paths = get_contig_fastq_file_list(contigspath, contigids)
    if writecontigfilelist:
        with open('contigfiles.lst', 'w') as fh:
            for f in contig_paths:
                fh.write(f + '\n')
    return contig_paths

def get_contig_fastq_file_list(basepath, names):
    '''
    Get list of paths for joining basepath and names.
    Report error on stdout for missing files and removes file from list
    Assumes contig files are named after being created by pathdiscov.parse_contigs
    '''
    # Get list of contif
    contig_fastqs = []
    for name in names:
        path = join(basepath, name + '.group.fq')
        if not exists(path):
            sys.stderr.write('{0} is missing\n'.format(path))
        else:
            contig_fastqs.append(path)
    return contig_fastqs

def get_seqrecs_from_files(filelist):
    '''
    Return list of parsed files
    '''
    parsed_fq = []
    for f in filelist:
        parsed_fq.append(SeqIO.parse(f, 'fastq'))
    return parsed_fq

def write_all_seqrecs(parsed_files, outputfile, overwrite):
    '''
    Chain all parsed files together and use biopython to write them
    into output file
    '''
    if exists(outputfile):
        if overwrite:
            os.unlink(outputfile)
        else:
            raise ValueError('{0} already exists'.format(outputfile))
    numwrote = 0
    with open(outputfile, 'a') as fh:
        for f in parsed_files:
            numwrote += SeqIO.write(f, fh, 'fastq')
    return numwrote

def get_qseqids_from_blastreport(blastreport, blastcol, blastval):
    '''
    Return all qseqids that match blastcol/blastval
    '''
    qseqids = []
    # returns generator of matching blast rows
    for row in make_summary.blast_results_for_(blastreport, blastcol, blastval):
        qseqids.append(row['qseqid'])
    return qseqids

def parse_args():
    parser = argparse.ArgumentParser(
        description='Get fastq of all reads that compose a requested blast result' \
    )

    parser.add_argument(
        'projectpath',
        help='Path to pathdiscov project'
    )

    parser.add_argument(
        'blastvalue',
        help='Blast value to search in --blastcol for to get reads'
    )

    parser.add_argument(
        '--blastcol',
        default='genus',
        help='Blast column to gather results from[Default: %(default)s]'
    )

    parser.add_argument(
        '--outputfile',
        default='blastresult.fastq',
        help='Output fastq file that will contain all results[Default: %(default)s]'
    )

    parser.add_argument(
        '--overwrite',
        default=False,
        action='store_true',
        help='Overwrite outputfile if exists[Default: %(default)s]'
    )

    parser.add_argument(
        '--write-used',
        default=False,
        action='store_true',
        help='Writes contigfiles.lst which contains contig files used and ' \
            'fastqids.lst which contains read ids ' \
            'that were used[Default: %(default)s]'
    )

    return parser.parse_args()

def main():
    args = parse_args()
    contig_files = get_fastq_records_for_contigs(
        args.projectpath, args.blastcol, args.blastvalue, args.write_used
    )
    contig_seqs = get_seqrecs_from_files(contig_files)
    read_seqs = get_fastq_records_for_unassembled(
        args.projectpath, args.blastcol, args.blastvalue, args.write_used
    )

    allseqs = contig_seqs + read_seqs
    numwrote = write_all_seqrecs(allseqs, args.outputfile, args.overwrite)
    sys.stdout.write("Wrote {0} records to {1}\n".format(numwrote, args.outputfile))

if __name__ == '__main__':
    main()
