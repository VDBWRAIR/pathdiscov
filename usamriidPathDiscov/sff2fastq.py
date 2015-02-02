import argparse
from Bio import SeqIO
import itertools
import gzip
from os.path import splitext

def sff_to_fastq(sfffiles, outfile):
    '''
    Convert list of sfffiles into fastq and write all records into outfile

    :param list sfffiles: list of sfffile paths
    :param file outfile: file handle or path to output file to write fastq records

    :rtype: int
    :return: The number of records written
    '''
    outfh = outfile
    if isinstance(outfile,str):
        outfh = open(outfile,'w')
    # Creates one continuous stream of seqrecords from all input sff
    all_records = itertools.chain.from_iterable([SeqIO.parse(*file_handle(sff,'rb')) for sff in sfffiles])
    # Write all records to outfh as fastq
    numwritten = SeqIO.write(all_records, outfh, 'fastq')
    return numwritten

def file_handle(filepath, mode):
    '''
    :param str filepath: path to file to open
    :param str mode: file mode to open with
    :return: (opened file handle, file extension)
    '''
    root, ext = splitext(filepath.replace('.gz',''))
    # Remove the period in extension
    ext = ext[1:]
    if filepath.endswith('.gz'):
        handle = gzip.open(filepath, mode)
    else:
        handle = open(filepath, mode)
    return (handle, ext)

def parse_args():
    parser = argparse.ArgumentParser('Convert sff to fastq')

    parser.add_argument(
        'sfffile',
        nargs='+',
        help='Sff file to convert to fastq'
    )

    parser.add_argument(
        'outfile',
        help='Output fastq file'
    )

    return parser.parse_args()

def main():
    args = parse_args()
    sff_to_fastq(args.sfffile, args.outfile)
