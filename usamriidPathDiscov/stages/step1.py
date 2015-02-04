import argparse
from os.path import exists, isdir, join, dirname, basename, abspath
import os

from Bio import SeqIO

from usamriidPathDiscov import util
from usamriidPathDiscov.sff2fastq import sff_to_fastq
from usamriidPathDiscov.linecount import linecount

def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        '--sample',
        required=True,
        help='Sample name'
    )

    parser.add_argument(
        '--outputdir',
        required=True,
        help='The output directory. All of the output will be ' \
            'written in this directory'
    )

    parser.add_argument(
        '--logs',
        required=True,
        help='The logs directory. stdout and stderr will be written ' \
            ' in this directory'
    )

    parser.add_argument(
        '--paramfile',
        required=True,
        help='The parameter file for the command'
    )

    parser.add_argument(
        '--R1',
        required=True,
        help='.fastq or .sff input file. If paired reads, this file represents' \
            'the R1 reads. Can be zipped or not.'
    )

    # It seems R1 and R2 should just be an nargs, but for backwards compatibility
    # they are separate
    parser.add_argument(
        '--R2',
        help='.fastq or .sff input file of the R2 reads. Can be zipped or not.' \
            'Required only for paired reads.'
    )

    # Why not just set --sample with this unique number?
    parser.add_argument(
        '--timestamp',
        required=True,
        help='A time stamp. This ensures your log files will be unique.'
    )

    parser.add_argument(
        '--example',
        help='Prints example parameter file.'
    )

    return parser.parse_args()

def step1(mates, outputdir, logs, samplename, timestamp, paramfile):
    '''
    Decompresses input files and rewrites them so they only contain sequential
    integers as their ids

    :param list mates: List of filepaths to R1 and R2 input files.
        Can be fastq, fastq.gz, sff, sff.gz
    :param str outputdir: Path to output directory for this stage
    :param str logs: log output directory
    :param str samplename: name for output logs
    :param str timestamp: unique timestamp for output logs
    '''
    if not isdir(outputdir):
        raise OSError('cannot chdir to {0}'.format(outputdir))
    if not isdir(logs):
        raise ValueError('{0} does not exist'.format(logs))

    logbase = join(logs, '{0}.{1}-out'.format(samplename,timestamp))
    stdoutlog = open(logbase+'.o','w')
    stderrlog = open(logbase+'.e','w')
    
    # Iter over mate list
    for matenum, mate in enumerate(mates, start=1):
        # If mate passed is 'none' then same as not having --R2
        if mate == "none":
            continue
        # Get easy reference to mate name
        matename = 'R{0}'.format(matenum)
        # Get easy reference to current mate paths (R1|R2).*
        matepath = join(outputdir, matename + '.fastq')
        countpath = join(outputdir, matename + '.count')
        idpath = join(outputdir, matename + '.id')
        step1mate = join(outputdir, 'step1.' + matename)
        # Get common file handle to either gzip/normal file
        mate_fh = util.get_file_handle(mate)
        # Create temp fastq file
        if '.sff' in mate:
            # Get tmp name
            tmpfq = join(outputdir, matename + '.sff' + '.fastq')
            # Touch file to ensure exists
            open(tmpfq,'w')
            # Open read/write tmp file
            tmp_fh = open(tmpfq, 'r+')
            # Write records
            # -- hack for an issue between gzip and SffIO.parse
            #    dealing with gzip streams
            import gzip
            gzip.READ = 'r'
            for rec in SeqIO.parse(mate_fh, 'sff'):
                SeqIO.write(rec, tmp_fh, 'fastq')
            # Set mate_fh to tmp opened file for later
            mate_fh.close()
            mate_fh = tmp_fh
            mate_fh.seek(0)
        # Create R1|R2.fastq with modified ids
        # also creates R1|R2.id
        util.change_fastq_id(mate_fh, matepath, idpath)
        mate_fh.seek(0)
        # Create R1|R2.count files
        linecount(mate_fh, 'rawfile', countpath, 'fastq', False)

        if not exists(step1mate):
            os.symlink(basename(matepath), step1mate)

    stdoutlog.close()
    stderrlog.close()

def main():
    args = parse_args()
    mates = [args.R1]
    if args.R2:
        mates.append(args.R2)
    step1(
        mates,
		args.outputdir,
		args.logs,
		args.sample,
		args.timestamp,
		args.paramfile
    )

if __name__ == '__main__':
    main()
