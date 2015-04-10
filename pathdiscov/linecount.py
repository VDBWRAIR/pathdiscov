from Bio import SeqIO
import argparse

def countlines(input,fileformat):
    '''
    Count lines in a file and return results
    fasta and fastq can be either single or multi-row formatted, but each sequence only
    counts as a single line

    :param file input: Open file handle to count or str of filepath to file
    :param str fileformat: same as Biopython's parse except you can specify None to signify to count all lines

    :rtype: int
    :return: Number of lines in file
    '''
    # Determine input type
    fh = input
    if isinstance(input,str):
        fh = open(input)

    # Start count at 0
    count = 0
    
    # How to count the file
    if fileformat is None:
        # Just count all lines
        for line in fh:
            count += 1
    else:
        # Count number of sequences
        for record in SeqIO.parse(fh,fileformat):
            count += 1
    return count

def linecount(inputfile, name, outfile, format, concat):
    '''
    Counts lines in a given file path based on the format and places
    the output into outfile as a tab separated file with name being the first column
    and the count from the file as the second file

    :param file inputfile: input path or file handle to count
    :param str name: name to use for first column of count file
    :param file outfile: path or filehandle to file to write results to
    :param str format: format of file being read. 1 - fastq, 2 - fasta, 0 - other(count all lines)
    :pram bool concat: 1|True - concat to outfile, 0|False - overwrite outfile
    '''
    # Determine input and output file handles and modes
    outfh = outfile
    if isinstance(outfh,str):
        if concat is True or concat == '1' or concat == 1:
            outfh = open(outfile,'a')
        else:
            outfh = open(outfile,'w')
    # Convert old style format 1, 2, 0 to fasta, fastq, None
    oldstyles = {'2':'fasta','1':'fastq','0':None,2:'fasta',1:'fastq',0:None}
    # Try to get oldstyle format from mapping or just use format
    format = oldstyles.get(format,format)
    # Count lines in file accordingly
    numlines = countlines(inputfile, format)
    # Write to the file
    outfh.write('{0}\t{1}\n'.format(name,numlines))

def parse_args():
    parser = argparse.ArgumentParser('Count lines in a file')

    parser.add_argument(
        'inputfile',
        help='Input file to count'
    )

    parser.add_argument(
        'name',
        help='Name to display in output'
    )

    parser.add_argument(
        'outfile',
        help='Path to output file to write to'
    )

    parser.add_argument(
        'format',
        help='Format of input file. 1 - fasta, 2 - fastq, 0 - other or any Biopython' \
            ' parse format to count sequence records'
    )

    parser.add_argument(
        'concat',
        help='Concat to file or write new file'
    )

    return parser.parse_args()

def main():
    args = parse_args()
    linecount(
        args.inputfile,
        args.name,
        args.outfile,
        args.format,
        args.concat
    )
