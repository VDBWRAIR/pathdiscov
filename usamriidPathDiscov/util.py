from Bio import SeqIO
import gzip
from os.path import splitext

def seqrec_id_to_num(input, start=1):
    '''
    Generator that converts input fastq sequence id's to sequential integers

    :param file input: fastq file to convert
    :param int start: What number to start the sequence ids with
    :rtype: generator
    :return: Generator of same input sequence records with id replaced with sequential integer
    '''
    count = start
    for rec in SeqIO.parse(input, 'fastq'):
        origid = rec.id
        # Set record id as a string
        rec.id = str(count)
        yield (origid, rec)
        count += 1

def change_fastq_id(input, output, idmapfile):
    '''
    Convert input fastq file such that:
    - Input sequence id's are converted to sequential integers starting with 1
    - All . in sequences replaced with N
    - Ensures + line in fastq record contains only + and nothing else

    Writes modified sequence records to output and writes mapping of::
        
        1\tid1
        2\tid2
        ...

    :param input: input fastq
    :type input: str or file
    :param output: output fastq
    :type output: str or file
    :param idmapfile: output new id integer\told id map
    :type idmapfile: str or file
    '''
    infh = input
    outfqfh = output
    outmapfh = idmapfile
    if isinstance(input,str):
        infh = open(input)
    if isinstance(output,str):
        outfqfh = open(output,'w')
    if isinstance(idmapfile,str):
        outmapfh = open(idmapfile,'w')

    for origid, record in seqrec_id_to_num(infh):
        # Write mapping
        outmapfh.write('{0}\t{1}\n'.format(record.id,record.description))
        # Remove description
        record.description = ''
        # Convert . -> N in sequence
        record.seq._data = record.seq._data.replace('.','N')
        # Write new sequence file
        SeqIO.write(record, outfqfh, 'fastq')

def get_common_uneven_files(input1, input2, inputformat, input1single, 
        input1paired, input2single, input2paired):
    '''
    Given two fasta/fastq files, find all paired and unpaired reads 
    for each file

    :param input1: fasta or fastq file
    :type input1: str or file
    :param input2: fasta or fastq file
    :type input2: str or file
    :param str inputformat: fasta or fastq
    :param input1single: input1's unpaired output
    :type input1single: str or file
    :param input1paired: input1's paired ouput
    :type input1paired: str or file
    :param input2single: input2's unpaired output
    :type input1single: str or file
    :param input2paired: input2's paired output
    :type input2paired: str or file
    '''
    i1fh = input1
    i2fh = input2
    o1sfh = input1single
    o1pfh = input1paired
    o2sfh = input2single
    o2pfh = input2paired
    if isinstance(input1,str):
        i1fh = open(input1)
    if isinstance(input2,str):
        i2fh = open(input2)
    if isinstance(input1single,str):
        o1sfh = open(input1single,'w')
    if isinstance(input1paired,str):
        o1pfh = open(input1paired,'w')
    if isinstance(input2single,str):
        o2sfh = open(input2single,'w')
    if isinstance(input2paired,str):
        o2pfh = open(input2paired,'w')

    # Retrieve paired, nunpaired ids
    p, s1, s2 = get_paired_reads(i1fh, i2fh, inputformat)
    # Seek to beginning of files as they need to be iterated again
    i1fh.seek(0)
    i2fh.seek(0)
    # Iterate over input1 and input2 and extract paired reads
    for p1r, p2r in zip(extract_reads(i1fh,inputformat,p),extract_reads(i2fh,inputformat,p)):
        SeqIO.write(p1r, o1pfh, inputformat)
        SeqIO.write(p2r, o2pfh, inputformat)
    # Write single records for input1
    i1fh.seek(0)
    for record in extract_reads(i1fh,inputformat,s1):
        SeqIO.write(record, o1sfh, inputformat)
    # Write single records for input2
    i2fh.seek(0)
    for record in extract_reads(i2fh,inputformat,s2):
        SeqIO.write(record, o2sfh, inputformat)

def get_paired_reads(input1, input2, inputformat):
    '''
    Find all paired and unpaired reads between two fastq or fasta files

    :param file input1: fasta or fastq input
    :param file input2: fasta or fastq input
    :param str inputformat: fasta or fastq
    :return tuple: paired ids, input1 single, input2 single, each is a set of read names
    '''
    id1 = set()
    paired = set()
    single1 = set()
    single2 = set()
    # Get all ids from first file
    for rec in SeqIO.parse(input1, inputformat):
        id1.add(rec.id)
    # Then determine paired, and singles
    for rec in SeqIO.parse(input2, inputformat):
        # In both file so paired
        if rec.id in id1:
            paired.add(rec.id)
        else:
            single2.add(rec.id)
    # Ids from input1 that are not paired
    single1 = id1 - paired
    return paired, single1, single2

def extract_reads(input, inputformat, idlist):
    '''
    Get all records from input that have id in idlist

    :param file input: input fasta or fastq file
    :param str inputformat: fasta or fastq
    :param set or dict idlist: any hashable sequence that can be used to look
        up ids
    :return generator: yields sequence records that are only in idlist
    '''
    for rec in SeqIO.parse(input,inputformat):
        if rec.id in idlist:
            yield rec

def get_file_handle(inputpath):
    '''
    Returns a normalize file handle for either gzipped or normal files.
    Allows easy opening of files as you don't need to know if gzipped
    or if normal.

    Inspects the file extension only(.gz == gzipped)

    :param str inputpath: input file path. If ends in .gz decrompress
    :rtype: file
    :return: opened file object representing decompressed data or normal data
    '''
    base, ext = splitext(inputpath)
    if ext == '.gz':
        fh = gzip.open(inputpath)
        # Fix a biopython bug
        if base.endswith('.sff'):
            gzip.READ = 'r'
            fh.mode = 'r'
        return fh
    else:
        return open(inputpath)
