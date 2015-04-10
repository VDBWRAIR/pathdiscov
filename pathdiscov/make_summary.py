#!/usr/bin/env python

import argparse
import os
from os.path import join, basename, dirname, exists, isdir, isfile
import sys
from glob import glob
import re
import csv
from collections import defaultdict
import itertools

# Exception class for when project files are missing
class MissingProjectFile( Exception ): pass

def get_n50(numlist):
    """
    :param list numlist: List of numbers

    Returns the N50 value of the passed list of numbers. 
    Code from http://seqanswers.com/forums/showpost.php?p=9414&postcount=7

    Based on the Broad Institute definition:
    https://www.broad.harvard.edu/crd/wiki/index.php/N50
    """
    numlist.sort()
    newlist = []
    for x in numlist :
        newlist += [x]*x
    # take the mean of the two middle elements if there are an even number
    # of elements.  otherwise, take the middle element
    if len(newlist) % 2 == 0:
        medianpos = len(newlist)/2  
        return float(newlist[medianpos] + newlist[medianpos-1]) /2
    else:
        medianpos = len(newlist)/2
        return newlist[medianpos]

def parse_tab_file( tabfile, fields=None ):
    delimiter = '\t'
    with open( tabfile ) as fh:
        if fields is not None:
            reader = csv.DictReader( fh, fieldnames=fields, delimiter=delimiter )
        else:
            reader = csv.reader( fh, delimiter=delimiter )
        for row in reader:
            yield row

def read_count( countfile ):
    '''
    Get counts for every line in count file
    returns a dictionary keyed by every item in the first column with value of the item in right column
    '''
    counts = {}
    for row in parse_tab_file( countfile, ['name', 'count'] ):
        v = row['count']
        try:
            counts[row['name']] = int( v )
        except ValueError as e:
            rv = round(float(v))
            sys.stderr.write( '{0} had a non integer {1} and was rounded to {2}\n'.format(countfile,v,rv) )
            counts[row['name']] = rv
    return counts

def total_reads( projdir ):
    ''' Number of reads total '''
    c = r1r2_count( join( projdir, 'results', 'step1' ) )
    return c['rawfile']

def non_host_num_reads( projdir ):
    ''' Number of non host reads from prinseq '''
    c = r1r2_count( join(projdir, 'results', 'quality_filter') )
    return c['prinseq']

def r1r2_count( steppath ):
    '''
    Merge together counts for R1 and R2 inside of steppath
    Raises MissingProjectFile if < 1 .count file found
    '''
    r1r2 = glob( join( steppath, '*.count' ) )
    if len(r1r2) == 0:
        raise MissingProjectFile( '{0} did not contain > 1 count files'.format(steppath) )
    counts = defaultdict( int )
    for f in r1r2:
        rc = read_count( f )
        for k,v in rc.items():
            counts[k] += int(v)
    return counts

def num_contig( projdir ):
    ''' Number of contig and blast contig from input and dc-blast '''
    contigdir = join( projdir, 'results', 'iterative_blast_phylo_1' )
    counts = r1r2_count( contigdir )
    return counts['input'], counts['dc-megablast']

def parse_blast_report( reportpath, filterf=None ):
    '''
    Return row by row results from a blast report in tab format

    @param filterf - Function that will filter the results being supplied one row
    and returns True for rows to keep and False for rows to discard
    '''
    with open( reportpath ) as fh:
        for row in csv.DictReader( fh, delimiter='\t' ):
            if filterf is None or filterf( row ):
                yield row

def blast_results_for_( blastfile, blastcol, filterval ):
    '''
    Return rows from blast output that match row[blastcol] == filterval
    '''
    def blastfilter(row):
        if filterval is None:
            return True
        elif row[blastcol] == filterval:
            return True
        else:
            return False

    for row in parse_blast_report( blastfile, blastfilter ):
        yield row

def contig_info( projdir ):
    '''
    Merge together results from contig.id and contig_numreads.txt file
    and return as a dictionary keyed by contigname
    '''
    d = join( projdir, 'results', 'ray2_assembly_1' )
    contiglenfile = join( d, 'contig_len.txt' )
    contignreadsfile = join( d, 'contig_numreads.txt' )
    info = {}
    if not (exists(contiglenfile) and exists(contignreadsfile)):
        raise MissingProjectFile( '{0} or {1} are missing from {2}'.format(contiglenfile,contignreadsfile,projdir) )
    cl = list( parse_tab_file( contiglenfile, ['contig', 'lenstr'] ) )
    cr = list( parse_tab_file( contignreadsfile, ['contig', 'nreads'] ) )
    # Sort on numeric portion of contig name
    sortkey = lambda x: int(x['contig'][1:])
    # Make sure sorted
    cl.sort( key=sortkey )
    cr.sort( key=sortkey )
    # Sometimes one list is not the same length as the other
    for l, nr in itertools.izip_longest( cl, cr, fillvalue={'contig':'','nreads':-1,'lenstr':'contig-c000000 -1 nucleotides'} ):
        contigname = l['contig']
        length = l['lenstr']
        nreads = nr['nreads']
        info[contigname] = (length,nreads)
    return info

def contigs_for( projdir, blastcol, blastval ):
    '''
    Merge together contig information from various sources and return
    information [{contigname, length, numreads, accession, family, genus, description},...]

    @param blastcol - Column in blast tab file that will be used to filter
    @param blastval - Value to filter the blastcol on
    '''
    # Contains contigname, length, numreads for every contig
    ci = contig_info( projdir )
    # Contains blast information(accession, family, genus, description)
    smallreport = glob( join( projdir, 'results', 'iterative_blast_phylo_1', 'reports', '*smallreport*.txt' ) )
    if len( smallreport ) != 1:
        bns = ' '.join( [basename(f) for f in smallreport] )
        raise MissingProjectFile( '{0} only has the following phylo_1 smallreport files: {1}'.format(projdir, bns) )
    smallreport = smallreport[0]
    bi = blast_results_for_( smallreport, blastcol, blastval )
    # Iterate over blast results that are filtered by blastcol and blastval
    for contig in bi:
        info = {}
        name = contig['qseqid']
        # Exception is here
        l, nr = ci[name]
        info['contigname'] = name
        info['length'] = int( l )
        info['numreads'] = int( nr )
        info['accession'] = contig['sseqid'].split('|')[3]
        info['family'] = contig['family']
        info['genus'] = contig['genus']
        info['description'] = contig['descrip']
        yield info

def unassembled_reads( projdir ):
    ''' Returns counts for total unassembled reads and num blast unassembled '''
    d = join( projdir, 'results', 'iterative_blast_phylo_2' )
    counts = r1r2_count( d )
    return counts['input'], counts['dc-megablast']

def group_blast_by_( blastfile, filtercol, filterval, groupbycol ):
    '''
    Filter blast results down by filtercol == filterval then
    group results by groupbycol
    '''
    # Will sort the results by groupbycol and put in list
    results = blast_results_for_( blastfile, filtercol, filterval )
    grouped = defaultdict( dict )
    for blastrow in results:
        groupval = blastrow[groupbycol]
        grouped[groupval].update( **blastrow )
        if 'count' not in grouped[groupval]:
            grouped[groupval]['count'] = 0
        grouped[groupval]['count'] += 1

    return grouped

def unassembled_report( projdir, kingdom, groupby='family' ):
    '''
    Returns the grouped blast results filtered first by kingdom and then
    grouped by groupby field
    Each key in the returned dictionary will be the groupby column uniq values
    Combines R1 & R2 results
    Adds a new key called accession with the parsed out accession
    '''
    smallreports = glob( join( projdir, 'results', 'iterative_blast_phylo_2', 'reports', 'R[12].*.top.smallreport.txt' ) )
    if len(smallreports) < 1:
        bns = ' '.join( [basename(f) for f in smallreports] )
        raise MissingProjectFile( '{0} only has the following phylo_2 smallreport files: {1}'.format(projdir, bns) )
    merged = {}
    for f in smallreports:
        # Begin with first small report
        #merged = group_blast_by_( f, 'superkingdom', kingdom, groupby )
        # Now update by adding any new keys and incrementing count for the ones that already exist
        for k, v in group_blast_by_( f, 'superkingdom', kingdom, groupby ).iteritems():
            # Already in merged so just increment count
            if k in merged:
                merged[k]['count'] += v['count']
            else:
                merged[k] = v
            # Make sure accession was parsed and put in
            if 'accession' not in merged[k]:
                accession = merged[k]['sseqid'].split('|')[3]
                merged[k]['accession'] = accession
    return merged

def summary( projdir, filtercol, filterval, groupby='family' ):
    summary = {}

    summary['numreads'] = total_reads( projdir )
    summary['nonhostreads'] = non_host_num_reads( projdir )
    summary['numcontig'], summary['numblastcontig'] = num_contig( projdir )
    summary['numreadsunassembled'], summary['numblastunassembled'] = unassembled_reads( projdir )
    summary['contigs'] = list( contigs_for( projdir, filtercol, filterval ) )
    summary['unassembled'] = dict( unassembled_report( projdir, filterval ) )
    contiglengths = [contiginfo['length'] for contiginfo in summary['contigs']]
    summary['n50'] = get_n50(contiglengths)
    summary['assemblylength'] = sum(contiglengths)

    return summary

def format_summary( summary ):
    '''
    Format summary into multiple rows based on longest of contig or unassembled reads
    Formats as rows of tab separated values
    '''
    rows = []
    import itertools
    # Iterate over longsest of the two and fill the other in with ''
    contigkeys = ('contigname','length','numreads','accession','family','genus','description')
    unasskeys = ('count','accession','family','genus','descrip')
    prefix = format_dict( summary, ('numreads','nonhostreads','numcontig','numblastcontig','n50','assemblylength') )
    unassembled = sorted( summary['unassembled'].items(), key=lambda x: x[1]['count'], reverse=True )
    for contig, unassembled in itertools.izip_longest( summary['contigs'], unassembled, fillvalue=None ):
        # Start a new row
        rows.append('\t'+prefix+'\t')
        # Put the contig information in
        if contig is not None:
            rows[-1] += format_dict( contig, contigkeys )
        else:
            # Insert blank cells when there is no more contig info but more unassembled
            rows[-1] += '\t'*len(contigkeys)

        # Then the number of unassembled reads
        if prefix[0] == '\t':
            rows[-1] += '\t'*2
        else:
            rows[-1] += '\t' + format_dict( summary, ('numreadsunassembled','numblastunassembled') )

        # Then the unassembled reads
        if unassembled is not None:
            rows[-1] += '\t' + format_dict( unassembled[1], unasskeys )
        else:
            rows[-1] += '\t'*4

        # Only the first time should prefix have values
        prefix = '\t'*5
    return rows

def format_dict( contig, keys ):
    format = []
    for k in keys:
        format.append( str(contig[k]) )
    return '\t'.join( format )

def main( ):
    args = parse_args()
    hdr = ['Sample Name']
    # These come from summary
    hdr += ['Num Reads', 'Non-Host Num reads', 'Num Ctg', 'Num blast0 Ctg', 'N50', 'Assembly Length']
    # These come from summary['contig']
    hdr += ['Ctg#', 'Ctg bp', 'numReads', 'Accession', 'Family', 'Genus', 'description']
    # These come from summary
    hdr += ['Num unassem', 'Num blast0 Unassem']
    # These come from summary['unassembled']
    hdr += ['num reads', 'Accession', 'Family', 'Virus Genus', 'descrip']
    print '\t'.join( hdr )
    for p in args.projdir:
        try:
            sys.stderr.write( p + '\n' )
            s = summary( p, args.filter_column, args.filter_value, args.group_by )
            rows = format_summary( s )
            samplename = basename(p)
            print samplename + ('\n'+samplename).join( rows )
        except Exception as e:
            import traceback
            traceback.print_exc()
            sys.stderr.write( 'Skipping {0} because of the above error\n'.format(p) )
        sys.stderr.write( '\n' )

def parse_args( args=sys.argv[1:] ):
    parser = argparse.ArgumentParser(
        description='Generate summary report for riidpipeline sample'
    )

    parser.add_argument(
        'projdir',
        nargs='+',
        help='Project directory path for riidpipeline project'
    )

    parser.add_argument(
        '--filter-column',
        dest='filter_column',
        default='superkingdom',
        help='What column in the blast reports to filter by[Default: %(default)s]'
    )

    parser.add_argument(
        '--filter-value',
        dest='filter_value',
        default=None,
        help='What value to filter on in the --filter-column[Default: %(default)s]'
    )

    parser.add_argument(
        '--group-by',
        dest='group_by',
        default='family',
        help='What column to group results by for consolidation[Default: %(default)s]'
    )

    return parser.parse_args( args )
