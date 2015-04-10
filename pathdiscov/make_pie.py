#!/usr/bin/env python

from glob import glob
from collections import defaultdict
import argparse
import sys
from os.path import *
import os

from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec

def main():
    args = parse_args()
    outdir = 'host_vector_pathogen'
    if not isdir( outdir ):
        os.makedirs( outdir )
    for proj in args.project_path:
        phylo_files = glob( join(proj, 'results/iterative*/*.top.blast.phylo') )
        hvp = host_vector_pathogen(phylo_files, args.hostclasses, args.vectorclasses, args.pathogenclasses)
        outpath = join( outdir, basename(normpath(proj)) + '.png' )
        create_image( proj, hvp=hvp, output_path=outpath )

def find_best_name(phylo_line):
    '''
    Find best name from a phylo_line
    Tries first to use species or column 8(index-0)
    If species is '-' then try columns to the left until non-dash
    is found and use that
    Fallback is to use description line(although I don't think there should ever
    be a case where that happens)

    :param list phylo_line: Line from blast.phylo file already split by tab
    '''
    start_col = 8
    # Work right to left over columns only 8 - 2
    for i in range(start_col, 1, -1):
        if phylo_line[i] != '-':
            return phylo_line[i]
    if phylo_line[9] == '-':
        raise ValueError('phylo line does not contain any names to use: {0}'.format(phylo_line))
    return phylo_line[9]

def host_vector_pathogen(phylo_files, hostclass, vectorclass, pathogenclass):
    '''
    Retrieves the sums of all hosts, vectors and pathogens by reading
    the supplied list of blast tab formatted files

    :param list phylo_files: list of .blast.phylo files
    :param str|list hostclass: name of class column to consider as host
    :param str|list vectorclass: name of class column to consider as vector
    :param str|list pathogenclass: name of class column to consider as pathogen

    :return: list of host counts, sumhosts, vector counts, sumvectors, pathogen counts, sumpathogens
    '''
    hosts = defaultdict(float)
    vectors = defaultdict(float)
    pathogens = defaultdict(float)
    hosts_ct, vectors_ct, pathogens_ct = 0,0,0
    for f in phylo_files:
        for line in open(f).read().splitlines()[1:]:
            parts = line.split('\t')
            count = float(parts[1])
            sk = parts[2]
            clss = parts[4]
            species = find_best_name(parts)
            if count < 1.0:
                continue
            if clss in hostclass:
                hosts[species] += count
                hosts_ct += count
            elif clss in vectorclass:
                vectors[species] += count
                vectors_ct += count
            elif sk in pathogenclass:
                pathogens[species] += count
                pathogens_ct += count

    return hosts, hosts_ct, vectors, vectors_ct, pathogens, pathogens_ct

def graph_( ax, counts, **kwargs ):
    ax.set_title( kwargs['title'] )
    ax.pie( [ct for l,ct in counts.items()], labels=[l for l,ct in counts.items()], autopct='%1.1f%%' )

def graph_vectors( ax, vector_counts ):
    graph_( ax, vector_counts, title='Vectors' )

def graph_hosts( ax, host_counts ):
    graph_( ax, host_counts, title='Hosts' )

def graph_pathogens( ax, pathogens_counts ):
    graph_( ax, pathogens_counts, title='Pathogens' )

def graph_all( ax, host_vector_pathogen, **kwargs ):
    hc,vc,pc = host_vector_pathogen[1::2]
    hvp = [('Mamalia',hc), ('Insecta',vc), ('Virus/Bacteria',pc)]
    for k,v in hvp:
        sys.stdout.write('{0}: {1}\n'.format(k,v))
    ax.set_title( kwargs['title'] )
    ax.pie( [ct for l,ct in hvp], labels=[l for l,ct in hvp], autopct='%1.1f%%' )

def create_image( project_path, **kwargs ):
    hvp = kwargs['hvp']
    h,hc,v,vc,p,pc = hvp

    fig = plt.figure()

    fig.suptitle( basename(kwargs['output_path']) )
    fig.set_size_inches( 50.0, 10.0 )

    gs = gridspec.GridSpec( 1, 4 )
    ax1 = plt.subplot( gs[0] )
    ax2 = plt.subplot( gs[1] )
    ax3 = plt.subplot( gs[2] )
    ax4 = plt.subplot( gs[3] )

    graph_vectors( ax1, v )
    graph_hosts( ax2, h )
    graph_pathogens( ax3, p )
    graph_all( ax4, hvp, title='Sample Composition' )

    gs.tight_layout(fig, w_pad=20)
    fig.savefig( kwargs['output_path'], dpi=100, pad_inches=0.1 )

def parse_args( args=sys.argv[1:] ):
    parser = argparse.ArgumentParser(
        description = 'Generate pie graph that shows vector, host and pathogen pie graphics'
    )

    parser.add_argument(
        'project_path',
        nargs='+',
        help='Path to pathdiscov project[s] path'
    )

    parser.add_argument(
        '--hostclasses',
        nargs='+',
        default='Mammalia',
        help='Names of blast classes to use to call as host [%(default)s]'
    )

    parser.add_argument(
        '--vectorclasses',
        nargs='+',
        default='Insecta',
        help='Names of blast classes to use to call as vector [%(default)s]'
    )

    parser.add_argument(
        '--pathogenclasses',
        nargs='+',
        default='Viruses Bacteria',
        help='Names of blast classes to use to call as pathogen [%(default)s]'
    )

    return parser.parse_args( args )

if __name__ == '__main__':
    main()
