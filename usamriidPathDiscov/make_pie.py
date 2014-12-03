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
        outpath = join( outdir, basename(proj) + '.png' )
        create_image( proj, output_path=outpath )

def host_vector_pathogen( phylo_files ):
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
            species = parts[8]
            if count < 1.0:
                continue
            if clss == 'Mammalia':
                hosts[species] += count
                hosts_ct += count
            elif clss == 'Insecta':
                vectors[species] += count
                vectors_ct += count
            elif sk in ('Viruses','Bacteria'):
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
    ax.set_title( kwargs['title'] )
    ax.pie( [ct for l,ct in hvp], labels=[l for l,ct in hvp], autopct='%1.1f%%' )

def create_image( project_path, **kwargs ):
    phylo_files = glob( join(project_path, 'results/iterative*/[12].contig.top.blast.phylo') )
    hvp = host_vector_pathogen( phylo_files )
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
        help='Path to riidpipeline project[s] path'
    )

    return parser.parse_args( args )
