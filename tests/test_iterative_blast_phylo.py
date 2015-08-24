import os
from os.path import *
import sys
import shutil
import inspect
import re

import unittest2 as unittest
import sh

from common import (
    TESTDATA, SCRATCH, PATHDISCOV, RIKKDB, TESTDIR,
    aexists, attr
)
import common

# Allows us to reference sh.quality_filter since quality_filter.pl
# cannot just work with sh import because of .pl
sh.iterative_blast_phylo = sh.Command('iterative_blast_phylo.pl')

blastnt = join(RIKKDB, 'rikkcdna')
blastnr = join(RIKKDB, 'rikkrna')
diamond = join(RIKKDB, 'rikknr')
taxnames = join(RIKKDB, 'names.dmp')
taxnodes = join(RIKKDB, 'nodes.dmp')

# Example contig fasta input file
fasta_input = join(TESTDIR, 'out.cap.fa')

# Param
paramtxt = '''
command iterative_blast_phylo
blast_db_list				{0},{0},{1}
blast_task_list				megablast,dc-megablast,blastx
blast_options_list			-evalue 1e-4 -word_size 28,-evalue 1e-4 -word_size 12,-evalue 1e-4 -word_size 12
ninst_list				    1,1,1
taxonomy_names				{2}
taxonomy_nodes				{3}
blast_pro_db                {4}
command orf_filter
getorf_options				-minsize 60 -find 0
'''.format(blastnt,blastnr,taxnames,taxnodes,blastnr)

expect_count_contig = [
    ('input', '86'),
    ('megablast', '76'),
    ('dc-megablast', '73'),
    ('orf_filter',73),
    ('blastx', 50)
]

expect_count_r1 = [
    ('input', '250'),
    ('megablast', '216'),
    ('dc-megablast', '214'),
    ('orf_filter','196'),
    ('blastx', '185')
]

expect_count_r2 = [
    ('input', '250'),
    ('megablast', '222'),
    ('dc-megablast', '218'),
    ('orf_filter', '201'),
    ('blastx', '183')
]

output_contig = 'iterative_blast_phylo_1.contig'
output_r1 = 'iterative_blast_phylo_1.R1'
output_r2 = 'iterative_blast_phylo_1.R2'
r1_noblast = '2.R1.noblast.fasta'
r2_noblast = '2.R2.noblast.fasta'
contig_noblast = '2.contig.noblast.fasta'
r1_smallreport = 'reports/R1.{0}.top.smallreport.txt'
r2_smallreport = 'reports/R2.{0}.top.smallreport.txt'
contig_smallreport = 'reports/contig.{0}.top.smallreport.txt'

# Options
# --contig use contig as name instead of r1/r2
# --run_iteration [Default: 1]
# --fasta same as --fastafile
# --fastafile
class TestIterativeBlastPhylo(common.StageTestBase):
    def test_iterative_blast_phylo_contig(self):
        param = self._write_param(paramtxt)
        outdir = inspect.stack()[0][3]
        logs = join(outdir,'logs')
        os.makedirs(logs)
        sh.iterative_blast_phylo(
            sample='testsample', paramfile=param,
            outputdir=outdir, logs=logs, timestamp='0',
            R1=fasta_input, fastafile='yes', contig='1'
        )
        aexists(join(outdir,output_contig))
        aexists(join(outdir,contig_noblast))
        aexists(join(outdir,contig_smallreport.format('testsample')))
        self._verify_countfile(expect_count_contig, join(outdir, 'contig.count'))

    def test_iterative_blast_phylo_r1r2(self):
        param = self._write_param(paramtxt)
        outdir = inspect.stack()[0][3]
        logs = join(outdir,'logs')
        os.makedirs(logs)
        sh.step1(
            sample='testsample', paramfile=param,
            outputdir=outdir, logs=logs, timestamp='0',
            R1=self.r1_fq, R2=self.r2_fq
        )
        r1 = join(outdir,'R1.fastq')
        r2 = join(outdir,'R2.fastq')
        sh.iterative_blast_phylo(
            sample='testsample', paramfile=param,
            outputdir=outdir, logs=logs, timestamp='0',
            R1=r1, R2=r2, fastafile='no', contig='0'
        )
        aexists(join(outdir,output_r1))
        aexists(join(outdir,output_r2))
        aexists(join(outdir,r1_noblast))
        aexists(join(outdir,r2_noblast))
        aexists(join(outdir,r1_smallreport.format('testsample')))
        aexists(join(outdir,r2_smallreport.format('testsample')))
        self._verify_countfile(expect_count_r1, join(outdir, 'R1.count'))
        self._verify_countfile(expect_count_r2, join(outdir, 'R2.count'))

    @attr('slow')
    def test_iterative_blast_getorf_diamond(self):
        paramtxt = '''
        command iterative_blast_phylo
        blast_db_list				{0},{0},{1}
        blast_task_list				megablast,dc-megablast,diamond
        blast_options_list			-evalue 1e-4 -word_size 28,-evalue 1e-4 -word_size 12,--compress 0 -p 1 -v -k 10 -w 28 --id 0.7 -c 4
        ninst_list				    1,1,1
        taxonomy_names				{2}
        taxonomy_nodes				{3}
        blast_pro_db                {4}
        command orf_filter
        getorf_options				-minsize 60 -find 0
        '''.format(blastnt,diamond,taxnames,taxnodes,blastnr)
        
        param = self._write_param(paramtxt)
        outdir = inspect.stack()[0][3]
        logs = join(outdir,'logs')
        os.makedirs(logs)
        sh.iterative_blast_phylo(
            sample='testsample', paramfile=param,
            outputdir=outdir, logs=logs, timestamp='0',
            R1=fasta_input, fastafile='yes', contig='1'
        )
        aexists(join(outdir,output_contig))
        aexists(join(outdir,contig_noblast))
        aexists(join(outdir,contig_smallreport.format('testsample')))
        orf_filter_dir = join(outdir,'tmp_contig_3','orf_filter')
        aexists(join(orf_filter_dir, 'orf_filter.contig'))
        aexists(join(orf_filter_dir, 'contig.orfout.fa'))
        # copy so we can modify and exclude orf_filter and blastx
        # counts as we will have our own
        econtig = expect_count_contig[:-1]
        econtig.append(('diamond','72'))
        self._verify_countfile(econtig, join(outdir, 'contig.count'))
