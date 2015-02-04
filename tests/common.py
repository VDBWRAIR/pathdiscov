from os.path import *
import tempfile
import shutil
import subprocess
import os
import sys
import glob

import unittest2 as unittest

# Testing directory
TESTDIR = dirname(abspath(__file__))

# Scratch directory
SCRATCH = join(TESTDIR, 'test_runs')

# testData directory
TESTDATA = join(dirname(TESTDIR), 'testData')

# All projects need these files regardless of -R1 -R2
# Paths are relative to --outdir
PROJECT_FILES = [
	'results/orf_filter/out.cap.fa',
	'results/orf_filter/logs/*-out.o',
	'results/orf_filter/logs/*-out.e',
	'results/orf_filter/orf_filter.R1',
	'results/quality_analysis/F_fastqc.zip',
	'results/quality_analysis/analysis_quality.log',
	'results/quality_analysis/R_fastqc.zip',
	'results/quality_analysis/F_fastqc.html',
	'results/quality_analysis/R_fastqc.html',
	'results/ray2_assembly_1/out.ray.fa.cap.concat',
	'results/ray2_assembly_1/results/FilePartition.txt',
	'results/ray2_assembly_1/results/NetworkTest.txt',
	'results/ray2_assembly_1/results/SequencePartition.txt',
	'results/ray2_assembly_1/results/RayPlatform_Version.txt',
	'results/ray2_assembly_1/results/ElapsedTime.txt',
	'results/ray2_assembly_1/results/NumberOfSequences.txt',
	'results/ray2_assembly_1/results/RaySmartCommand.txt',
	'results/ray2_assembly_1/results/RayCommand.txt',
	'results/ray2_assembly_1/results/RayVersion.txt',
	'results/ray2_assembly_1/out.cap.fa',
	'results/ray2_assembly_1/contig_numreads.txt',
	'results/ray2_assembly_1/head.1.R1.unmap.fastq',
	'results/ray2_assembly_1/assembly.count',
	'results/ray2_assembly_1/out.ray.fa',
	'results/ray2_assembly_1/R2.single.fastq',
	'results/ray2_assembly_1/R2.paired.fastq',
	'results/ray2_assembly_1/logs_assembly/assembly.o',
	'results/ray2_assembly_1/logs_assembly/assembly.e',
	'results/ray2_assembly_1/ray2_assembly_1.fasta',
	'results/ray2_assembly_1/bowtie2_mapping/R2.unmap.id',
	'results/ray2_assembly_1/bowtie2_mapping/out.bam',
	'results/ray2_assembly_1/bowtie2_mapping/R1.map.id',
	'results/ray2_assembly_1/bowtie2_mapping/R1.unmap.id',
	'results/ray2_assembly_1/bowtie2_mapping/R2.map.id',
	'results/ray2_assembly_1/cap3.out',
	'results/ray2_assembly_1/head.1.R2.unmap.fastq',
	'results/ray2_assembly_1/logs/*-out.o',
	'results/ray2_assembly_1/logs/*-out.e',
	'results/ray2_assembly_1/R1.single.fastq',
	'results/ray2_assembly_1/R1.paired.fastq',
	'results/ray2_assembly_1/contig_len.txt',
	'results/output/{projectdir}.aug.report.txt',
	'results/output/2.contig.noblast.fasta',
	'results/output/out.cap.fa',
	'results/output/R1.{projectdir}.top.phylo.txt',
	'results/output/R2.{projectdir}.top.phylo.txt',
	'results/output/{projectdir}.R2.count.txt',
	'results/output/{projectdir}.R1.count.txt',
	'results/output/contig.{projectdir}.top.report.txt',
	'results/output/contig.{projectdir}.top.phylo.txt',
	'results/output/stats.txt',
	'results/output/reads.txt',
	'results/output/qstats.txt',
	'results/iterative_blast_phylo_1/2.contig.top.blast',
	'results/iterative_blast_phylo_1/2.contig.noblast.fasta',
	'results/iterative_blast_phylo_1/1.contig.top.blast',
	'results/iterative_blast_phylo_1/reports/contig.{projectdir}.phylo.txt',
	'results/iterative_blast_phylo_1/reports/contig.{projectdir}.top.report.txt',
	'results/iterative_blast_phylo_1/reports/contig.{projectdir}.top.phylo.txt',
	'results/iterative_blast_phylo_1/reports/contig.{projectdir}.top.smallreport.txt',
	'results/iterative_blast_phylo_1/contig.count.superclass',
	'results/iterative_blast_phylo_1/1.contig.fasta',
	'results/iterative_blast_phylo_1/contig.count',
	'results/iterative_blast_phylo_1/logs/*-out.o',
	'results/iterative_blast_phylo_1/logs/*-out.e',
	'results/iterative_blast_phylo_1/3.contig.fasta',
	'results/iterative_blast_phylo_1/iterative_blast_phylo_1.contig',
	'results/iterative_blast_phylo_1/contig.top.count.superclass',
	'results/iterative_blast_phylo_1/2.contig.fasta',
	'results/iterative_blast_phylo_1/1.contig.noblast.fasta',
	'results/step1/R2.fastq',
	'results/step1/R1.count',
	'results/step1/R1.id',
	'results/step1/R2.count',
	'results/step1/R2.id',
	'results/step1/step1.R2',
	'results/step1/R1.fastq',
	'results/step1/logs/*-out.o',
	'results/step1/logs/*-out.e',
	'results/step1/step1.R1',
	'results/analysis.log',
	'results/iterative_blast_phylo_2/reports/{projectdir}.joinR1R2.smallreport.txt',
	'results/iterative_blast_phylo_2/logs/*-out.o',
	'results/iterative_blast_phylo_2/logs/*-out.e',
	'results/quality_filter/logs/*-out.o',
	'results/quality_filter/logs/*-out.e',
	'results/logs/{projectdir}.*_param.txt',
	'results/host_map_1/host_map_1.R1',
	'results/host_map_1/R1.discard',
	'results/host_map_1/map_1/R2.unmap.id',
	'results/host_map_1/map_1/out.bam',
	'results/host_map_1/map_1/R1.map.id',
	'results/host_map_1/map_1/R1.unmap.id',
	'results/host_map_1/map_1/R2.unmap.fastq',
	'results/host_map_1/map_1/R1.unmap.fastq',
	'results/host_map_1/map_1/R2.map.id',
	'results/host_map_1/map_2/out.bam',
	'results/host_map_1/map_2/singleton.map.id',
	'results/host_map_1/map_2/R1.unmap.fastq',
	'results/host_map_1/map_2/singleton.unmap.id',
	'results/host_map_1/R1.count',
	'results/host_map_1/R2.discard',
	'results/host_map_1/R2.count',
	'results/host_map_1/logs/*-out.o',
	'results/host_map_1/logs/*-out.e',
	'results/host_map_1/host_map_1.R2',
	'quality_analysis',
	'output',
	'R1.count',
	'contig_reports',
	'R2.count',
	'analysis.log',
	'unassembled_read_reports',
	'input/param.txt',
]

def verify_project(projpath, expectedfiles, r2=True):
    '''
    Just makes sure all files are present in project
    if r2 is false then ignore files that have R2 in them
    expectedfiles should contain paths to files that are relative to
    projpath(like [join(projpath, f) for f in expectedfiles])
    '''
    # Will use this to replace {projdir}
    projname = basename(projpath)
    # Keep track of missing files
    missing = []
    for f in expectedfiles:
        # Skip R2 file checks if specified
        if not r2 and ('R2' in f or 'R_' in f):
            continue
        p = join(projpath, f)
        try:
            # Attempt to replace projectname if exists
            p = p.format(projectdir=projname)
        except:
            pass
        # Expand glob expressions if there
        files = glob.glob(p)
        if not files:
            # Append since glob didn't match anything
            missing.append(p)
        else:
            for f in files:
                if not exists(f):
                    missing.append(f)
    return missing

def print_list(lst):
    '''
    Print list each item on own line
    '''
    print '\n'.join(lst)

def run_path_discov(args):
    '''
    Run usamriidPathDiscov_cli with args
    '''
    cmd = ['usamriidPathDiscov_cli'] + args
    print "Running {0}".format(' '.join(cmd))
    try:
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        # Return stdout, stderr
        sout, serr = p.communicate()
        ret = p.returncode
    except OSError as e:
        print "error running {0}, maybe command not in path?".format(cmd)
        sout, serr = ('',str(e))
        ret = -1
    return sout, serr, ret

class TempDir(unittest.TestCase):
    def setUp(self):
        self.keep_temp_dir = False
        if not isdir(SCRATCH):
            os.makedirs(SCRATCH)
        # Create temporary directory to run tests in inside scratch dir
        self.testdir = tempfile.mkdtemp(suffix=__name__, dir=SCRATCH)
        # Change directory into the testdir so accidental files will be
        # put into it hopefully instead of somwhere random
        os.chdir(self.testdir)

    def tearDown(self):
        os.chdir(TESTDIR)
        if self.keep_temp_dir:
            print "Not removing {0} because keep_temp_dir set".format(self.testdir)
        else:
            shutil.rmtree(self.testdir)
