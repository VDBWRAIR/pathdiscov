from os.path import *
import fileinput
import sys
import os
import re

import unittest2 as unittest
from nose.plugins.attrib import attr
import mock

import common

DBDIR = join(common.TESTDIR, 'rikkcdna')
RIKKNT = join(DBDIR, 'rikkcdna')
RIKKNR = join(DBDIR, 'rikkrna')
DIAMOND = join(DBDIR, 'rikknr')
SNAPDNA = join(DBDIR, 'rikkdna')
SNAPRNA = join(DBDIR, 'rikkrna')
BOWTIEDNA = join(DBDIR, 'rikkdna')
BOWTIERNA = join(DBDIR, 'rikkrna')

@attr('slow')
class TestInputOutputArguments(common.TempDir):
    def setUp(self):
        super(TestInputOutputArguments, self).setUp()
        self.keep_temp_dir = True

    def get_abspath_symlinks(self, rootpath):
        '''
        Checks to make sure project does not contain any absolute symlinks
        '''
        abspathsymlinks = []
        for root, dirs, files in os.walk(rootpath):
            for file in files:
                p = join(root,file)
                if islink(p) and isabs(os.readlink(p)):
                    abspathsymlinks.append(p)
        return abspathsymlinks

    def make_param(self, args):
        '''
        run with --param so that param.txt can be modified and then run
        extract --outdir from args
        '''
        # outdir is argument after --outdir
        outdir = args[args.index('--outdir') + 1]
        # Create param.txt
        o,e,r = common.run_path_discov(args + ['--param'])
        paramtxt = join(outdir, 'input', 'param.txt')
        return paramtxt

    def run_with_rikkcdna(self, args):
        ''' run --param and change param.txt to rikkcdna to run faster '''
        # Generate param.txt base project
        paramtxt = self.make_param(args)
        param = None
        for line in fileinput.input(paramtxt, inplace=True):
            nocommentline = line.split('#')[0].rstrip()
            if 'blast_db_list' in line:
                sys.stdout.write('blast_db_list {0},{0}\n'.format(RIKKNT))
            elif 'blast_task_list' in line:
                sys.stdout.write('blast_task_list megablast,dc-megablast\n')
            elif 'blast_options_list' in line:
                sys.stdout.write('blast_options_list -evalue 1e-4 -word_size 28,-evalue 1e-4 -word_size 12\n')
            elif 'ninst_list' in line:
                ninst = nocommentline.replace('ninst_list','').strip().split(',')[0]
                sys.stdout.write('ninst_list {0},{0}\n'.format(ninst))
            else:
                sys.stdout.write(line)
        # Now run with --noparam to use modified param.txt
        return common.run_path_discov(args + ['--noparam'])

    def edit_hostmap(self, paramtxt, mappers, dbs, names, options):
        for line in fileinput.input(paramtxt, inplace=True):
            nocommentline = line.split('#')[0].rstrip()
            p = nocommentline.split(None, 1)
            opname = ''
            if len(p) == 2:
                opname,opvalue = p
            if opname == 'mapper_program_list':
                sys.stdout.write('mapper_program_list {0}\n'.format(mappers))
            elif opname == 'mapper_db_list':
                sys.stdout.write('mapper_db_list {0}\n'.format(dbs))
            elif opname == 'mapper_name_list':
                sys.stdout.write('mapper_name_list {0}\n'.format(names))
            elif opname == 'mapper_options_list':
                sys.stdout.write('mapper_options_list {0}\n'.format(options))
            else:
                sys.stdout.write(nocommentline + '\n')

    def run_with_snap(self, args):
        paramtxt = self.make_param(args)
        self.edit_hostmap(
            paramtxt,
            'snap,snap',
            ','.join([SNAPDNA,SNAPRNA]),
            'snap_dna,snap_rna',
            '-t 1,-t 1', # Use 1 thread per run
        )
        return common.run_path_discov(args + ['--noparam'])

    def run_with_diamond(self, args):
        paramtxt = self.make_param(args)
        # Add diamond to db list
        for line in fileinput.input(paramtxt, inplace=True):
            nocommentline = line.split('#')[0].rstrip()
            splitlist = nocommentline.split(',')
            if 'blast_db_list' in line:
                sys.stdout.write('blast_db_list {0},{1}\n'.format(
                    RIKKNT, DIAMOND
                ))
            elif 'blast_task_list' in line:
                sys.stdout.write('blast_task_list megablast,diamond\n')
            elif line.startswith('blast_options_list'):
                tmpdir = join(self.outdir,'temp')
                tmpdir = ''
                sys.stdout.write(
                    splitlist[0] + ',--compress 0 -p 0 -v -k 10 -w 28 ' \
                    '--id 0.7 -c 4 -t {0}\n'.format(tmpdir)
                )
            elif 'ninst_list' in line:
                ninst = nocommentline.replace('ninst_list','').strip().split(',')[0]
                sys.stdout.write('ninst_list {0},{0}\n'.format(ninst))
                sys.stdout.write('blast_pro_db {0}\n'.format(RIKKNR))
            else:
                sys.stdout.write(line)
        # Now run with --noparam to use modified param.txt
        return common.run_path_discov(args + ['--noparam'])

    def test_diamond_db(self):
        self.outdir = 'diamond_test'
        f = self.f_fastq
        r = self.r_fastq
        args = ['-R1', f, '-R2', r, '--outdir', self.outdir]
        o,e,r = self.run_with_diamond(args)
        missingfiles = common.verify_project(
            self.outdir, common.PROJECT_FILES, r1r2='R1R2'
        )
        common.print_list(missingfiles)
        print o
        print e
        self.assertEqual([], missingfiles, 'Required files missing from project')
        self.assertEqual(0, r, 'Return code was not 0')
        self.assertEqual([], self.get_abspath_symlinks(self.outdir))

    def test_r1only_abspath(self):
        # relative path outdir
        self.outdir = 'r1_abspath_outdir_relpath'
        # abolute path read
        f = self.f_fastq
        args = ['-R1', f, '--outdir', self.outdir]
        o,e,r = self.run_with_rikkcdna(args)
        missingfiles = common.verify_project(
            self.outdir, common.PROJECT_FILES, r1r2='R1'
        )
        common.print_list(missingfiles)
        print o
        print e
        self.assertEqual([], missingfiles, 'Required files missing from project')
        self.assertEqual(0, r, 'Return code was not 0')
        self.assertEqual([], self.get_abspath_symlinks(self.outdir))

    def test_r1r2_outdir_abspath(self):
        # abspath outdir
        self.outdir = join(self.testdir, 'r1r2_relpath_outdir_abspath')
        # relative path r1r2
        f = relpath(self.f_fastq)
        r = relpath(self.r_fastq)
        args = ['-R1', f, '-R2', r, '--outdir', self.outdir]
        o,e,r = self.run_with_rikkcdna(args)
        missingfiles = common.verify_project(
            self.outdir, common.PROJECT_FILES, r1r2='R1R2'
        )
        common.print_list(missingfiles)
        print o
        print e
        self.assertEqual([], missingfiles, 'Required files missing from project')
        self.assertEqual(0, r, 'Return code was not 0')
        self.assertEqual([], self.get_abspath_symlinks(self.outdir))

    @attr('current')
    def test_r1gzip(self):
        self.outdir = join(self.testdir, 'r1gzip')
        args = ['-R1', self.f_fastq_gz, '--outdir', self.outdir]
        o,e,r = self.run_with_rikkcdna(args)
        missingfiles = common.verify_project(
            self.outdir, common.PROJECT_FILES, r1r2='R1'
        )
        common.print_list(missingfiles)
        print o
        print e
        self.assertEqual([], missingfiles, 'Required files missing from project')
        self.assertEqual(0, r, 'Return code was not 0')
        self.assertEqual([], self.get_abspath_symlinks(self.outdir))

    def test_r1sff(self):
        skip = [
            'results/quality_analysis/F_fastqc.zip',
            'results/quality_analysis/F_fastqc.html',
            'results/orf_filter/R1.unmap.fastq', # There will be no unmapped reads so these will not be generated
            'results/orf_filter/R1.orfout.fa'    # so we will just skip their checks
        ]
        self.outdir = join(self.testdir, 'r1sff')
        args = ['-R1', self.f_sff, '--outdir', self.outdir]
        o,e,r = self.run_with_rikkcdna(args)
        missingfiles = common.verify_project(
            self.outdir, common.PROJECT_FILES, r1r2='R1', skiplist=skip
        )
        common.print_list(missingfiles)
        print o
        print e
        self.assertEqual([], missingfiles, 'Required files missing from project')
        self.assertEqual(0, r, 'Return code was not 0')
        self.assertEqual([], self.get_abspath_symlinks(self.outdir))

    def test_snap_hostmap(self):
        self.outdir = 'snap_test'
        f = self.f_fastq
        r = self.r_fastq
        args = ['-R1', f, '-R2', r, '--outdir', self.outdir]
        o,e,r = self.run_with_snap(args)
        missingfiles = common.verify_project(
            self.outdir, common.PROJECT_FILES, r1r2='R1R2'
        )
        common.print_list(missingfiles)
        print o
        print e
        self.assertEqual([], missingfiles, 'Required files missing from project')
        self.assertEqual(0, r, 'Return code was not 0')
        self.assertEqual([], self.get_abspath_symlinks(self.outdir))
