from os.path import *
import fileinput
import sys
import os

import unittest2 as unittest
from nose.plugins.attrib import attr
import mock

import common

@attr('slow')
class TestInputOutputArguments(common.TempDir):
    def setUp(self):
        super(TestInputOutputArguments, self).setUp()
        self.keep_temp_dir = True

    def run_with_rikkcdna(self, args):
        ''' run --param and change param.txt to rikkcdna to run faster '''
        # Some day will do something special, but just run it for now as
        # the param.txt is
        #return common.run_path_discov(args)
        # outdir is argument after --outdir
        outdir = args[args.index('--outdir') + 1]
        # Create param.txt
        o,e,r = common.run_path_discov(args + ['--param'])
        paramtxt = join(outdir, 'input', 'param.txt')
        # modify param.txt to use rikkcdna
        rikkcdnapath = join(common.TESTDIR, 'rikkcdna','rikkcdna')
        param = None
        for line in fileinput.input(paramtxt, inplace=True):
            if 'blast_db_list' in line:
                sys.stdout.write('blast_db_list {0},{0}\n'.format(rikkcdnapath))
            else:
                sys.stdout.write(line)
        # Now run with --noparam to use modified param.txt
        return common.run_path_discov(args + ['--noparam'])

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

    def test_r1sff(self):
        skip = [
            'results/quality_analysis/F_fastqc.zip',
            'results/quality_analysis/F_fastqc.zip',
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
