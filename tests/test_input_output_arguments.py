from os.path import *

import unittest2 as unittest
from nose.plugins.attrib import attr

import common

@attr('slow')
class TestInputOutputArguments(common.TempDir):
    def setUp(self):
        super(TestInputOutputArguments, self).setUp()
        self.f_fastq = join(common.TESTDATA, 'F.fastq')
        self.f_fastq_gz = join(common.TESTDATA, 'F.fastq.gz')
        self.r_fastq = join(common.TESTDATA, 'R.fastq')
        self.r_fastq_gz = join(common.TESTDATA, 'R.fastq.gz')
        self.f_sff = join(common.TESTDATA, '454Reads.sff')
        self.f_sff_gz = join(common.TESTDATA, '454Reads.sff.gz')

    def test_r1only_abspath(self):
        # relative path outdir
        self.outdir = 'r1_abspath_outdir_relpath'
        # abolute path read
        f = self.f_fastq
        args = ['-R1', f, '--outdir', self.outdir]
        o,e,r = common.run_path_discov(args)
        missingfiles = common.verify_project(
            self.outdir, common.PROJECT_FILES, r1r2='R1'
        )
        common.print_list(missingfiles)
        print r
        self.assertEqual([], missingfiles, 'Required files missing from project')
        self.assertEqual(0, r, 'Return code was not 0')

    def test_r1r2_outdir_abspath(self):
        # abspath outdir
        self.outdir = join(self.testdir, 'r1r2_relpath_outdir_abspath')
        # relative path r1r2
        f = relpath(self.f_fastq)
        r = relpath(self.r_fastq)
        args = ['-R1', f, '-R2', r, '--outdir', self.outdir]
        o,e,r = common.run_path_discov(args)
        missingfiles = common.verify_project(
            self.outdir, common.PROJECT_FILES, r1r2='R1R2'
        )
        common.print_list(missingfiles)
        print r
        self.assertEqual([], missingfiles, 'Required files missing from project')
        self.assertEqual(0, r, 'Return code was not 0')

    def test_r1gzip(self):
        self.outdir = join(self.testdir, 'r1gzip')
        args = ['-R1', self.f_fastq_gz, '--outdir', self.outdir]
        o,e,r = common.run_path_discov(args)
        missingfiles = common.verify_project(
            self.outdir, common.PROJECT_FILES, r1r2='R1'
        )
        common.print_list(missingfiles)
        print r
        self.assertEqual([], missingfiles, 'Required files missing from project')
        self.assertEqual(0, r, 'Return code was not 0')

    def test_r1sff(self):
        skip = [
            'results/quality_analysis/F_fastqc.zip',
            'results/quality_analysis/F_fastqc.zip',
        ]
        self.outdir = join(self.testdir, 'r1sff')
        args = ['-R1', self.f_sff, '--outdir', self.outdir]
        o,e,r = common.run_path_discov(args)
        missingfiles = common.verify_project(
            self.outdir, common.PROJECT_FILES, r1r2='R1', skiplist=skip
        )
        common.print_list(missingfiles)
        print r
        self.assertEqual([], missingfiles, 'Required files missing from project')
        self.assertEqual(0, r, 'Return code was not 0')
