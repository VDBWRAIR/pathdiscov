from os.path import *

import unittest2 as unittest
from nose.plugins.attrib import attr

import common

class TestR1Only(common.TempDir):
    def setUp(self):
        super(TestR1Only, self).setUp()
        self.f_fastq = join(common.TESTDATA, 'F.fastq')
        self.f_fastq_gz = join(common.TESTDATA, 'F.fastq.gz')
        self.f_sff = join(common.TESTDATA, '454Reads.sff')
        self.f_sff_gz = join(common.TESTDATA, '454Reads.sff.gz')

    @attr('current')
    def test_run_r1only_abspath(self):
        self.outdir = join(self.testdir, 'r1only_abspath')
        f = self.f_fastq
        args = ['-R1', f, '--outdir', self.outdir]
        r = common.run_path_discov(args)
        missingfiles = common.verify_project(self.outdir, common.PROJECT_FILES, False)
        common.print_list(missingfiles)
        print r
        self.assertEqual([], missingfiles)

    def test_run_r1only_relpath(self):
        self.outdir = join(self.testdir, 'r1only_relpath')
        f = relpath(self.f_fastq)
        args = ['-R1', f, '--outdir', self.outdir]
        r = common.run_path_discov(args)
        missingfiles = common.verify_project(self.outdir, common.PROJECT_FILES, False)
        common.print_list(missingfiles)
        print r
        self.assertEqual([], missingfiles)

    def test_run_r1only_sff(self):
        self.outdir = join(self.testdir, 'r1only_sff')
        f = self.f_sff
        args = ['-R1', f, '--outdir', self.outdir]
        r = common.run_path_discov(args)
        missingfiles = common.verify_project(self.outdir, common.PROJECT_FILES, False)
        common.print_list(missingfiles)
        print r
        self.assertEqual([], missingfiles)

    def test_run_r1only_fastq_gzip(self):
        self.outdir = join(self.testdir, 'r1only_fastq_gzip')
        f = self.f_fastq_gz
        args = ['-R1', f, '--outdir', self.outdir]
        r = common.run_path_discov(args)
        missingfiles = common.verify_project(self.outdir, common.PROJECT_FILES, False)
        common.print_list(missingfiles)
        print r
        self.assertEqual([], missingfiles)
