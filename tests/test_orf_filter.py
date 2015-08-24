import os
from os.path import *
import sys
import shutil
import inspect

import unittest2 as unittest
import sh

from common import TESTDATA, SCRATCH, PATHDISCOV, RIKKDB, TESTDIR
import common

sh.orf_filter = sh.Command('orf_filter.pl')

# Example contig fasta input file
fasta_input = join(TESTDIR, 'out.cap.fa')

# Param
paramtxt = '''
command orf_filter
getorf_options -minsize 60 -find 0
'''

expect_count_contig = [
    ('input', '86'),
    ('orf_filter', '86'),
]

expect_count_r1 = [
    ('input', '250'),
    ('orf_filter', '232'),
]

expect_count_r2 = [
    ('input', '250'),
    ('orf_filter', '230'),
]

class TestOrfFilter(common.StageTestBase):
    def test_can_supply_only_r2(self):
        param = self._write_param(paramtxt)
        outdir = inspect.stack()[0][3]
        logs = join(outdir,'logs')
        os.makedirs(logs)
        sh.orf_filter(
            sample='testsample', paramfile=param,
            outputdir=outdir, logs=logs, timestamp='0',
            R2=self.r2_fq
        )
        self.assertTrue(exists(join(outdir,'orf_filter.R2')))
        self.assertTrue(exists(join(outdir,'R2.orfout.fa')))
        self._verify_countfile(expect_count_r2, join(outdir, 'R2.count'))

    def test_orf_filter_fastq_r1r2(self):
        param = self._write_param(paramtxt)
        outdir = inspect.stack()[0][3]
        logs = join(outdir,'logs')
        os.makedirs(logs)
        sh.orf_filter(
            sample='testsample', paramfile=param,
            outputdir=outdir, logs=logs, timestamp='0',
            R1=self.r1_fq, R2=self.r2_fq
        )

        self.assertNotEmpty(join(outdir,'orf_filter.R1'))
        self.assertNotEmpty(join(outdir,'orf_filter.R2'))
        self.assertTrue(exists(join(outdir,'R1.orfout.fa')))
        self.assertTrue(exists(join(outdir,'R2.orfout.fa')))
        self._verify_countfile(expect_count_r1, join(outdir, 'R1.count'))
        self._verify_countfile(expect_count_r2, join(outdir, 'R2.count'))

    def test_orf_filter_fasta_contig(self):
        param = self._write_param(paramtxt)
        outdir = inspect.stack()[0][3]
        logs = join(outdir,'logs')
        os.makedirs(logs)
        sh.orf_filter(
            sample='testsample', paramfile=param,
            outputdir=outdir, logs=logs, timestamp='0',
            R1=fasta_input, contig='1'
        )

        self.assertTrue(exists(join(outdir,'orf_filter.contig')))
        self.assertTrue(exists(join(outdir,'contig.orfout.fa')))
        self._verify_countfile(expect_count_contig, join(outdir, 'contig.count'))
