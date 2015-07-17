import os
from os.path import *
import sys
import shutil

import unittest2 as unittest
import sh

from common import TESTDATA, SCRATCH, PATHDISCOV
import common

# Allows us to reference sh.quality_filter since quality_filter.pl
# cannot just work with sh import because of .pl
sh.quality_filter = sh.Command(join(
    PATHDISCOV, 'quality_filter', 'quality_filter.pl'
))

# Param for quality_filter
paramtxt = '''
command quality_filter
cutadapt_options_R1			-a GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -g GCCGGAGCTCTGCAGATATC -a GATATCTGCAGAGCTCCGGC -m 50 --match-read-wildcards
cutadapt_options_R2			-a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -g GCCGGAGCTCTGCAGATATC -a GATATCTGCAGAGCTCCGGC -m 50 --match-read-wildcards
prinseq_options				-min_len 50 -derep 14 -lc_method dust -lc_threshold 3 -trim_ns_left 1 -trim_ns_right 1 -trim_qual_right 15
'''

expect_count_r1 = [
    ('input', '250'),
    ('cut_adapt', '236'),
    ('prinseq', '159')
]

expect_count_r2 = [
    ('input', '250'),
    ('cut_adapt', '232'),
    ('prinseq', '158')
]

class TestQualityFilter(common.StageTestBase):
    def test_runs_on_r1(self):
        param = self._write_param(paramtxt)
        outdir = 'r1'
        logs = join(outdir,'logs')
        os.makedirs(logs)
        sh.quality_filter(
            sample="testsample", paramfile=param, outputdir=outdir,
            logs=logs, timestamp='0',
            R1=self.r1_fq
        )
        self.assertTrue(join(outdir,'quality_filter.R1'))
        self._verify_countfile(expect_count_r1, join(outdir, 'R1.count'))

    def test_runs_on_r1r2(self):
        param = self._write_param(paramtxt)
        outdir = 'r1r2'
        logs = join(outdir,'logs')
        os.makedirs(logs)
        sh.quality_filter(
            sample="testsample", paramfile=param, outputdir=outdir,
            logs=logs, timestamp='0',
            R1=self.r1_fq, R2=self.r2_fq
        )
        self.assertTrue(join(outdir,'quality_filter.R1'))
        self.assertTrue(join(outdir,'quality_filter.R2'))
        print "Check R1.count"
        self._verify_countfile(expect_count_r1, join(outdir, 'R1.count'))
        print "Check R2.count"
        self._verify_countfile(expect_count_r2, join(outdir, 'R2.count'))
