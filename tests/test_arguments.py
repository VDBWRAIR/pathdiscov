from os.path import *
import os
import re

import unittest2 as unittest

import common

class TestArguments(common.TempDir):
    def setUp(self):
        super(TestArguments,self).setUp()
        self.r1 = join(common.TESTDATA, 'F.fastq')

    def test_create_param(self):
        outdir = 'create_param'
        args = ['--param', '--outdir', outdir, '-R1', self.r1]
        o,e,r = common.run_path_discov(args)
        print o
        print
        print e
        self.assertTrue(
            exists(join(outdir, 'input', 'param.txt')),
            'input/param.txt was not created'
        )
        self.assertEqual(0, r, 'Return code was not 0')

    def test_cpu_option_sets_param_txt_ninst(self):
        outdir = 'cpu_param'
        args = ['--param', '--cpuNum', '99', '--outdir', outdir, '-R1', self.r1]
        o,e,r = common.run_path_discov(args)
        paramtxt = join(outdir, 'input', 'param.txt')
        param_contents = None
        fh = open(paramtxt)
        param_contents = fh.read()
        fh.close()

        ninst = re.findall('ninst(?:_list){0,1}\s+(?:(\d+)(?:,(\d+))*)', param_contents)
        for l,r in ninst:
            self.assertEqual('99', l, 'Did not set ninst inside param.txt')
            self.assertIn(r, ('99',''), 'Did not set ninst_list inside param.txt')

    def test_missing_outdir(self):
        o,e,r = common.run_path_discov(
            ['-R1', self.r1]
        )
        self.assertEqual(2, r, 'Did not ensure -R1 was specified')

    def test_missing_r1(self):
        o,e,r = common.run_path_discov(
            ['--outdir', 'outdir']
        )
        self.assertEqual(2, r, 'Did not ensure --outdir was specified')
