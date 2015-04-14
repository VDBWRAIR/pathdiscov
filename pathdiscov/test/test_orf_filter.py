import unittest2 as unittest

from os.path import *
import sys
import subprocess
import os
import tempfile
import shutil

import common

SCRIPTS = dirname(join(common.THIS))
ORFFILTER = join(SCRIPTS,'orf_filter')

class TestOrfFilter(unittest.TestCase):
    def setUp(self):
        self.tdir = tempfile.mkdtemp(suffix='orffiltertest')
        self.f = join(
            common.FIXTURES,
            'example_data',
            'F.fastq'
        )
        self.r = join(
            common.FIXTURES,
            'example_data',
            'R.fastq'
        )
        self.contig = join(
            common.FIXTURES,
            'orf_filter_input.contig'
        )
        self.sample = 'foo'
        self.paramfile = join(
            common.FIXTURES,
            'orf_filter.param.txt'
        ) 
        self.timestamp = 1
        self.outputdir = self.tdir
        self.logs = join(self.outputdir,'logs')
        self.r1outputs = [
            join(self.outputdir, 'R1.orfout.fa'),
            join(self.outputdir, 'orf_filter.R1'),
            join(self.outputdir, 'R1.count'),
        ]

    def tearDown(self):
        return
        shutil.rmtree(self.tdir)

    def run_orf_filter(self, **kwargs):
        os.makedirs(kwargs['logs'])
        script = join(ORFFILTER,'orf_filter.pl')
        cmd = script
        cmd += ' --sample {sample} ' \
                '--paramfile {paramfile} ' \
                '--outputdir {outputdir} ' \
                '--logs {logs} ' \
                '--timestamp {timestamp} ' \
                '--R1 {R1}'.format(
                    **kwargs
                )
        if 'R2' in kwargs:
            cmd += ' --R2 {0}'.format(kwargs['R2'])
        if 'contig' in kwargs:
            cmd += ' --contig {0}'.format(kwargs['contig'])
        p = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            shell=True
        )
        sout,serr = p.communicate()
        r = p.returncode
        return r,sout,serr

    def notEmpty(self, path):
        self.assertTrue(
            exists(path),
            '{0} does not exist'.format(path)
        )
        self.assertTrue(
            os.stat(path).st_size != 0,
            '{0} size is 0'.format(path)
        )

    def test_produces_contig_r1name(self):
        rc,sout,serr = self.run_orf_filter(
            sample=self.sample,
            paramfile=self.paramfile,
            outputdir=self.outputdir,
            logs=self.logs,
            timestamp=self.timestamp,
            R1=self.contig,
            contig=1
        )
        self.assertEqual(0,rc)
        for f in self.r1outputs:
            self.notEmpty(f.replace('R1','contig'))

    def test_runs_on_only_r1(self):
        rc,sout,serr = self.run_orf_filter(
            sample=self.sample,
            paramfile=self.paramfile,
            outputdir=self.outputdir,
            logs=self.logs,
            timestamp=self.timestamp,
            R1=self.f
        )
        self.assertEqual(0,rc)
        for f in self.r1outputs:
            self.notEmpty(f)

    def test_runs_on_r1r2(self):
        rc,sout,serr = self.run_orf_filter(
            sample=self.sample,
            paramfile=self.paramfile,
            outputdir=self.outputdir,
            logs=self.logs,
            timestamp=self.timestamp,
            R1=self.f,
            R2=self.r
        )
        self.assertEqual(0,rc)
        for f in self.r1outputs:
            self.notEmpty(f)
            self.notEmpty(f.replace('R1','R2'))
