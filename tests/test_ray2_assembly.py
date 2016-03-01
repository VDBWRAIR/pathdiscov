import os
from os.path import *
import sys
import shutil
import inspect
import re

import unittest2 as unittest
import sh

from common import TESTDATA, SCRATCH, PATHDISCOV, RIKKDB
import common

# Allows us to reference sh.quality_filter since quality_filter.pl
# cannot just work with sh import because of .pl
sh.ray2_assembly = sh.Command('ray2_assembly.pl')

# Param
paramtxt = '''
command ray2_assembly
kmer					25
ninst   				1
cap   					1
map2contigs				yes
bowtie2_options			--local -p 1
parse_contigs           yes
'''

expect_count_assembly = [
    ('ray_contigs', '89'),
    ('cap_contigs', '89'),
]

expect_count_r1 = [
    ('input', '250'),
    ('unassembled_reads', '162'),
]

expect_count_r2 = [
    ('input', '250'),
    ('unassembled_reads', '162'),
]

r1unmap_sym = '1.R1.unmap.fastq'
r2unmap_sym = '1.R2.unmap.fastq'
r1unmap_fq = 'bowtie2_mapping/R1.unmap.fastq'
r2unmap_fq = 'bowtie2_mapping/R2.unmap.fastq'
outray_fa = 'out.ray.fa'
outcap_fa = 'out.cap.fa'
ray2assembly_fa = 'ray2_assembly_1.fasta'
readsbycontig = 'reads_by_contig'

class TestRay2Assembly(common.StageTestBase):
    def test_mpiexec_missing_is_resolved(self):
        mpiexecpath = self._get_mpiexec_path()
        # Ensure mpiexec not in path so we can ensure that 
        # stage fixes that issue if it exists
        if mpiexecpath:
            os.environ['PATH'] = os.environ['PATH'].replace(mpiexecpath+':', '')
        param = self._write_param(paramtxt)
        outdir = inspect.stack()[0][3]
        logs = join(outdir,'logs')
        os.makedirs(logs)
        sh.ray2_assembly(
            sample='testsample', paramfile=param,
            outputdir=outdir, logs=logs, timestamp='0',
            R1=self.r1_fq, R2=self.r2_fq
        )
        self.assertTrue(exists(join(outdir,r1unmap_sym)))
        self.assertTrue(exists(join(outdir,r2unmap_sym)))
        self.assertTrue(exists(join(outdir,outray_fa)))
        self.assertTrue(exists(join(outdir,outcap_fa)))
        self.assertTrue(exists(join(outdir,ray2assembly_fa)))
        self.assertTrue(exists(join(outdir,readsbycontig)))
        self.assertTrue(89, len(os.listdir(join(outdir,readsbycontig))))
        self._verify_countfile(expect_count_r1, join(outdir,'R1.count'))
        self._verify_countfile(expect_count_r1, join(outdir,'R2.count'))
        self._verify_countfile(expect_count_assembly, join(outdir,'assembly.count'))

    def _get_mpiexec_path(self):
        out = sh.which('mpiexec')
        mpiexec_dir = ''
        # mpiexec was found
        if out is not None:
            mpiexec_dir = dirname(out)
        return mpiexec_dir

    def test_creates_contigs_unmapped_parse_contigs(self):
        param = self._write_param(paramtxt)
        outdir = inspect.stack()[0][3]
        logs = join(outdir,'logs')
        os.makedirs(logs)
        sh.ray2_assembly(
            sample='testsample', paramfile=param,
            outputdir=outdir, logs=logs, timestamp='0',
            R1=self.r1_fq, R2=self.r2_fq
        )
        self.assertTrue(exists(join(outdir,r1unmap_sym)))
        self.assertTrue(exists(join(outdir,r2unmap_sym)))
        self.assertTrue(exists(join(outdir,outray_fa)))
        self.assertTrue(exists(join(outdir,outcap_fa)))
        self.assertTrue(exists(join(outdir,ray2assembly_fa)))
        self.assertTrue(exists(join(outdir,readsbycontig)))
        self.assertTrue(89, len(os.listdir(join(outdir,readsbycontig))))
        self._verify_countfile(expect_count_r1, join(outdir,'R1.count'))
        self._verify_countfile(expect_count_r1, join(outdir,'R2.count'))
        self._verify_countfile(expect_count_assembly, join(outdir,'assembly.count'))

    def test_no_cap(self):
        txt = re.sub('cap\s+1', 'cap 0', paramtxt)
        param = self._write_param(txt, 'nocapparam.txt')
        outdir = inspect.stack()[0][3]
        logs = join(outdir,'logs')
        os.makedirs(logs)
        sh.ray2_assembly(
            sample='testsample', paramfile=param,
            outputdir=outdir, logs=logs, timestamp='0',
            R1=self.r1_fq, R2=self.r2_fq
        )
        self.assertTrue(exists(join(outdir,r1unmap_sym)))
        self.assertTrue(exists(join(outdir,r2unmap_sym)))
        self.assertTrue(exists(join(outdir,outray_fa)))
        self.assertFalse(exists(join(outdir,outcap_fa)))
        self.assertTrue(exists(join(outdir,readsbycontig)))
        self.assertTrue(exists(join(outdir,ray2assembly_fa)))
        self.assertTrue(89, len(os.listdir(join(outdir,readsbycontig))))
        self._verify_countfile(expect_count_r1, join(outdir,'R1.count'))
        self._verify_countfile(expect_count_r2, join(outdir,'R2.count'))
        eca = [(name,count) for name,count in expect_count_assembly if name != 'cap_contigs']
        self._verify_countfile(expect_count_assembly, join(outdir,'assembly.count'))

    def test_no_map2contig(self):
        txt = re.sub('map2contigs\s+yes', 'map2contigs no', paramtxt)
        param = self._write_param(txt, 'nomap2contigparam.txt')
        outdir = inspect.stack()[0][3]
        logs = join(outdir,'logs')
        os.makedirs(logs)
        sh.ray2_assembly(
            sample='testsample', paramfile=param,
            outputdir=outdir, logs=logs, timestamp='0',
            R1=self.r1_fq, R2=self.r2_fq
        )
        self.assertFalse(exists(join(outdir,r1unmap_sym)))
        self.assertFalse(exists(join(outdir,r2unmap_sym)))
        self.assertTrue(exists(join(outdir,outray_fa)))
        self.assertTrue(exists(join(outdir,outcap_fa)))
        self.assertFalse(exists(join(outdir,readsbycontig)))
        self.assertTrue(exists(join(outdir,ray2assembly_fa)))
        self._verify_countfile(expect_count_r1, join(outdir,'R1.count'))
        self._verify_countfile(expect_count_r1, join(outdir,'R2.count'))
        eca = [(name,count) for name,count in expect_count_r1 if name != 'unassembled_reads']
        eca = [(name,count) for name,count in expect_count_r2 if name != 'unassembled_reads']
        self._verify_countfile(expect_count_r1, join(outdir,'R1.count'))
        self._verify_countfile(expect_count_r2, join(outdir,'R2.count'))
        self._verify_countfile(expect_count_assembly, join(outdir,'assembly.count'))
