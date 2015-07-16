import os
from os.path import *
import sys
import shutil
import inspect

import unittest2 as unittest
import sh

from common import TESTDATA, SCRATCH, PATHDISCOV, RIKKDB
import common

# Allows us to reference sh.quality_filter since quality_filter.pl
# cannot just work with sh import because of .pl
sh.host_map = sh.Command(join(
    PATHDISCOV, 'host_map', 'host_map.pl'
))

dna_path = join(RIKKDB, 'rikkdna')
dna_path = join(RIKKDB, 'host_map')
rna_path = join(RIKKDB, 'rikkrna')

# Param
paramtxt = '''
command host_map
mapper_program_list		bowtie2,bowtie2
mapper_db_list			{0},{1}
mapper_name_list		bowtie2_1,bowtie2_2
mapper_options_list		--local -p 1,--local -p 1
'''.format(dna_path,rna_path)

expect_count_r1 = [
    ('input', '250'),
    ('bowtie2_1', '52'),
    ('bowtie2_2', '52'),
]

expect_count_r2 = [
    ('input', '250'),
    ('bowtie2_1', '53'),
    ('bowtie2_2', '53'),
]

class TestHostMap(common.StageTestBase):
    def test_removes_host(self):
        param = self._write_param(paramtxt)
        outdir = inspect.stack()[0][3]
        logs = join(outdir,'logs')
        os.makedirs(logs)
        sh.host_map(
            sample='testsample', paramfile=param,
            outputdir=outdir, logs=logs, timestamp='0',
            R1=self.r1_fq, R2=self.r2_fq
        )

        self.assertTrue(exists(join(outdir,'host_map_1.R1')))
        self.assertTrue(exists(join(outdir,'host_map_1.R2')))
        self._verify_countfile(expect_count_r1, join(outdir, 'R1.count'))
        self._verify_countfile(expect_count_r2, join(outdir, 'R2.count'))
