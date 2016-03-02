import os
from os.path import *
import shutil
import sys

from nose.plugins.attrib import attr
import mock

import common

@attr('fast')
class TestEnvironSetup(common.TempDir):
    def setUp(self):
        super(TestEnvironSetup,self).setUp()
        self.mock_args = mock.Mock()
        self.mock_args.R1 = self.f_fastq
        self.mock_args.R2 = None
        self.mock_args.outdir = 'outdir'
        self.upd = join(dirname(common.TESTDIR), 'pathdiscov')

#    def test_sets_up_environment(self):
#        # Get globals
#        g = common.exec_main(self.mock_args)
#        # Ensure PATH contains correct locations
#        paths = [
#            self.upd,
#            join(self.upd, 'bin'),
#            join(self.upd, 'scripts'),
#            join(self.upd, 'step1'),
#            '/usr/lib64/openmpi/bin'
#        ]
#        result_paths = os.environ['PATH'].split(os.pathsep)
#        for p in paths:
#            self.assertIn(p, result_paths)
#
#        # Make sure openmpi lib is set in LDD_LIBRARY_PATH
#        result_paths = os.environ['LD_LIBRARY_PATH'].split(os.pathsep)
#        self.assertIn('/usr/lib64/openmpi/lib', result_paths)
#
#        # Make sure PERL5LIB is set
#        self.assertEqual(join(self.upd,'Local_Module'), os.environ['PERL5LIB'])
#
#        # Make sure R_LIBS are set
#        self.assertEqual(join(self.upd,'scripts'), os.environ['R_LIBS'])
#
#        # Make sure there are 9 INNO_ keys
#        self.assertEqual(9, sum([1 for k in os.environ if 'INNO' in k]))
