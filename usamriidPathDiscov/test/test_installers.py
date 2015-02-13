import tempfile
import shutil
from os.path import *
import os
import sys

import unittest2 as unittest

import mock

TESTDIR = dirname(__file__)
# Path to downloads
DOWNLOADDIR = join(dirname(TESTDIR),'download')

from .. import installers

class TestInstallDiamond(unittest.TestCase):
    def setUp(self):
        self.tdir = tempfile.mkdtemp()
        os.chdir(self.tdir)
        self.addCleanup(shutil.rmtree, self.tdir)

        diamond_src_path = join(DOWNLOADDIR, 'diamond-linux64.tar.gz')
        
        self.options = mock.Mock()
        self.options.diamond.src = diamond_src_path
        # Should contain bin/diamond
        self.options.diamond.install_to = join(self.tdir, 'diamond_install_to')

        self.expected_bin_path = join(self.options.diamond.install_to,'bin','diamond')

    def assertExists(self, path, msg=None):
        if msg is None:
            msg = '{0} does not exist'.format(path)
        self.assertTrue(exists(path), msg)

    def assertExecutable(self, path, msg=None):
        if msg is None:
            msg = '{0} is not executable'.format(path)
        self.assertExists(path)
        self.assertTrue(os.access(path,os.X_OK),msg)

    def test_installs_diamond(self):
        r = installers.install_diamond(self.options)
        self.assertExecutable(self.expected_bin_path)

    def test_does_not_reinstall_diamond(self):
        self.expected_bin_path = join(self.options.diamond.install_to,'bin','diamond')
        r = installers.install_diamond(self.options)
        _stat = os.stat(self.expected_bin_path)
        r = installers.install_diamond(self.options)
        _stat2 = os.stat(self.expected_bin_path)
        self.assertEqual(_stat.st_atime, _stat2.st_atime)
        self.assertEqual(_stat.st_mtime, _stat2.st_mtime)
        self.assertEqual(_stat.st_ctime, _stat2.st_ctime)
        self.assertExecutable(self.expected_bin_path)
