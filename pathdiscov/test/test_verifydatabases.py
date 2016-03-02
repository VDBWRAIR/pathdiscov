from os.path import *
import os
import sys
from glob import glob
import shutil

import unittest2 as unittest

from nose.tools import *
from nose.plugins.attrib import attr
from mock import *

from pathdiscov import verifydatabases

def configyaml():
    return {
        'databases': '/tmp/databases',
        'host_dna': 'humandna/hg38',
        'host_rna': 'humanrna/hg38_mrna',
        'nt_db': 'ncbi/blast/nt/nt',
        'tax_nodes': 'ncbi/taxonomy/nodes.dmp',
        'tax_names': 'ncbi/taxonomy/names.dmp',
    }

@patch('pathdiscov.verifydatabases.glob')
class TestDatabaseFilesHaveSameNumExtensions(unittest.TestCase):
    def test_has_no_files(self, mglob):
        mglob.glob.return_value = []
        self.assertRaises(
            ValueError,
            verifydatabases.all_db_files_same_prefix, 'db'
        )

    def test_missing_one_extension(self, mglob):
        mglob.glob.return_value = [
            'db.00.1', 'db.00.2', 'db.00.3',
            'db.01.1', 'db.01.3'
        ]
        r = verifydatabases.all_db_files_same_prefix('db')
        eq_(['.2'], r)

    def test_single_db(self, mglob):
        mglob.glob.return_value = [
            'db.00.1', 'db.00.2',
        ]
        r = verifydatabases.all_db_files_same_prefix('db')
        eq_([], r)

    def test_has_no_extensions(self, mglob):
        mglob.glob.return_value = [
            'db.00', 'db.01',
        ]
        r = verifydatabases.all_db_files_same_prefix('db')
        eq_([], r)

    def test_has_all_extensions(self, mglob):
        mglob.glob.return_value = [
            '/path/db.00.1', '/path/db.00.2', '/path/db.00.3',
            '/path/db.01.1', '/path/db.01.2', '/path/db.01.3'
        ]
        r = verifydatabases.all_db_files_same_prefix('/path/db')
        eq_([], r)

    def test_one_extension_for_one_file_different_same_count(self, mglob):
        mglob.glob.return_value = [
            'db.00.1', 'db.00.2', 'db.00.3',
            'db.01.1', 'db.01.3', 'db.01.4'
        ]
        r = verifydatabases.all_db_files_same_prefix('db')
        eq_(['.2', '.4'], r)

    def test_excludes_extensions(self, mglob):
        mglob.glob.return_value = [
            'db.00.1', 'db.00.2', 'db.00.3',
            'db.01.1', 'db.01.3', 'db.01.4',
            'db.00.fa', 'db.00.gz'
        ]
        r = verifydatabases.all_db_files_same_prefix('db', ['.fa', '.gz'])
        eq_(['.2', '.4'], r)

@attr('current')
class TestVerifyDatabases(object):
    def setUp(self):
        self.config = configyaml()
        self.nt = ['nt.00.1', 'nt.00.2', 'nt.01.1', 'nt.01.2', 'nt.nal']
        self.hdna = ['hdna.1.bt2', 'hdna.2.bt2', 'hdna.rev.1.bt2', 'hdna.rev.2.bt2']
        self.hrna = ['hrna.1.bt2', 'hrna.2.bt2', 'hrna.rev.1.bt2', 'hrna.rev.2.bt2']
        self.taxnames = ['names.dmp']
        self.taxnodes = ['names.dmp']
        self.dbs = [
            self.hdna, self.hrna, self.nt,
            self.taxnames, self.taxnodes
        ]

    @patch('pathdiscov.verifydatabases.glob')
    def test_contains_all_databases(self, mglob):
        self.hdna.append('hdna_all.fa')
        self.hdna.append('hdna.chromFa.tar.gz')
        mglob.glob.side_effect = self.dbs
        r = verifydatabases.verifydatabases(self.config)        
        assert_true(r)

    def missing_entire_db(self, mglob, i):
        self.dbs[i] = []
        mglob.glob.side_effect = self.dbs
        r = verifydatabases.verifydatabases(self.config)
        assert_false(r)
        ok_(r is not None)

    def test_database_does_not_exist(self):
        self.setUp()
        with patch('pathdiscov.verifydatabases.glob') as mglob:
            for i in range(len(self.dbs)):
                yield self.missing_entire_db, mglob, i

    @patch('pathdiscov.verifydatabases.glob')
    def test_nt_missing_file(self, mglob):
        self.missing_db_file(mglob, 2)

    def missing_db_file(self, mglob, i):
        del self.dbs[i][0]
        mglob.glob.side_effect = self.dbs
        r = verifydatabases.verifydatabases(self.config)        
        assert_false(r)
        ok_(r is not None)

    @patch('pathdiscov.verifydatabases.glob')
    def test_converts_environmental_vars(self, mglob):
        self.config['databases'] = '$HOME/databases'
        mglob.glob.side_effect = self.dbs
        r = verifydatabases.verifydatabases(self.config)        
        assert_true(r)

    @patch('pathdiscov.verifydatabases.glob')
    def test_converts_tilda_paths(self, mglob):
        self.config['databases'] = '~/databases'
        mglob.glob.side_effect = self.dbs
        r = verifydatabases.verifydatabases(self.config)        
        assert_true(r)
