import os
import sys
from os.path import *
import subprocess
import re
from glob import glob
import shutil
import tempfile
import subprocess

from nose.tools import eq_, ok_, raises, assert_raises
from nose.plugins.attrib import attr
from nose.plugins.skip import SkipTest
from mock import Mock, MagicMock, patch, call, mock_open as _mock_open

THIS = dirname(abspath(__file__))
FIXTURES = join(THIS,'fixtures')

if 'check_output' not in dir(subprocess):
    def check_output(*popenargs, **kwargs):
        r"""Run command with arguments and return its output as a byte string.
        Backported from Python 2.7 as it's implemented as pure python on stdlib.
        >>> check_output(['/usr/bin/python', '--version'])
        Python 2.6.2
        """
        process = subprocess.Popen(stdout=subprocess.PIPE, *popenargs, **kwargs)
        output, unused_err = process.communicate()
        retcode = process.poll()
        if retcode:
            cmd = kwargs.get("args")
            if cmd is None:
                cmd = popenargs[0]
            error = subprocess.CalledProcessError(retcode, cmd)
            error.output = output
            raise error
        return output 
    subprocess.check_output = check_output

class MockFile(file):
    '''
    I wish mock_open would just do this instead

    Essentially mock file object(open call) but also be iterable
    '''
    def __init__(self, read_data=None):
        self.pos = 0
        self.lines = []
        self.data = read_data
        if read_data:
            self.lines = self.data.splitlines(True)
        self.curlineindex = 0

    def __call__(self,filepath,mode='r'):
        return self

    def __iter__(self):
        return iter(self.lines)

    def read(self,bytes=None):
        self.pos = len(self.data)
        return '\n'.join(self.lines)

    def next(self):
        if self.curlineindex == len(self.lines):
            raise StopIteration
        self.curlineindex += 1
        line = self.lines[self.curlineindex-1]
        self.pos += len(line)
        return line

    def readline(self):
        try:
            return self.next()
        except StopIteration:
            return ''

    def readlines(self):
        self.pos = len(self.data)
        return [line.rstrip() for line in self.lines]

    def seek(self, pos):
        self.pos = pos
        if self.data:
            self.lines = self.data.splitlines(True)
        self.curlineindex = 0

    def tell(self):
        return self.pos

def mock_open(mock=None, read_data=None):
    '''
    Add iteration for read_data
    '''
    mo = MockFile(read_data)
    return mo

class Base(object):
    def setUp(self):
        self.sff_file = join(FIXTURES, 'reads.sff')
        self.sffgz_file = join(FIXTURES, 'reads.sff.gz')
        self.fasta_file = join(FIXTURES, 'reads.fasta')
        self.fastagz_file = join(FIXTURES, 'reads.fasta.gz')
        self.fastq_file = join(FIXTURES, 'reads.fastq')
        self.fastqgz_file = join(FIXTURES, 'reads.fastq.gz')
        # Mock fasta file with one sequence having multiple lines
        self.mock_fasta_handle = mock_open(
            read_data = '\n'.join([
                '>id1',
                'ATGC',
                '>id2',
                'ATGC',
                'CGTA',
                '>id3',
                'ATGCA'
            ])
        )
        # Mock fastq file
        self.mock_fastq_handle = mock_open(
            read_data = '\n'.join([
                '@id1',
                'ATGC',
                '+',
                '!!!!',
                '@id2',
                'ATGCCGTA',
                '+',
                '!!!!!!!!',
                '@id3',
                'ATGCA',
                '+',
                '!!!!!'
            ])
        )

    def tearDown(self):
        pass

class BaseTester(Base):
    def setUp(self):
        super(BaseTester,self).setUp()

    def _C( self, *args, **kwargs):
        '''
        Set modulepath as instance variable to be the name of the module
        Set functionname as instance variable to automagically run that function
        with self._C
        '''
        m = __import__( self.modulepath, fromlist=[self.functionname] )
        return getattr(m,self.functionname)( *args, **kwargs ) 

class BaseTempDir(Base):
    def setUp( self ):
        super(BaseTempDir,self).setUp()
        self.tdir = tempfile.mkdtemp( suffix='riidtest' )
        os.chdir( self.tdir )

    def tearDown( self ):
        os.chdir( '/tmp' )
        shutil.rmtree( self.tdir )
