from common import *
from mock import mock_open as _mock_open

from pathdiscov.linecount import countlines

class Base(BaseTester):
    modulepath = 'pathdiscov.sff2fastq'

@patch('pathdiscov.sff2fastq.argparse')
class TestMain(Base, BaseTempDir):
    functionname = 'main'

    def setUp(self):
        super(TestMain,self).setUp()
        self.parse_args = Mock()
        self.parser = Mock()
        self.parser.return_value.parse_args.return_value = self.parse_args
        self.in_sffs = [self.sff_file, self.sffgz_file]

    def test_accepts_multiple_sff(self, mock_argparse):
        mock_argparse.ArgumentParser = self.parser
        self.parse_args.sfffile = self.in_sffs
        self.parse_args.outfile = 'output.fastq'

        self._C()

        numseq = countlines('output.fastq','fastq')
        eq_(20,numseq)

class TestSffToFastq(Base):
    functionname = 'sff_to_fastq'

    def setUp(self):
        super(TestSffToFastq,self).setUp()
        self.outfh = MagicMock(file)()
        self.in_sffs = [self.sff_file, self.sffgz_file]

    def test_empty_sff_list(self):
        r = self._C([], self.outfh)
        eq_(0, r)
        # Should not have written at all
        eq_(0,self.outfh.write.call_count)

    def test_multiple_sff(self):
        r = self._C(self.in_sffs, self.outfh)
        eq_(20,self.outfh.write.call_count)
        eq_(20,r)

import unittest2 as unittest
import mock

from .. import sff2fastq

class TestFileHandle(unittest.TestCase):
    def test_opens_gzip(self):
        with mock.patch.object(sff2fastq, 'gzip') as mgzip:
            r = sff2fastq.file_handle('/path/to/foo.sff.gz', 'rb')
            self.assertEqual((mgzip.open.return_value, 'sff'), r)
            mgzip.open.assert_called_once_with('/path/to/foo.sff.gz', 'rb')

    def test_opens_normal(self):
        with mock.patch('__builtin__.open') as mopen:
            r = sff2fastq.file_handle('/path/to/foo.sff', 'rb')
            self.assertEqual((mopen.return_value, 'sff'), r)
            mopen.assert_called_once_with('/path/to/foo.sff', 'rb')
