from common import *
from mock import mock_open as _mock_open


class Base(BaseTester):
    modulepath = 'pathdiscov.linecount'

    def setUp(self):
        super(Base,self).setUp()
        # 3 Lines with content 1,2,3
        self.mock_other_handle = mock_open(
            read_data = '\n'.join([str(i) for i in range(1,4)]) + '\n'
        )
        self.fasta_fh = self.mock_fasta_handle
        self.fastq_fh = self.mock_fastq_handle
        self.other_fh = self.mock_other_handle

class TestCountLines(Base):
    functionname = 'countlines'

    def test_input_is_string(self):
        with patch('__builtin__.open',self.mock_other_handle):
            r = self._C('/count/this/file',None)
            eq_(3,r)

    def test_counts_fasta(self):
        r = self._C(self.fasta_fh,'fasta') 
        eq_(3,r)

    def test_counts_fastq(self):
        r = self._C(self.fastq_fh,'fastq')
        eq_(3,r)

    def test_counts_other(self):
        r = self._C(self.other_fh,None)
        eq_(3,r)

    def test_counts_zero_filesize(self):
        mock_fh = mock_open(read_data='')
        r = self._C(mock_fh,None)
        eq_(0,r)

    @raises(ValueError)
    def test_unknown_fileformat_raises_exception(self):
        r = self._C(self.fasta_fh,'Unknown')

class TestLineCount(Base):
    functionname = 'linecount'

    def test_counts_fasta(self):
        outfh = MagicMock(file)()
        self._C(self.fasta_fh, 'test', outfh, 'fasta', False)
        outfh.write.assert_called_once_with('test\t3\n')

    def test_counts_fastq(self):
        outfh = MagicMock(file)()
        self._C(self.fastq_fh, 'test', outfh, 'fastq', False)
        outfh.write.assert_called_once_with('test\t3\n')

    def test_counts_other(self):
        outfh = MagicMock(file)()
        self._C(self.other_fh, 'test', outfh, None, False)
        outfh.write.assert_called_once_with('test\t3\n')

    def test_output_is_concatted_to(self):
        outfh = MagicMock(file)()
        self._C(self.fasta_fh, 'test', outfh, 'fasta', True)
        outfh.write.assert_called_once_with('test\t3\n')

    def test_input_is_string(self):
        outfh = MagicMock(file)()
        with patch('__builtin__.open', _mock_open()) as mck_open:
            self._C('/path/to/file.txt', 'test', outfh, None, True)
            mck_open.assert_called_once_with('/path/to/file.txt')

    def test_output_is_string_no_concat(self):
        m = _mock_open()
        with patch('__builtin__.open', m, create=True) as mck_open:
            self._C(self.fasta_fh, 'test', '/path/to/out.txt', 'fasta', False)
            mck_open.assert_called_once_with('/path/to/out.txt', 'w')
            mck_open.return_value.write.assert_called_once_with('test\t3\n')

    def test_output_is_string_no_concat_arg_str(self):
        m = _mock_open()
        with patch('__builtin__.open', m, create=True) as mck_open:
            self._C(self.fasta_fh, 'test', '/path/to/out.txt', 'fasta', '0')
            mck_open.assert_called_once_with('/path/to/out.txt', 'w')
            mck_open.return_value.write.assert_called_once_with('test\t3\n')

    def test_output_is_string_concat(self):
        m = _mock_open()
        with patch('__builtin__.open', m, create=True) as mck_open:
            self._C(self.fasta_fh, 'test', '/path/to/out.txt', 'fasta', True)
            mck_open.assert_called_once_with('/path/to/out.txt', 'a')
            mck_open.return_value.write.assert_called_once_with('test\t3\n')

    def test_concat_arg_can_be_str_int_bool(self):
        outfh = MagicMock(file)()
        self._C(self.fasta_fh, 'test', outfh, 'fasta', '1')
        self._C(self.fasta_fh, 'test', outfh, 'fasta', 1)
        self._C(self.fasta_fh, 'test', outfh, 'fasta', True)

    def test_format_arg_oldstyle_fasta(self):
        outfh = MagicMock(file)()
        self._C(self.fastq_fh, 'test', outfh, '1', False)
        outfh.write.assert_called_once_with('test\t3\n')

    def test_format_arg_oldstyle_fastq(self):
        outfh = MagicMock(file)()
        self._C(self.fasta_fh, 'test', outfh, '2', False)
        outfh.write.assert_called_once_with('test\t3\n')

    def test_format_arg_oldstyle_other(self):
        outfh = MagicMock(file)()
        self._C(self.other_fh, 'test', outfh, '0', False)
        outfh.write.assert_called_once_with('test\t3\n')
