from common import *
from mock import mock_open as _mock_open

from usamriidPathDiscov.linecount import countlines

class Base(BaseTester):
    modulepath = 'usamriidPathDiscov.sff2fastq'

@patch('usamriidPathDiscov.sff2fastq.argparse')
class TestMain(Base, BaseTempDir):
    functionname = 'main'

    def setUp(self):
        super(TestMain,self).setUp()
        self.parse_args = Mock()
        self.parser = Mock()
        self.parser.return_value.parse_args.return_value = self.parse_args
        self.in_sffs = [self.sff_file]

    def test_accepts_multiple_sff(self, mock_argparse):
        mock_argparse.ArgumentParser = self.parser
        self.parse_args.sfffile = self.in_sffs*10
        self.parse_args.outfile = 'output.fastq'

        self._C()

        numseq = countlines('output.fastq','fastq')
        eq_(100,numseq)

    def test_accepts_multiple_single_sff(self, mock_argparse):
        mock_argparse.ArgumentParser = self.parser
        self.parse_args.sfffile = self.in_sffs
        self.parse_args.outfile = 'output.fastq'

        self._C()
        
        numseq = countlines('output.fastq','fastq')
        eq_(10,numseq)
    
class TestSffToFastq(Base):
    functionname = 'sff_to_fastq'

    def setUp(self):
        super(TestSffToFastq,self).setUp()
        self.outfh = MagicMock(file)()
        self.in_sffs = [self.sff_file]

    def test_empty_sff_list(self):
        r = self._C([], self.outfh)
        eq_(0, r)
        # Should not have written at all
        eq_(0,self.outfh.write.call_count)

    def test_single_sff(self):
        r = self._C(self.in_sffs, self.outfh)
        eq_(10,r)
        eq_(10,self.outfh.write.call_count)

    def test_multiple_sff(self):
        r = self._C(self.in_sffs*10, self.outfh)
        eq_(10*10,r)
        eq_(10*10,self.outfh.write.call_count)

    def test_output_arg_string(self):
        m = MagicMock(file)
        outfh = MagicMock(file)
        sfffh = open(self.sff_file,'rb')
        m.side_effect = [outfh, sfffh]
        with patch('__builtin__.open', m, create=True) as m:
            r = self._C(self.in_sffs, '/path/to/out.fastq')
            m.assert_has_calls([
                call('/path/to/out.fastq','w'),
                call(self.in_sffs[0],'rb')
            ])
            eq_(10,r)
            eq_(10,outfh.write.call_count)
