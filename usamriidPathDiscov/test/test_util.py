from common import *

class Base(BaseTester):
    modulepath = 'usamriidPathDiscov.util'

    def setUp(self):
        super(Base, self).setUp()
        self.fastq_data = '@1\nATGC\n+\n!!!!\n' \
            '@2\nATGC\n+\n!!!!\n' \
            '@3\nATGC\n+\n!!!!\n'
        self.fasta_data = '>1\nATGC\n>2\nATGC\n>3\nATGC\n'

class TestSeqrecIdToNum(Base):
    functionname = 'seqrec_id_to_num'

    def test_converts_to_sequential_integers_start_1(self):
        r = self._C(self.mock_fastq_handle)
        for i, item in enumerate(r,start=1):
            origid, rec = item
            eq_(str(i), rec.id)
            eq_('id{0}'.format(i), origid)

    def test_converts_to_sequential_integers_start_not_1(self):
        r = self._C(self.mock_fastq_handle, 10)
        for i, item in enumerate(r,start=10):
            origid, rec = item
            eq_(str(i), rec.id)
            eq_('id{0}'.format(i-9), origid)

class TestChangeFastqId(Base):
    functionname = 'change_fastq_id'
    
    def setUp(self):
        super(TestChangeFastqId,self).setUp()
        self.outfqfh = MagicMock(file)
        self.outmapfh = MagicMock(file)

    def test_writes_files(self):
        self._C(self.mock_fastq_handle, self.outfqfh, self.outmapfh)
        eq_(3, self.outfqfh.write.call_count)
        self.outmapfh.write.assert_has_calls([
            call('1\tid1\n'),
            call('2\tid2\n'),
            call('3\tid3\n'),
        ])

    def test_trims_plusline(self):
        # Add line with that extra stuff on + line
        self.mock_fastq_handle.lines = [
            '@id1\n',
            'ATGC\n',
            '+id1\n',
            '!!!!\n'
        ]
        self._C(self.mock_fastq_handle, self.outfqfh, self.outmapfh)
        self.outmapfh.write.assert_has_calls(call('1\tid1\n'))
        self.outfqfh.write.assert_has_calls(call(
            '@1\nATGC\n+\n!!!!\n'
        ))

    def test_converts_dot_to_N(self):
        # Add line with that extra stuff on + line
        self.mock_fastq_handle.lines = [
            '@id1\n',
            'AT.C\n',
            '+\n',
            '!!!!\n'
        ]
        self._C(self.mock_fastq_handle, self.outfqfh, self.outmapfh)
        self.outmapfh.write.assert_has_calls(call('1\tid1\n'))
        self.outfqfh.write.assert_has_calls(call(
            '@1\nATNC\n+\n!!!!\n'
        ))

    def test_file_args_are_strings(self):
        mock_open = MagicMock(file)
        mock_open.side_effect = [self.mock_fastq_handle, self.outfqfh, self.outmapfh]
        with patch('__builtin__.open', mock_open, create=True):
            self._C('input.fastq', 'output.fastq', 'map.txt')
            eq_(3, self.outfqfh.write.call_count)
            self.outmapfh.write.assert_has_calls([
                call('1\tid1\n'),
                call('2\tid2\n'),
                call('3\tid3\n'),
            ])
            mock_open.assert_has_calls([
                call('input.fastq'),
                call('output.fastq','w'),
                call('map.txt','w'),
            ])

class TestGetPairedReads(Base):
    functionname = 'get_paired_reads'

    def setUp(self):
        super(TestGetPairedReads,self).setUp()
        self.data = self.fastq_data

    def test_has_no_paired(self):
        input1 = mock_open(read_data=self.data)
        data2 = self.data.replace('@1','@4').replace('@2','@5').replace('@3','@6')
        input2 = mock_open(read_data=data2)
        r = self._C(input1, input2, 'fastq')
        p, s1, s2 = r
        eq_(set([]), p)
        eq_(set(['1','2','3']), s1)
        eq_(set(['4','5','6']), s2)

    def test_has_no_single(self):
        eset = set(['1','2','3'])
        input1 = mock_open(read_data=self.data)
        input2 = mock_open(read_data=self.data)
        r = self._C(input1, input2, 'fastq')
        p, s1, s2 = r
        eq_(eset, p)
        eq_(set([]), s1)
        eq_(set([]), s2)

    def test_has_paired_and_single(self):
        input1 = mock_open(read_data=self.data)
        data2 = self.data.replace('@2','@4')
        input2 = mock_open(read_data=data2)
        r = self._C(input1, input2, 'fastq')
        p, s1, s2 = r
        eq_(set(['1','3']), p)
        eq_(set(['2']), s1)
        eq_(set(['4']), s2)

    def test_fasta_file(self):
        data2 = self.fasta_data.replace('>2','>4')
        input1 = mock_open(read_data=self.fasta_data)
        input2 = mock_open(read_data=data2)
        r = self._C(input1, input2, 'fasta')
        p, s1, s2 = r
        eq_(set(['1','3']), p)
        eq_(set(['2']), s1)
        eq_(set(['4']), s2)

class TestExtractReads(Base):
    functionname = 'extract_reads'

    def test_returns_generator(self):
        from types import GeneratorType
        r = self._C(self.mock_fastq_handle, 'fastq', set(['id1','id2','id3']))
        ok_(isinstance(r, GeneratorType), 'Did not return generator')

    def test_all_can_be_returned(self):
        ids = ['id1','id2','id3']
        r = self._C(self.mock_fastq_handle, 'fastq', set(ids))
        for eid, rrecord in zip(ids, list(r)):
            eq_(eid, rrecord.id)

    def test_empty_idlist_returns_empty(self):
        r = self._C(self.mock_fastq_handle, 'fastq', set([]))
        eq_([], list(r))

    def test_gets_only_from_list_fastq(self):
        ids = ['id1','id3']
        r = self._C(self.mock_fastq_handle, 'fastq', set(['id1','id3']))
        for eid, rrecord in zip(ids, list(r)):
            eq_(eid, rrecord.id)

    def test_gets_only_from_list_fasta(self):
        ids = ['id1','id3']
        r = self._C(self.mock_fasta_handle, 'fasta', set(['id1','id3']))
        for eid, rrecord in zip(ids, list(r)):
            eq_(eid, rrecord.id)

class TestGetCommonUnevenFiles(Base):
    functionname = 'get_common_uneven_files'

    def setUp(self):
        super(TestGetCommonUnevenFiles,self).setUp()
        # Output1
        self.o1s = MagicMock(file)
        self.o1p = MagicMock(file)
        # output2
        self.o2s = MagicMock(file)
        self.o2p = MagicMock(file)

    def test_writes_paired_unpaired_files_fastq(self):
        input1 = mock_open(read_data=self.fastq_data)
        input2 = mock_open(read_data=self.fastq_data.replace('@2','@4'))

        self._C(input1,input2,'fastq',self.o1s,self.o1p,self.o2s,self.o2p)

        self.o1s.write.assert_called_once_with('@2\nATGC\n+\n!!!!\n')
        self.o2s.write.assert_called_once_with('@4\nATGC\n+\n!!!!\n')

        pairedcalls = [call('@1\nATGC\n+\n!!!!\n'),call('@3\nATGC\n+\n!!!!\n')]
        self.o1p.write.assert_has_calls(pairedcalls)
        self.o2p.write.assert_has_calls(pairedcalls)

    def test_writes_paired_unpaired_files_fasta(self):
        input1 = mock_open(read_data=self.fasta_data)
        input2 = mock_open(read_data=self.fasta_data.replace('>2','>4'))

        self._C(input1,input2,'fasta',self.o1s,self.o1p,self.o2s,self.o2p)

        self.o1s.write.assert_has_calls([call('>2\n'),call('ATGC\n')])
        self.o2s.write.assert_has_calls([call('>4\n'),call('ATGC\n')])

        pairedcalls = [call('>1\n'),call('ATGC\n'),call('>3\n'),call('ATGC\n')]
        self.o1p.write.assert_has_calls(pairedcalls)
        self.o2p.write.assert_has_calls(pairedcalls)

    def test_file_args_are_stringpaths(self):
        input1 = mock_open(read_data=self.fastq_data)
        input2 = mock_open(read_data=self.fastq_data.replace('@2','@4'))

        mopen = MagicMock(file)
        mocked_files = [input1,input2,self.o1s,self.o1p,self.o2s,self.o2p]
        mocked_names = ['input1.fastq','input2.fastq','fastq','o1s.fastq','o1p.fastq','o2s.fastq','o2p.fastq']
        mopen.side_effect = mocked_files

        with patch('__builtin__.open', mopen, create=True):
            self._C(*mocked_names)

        print mopen.mock_calls
        calls = [
            call('input1.fastq'), call('input2.fastq'),
            call('o1s.fastq','w'), call('o1p.fastq','w'),
            call('o2s.fastq','w'), call('o2p.fastq','w'),
        ]
        mopen.assert_has_calls(calls)

class TestGetFileHandle(Base):
    functionname = 'get_file_handle'

    def test_inputpath_gzipped(self):
        with patch('usamriidPathDiscov.util.gzip') as mock_gzip:
            mock_gzip.open.return_value.mode = 1
            r = self._C('/path/to/file.gz')
            mock_gzip.open.assert_called_with('/path/to/file.gz')
            eq_(mock_gzip.open.return_value, r)
            eq_(mock_gzip.open.return_value.mode, 1)

    def test_inputpath_gzipped_sff(self):
        with patch('usamriidPathDiscov.util.gzip') as mock_gzip:
            mock_gzip.open.return_value.mode = 1
            r = self._C('/path/to/file.sff.gz')
            mock_gzip.open.assert_called_with('/path/to/file.sff.gz')
            eq_(mock_gzip.open.return_value.mode, 'r')

    def test_inputpath_noncompressed(self):
        with patch('__builtin__.open') as mock_open:
            r = self._C('/path/to/file')
            mock_open.assert_called_with('/path/to/file')
            eq_(mock_open.return_value, r)
