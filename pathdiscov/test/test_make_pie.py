import unittest2 as unittest
import mock

from StringIO import StringIO

from .. import make_pie

PHYLO = '''taxid\tcount\tsuperkingdom\tkingdom\tclass\torder\tfamily\tgenus\tspecies\tdescrip
1\t1.0\tBacteria\tKingdom\tClass\tOrder\tFamily\tGenus\tBact1\tDescription
1\t1.0\tBacteria\tBact2\t-\t-\t-\t-\t-\tDescription
1\t1.0\tViruses\tKingdom\tClass\tOrder\tFamily\tGenus\tVirus1\tDescription
1\t1.0\tViruses\t-\t-\t-\t-\t-\t-\tVirus2
1\t1.0\tWhatever\tKingdom\tMammalia\tOrder\tFamily\tGenus\tMammal1\tDescription
1\t1.0\tWhatever\tKingdom\tInsecta\tOrder\tFamily\tGenus\tInsecta1\tDescription
1\t1.0\tWhatever\tKingdom\tFoo\tOrder\tFamily\tGenus\tBar1\tDescription
1\t1.0\tWhatever\tKingdom\tFoo2\tOrder\tFamily\tGenus\tBar2\tDescription
1\t1.0\tPathogens1\tKingdom\tFoo3\tOrder\tFamily\tGenus\tBar3\tDescription
1\t1.0\tPathogens2\tKingdom\tFoo4\tOrder\tFamily\tGenus\tBar4\tDescription
'''

class TestFindBestName(unittest.TestCase):
    def setUp(self):
        self.phylo_headers = [
            'species', 'genus', 'family',
            'order', 'class', 'kingdom',
            'superkingdom'
        ]
        self.row = {
            'count': '1.0',
            'superkingdom': 'superkingdom',
            'species': 'species',
            'family': 'family',
            'taxid': '1',
            'class': 'class',
            'descrip': 'Description',
            'kingdom': 'kingdom',
            'genus': 'genus',
            'order': 'order',
        }

    def test_uses_species(self):
        cols = self.phylo_headers
        r = make_pie.find_best_name(self.row, self.phylo_headers)
        self.assertEqual('species', r)

    def test_works_backwords_and_finds_non_dash(self):
        # Make all -
        row = {}
        for k,v in self.row.items():
            row[k] = '-' 
        row['superkingdom'] = 'foo'
        r = make_pie.find_best_name(row, self.phylo_headers)
        self.assertEqual('foo', r)

    def test_falls_back_on_description(self):
        # Make all -
        row = {}
        for k,v in self.row.items():
            row[k] = '-' 
        row['descrip'] = 'foo'
        r = make_pie.find_best_name(row, self.phylo_headers)
        self.assertEqual('foo', r)

    def test_raises_valueerror_if_none_found(self):
        # Make all -
        row = {}
        for k,v in self.row.items():
            row[k] = '-' 
        row['descrip'] = '-'
        self.assertRaises(
            ValueError,
            make_pie.find_best_name, row, self.phylo_headers
        )

class TestHostVectorPathogen(unittest.TestCase):
    def setUp(self):
        self.phylo = PHYLO
        self.phylo_file = mock.MagicMock(
            __enter__ = mock.MagicMock(return_value=StringIO(self.phylo))
        )

    @mock.patch('__builtin__.open')
    def test_correct_calls(self, mopen):
        #mopen.return_value.read.return_value = self.phylo
        mopen.return_value = self.phylo_file
        r = make_pie.host_vector_pathogen(
            ['foo.phylo','bar.phylo'],
            ['Mammalia'], ['Insecta'], ['Viruses','Bacteria']
        )
        mopen.assert_has_calls([
            mock.call('foo.phylo'),
            mock.call().__enter__(),
            mock.call().__exit__(None, None, None),
            mock.call('bar.phylo'),
            mock.call().__enter__(),
            mock.call().__exit__(None, None, None),
        ])

    @mock.patch('__builtin__.open')
    def test_correct_defaults(self, mopen):
        #mopen.return_value.read.return_value = self.phylo
        mopen.return_value = self.phylo_file
        r = make_pie.host_vector_pathogen(
            ['foo.phylo'],
            ['Mammalia'], ['Insecta'], ['Viruses','Bacteria']
        )
        self.assertEqual(1.0, r[0]['Mammal1'])
        self.assertEqual(1, r[1])
        self.assertEqual(1.0, r[2]['Insecta1'])
        self.assertEqual(1, r[3])
        self.assertEqual(1.0, r[4]['Bact1'])
        self.assertEqual(1.0, r[4]['Bact2'])
        self.assertEqual(1.0, r[4]['Virus1'])
        self.assertEqual(1.0, r[4]['Viruses'])
        self.assertEqual(4, r[5])

    @mock.patch('__builtin__.open')
    def test_can_change_hostvectorpath_classnames(self, mopen):
        #mopen.return_value.read.return_value = self.phylo
        mopen.return_value = self.phylo_file
        r = make_pie.host_vector_pathogen(
            ['foo.phylo'],
            ['Foo'], ['Foo2'], ['Pathogens1','Pathogens2']
        )
        self.assertEqual(1.0, r[0]['Bar1'])
        self.assertEqual(1, r[1])
        self.assertEqual(1.0, r[2]['Bar2'])
        self.assertEqual(1, r[3])
        self.assertEqual(1.0, r[4]['Bar3'])
        self.assertEqual(1.0, r[4]['Bar4'])
        self.assertEqual(2, r[5])
