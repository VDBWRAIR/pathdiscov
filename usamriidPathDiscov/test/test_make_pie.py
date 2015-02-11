import unittest2 as unittest
import mock

from .. import make_pie

PHYLO = '''taxid	count	superkingdom	kingdom	class	order	family	genus	species	descrip
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
            'taxid', 'count', 'superkingdom',
            'kingdom', 'class', 'order',
            'family', 'genus', 'species',
            'descrip'
        ]

    def test_uses_species(self):
        cols = self.phylo_headers
        r = make_pie.find_best_name(cols)
        self.assertEqual(
            cols[8],
            r
        )

    def test_works_backwords_and_finds_non_dash(self):
        # Make all -
        cols = ['-' for i in range(len(self.phylo_headers))]
        cols[2] = 'Superkingdom'
        r = make_pie.find_best_name(cols)
        self.assertEqual(
            cols[2],
            r
        )

    def test_falls_back_on_description(self):
        # Make all -
        cols = ['-' for i in range(len(self.phylo_headers))]
        cols[-1] = 'description'
        r = make_pie.find_best_name(cols)
        self.assertEqual(
            cols[-1],
            r
        )

    def test_raises_valueerror_if_none_found(self):
        # Make all -
        cols = ['-' for i in range(len(self.phylo_headers))]
        self.assertRaises(
            ValueError,
            make_pie.find_best_name, cols
        )

class TestHostVectorPathogen(unittest.TestCase):
    def setUp(self):
        self.phylo = PHYLO

    @mock.patch('__builtin__.open')
    def test_correct_calls(self, mopen):
        mopen.return_value.read.return_value = self.phylo
        r = make_pie.host_vector_pathogen(
            ['foo.phylo','bar.phylo'],
            ['Mammalia'], ['Insecta'], ['Viruses','Bacteria']
        )
        mopen.assert_has_calls([
            mock.call('foo.phylo'),
            mock.call().read(),
            mock.call('bar.phylo'),
            mock.call().read(),
        ])

    @mock.patch('__builtin__.open')
    def test_correct_defaults(self, mopen):
        mopen.return_value.read.return_value = self.phylo
        r = make_pie.host_vector_pathogen(
            ['foo.phylo','bar.phylo'],
            ['Mammalia'], ['Insecta'], ['Viruses','Bacteria']
        )
        self.assertEqual(2.0, r[0]['Mammal1'])
        self.assertEqual(2, r[1])
        self.assertEqual(2.0, r[2]['Insecta1'])
        self.assertEqual(2, r[3])
        self.assertEqual(2.0, r[4]['Bact1'])
        self.assertEqual(2.0, r[4]['Bact2'])
        self.assertEqual(2.0, r[4]['Virus1'])
        self.assertEqual(2.0, r[4]['Viruses'])
        self.assertEqual(8, r[5])

    @mock.patch('__builtin__.open')
    def test_can_change_hostvectorpath_classnames(self, mopen):
        mopen.return_value.read.return_value = self.phylo
        r = make_pie.host_vector_pathogen(
            ['foo.phylo','bar.phylo'],
            ['Foo'], ['Foo2'], ['Pathogens1','Pathogens2']
        )
        self.assertEqual(2.0, r[0]['Bar1'])
        self.assertEqual(2, r[1])
        self.assertEqual(2.0, r[2]['Bar2'])
        self.assertEqual(2, r[3])
        self.assertEqual(2.0, r[4]['Bar3'])
        self.assertEqual(2.0, r[4]['Bar4'])
        self.assertEqual(4, r[5])
