import os
import sys
from os.path import *
import subprocess
import re
from glob import glob
import shutil
import tempfile

from nose.tools import eq_, ok_, raises, assert_raises
from nose.plugins.attrib import attr

# I guess globl import will be ok this time
from usamriidPathDiscov.make_summary import MissingProjectFile

class Base(object):
    def setUp( self ):
        self.tdir = tempfile.mkdtemp( suffix='riidtest' )
        os.chdir( self.tdir )

    def tearDown( self ):
        os.chdir( '/tmp' )
        shutil.rmtree( self.tdir )

class BaseTest(Base):
    def setUp( self ):
        super( BaseTest, self ).setUp()
        self.mockproj = join( dirname(__file__), 'fixtures', 'mock_project' )
        self.setup_dirs()
    
    def setup_dirs( self ):
        self.results = join( self.mockproj, 'results' )
        self.step1 = join( self.results, 'step1' )
        self.phylo1 = join( self.results, 'iterative_blast_phylo_1' )
        self.phylo2 = join( self.results, 'iterative_blast_phylo_2' )
        self.ray2 = join( self.results, 'ray2_assembly_1' )

    def make_tmp_proj( self ):
        mockproj = join( self.tdir, 'mockproj' )
        shutil.copytree( self.mockproj, mockproj )
        self.mockproj = mockproj
        self.setup_dirs()

    def mock_summary( self, nr=1, nhr=2, nc=3, nbc=4, nru=5, nbu=6, n50=7, alen=8, contigs=[], unassembled={} ):
        return {
            'numreads': nr,
            'nonhostreads': nhr,
            'numcontig': nc,
            'numblastcontig': nbc,
            'numreadsunassembled': nru,
            'numblastunassembled': nbu,
            'contigs': contigs,
            'unassembled': unassembled,
            'n50': n50,
            'assemblylength': alen
        }

    def mock_contig( self, *args, **kwargs ):
        mock = {
            'accession': 'accession',
            'contigname': 'contigname',
            'description': 'description',
            'family': 'family',
            'genus': 'genus',
            'length': 1,
            'numreads': 1
        }
        mock.update( **kwargs )
        return mock

    def mock_unassembled( self, *args, **blastcols ):
        mock = {
            'bitscore': '80.5',
            'blast_alg': 'megablast',
            'count': 2,
            'descrip': 'Choristoneura_occidentalis_granulovirus,_complete_genome',
            'evalue': '3e-12',
            'family': 'Baculoviridae',
            'gapopen': '3',
            'genus': 'Betabaculovirus',
            'length': '82',
            'mismatch': '7',
            'order': '-',
            'pident': '85.37',
            'qend': '150',
            'qlen': '150',
            'qseqid': '1380743',
            'qstart': '69',
            'send': '25976',
            'slen': '104710',
            'sseqid': 'gi|84683224|gb|DQ333351.1|',
            'sstart': '25900',
            'superkingdom': 'Bacteria',
            'taxid': '364745',
            'accession': 'DQ333351.1'
        }
        mock.update( **blastcols )
        return mock

class TestN50(BaseTest):
    def _C(self, *args, **kwargs):
        from usamriidPathDiscov.make_summary import get_n50
        return get_n50(*args, **kwargs)

    def test_generates_correct_n50_from_list(self):
        r = self._C([2,2,2,3,3,4,8,8])
        eq_(6, r)

    def test_mockproj_contig_n50(self):
        from usamriidPathDiscov.make_summary import parse_tab_file
        print self.mockproj
        contiglenfile = join(self.ray2,'contig_len.txt')
        contiglens = parse_tab_file(contiglenfile)
        lens = []
        for name, len in contiglens:
            lens.append(int(len))
        r = self._C(lens)
        eq_(175, r)

class TestParseTabFile( BaseTest ):
    def setUp( self ):
        super( TestParseTabFile, self ).setUp()
        self.tfile = join( self.step1, 'R1.count' )

    def _C( self, *args, **kwargs ):
        from usamriidPathDiscov.make_summary import parse_tab_file
        return parse_tab_file( *args, **kwargs )

    def test_parses_nofields( self ):
        r = self._C( self.tfile )
        sum = 0
        for _, count in r:
            sum += int( count )
        eq_( 658336, sum )

    def test_parses_fields( self ):
        r = self._C( self.tfile, ['name', 'count'] )
        sum = 0
        for row in r:
            sum += int( row['count'] )
        eq_( 658336, sum )

    def test_inputfile_string( self ):
        r = self._C( self.tfile, ['name', 'count'] )
        sum = 0
        for row in r:
            sum += int( row['count'] )
        eq_( 658336, sum )

class TestReadCount( BaseTest ):
    def setUp( self ):
        super( TestReadCount, self ).setUp()
        self.iter1r1r2 = sorted( glob( join( self.phylo1, '*.count' ) ) )
        self.iter2r1r2 = sorted( glob( join( self.phylo2, '*.count' ) ) )
        self.step1r1r2 = sorted( glob( join( self.step1, '*.count' ) ) )

    def _C( self, *args, **kwargs ):
        from usamriidPathDiscov.make_summary import read_count
        return read_count( *args, **kwargs )

    def test_correct_count( self ):
        r = self._C( self.iter1r1r2[0] )
        eq_( 430, r['input'] )
        eq_( 365, r['megablast'] )
        eq_( 310, r['dc-megablast'] )

        r = self._C( self.iter2r1r2[0] )
        eq_( 101898, r['input'] )
        eq_( 6041, r['megablast'] )
        eq_( 5127, r['dc-megablast'] )

        for f in self.step1r1r2:
            r = self._C( f )
            eq_( 658336, r['rawfile'] )

    def test_nonint( self ):
        with open( 'R1.count', 'w' ) as fh:
            fh.write( 'intval\t1\n' )
            fh.write( 'floatval\t0.5\n' )
        r = self._C( 'R1.count' )
        eq_( 1, r['intval'] )
        eq_( 1, r['floatval'] )

class TestTotalReads( BaseTest ):
    def _C( self, *args, **kwargs ):
        from usamriidPathDiscov.make_summary import total_reads
        return total_reads( *args, **kwargs )

    def test_correct_count( self ):
        r = self._C( self.mockproj )
        eq_( 658336 + 658336, r )

    @raises(MissingProjectFile)
    def test_missing_files( self ):
        self.make_tmp_proj()
        countfiles = glob( join(self.mockproj, 'results', 'step1', '*.count') )
        for f in countfiles:
            os.unlink(f)
        self._C( self.mockproj )

class TestNonHostNumReads( BaseTest ):
    def _C( self, *args, **kwargs ):
        from usamriidPathDiscov.make_summary import non_host_num_reads
        return non_host_num_reads( *args, **kwargs )

    def test_correct_count( self ):
        r = self._C( self.mockproj )
        eq_( 189710 + 290257, r )

    @raises(MissingProjectFile)
    def test_missing_files( self ):
        self.make_tmp_proj()
        countfiles = glob( join(self.mockproj, 'results', 'quality_filter', '*.count') )
        for f in countfiles:
            os.unlink(f)
        self._C( self.mockproj )

class TestNumContig( BaseTest ):
    def _C( self, *args, **kwargs ):
        from usamriidPathDiscov.make_summary import num_contig
        return num_contig( *args, **kwargs )

    def test_correct_result( self ):
        r = self._C( self.mockproj )
        eq_( 430, r[0] )
        eq_( 310, r[1] )

    @raises(MissingProjectFile)
    def test_missing_files( self ):
        self.make_tmp_proj()
        countfiles = glob( join(self.mockproj, 'results', 'iterative_blast_phylo_1', '*.count') )
        for f in countfiles:
            os.unlink( f )
        self._C( self.mockproj )

class TestR1R2Count( BaseTest ):
    def _C( self, *args, **kwargs ):
        from usamriidPathDiscov.make_summary import r1r2_count
        return r1r2_count( *args, **kwargs )

    def test_correct_count_phylo1( self ):
        r = self._C( self.phylo1 )
        eq_( 430, r['input'] )
        eq_( 365, r['megablast'] )
        eq_( 310, r['dc-megablast'] )

    def test_correct_count_phylo2( self ):
        r = self._C( self.phylo2 )
        eq_( 101898, r['input'] )
        eq_( 6041, r['megablast'] )
        eq_( 5127, r['dc-megablast'] )

    @raises(MissingProjectFile)
    def test_missing_files( self ):
        r = self._C( self.tdir )

class TestParseBlastReport( BaseTest ):
    def setUp( self ):
        super( TestParseBlastReport, self ).setUp()
        self.blastfile = join( self.phylo1, 'reports', 'contig.mock_project.top.smallreport.txt' )

    def _C( self, *args, **kwargs ):
        from usamriidPathDiscov.make_summary import parse_blast_report
        return parse_blast_report( *args, **kwargs )

    def test_parses_nofilter( self ):
        r = self._C( self.blastfile )
        r = list( r )
        eq_( 8, len( r ), 'Not correct number of rows' )
        eq_( 'c1', r[0]['qseqid'] )
        eq_( 'gi|470488608|ref|NR_076881.1|', r[0]['sseqid'] )
        eq_( 'Bacteria', r[0]['superkingdom'] )
        eq_( 'Burkholderiales', r[0]['order'] )
        eq_( 'Comamonadaceae', r[0]['family'] )
        eq_( 'Delftia', r[0]['genus'] )

    def test_parses_filter( self ):
        def filt( blastrow ):
            return blastrow['superkingdom'] == 'Bacteria'

        r = self._C( self.blastfile, filt )
        r = list( r )
        eq_( 6, len(r), 'Did not filter correctly. Returned {0} rows instead of {1}'.format(len(r),5) )
        eq_( 'c1', r[0]['qseqid'] )

class TestBlastResultsFor( BaseTest ):
    def setUp( self ):
        super( TestBlastResultsFor, self ).setUp()
        self.blastfile = join( self.phylo1, 'reports', 'contig.mock_project.top.smallreport.txt' )

    def _C( self, *args, **kwargs ):
        from usamriidPathDiscov.make_summary import blast_results_for_
        return blast_results_for_( *args, **kwargs )

    def test_filters_on_column( self ):
        r = self._C( self.blastfile, 'superkingdom', 'Bacteria' )
        r = list( r )
        eq_( 6, len(r), 'Did not filter correctly' )
        eq_( 'c1', r[0]['qseqid'] )

class TestContigInfo( BaseTest ):
    def _C( self, *args, **kwargs ):
        from usamriidPathDiscov.make_summary import contig_info
        return contig_info( *args, **kwargs )

    def test_gets_correct_info( self ):
        r = self._C( self.mockproj )
        # Just make sure that all the contigs made it into dict
        # and have length and num reads
        for i in range( 1, 9 ):
            l, nr = r['c{0}'.format(i)]
            int(l)
            int(nr)

    def test_missing_files( self ):
        self.make_tmp_proj()
        files = [
            join(self.mockproj, 'results', 'ray2_assembly_1', 'contig_numreads.txt'),
            join(self.mockproj, 'results', 'ray2_assembly_1', 'contig_len.txt'),
        ]
        ok_( len(files) > 0, files )
        for f in files:
            print "Removing " + f
            os.unlink(f)
            assert_raises( MissingProjectFile, self._C, self.mockproj )
            open(f,'w').close()

    def test_contig_files_not_same_size( self ):
        from contextlib import nested
        # Sometimes contig_numreads has less contigs than contig_len
        self.make_tmp_proj()
        clenf  = join(self.mockproj, 'results', 'ray2_assembly_1', 'contig_len.txt' )
        nreadsf = join(self.mockproj, 'results', 'ray2_assembly_1', 'contig_numreads.txt' )
        # Create a mocked smaller version with 10 contigs
        contigs = [('c'+str(i),i,i) for i in range(1,11)]
        with nested(open(nreadsf,'w'), open(clenf,'w')) as (nrfh,clfh):
            for cn, l, nr in contigs:
                # Contig number
                cnum = int(cn[1:])
                if cnum <= 8:
                    # Write numreads for contigs <= 8 so c9,c10 do not exist
                    nrfh.write( '{0}\t1{1}\n'.format(cn,nr) )
                # Write all contig.id
                clfh.write( '{0}\t{2}\n'.format(cn,cnum-1,l) )
        r = self._C( self.mockproj )
        # Assert that 1-8 exist
        ok_( [r['c'+str(i)] for i in range(1,9)] )
        # Assert that 9 & 10 are missing
        # and were set to -1
        eq_( ('9',-1), r.get( 'c9', False ) )
        eq_( ('10',-1), r.get( 'c10', False ) )

class TestContigsFor( BaseTest ):
    def _C( self, *args, **kwargs ):
        from usamriidPathDiscov.make_summary import contigs_for
        return contigs_for( *args, **kwargs )

    def test_gets_correct_info( self ):
        r = self._C( self.mockproj, 'superkingdom', 'Bacteria' )
        r = list( r )
        eq_( 6, len(r) )
        # Check the first item
        v = r[0]
        eq_( 'c1', v['contigname'] )
        eq_( 257, v['length'] )
        eq_( 2, v['numreads'] )
        eq_( 'NR_076881.1', v['accession'] )
        eq_( 'Comamonadaceae', v['family'] )
        eq_( 'Delftia', v['genus'] )
        # Check the last item too
        v = r[-1]
        eq_( 'c5', v['contigname'] )

    @raises(MissingProjectFile)
    def test_missing_files( self ):
        self.make_tmp_proj()
        files = glob( join( self.mockproj, 'results', 'iterative_blast_phylo_1', 'reports', '*smallreport*.txt' ) )
        files += [
            join(self.mockproj, 'results', 'ray2_assembly_1', 'contig_numreads.txt'),
            join(self.mockproj, 'results', 'ray2_assembly_1', 'contig_len.txt'),
        ]
        for f in files:
            os.unlink( f )
            r = self._C( self.mockproj, 'superkingdom', 'Bacteria' )
            next(r)
            open(f,'w').close()

    def test_missing_contig( self ):
        self.make_tmp_proj()
        lenfile = join( self.ray2, 'contig.id' )
        readsfile = join( self.ray2, 'contig_numreads.txt' )
        # Truncate numreads file so no contigs exist
        open(readsfile, 'w' ).close()
        r = list( self._C( self.mockproj, 'superkingdom', 'Bacteria' ) )
        for c in r:
            eq_( c['numreads'], -1 )

class TestUnassembledReads( BaseTest ):
    def _C( self, *args, **kwargs ):
        from usamriidPathDiscov.make_summary import unassembled_reads
        return unassembled_reads( *args, **kwargs )

    def test_correct_value( self ):
        r = self._C( self.mockproj )
        eq_( 101898, r[0] )
        eq_( 5127, r[1] )

class TestGroupBlastBy( BaseTest ):
    def setUp( self ):
        super( TestGroupBlastBy, self ).setUp()
        self.contigblast = join(self.phylo2,'reports','contig.mock_project.top.smallreport.txt')

    def _C( self, *args, **kwargs ):
        from usamriidPathDiscov.make_summary import group_blast_by_
        return group_blast_by_( *args, **kwargs )

    def test_groupby_order( self ):
        r = self._C(self.contigblast, 'superkingdom', None, 'order')
        eq_( 3, len(r) )
        eq_( 5, r['Diptera']['count'] )
        eq_( 2, r['Rickettsiales']['count'] )
        eq_( 1, r['Polyopisthocotylea']['count'] )

class TestUnassembledReport( BaseTest ):
    def _C( self, *args, **kwargs ):
        from usamriidPathDiscov.make_summary import unassembled_report
        return unassembled_report( *args, **kwargs )

    def test_correctly_parsed_accession( self ):
        r = self._C( self.mockproj, 'Eukaryota' )
        eq_( 'LM975171.1', r['Polystomatidae']['accession'] )

    def test_correct_report( self ):
        r = self._C( self.mockproj, 'Bacteria' )
        print r
        eq_( 1, len(r) )
        eq_( 2, r['Rickettsiaceae']['count'] )

    def test_parses_accession( self ):
        r = self._C( self.mockproj, 'Bacteria' )
        print r
        for sk, values in r.items():
            ok_( 'accession' in values, 'Accession was not parsed and put int report' )

    def test_missing_smallreports( self ):
        self.make_tmp_proj()
        smreport = glob( join( self.mockproj, 'results', 'iterative_blast_phylo_2', 'reports', 'contig.*.top.smallreport.txt' ) )
        for f in smreport:
            os.unlink( f )
        assert_raises( MissingProjectFile, self._C, self.mockproj, 'Bacteria' )

class TestSummary( BaseTest ):
    def setUp( self ):
        super( TestSummary, self ).setUp()

        self.make_tmp_proj()
        self.phylo1files = glob( join( self.mockproj, 'results', 'iterative_blast_phylo_1', 'reports', 'contig.*smallreport*.txt' ) )
        self.phylo2files = glob( join( self.mockproj, 'results', 'iterative_blast_phylo_2', 'reports', 'contig.*.top.smallreport.txt' ) )
        self.step1files = glob( join(self.mockproj, 'results', 'step1', '*.count') )
        self.qualfiles = glob( join(self.mockproj, 'results', 'quality_filter', '*.count') )
        self.contigcountfiles = glob( join(self.mockproj, 'results', 'iterative_blast_phylo_1', '*.count' ) )

    def _C( self, *args, **kwargs ):
        from usamriidPathDiscov.make_summary import summary
        return summary( *args, **kwargs )

    def test_correct_summary( self ):
        r = self._C( self.mockproj, 'superkingdom', 'Eukaryota' )
        eq_( 1316672, r['numreads'] )
        eq_( 479967, r['nonhostreads'] )
        eq_( 430, r['numcontig'] )
        eq_( 310, r['numblastcontig'] )
        eq_( 101898, r['numreadsunassembled'] )
        eq_( 5127, r['numblastunassembled'] )
        eq_( 144, r['n50'] )
        eq_( 250, r['assemblylength'] )

        contigs = list( r['contigs'] )
        contigs[0]['contigname'] = 'c1'
        contigs[-1]['contigname'] = 'c8'

        unreads = r['unassembled']
        eq_(1, len(unreads) )
        eq_(1, unreads['Polystomatidae']['count'])

    def test_missing_files( self ):
        for files in (self.phylo1files,self.phylo2files,self.step1files,self.qualfiles,self.contigcountfiles):
            # Temp move the files for this step
            for i,f in enumerate(files):
                shutil.copy( f, 'tmp.'+str(i) )
                os.unlink(f)
            try:
                self._C( self.mockproj, 'superkingdom', 'Bacteria' )
                ok_(False, "Did not raise MissingProjectFile")
            except MissingProjectFile as e:
                ok_(True)

            # Restore the files
            for i,f in enumerate(files):
                shutil.copy( 'tmp.'+str(i), f )
                os.unlink( 'tmp.'+str(i) )

    def test_only_r1( self ):
        allfiles = self.phylo1files + self.phylo2files + self.step1files + self.qualfiles + self.contigcountfiles
        # Remove all R2 files
        for f in allfiles:
            print f
            if 'R2' in f:
                print "Removing " + f
                os.unlink( f )
        self._C( self.mockproj, 'superkingdom', 'Bacteria' )

class TestFormatSummary( BaseTest ):
    def _C( self, *args, **kwargs ):
        from usamriidPathDiscov.make_summary import format_summary
        return format_summary( *args, **kwargs )

    def test_unassem_gt_contigs( self ):
        una = {}
        for i in range(1,4):
            una['Virus'+str(i)] = self.mock_unassembled(count=i)
        contig = [self.mock_contig(length=i,numreads=i,contigname='c'+str(i)) for i in range(1,3)]
        summary = self.mock_summary( contigs=contig, unassembled=una )
        r = self._C( summary )
        eq_( 3, len(r) )

    def test_unassem_lt_contigs( self ):
        una = {}
        for i in range(1,3):
            una['Virus'+str(i)] = self.mock_unassembled(count=i)
        contig = [self.mock_contig(length=i,numreads=i,contigname='c'+str(i)) for i in range(1,4)]
        summary = self.mock_summary( contigs=contig, unassembled=una )
        r = self._C( summary )
        eq_( 3, len(r) )

    def test_check_csv( self ):
        una = {
            'Virus1':
                self.mock_unassembled(count=1,accession='acc1',family='family1',genus='genus1',descrip='descrip1'),
            'Virus2':
                self.mock_unassembled(count=2,accession='acc2',family='family2',genus='genus2',descrip='descrip2')
        }
        contig = [
            self.mock_contig(length=1,numreads=2,contigname='c1',accession='ca',family='cfam',genus='cgen',description='cdesc'),
            self.mock_contig(length=10,numreads=20,contigname='c2',accession='cb',family='cfam',genus='cgen',description='cdesc')
        ]
        summary = self.mock_summary(contigs=contig, unassembled=una)
        # Check summary line with contig and unassembled read
        r = self._C( summary )
        e = ['','1','2','3','4','7','8','c1','1','2','ca','cfam','cgen','cdesc','5','6','2','acc2','family2','genus2','descrip2']
        print e
        print r[0].split('\t')
        print '---------'
        eq_(e, r[0].split('\t'))
        # Check summary line with only unassembled read
        e = ['','','','','','','','c2','10','20','cb','cfam','cgen','cdesc','','','1','acc1','family1','genus1','descrip1']
        print e
        print r[1].split('\t')
        eq_( e, r[1].split('\t') )

    def test_sorted_unassembled( self ):
        una = {}
        for i in range(1,11):
            mu = self.mock_unassembled(count=i,accession='acc1',family='family1',genus='genus1',descrip='descrip1')
            una['Virus'+str(i)] = mu
        contig = []
        summary = self.mock_summary(nr=10,nhr=10,nc=0,nbc=0,nru=10,nbu=10,n50=0,contigs=[], unassembled=una)
        rows = self._C( summary )
        # First item should be count == 10
        for i, r in enumerate( rows, 1 ):
            e = 11 - i
            print r
            line = r.split('\t')
            print line
            count = line[17]
            eq_( e, int(count), 'Count should be {0} but got {1} at index {2}'.format(e,count,i) )

class TestFormatDic( BaseTest ):
    def _C( self, *args, **kwargs ):
        from usamriidPathDiscov.make_summary import format_dict
        return format_dict( *args, **kwargs )

    def test_formats_contig( self ):
        from usamriidPathDiscov.make_summary import contigs_for
        contigs = list( contigs_for( self.mockproj, 'superkingdom', 'Bacteria' ) )
        r = self._C( contigs[0], ('contigname','length','numreads') )
        eq_( 'c1\t257\t2', r )

    def test_formats_unassem( self ):
        from usamriidPathDiscov.make_summary import unassembled_report
        ur = unassembled_report( self.mockproj, 'Eukaryota' )
        f = 'Polystomatidae'
        r = self._C( ur[f], ('count','accession','family') )
        eq_( '1\tLM975171.1\tPolystomatidae', r )
