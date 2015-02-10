from common import *
from Bio import SeqIO

from usamriidPathDiscov.util import get_file_handle

class Base(BaseTester):
    modulepath = 'usamriidPathDiscov.stages.step1'

    def _check_id_fastq(self, filepath):
        ''' Check filepath to make sure sequential id renamed '''
        print '_check_id_fastq'
        for i, rec in enumerate(SeqIO.parse(get_file_handle(filepath),'fastq'),start=1):
            eq_(str(i), rec.id)

    def _count_reads(self, readpath, inputformat):
        print '_count_reads'
        count = 0
        for rec in SeqIO.parse(get_file_handle(readpath),inputformat):
            count += 1
        return count

    def _check_countfile(self, readpath, countpath, inputformat):
        print '_check_countfile'
        print readpath
        print countpath
        print inputformat
        ereads = self._count_reads(readpath,inputformat)
        rcount = None
        with open(countpath) as fh:
            line = fh.readline().rstrip().split()
            assert line[0] == 'rawfile', 'countfile did not have \'rawfile\''
            rcount = int(line[1])
        eq_(ereads, rcount)

    def _id_map(self, idmapfile):
        ''' Return [(,),] of renamed id: orig id  from idmapfile '''
        print '_id_map'
        idmap = []
        with open(idmapfile) as fh:
            for line in fh:
                line = line.rstrip().split()
                idmap.append(line)
        return idmap

    def _check_id_file(self, readfile, inputformat, idmapfile):
        print '_check_id_file'
        idmap = self._id_map(idmapfile)
        fh = get_file_handle(readfile)
        recs = SeqIO.parse(fh, inputformat)
        for rec, mapping in zip(recs, idmap):
            newid, oldid = mapping
            eq_(oldid, rec.id)

    def _check_logfilenames(self, logpath, samplename, timestamp):
        print '_check_logfilenames'
        logname = "{0}.{1}-out".format(samplename, timestamp)
        soutlog = join(logpath, logname + '.o')
        serrlog = join(logpath, logname + '.e')
        ok_(exists(soutlog), 'stdout log {0} is missing'.format(soutlog))
        ok_(exists(serrlog), 'stderr log {0} is missing'.format(serrlog))

class TestStep1(Base, BaseTempDir):
    def setUp(self):
        super(TestStep1,self).setUp()
        # Script name
        THIS = abspath(__file__) # this file abspath
        THIS = dirname(THIS) # test directory
        THIS = dirname(THIS) # usamriidPathDiscov directory
        script = join(THIS, 'step1','step1.pl')
        self.scriptpath = ['perl', script]
        script = join(THIS, 'stages', 'step1.py')
        self.scriptpath = ['python', script]
        #self.scriptpath = ['step1']
        # Just defaults
        self.kwargs = {
            '--sample': 'samplename',
            '--outputdir': self.tdir,
            '--paramfile': 'param.txt',
            '--R1': self.fastq_file,
            '--R2': self.fastq_file,
            '--timestamp': '1',
        }
        self.kwargs['--logs'] = join(self.kwargs['--outputdir'], 'logs')

    def _print_dir(self, dirpath):
        for root, dirs, files in os.walk(dirpath):
            for f in files:
                print join(root, f)

    def _check_step1(self):
        '''
        Using self.kwargs check everything is good
        '''
        projdir = self.kwargs['--outputdir']
        print subprocess.check_output('find {0}'.format(projdir), shell=True)
        self._print_dir(projdir)
        r1 = self.kwargs['--R1']
        r2 = self.kwargs.get('--R2','')
        # These two have to be valid symlinks
        # so instead of checking 2 things
        # just use them instead of R1.fastq/R2.fastq
        r1path = join(projdir, 'step1.R1')
        r2path = join(projdir, 'step1.R2')
        r1idpath = join(projdir, 'R1.id')
        r2idpath = join(projdir, 'R2.id')
        r1countpath = join(projdir, 'R1.count')
        r2countpath = join(projdir, 'R2.count')
        # Should get us sff or fastq
        inputformat = splitext(r1)
        if inputformat[1] == '.gz':
            inputformat = splitext(inputformat[0])
        inputformat = inputformat[1][1:]

        self._check_logfilenames(
            self.kwargs['--logs'],
            self.kwargs['--sample'],
            self.kwargs['--timestamp']
        )
        print '----'
        print 'R1: {0}'.format(r1)
        print 'R2: {0}'.format(r2)
        print r1path
        print os.readlink(r1path)
        self._check_id_fastq(r1path)
        self._check_id_file(r1, inputformat, r1idpath)
        self._check_countfile(r1, r1countpath, inputformat)
        r2 = self.kwargs.get('--R2',False)
        if r2 and exists(r2):
            self._check_id_fastq(r2path)
            self._check_id_file(r2, inputformat, r2idpath)
            self._check_countfile(r2, r2countpath, inputformat)

    def _make_dirs(self):
        ''' Make outputdir and logs '''
        os.makedirs(self.kwargs['--logs'])

    def run_script(self, script, **kwargs):
        '''
        Runs a script with kwargs created from arg_parse like paramters
        Script should be a list
        '''
        cmd = script
        for kwarg, arg in kwargs.iteritems():
            # Just append flag args without a value
            if arg is True:
                cmd.append(arg)
            else:
                cmd.append(kwarg)
                cmd.append(arg)
            
        return subprocess.check_output(cmd, stderr=subprocess.STDOUT)

    def test_missing_logdir(self):
        self.kwargs['--logs'] = join(self.tdir, 'missingdir')
        try:
            output = self.run_script(self.scriptpath,**self.kwargs)
            ok_(False,'Did not fail on missing logdir')
        except subprocess.CalledProcessError as e:
            ok_('missingdir does not exist' in e.output)

    def test_missing_outputdir(self):
        self.kwargs['--outputdir'] = join(self.tdir, 'missingdir')
        try:
            output = self.run_script(self.scriptpath,**self.kwargs)
            ok_(False, 'Did not fail on missing outputdir')
        except subprocess.CalledProcessError as e:
            ok_('cannot chdir to' in e.output)

    def test_only_r1(self):
        self._make_dirs()
        del self.kwargs['--R2']
        print 'Script output'
        print self.run_script(self.scriptpath,**self.kwargs)
        print '---------'
        self._check_step1()

    def test_r2_is_none_str(self):
        self._make_dirs()
        self.kwargs['--R2'] = 'none'
        print 'Script output'
        print self.run_script(self.scriptpath,**self.kwargs)
        print '---------'
        self._check_step1()

    def test_mate_not_exist(self):
        self._make_dirs()
        self.kwargs['--R2'] = '/path/to/none'
        print 'Script output'
        print self.run_script(self.scriptpath,**self.kwargs)
        print '---------'
        self._check_step1()

    def test_input_is_fastq(self):
        self._make_dirs()
        print self.run_script(self.scriptpath,**self.kwargs)
        self._check_step1()

    def test_input_is_fastq_gzip(self):
        self._make_dirs()
        self.kwargs['--R1'] = self.fastqgz_file
        self.kwargs['--R2'] = self.fastqgz_file
        print self.run_script(self.scriptpath,**self.kwargs)
        self._check_step1()

    def test_input_is_sff(self):
        self._make_dirs()
        self.kwargs['--R1'] = self.sff_file
        self.kwargs['--R2'] = self.sff_file
        print self.run_script(self.scriptpath,**self.kwargs)
        #with open(join(self.tdir,'R1.fastq')) as fh:
        #    print fh.read()
        self._check_step1()

    def test_input_is_sff_gzip(self):
        #raise SkipTest('Biopython is buggy with gzip right now')
        self._make_dirs()
        self.kwargs['--R1'] = self.sffgz_file
        self.kwargs['--R2'] = self.sffgz_file
        print self.run_script(self.scriptpath,**self.kwargs)
        self._check_step1()

    def test_symlink_already_exist(self):
        self._make_dirs()
        print self.run_script(self.scriptpath, **self.kwargs)
        try:
            print self.run_script(self.scriptpath, **self.kwargs)
        except subprocess.CalledProcessError as e:
            print e.output
            ok_(False, 'Failed')
