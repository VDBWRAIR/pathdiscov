from StringIO import StringIO
from os.path import *
import os

from nose.tools import eq_, ok_, raises
from nose.plugins.attrib import attr
from mock import Mock, MagicMock, patch, call

from pathdiscov import verifyproject

class Base(object):
    modulepath = 'pathdiscov.verifyproject'

    def _C( self, *args, **kwargs):
        '''
        Set modulepath as instance variable to be the name of the module
        Set functionname as instance variable to automagically run that function
        with self._C
        '''
        m = __import__( self.modulepath, fromlist=[self.functionname] )
        return getattr(m,self.functionname)( *args, **kwargs )

    def setUp(self):
        self.filestemplate = '''${projpath}
${projpath}/file1
${projpath}/${projname}
${projpath}/${projname}/file2
'''
        self.projpath = '/path/to/project'
        self.projname = 'project'

        self.mock_open = Mock(return_value=MagicMock(spec=file))
        self.mock_exists = Mock(return_value=True)
        self.mock_stat = Mock(return_value=Mock(st_size=100))
        verifyproject.open = self.mock_open
        verifyproject.exists = self.mock_exists
        verifyproject.os.stat = self.mock_stat

    def reset_imports(self):
        reload(os.path)
        reload(os)
        verifyproject.open = open
        verifyproject.exists = exists
        verifyproject.os.stat = os.stat
        reload(verifyproject)

    def with_filestemplate_to_listing_mock(self):
        self.ftl_mock = Mock(
            return_value=verifyproject.filestemplate_to_listing(
                self.projpath, self.projname, self.filestemplate
            )
        )
        verifyproject.filestemplate_to_listing = self.ftl_mock

class TestVerifyProject(Base):
    functionname = 'verify_project'
    
    def setUp(self):
        super(TestVerifyProject,self).setUp()
        self.with_filestemplate_to_listing_mock()
        self.templates = [
            'template1.lst',
            'template2.lst',
        ]

    def test_returns_compiled_list(self):
        self.mock_exists.return_value = False
        r = self._C(self.projpath, self.projname, self.templates)
        eq_(8, len(r))

class TestFilestemplateToListing(Base):
    functionname = 'filestemplate_to_listing'

    def test_returns_list(self):
        r = self._C(self.projpath,self.projname,self.filestemplate)
        ok_(isinstance(r,list), 'Did not return list')
        for item in r:
            ok_(item.startswith(self.projpath),'Did not have projpath in list items')
        eq_(4, len(r))

class TestVerifyFiles(Base):
    functionname = 'verify_files'

    def setUp(self):
        super(TestVerifyFiles,self).setUp()
        self.with_filestemplate_to_listing_mock()

    def test_has_all_files_dirs(self):
        r = self._C(self.projpath, '')
        eq_([], r)

    @raises(ValueError)
    def test_missing_template_file_raises_exception(self):
        self.mock_open.side_effect = IOError
        self._C('','')

    def test_missing_filedirectory(self):
        self.mock_exists.side_effect = [True,True,False,False]
        r = self._C(self.projpath, '')
        expected = [
            (
                join(self.projpath,basename(self.projpath)),
                'Missing'
            ),
            (
                join(self.projpath,self.projname,'file2'),
                'Missing'
            ),
        ]
        eq_(expected, r)

    def test_file_size_0(self):
        self.mock_stat.return_value.st_size = 0
        r = self._C(self.projpath, '')
        expected = [
            (join(self.projpath), 'Size zero'),
            (join(self.projpath,'file1'), 'Size zero'),
            (join(self.projpath,self.projname), 'Size zero'),
            (join(self.projpath,self.projname,'file2'), 'Size zero'),
        ]
        eq_(expected, r)

class TestVerifyStandardStagesFiles(Base):
    functionname = 'verify_standard_stages_files'

    def setUp(self):
        super(TestVerifyStandardStagesFiles,self).setUp()
        this = dirname(abspath(__file__))
        self.mockproj = join(this, 'fixtures', 'mock_project')
        self.templatedir = join(dirname(this), 'output_files_templates')
        self.reset_imports()

    def test_runs_on_mock_proj(self):
        r = self._C(self.mockproj, self.templatedir)
        for item in r:
            if item[1] != 'Size zero':
                ok_(False, 'Found a missing file {0}'.format(item))
