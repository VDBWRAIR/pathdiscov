import unittest2 as unittest
import mock

from .. import helpers

@mock.patch.object(helpers, 'exists')
@mock.patch.object(helpers, 'os')
class TestSymlink(unittest.TestCase):
    def setUp(self):
        self.src = '/path/to/dir/foo'
        self.dst = '/path/to/bar'

    def test_unlinks_existing_destination(self, mos, mexists):
        mexists.side_effect = [True, True]
        helpers.symlink(self.src, self.dst)
        mos.symlink.assert_called_once_with(
            'dir/foo', self.dst
        )
        mos.unlink.assert_called_once_with(
            self.dst
        )

    def test_src_does_not_exist_so_does_nothing(self, mos, mexists):
        mexists.side_effect = [False, True]
        helpers.symlink(self.src, self.dst)
        self.assertFalse(mos.symlink.return_value.called)
        self.assertFalse(mos.unlink.return_value.called)

    def test_symlinks_correctly(self, mos, mexists):
        mexists.side_effect = [False, True]
        helpers.symlink(self.src, self.dst)
        mos.symlink.assert_called_once_with(
            'dir/foo', self.dst
        )
        self.assertFalse(mos.unlink.return_value.called)

@mock.patch.object(helpers, 'exists')
@mock.patch.object(helpers, 'os')
class TestCreateNewProject(unittest.TestCase):
    def setUp(self):
        self.projdir = '/path/to/projdir'

    def test_copies_existing_project_if_exists(self, mos, mexists):
        mexists.side_effect = [True, False]
        helpers.create_new_project(self.projdir)
        mos.makedirs.assert_called_once_with(self.projdir)
        mos.rename.assert_called_once_with(self.projdir, self.projdir + '.bk')

    def test_removes_backup_directory_if_exists(self, mos, mexists):
        mexists.side_effect = [True, True]
        with mock.patch.object(helpers, 'shutil') as mshutil:
            helpers.create_new_project(self.projdir)
            mshutil.rmtree.assert_called_once_with(self.projdir + '.bk')
        mos.makedirs.assert_called_once_with(self.projdir)
        mos.rename.assert_called_once_with(self.projdir, self.projdir + '.bk')

    def test_creates_new_project_if_notexists(self, mos, mexists):
        mexists.side_effect = [False, False]
        helpers.create_new_project(self.projdir)
        mos.makedirs.assert_called_once_with(self.projdir)

class TestRunCommand(unittest.TestCase):
    def setUp(self):
        self.patch_popen = mock.patch('pathdiscov.helpers.subprocess.Popen')
        self.patch_stderr = mock.patch('pathdiscov.helpers.sys.stderr')
        self.patch_stdout = mock.patch('pathdiscov.helpers.sys.stdout')
        self.mock_stderr = self.patch_stderr.start()
        self.mock_stdout = self.patch_stdout.start()
        self.mock_popen = self.patch_popen.start()
        self.addCleanup(self.patch_popen.stop)
        self.addCleanup(self.patch_stderr.stop)
        self.addCleanup(self.patch_stdout.stop)
        self.mock_popen.return_value.communicate.return_value = (None,None)

    def test_runs_command_specified(self):
        import subprocess
        cmd = 'foo bar baz'
        helpers.runCommand(cmd, True)
        self.mock_popen.assert_called_once_with(
            'foo bar baz', shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        self.mock_popen.return_value.communicate.assert_called()

    def test_writes_cmd_to_console_when_specified_to_do_so(self):
        helpers.runCommand('foo bar baz', True)
        self.assertEqual(
            [mock.call('foo bar baz',), mock.call('\n',)],
            self.mock_stdout.write.call_args_list
        )

    def test_just_runs_command_with_no_print_cmd_option(self):
        helpers.runCommand('foo bar baz', False)
        self.assertEqual(
            [],
            self.mock_stdout.write.call_args_list
        )
