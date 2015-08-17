import tempfile
import os
try:
    import unittest2 as unittest
except ImportError:
    import unittest
from StringIO import StringIO

import mock

from pathdiscov import parallel_blast

#http://talby.rcs.manchester.ac.uk/~ri/_notes_sge/par_envs_and_integration.html
# nodename cpucount queue processorrange
#$PE_HOSTFILE
SGE_HOSTFILE = '''node1.localhost 1 queue UNDEFINED
node2.localhost 2 queue UNDEFINED
node3.localhost 3 queue UNDEFINED
'''

#$PBS_NODEFILE
PBS_MACHINEFILE = '''node1.localhost
node2.localhost
node2.localhost
node3.localhost
node3.localhost
node3.localhost
'''

hosts = [
    ('node1.localhost', 1),
    ('node2.localhost', 2),
    ('node3.localhost', 3),
]
sshlogins = []
for x in hosts:
    sshlogins += ['--sshlogin', '{1}/{0}'.format(*x)]

class TestParseHostfile(unittest.TestCase):
    def test_parses_pbs_hostfile(self):
        r = parallel_blast.parse_hostfile(
            StringIO(PBS_MACHINEFILE)
        )
        self.assertListEqual(hosts, r)

    def test_parses_sge_hostfile(self):
        r = parallel_blast.parse_hostfile(
            StringIO(SGE_HOSTFILE)
        )
        self.assertListEqual(hosts, r)

    def test_unknown_hostfile(self):
        self.assertRaises(
            ValueError,
            parallel_blast.parse_hostfile,
            StringIO('node1.localhost something something something')
        )

class TestGetHostfile(unittest.TestCase):
    @mock.patch.dict('pathdiscov.parallel_blast.os.environ', {'PBS_NODEFILE': 'foo'})
    def test_detects_pbs(self):
        r = parallel_blast.get_hostfile()
        self.assertEqual('foo', r)
    
    @mock.patch.dict('pathdiscov.parallel_blast.os.environ', {'PE_HOSTFILE': 'foo'})
    def test_detects_sge(self):
        r = parallel_blast.get_hostfile()
        self.assertEqual('foo', r)

    def test_detects_local(self):
        r = parallel_blast.get_hostfile()
        self.assertEqual('', r)

class TestGenerateSSHLogins(unittest.TestCase):
    def setUp(self):
        _, self.hostfile = tempfile.mkstemp()
        self.addCleanup(os.unlink, self.hostfile)
        
    def test_pbs_sshlogins(self):
        with open(self.hostfile, 'w') as fh:
            fh.write(PBS_MACHINEFILE)
        with mock.patch.dict('pathdiscov.parallel_blast.os.environ', {'PBS_NODEFILE': self.hostfile}):
            r = parallel_blast.generate_sshlogins()
            self.assertListEqual(sshlogins, r)

    def test_sge_sshlogins(self):
        with open(self.hostfile, 'w') as fh:
            fh.write(SGE_HOSTFILE)
        with mock.patch.dict('pathdiscov.parallel_blast.os.environ', {'PE_HOSTFILE': self.hostfile}):
            r = parallel_blast.generate_sshlogins()
            self.assertListEqual(sshlogins, r)

    def test_local_sshlogins(self):
        r = parallel_blast.generate_sshlogins(2)
        self.assertListEqual(['--sshlogin', '2/:'], r)

    def test_local_sshlogins_noninst(self):
        r = parallel_blast.generate_sshlogins()
        self.assertListEqual(['--sshlogin', ':'], r)

class TestParallelBlast(unittest.TestCase):
    def setUp(self):
        _, self.hostfile = tempfile.mkstemp()
        self.patch_sh_cmd = mock.patch('pathdiscov.parallel_blast.sh.Command')
        self.patch_sh_which = mock.patch('pathdiscov.parallel_blast.sh.which')
        self.mock_sh_which = self.patch_sh_which.start()
        self.mock_sh_cmd = self.patch_sh_cmd.start()
        self.addCleanup(self.patch_sh_cmd.stop)
        self.addCleanup(self.patch_sh_which.stop)
        self.infile = StringIO('>test\nACGT\n')
        self.outfile = StringIO()

    def test_command_string_is_correct(self):
        self.mock_sh_which.return_value = '/path/to/foon'
        parallel_blast.parallel_blast(
            self.infile, self.outfile, 5, '/path/db/nt', 'foon', 'barn',
            '-evalue 0.01 -otherblast arg'
        )
        self.mock_sh_cmd.assert_called_once_with('parallel')
        r = self.mock_sh_cmd.return_value.call_args
        blastcmd = r[0][8]
        self.assertIn('-task barn', blastcmd)
        self.assertIn('-db /path/db/nt', blastcmd)
        self.assertIn('-otherblast arg', blastcmd)
        self.assertIn(
            '-max_target_seqs {0}'.format(parallel_blast.MAX_TARGET_SEQS),
            blastcmd
        )
        self.assertIn('-outfmt "{0}"'.format(parallel_blast.BLAST_FORMAT), blastcmd)

    def test_localhost(self):
        self.mock_sh_which.return_value = '/path/to/foon'
        parallel_blast.parallel_blast(
            self.infile, self.outfile, 5, '/path/db/nt', 'foon', 'barn',
            '-evalue 0.01 -otherblast arg'
        )
        self.mock_sh_cmd.assert_called_once_with('parallel')
        r = self.mock_sh_cmd.return_value.call_args[0]
        self.assertIn('--sshlogin', r)
        self.assertIn('5/:', r)

    def test_remote_hosts(self):
        with open(self.hostfile, 'w') as fh:
            fh.write(PBS_MACHINEFILE)
        with mock.patch.dict('pathdiscov.parallel_blast.os.environ', {'PBS_NODEFILE': self.hostfile}):
            self.mock_sh_which.return_value = '/path/to/foon'
            parallel_blast.parallel_blast(
                self.infile, self.outfile, 5, '/path/db/nt', 'foon', 'barn',
                '-evalue 0.01 -otherblast arg'
            )
            self.mock_sh_cmd.assert_called_once_with('parallel')
            r = self.mock_sh_cmd.return_value.call_args[0]
            self.assertEqual(3, r.count('--sshlogin'))
            self.assertIn('1/node1.localhost', r)
            self.assertIn('2/node2.localhost', r)
            self.assertIn('3/node3.localhost', r)
