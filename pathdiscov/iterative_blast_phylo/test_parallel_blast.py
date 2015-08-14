import tempfile
import os
try:
    import unittest2 as unittest
except ImportError:
    import unittest
from StringIO import StringIO

import mock

import parallel_blast

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
sshlogins = ['--sshlogin {1}/{0}'.format(*x) for x in hosts]

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
    @mock.patch.dict('parallel_blast.os.environ', {'PBS_NODEFILE': 'foo'})
    def test_detects_pbs(self):
        #with mock.patch.dict('parallel_blast.os.environ', {'PBS_NODEFILE': 'foo'}):
        r = parallel_blast.get_hostfile()
        self.assertEqual('foo', r)
    
    @mock.patch.dict('parallel_blast.os.environ', {'PE_HOSTFILE': 'foo'})
    def test_detects_sge(self):
        #with mock.patch.dict('parallel_blast.os.environ', {'PE_HOSTFILE': 'foo'}):
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
        with mock.patch.dict('parallel_blast.os.environ', {'PBS_NODEFILE': self.hostfile}):
            r = parallel_blast.generate_sshlogins()
            self.assertListEqual(sshlogins, r)

    def test_sge_sshlogins(self):
        with open(self.hostfile, 'w') as fh:
            fh.write(SGE_HOSTFILE)
        with mock.patch.dict('parallel_blast.os.environ', {'PE_HOSTFILE': self.hostfile}):
            r = parallel_blast.generate_sshlogins()
            self.assertListEqual(sshlogins, r)

    def test_local_sshlogins(self):
        r = parallel_blast.generate_sshlogins(2)
        self.assertListEqual(['--sshlogin 2/:'], r)

    def test_local_sshlogins_noninst(self):
        r = parallel_blast.generate_sshlogins()
        self.assertListEqual(['--sshlogin :'], r)

class TestParallelBlast(unittest.TestCase):
    def test_test(self):
        parallel_blast.parallel_blast(
            'test.fasta', 'out.blast', 5, '/path/db/nt', 'blastn', 'blastn',
            '-evalue 0.01'
        )
        self.assertFalse(True)
