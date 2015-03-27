import unittest2 as unittest
import mock

from .. import tasks
    
@mock.patch('usamriidPathDiscov.tasks.runCommand')
class TestCreateQuality(unittest.TestCase):
    def test_calls_runcommand_correctly(self, mock_runcommand):
        cmd = 'fastqc foo.fasta -o somewhere 2>&1 | ' \
                'tee -a somewhere/analysis_quality.log'
        tasks.createQuality('foo.fasta', 'somewhere')
        mock_runcommand.assert_called_once_with(
            cmd, True
        )
