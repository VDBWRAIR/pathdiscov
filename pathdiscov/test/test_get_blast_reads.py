from os.path import *
import os

import unittest2 as unittest
import mock
import sh
import tempdir

from common import FIXTURES, THIS

FIXTURE = join(FIXTURES, 'get_blast_reads')

class TestGetBlastReads(unittest.TestCase):
    def setUp(self):
        self.tdir = tempdir.TempDir()
        os.chdir(self.tdir.name)
        self.FILES = [
            'contigfiles.lst',
            'fastqids.lst',
            'blastresult.fastq'
        ]
        self.r1_reads = [2,7,8,9]
        self.r2_reads = [2,3,4,7,9]
        # Each contig fastq is only 1 read
        self.contigs = [2,7,12,19,20,22,25,29,32,34,44,46,47,48,61,65,77,82]
        self.addCleanup(os.chdir, '/')

    def _check_fqids_file(self, read_count=None):
        if read_count is None:
            read_count = len(self.r1_reads+self.r2_reads)
        out = sh.wc(self.FILES[1], l=True)
        self.assertEqual(
            read_count, int(out.split()[0]),
            '{0} had incorrect number of lines'.format(self.FILES[1])
        )

    def _check_output_file_numreads(self, outputpath,
        contig_count=None, read_count=None):
        if contig_count is None:
            contig_count = len(self.contigs)
        if read_count is None:
            read_count = len(self.r1_reads+self.r2_reads)
        _contig_count = 0
        _read_count = 0
        with open(outputpath) as fh:
            for line in fh:
                if 'rickettsia' in line.lower():
                    _contig_count += 1
                elif line.startswith('@'):
                    _read_count += 1
        self.assertEqual(
            contig_count, _contig_count,
            'Contig Count was incorrect'
        )
        self.assertEqual(
            read_count, _read_count,
            'Read count was incorrect'
        )

    def test_combines_for_rickettsia(self):
        print sh.get_blast_reads(
            FIXTURE, 'Rickettsia', '--write-used', blastcol='genus',
            outputfile='output.fastq'
        )
        self.FILES[2] = 'output.fastq'
        for f in self.FILES:
            self.assertTrue(exists(f), '{0} did not exist'.format(f))
        self._check_output_file_numreads(self.FILES[2])
        self._check_fqids_file()

    def test_overwrites_existing_output(self):
        for f in self.FILES:
            open(f, 'w').close()
        print sh.get_blast_reads(
            FIXTURE, 'Rickettsia', '--write-used', '--overwrite'
        )
        for f in self.FILES:
            self.assertTrue(exists(f), '{0} did not exist'.format(f))
        self._check_output_file_numreads(self.FILES[2])
        self._check_fqids_file()

    def test_does_not_overwrite_existing_output(self):
        for f in self.FILES:
            open(f, 'w').close()
        out = sh.get_blast_reads(
            FIXTURE, 'Rickettsia', '--write-used', _ok_code=[1],
            _err_to_out=True
        )
        self.assertIn(
            'ValueError: {0} already exists'.format(self.FILES[2]),
            out
        )
        for f in self.FILES:
            self.assertTrue(exists(f), '{0} did not exist'.format(f))
        self._check_output_file_numreads(self.FILES[2], 0, 0)
        self._check_fqids_file()

    def test_does_not_write_used_if_not_specified(self):
        print sh.get_blast_reads(
            FIXTURE, 'Rickettsia'
        )
        print os.listdir('.')
        outfile = self.FILES[2]
        del self.FILES[2]
        for f in self.FILES:
            self.assertFalse(exists(f), '{0} exists'.format(f))
        self._check_output_file_numreads(outfile)

    def test_results_for_other_value(self):
        print sh.get_blast_reads(
            FIXTURE, 'Rickettsia_rickettsii_str._Morgan,_complete_genome',
            '--write-used', blastcol='descrip'
        )
        for f in self.FILES:
            self.assertTrue(exists(f), '{0} exists'.format(f))
        self._check_output_file_numreads(self.FILES[2], 0, 2)
        self._check_fqids_file(2)
