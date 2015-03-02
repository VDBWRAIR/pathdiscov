import glob
from collections import defaultdict
from os.path import (
    join, expandvars, expanduser, splitext,
    basename, dirname
)
import sys

import argparse

import yaml

def main():
    args = parse_args()
    with open(args.configfile) as fh:
        if verifydatabases(yaml.load(fh)):
            sys.stdout.write('All databases look to be ok\n')

def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        'configfile',
        help='Path to config.yaml to check databases paths'
    )

    return parser.parse_args()

def all_db_files_same_prefix(dbpath, excludeext=[]):
    '''
    Ensure that all db files contain the same amount of extensions
    For example, bowtie indexed dbs should all have 1, .bwt2
    NCBI has 10

    :param str dbpath: Path to prefix of db such as /path/to/nt, where nt is the
                        basename of all the database files
    '''
    allfiles = glob.glob(dbpath + '*')
    if not allfiles:
        raise ValueError('No database for {0}'.format(dbpath))
    # for each db file, get a set of extensions
    dbfiles = defaultdict(set)
    for f in allfiles:
        pth, ext = splitext(f)
        if ext in excludeext:
            continue
        dbfile = basename(pth)
        dbfiles[dbfile].add(ext)
    # Also will need to know max key(longest list of extensions)
    mlen = 0
    mkey = None
    for k,v in dbfiles.iteritems():
        if len(v) > mlen:
            mkey = k
    # Now get all missing extensions by comparing maxlen to rest
    missing = set()
    for k,v in dbfiles.iteritems():
        _missing = dbfiles[mkey].symmetric_difference(v)
        missing.update(_missing)
    return list(missing)

def verifydatabases(config):
    '''
    Given a config from config.yaml ensure that databases exist

    :param dict config: Config that contains databases setup
    :return: True or False if the databases are correct
    :rtype: bool
    '''
    dbpath = expandvars(expanduser(config['databases']))
    dbs = (
        'human_dna', 'human_rna', 'nt_db',
        'tax_nodes', 'tax_names','diamond_db', 'nr_db'
    )
    ntdb = join(dbpath, config['nt_db'])

    ok = True
    for db in dbs:
        pth = join(dbpath, config[db])
        try:
            missing = all_db_files_same_prefix(pth, ['.fa', '.gz', '.nal', 'dmnd'])
            if missing:
                ok = False
                sys.stderr.write('Please check {0} as the following database extensions' \
                    'are missing from some database files {1}\n'.format(pth,missing))
        except ValueError as e:
            sys.stderr.write(str(e) + '\n')
            ok = False
    return ok
