from os.path import *

def test_tst_databases_are_built():
    assert exists('rikkcdna.nin'), "You did not build databases. See rikkcdna/Readme.rst"
