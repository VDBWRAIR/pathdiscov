# -*- coding: utf-8 -*-
from __future__ import print_function
from ez_setup import use_setuptools
use_setuptools()

import os
import sys
import imp
import subprocess
from glob import glob
#from paver.easy import  *

# Python 2.6 subprocess.check_output compatibility. Thanks Greg Hewgill!
if 'check_output' not in dir(subprocess):
    def check_output(cmd_args, *args, **kwargs):
        proc = subprocess.Popen(
            cmd_args, *args,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, **kwargs)
        out, err = proc.communicate()
        if proc.returncode != 0:
            raise subprocess.CalledProcessError(args)
        return out
    subprocess.check_output = check_output

from setuptools import setup, find_packages
from setuptools.command.test import test as TestCommand

try:
    import colorama
    colorama.init()  # Initialize colorama on Windows
except ImportError:
    # Don't require colorama just for running paver tasks. This allows us to
    # run `paver install' without requiring the user to first have colorama
    # installed.
    pass

# Add the current directory to the module search path.
sys.path.append('.')

# Constants
CODE_DIRECTORY = 'pathdiscov'
DOCS_DIRECTORY = 'docs'
TESTS_DIRECTORY = 'tests'
PYTEST_FLAGS = ['--doctest-modules']
# Import metadata. Normally this would just be:
#
#     from pathdiscov import metadata
#
# However, when we do this, we also import `pathdiscov/__init__.py'. If this
# imports names from some other modules and these modules have third-party
# dependencies that need installing (which happens after this file is run), the
# script will crash. What we do instead is to load the metadata module by path
# instead, effectively side-stepping the dependency problem. Please make sure
# metadata has no dependencies, otherwise they will need to be added to
# the setup_requires keyword.
metadata = imp.load_source(
    'metadata', os.path.join(CODE_DIRECTORY, 'metadata.py'))


# Miscellaneous helper functions

def get_project_files():
    """Retrieve a list of project files, ignoring hidden files.

    :return: sorted list of project files
    :rtype: :class:`list`
    """
    if is_git_project():
        return get_git_project_files()

    project_files = []
    for top, subdirs, files in os.walk('.'):
        for subdir in subdirs:
            if subdir.startswith('.'):
                subdirs.remove(subdir)

        for f in files:
            if f.startswith('.'):
                continue
            project_files.append(os.path.join(top, f))

    return project_files


def is_git_project():
    return os.path.isdir('.git')


def get_git_project_files():
    """Retrieve a list of all non-ignored files, including untracked files,
    excluding deleted files.

    :return: sorted list of git project files
    :rtype: :class:`list`
    """
    cached_and_untracked_files = git_ls_files(
        '--cached',  # All files cached in the index
        '--others',  # Untracked files
        # Exclude untracked files that would be excluded by .gitignore, etc.
        '--exclude-standard')
    uncommitted_deleted_files = git_ls_files('--deleted')

    # Since sorting of files in a set is arbitrary, return a sorted list to
    # provide a well-defined order to tools like flake8, etc.
    return sorted(cached_and_untracked_files - uncommitted_deleted_files)


def git_ls_files(*cmd_args):
    """Run ``git ls-files`` in the top-level project directory. Arguments go
    directly to execution call.

    :return: set of file names
    :rtype: :class:`set`
    """
    cmd = ['git', 'ls-files']
    cmd.extend(cmd_args)
    return set(subprocess.check_output(cmd).splitlines())


def print_success_message(message):
    """Print a message indicating success in green color to STDOUT.

    :param message: the message to print
    :type message: :class:`str`
    """
    try:
        import colorama
        print(colorama.Fore.GREEN + message + colorama.Fore.RESET)
    except ImportError:
        print(message)


def print_failure_message(message):
    """Print a message indicating failure in red color to STDERR.

    :param message: the message to print
    :type message: :class:`str`
    """
    try:
        import colorama
        print(colorama.Fore.RED + message + colorama.Fore.RESET,
              file=sys.stderr)
    except ImportError:
        print(message, file=sys.stderr)


def read(filename):
    """Return the contents of a file.

    :param filename: file path
    :type filename: :class:`str`
    :return: the file's content
    :rtype: :class:`str`
    """
    with open(os.path.join(os.path.dirname(__file__), filename)) as f:
        return f.read()


def _lint():
    """Run lint and return an exit code."""
    # Flake8 doesn't have an easy way to run checks using a Python function, so
    # just fork off another process to do it.

    # Python 3 compat:
    # - The result of subprocess call outputs are byte strings, meaning we need
    #   to pass a byte string to endswith.
    project_python_files = [filename for filename in get_project_files()
                            if filename.endswith(b'.py')]
    retcode = subprocess.call(
        ['flake8', '--max-complexity=10'] + project_python_files)
    if retcode == 0:
        print_success_message('No style errors')
    return retcode


def _test():
    """Run the unit tests.

    :return: exit code
    """
    # Make sure to import pytest in this function. For the reason, see here:
    # <http://pytest.org/latest/goodpractises.html#integration-with-setuptools-test-commands>  # NOPEP8
    import pytest
    # This runs the unit tests.
    # It also runs doctest, but only on the modules in TESTS_DIRECTORY.
    return pytest.main(PYTEST_FLAGS + [TESTS_DIRECTORY])


def _test_all():
    """Run lint and tests.

    :return: exit code
    """
    # return _lint() + _test()
    return _test()


# The following code is to allow tests to be run with `python setup.py test'.
# The main reason to make this possible is to allow tests to be run as part of
# Setuptools' automatic run of 2to3 on the source code. The recommended way to
# run tests is still `paver test_all'.
# See <http://pythonhosted.org/setuptools/python3.html>
# Code based on <http://pytest.org/latest/goodpractises.html#integration-with-setuptools-test-commands>  # NOPEP8
class TestAllCommand(TestCommand):

    def finalize_options(self):
        TestCommand.finalize_options(self)
        # These are fake, and just set to appease distutils and setuptools.
        self.test_suite = True
        self.test_args = []

    def run_tests(self):
        raise SystemExit(_test_all())


# define install_requires for specific Python versions
python_version_specific_requires = []

# as of Python >= 2.7 and >= 3.2, the argparse module is maintained within
# the Python standard library, otherwise we install it as a separate package
if sys.version_info < (2, 7) or (3, 0) <= sys.version_info < (3, 3):
    python_version_specific_requires.append('argparse')


# See here for more options:
# <http://pythonhosted.org/setuptools/setuptools.html>
setup_dict = dict(
    name=metadata.package,
    version=metadata.version,
    author=metadata.authors[0],
    author_email=metadata.emails[0],
    maintainer=metadata.authors[0],
    maintainer_email=metadata.emails[0],
    url=metadata.url,
    description=metadata.description,
    long_description=read('README.rst'),
    # Find a list of classifiers here:
    # <http://pypi.python.org/pypi?%3Aaction=list_classifiers>
    classifiers=[
        'Development Status :: Production',
        'Environment :: Console',
        'Intended Audience :: General User - Biologist',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: Implementation :: PyPy',
        'Topic :: Documentation',
        'Topic :: Software Development :: Libraries :: Python Modules',
        'Topic :: System :: Installation/Setup',
        'Topic :: System :: Software Distribution',
    ],
    # bwa=Bunch(
    #sdir = path('pathdiscov/download/bwa'),
    #bindir =path('pathdiscov/bin'),
    #),
    # minilib = Bunch(
    # extra_files=['doctools','virtual']
    #),
    packages=find_packages(exclude=(TESTS_DIRECTORY,)),
    install_requires=[
        # your module dependencies
    ] + python_version_specific_requires,
    # Allow tests to be run with `python setup.py test'.
    tests_require=[
        'pytest==2.5.1',
        'mock==1.0.1',
        'flake8==2.1.0',
    ],
    cmdclass={'test': TestAllCommand},
    zip_safe=False,  # don't use eggs
    entry_points={
        'console_scripts': [
            'pathdiscov_cli = pathdiscov.main:main',
            'verifyproject = pathdiscov.verifyproject:main',
            'make_summary = pathdiscov.make_summary:main',
            'make_pie = pathdiscov.make_pie:main',
            'verifydatabases = pathdiscov.verifydatabases:main',
            'linecount = pathdiscov.linecount:main',
            'sff2fastq = pathdiscov.sff2fastq:main',
            'step1 = pathdiscov.stages.step1:main',
            'get_blast_reads = pathdiscov.get_blast_reads:main',
        ],
    },
    # These all get copied to our installation's bin folder for us
    scripts = [
        'pathdiscov/download/bwa/bwa',
        'pathdiscov/download/samtools/samtools',
        'pathdiscov/download/EMBOSS-6.6.0/emboss/getorf',
        'pathdiscov/download/CAP3/cap3',
        'pathdiscov/download/wkhtmltopdf',
        'pathdiscov/download/ray/Ray',
        'pathdiscov/download/Ray2',
        'pathdiscov/download/bowtie2/bowtie2',
        'pathdiscov/download/diamond64/diamond',
        'pathdiscov/download/snap/snap',
        'pathdiscov/quality_filter/quality_filter.pl',
        'pathdiscov/host_map/host_map.pl',
        'pathdiscov/host_map/bwa_filter_host.sh',
        'pathdiscov/ray2_assembly/count_mapped.pl',
        'pathdiscov/ray2_assembly/joinlines.sh',
        'pathdiscov/ray2_assembly/ray2_assembly.pl',
		'pathdiscov/iterative_blast_phylo/annotate_blast.sh',
		'pathdiscov/iterative_blast_phylo/blast_wrapper.pl',
		'pathdiscov/iterative_blast_phylo/block_blast.sh',
		'pathdiscov/iterative_blast_phylo/fastq2fasta.awk',
		'pathdiscov/iterative_blast_phylo/format_iterative_blast_phylo.pl',
		'pathdiscov/iterative_blast_phylo/get_unblast_reads.pl',
		'pathdiscov/iterative_blast_phylo/iterative_blast_phylo.pl',
		'pathdiscov/iterative_blast_phylo/make_phylogeny_pipe.pl',
		'pathdiscov/iterative_blast_phylo/make_phylogeny.pl',
		'pathdiscov/iterative_blast_phylo/par_block_blast.pl',
		'pathdiscov/iterative_blast_phylo/phylogeny_wrapper.sh',
		'pathdiscov/iterative_blast_phylo/plot_phylo_percents.r',
		'pathdiscov/iterative_blast_phylo/plot_pie.r',
		'pathdiscov/iterative_blast_phylo/prepare_plot_phylo_percents.sh',
		'pathdiscov/iterative_blast_phylo/tableconcatlines',
		'pathdiscov/iterative_blast_phylo/taxid2queryid.pl',
		'pathdiscov/iterative_blast_phylo/weighted_count.pl',
        'pathdiscov/orf_filter/orf_filter.pl',
    ] + glob('pathdiscov/download/bowtie2/bowtie2-*') +
        glob('pathdiscov/download/blast-2.2.28/bin/*') +
        glob('pathdiscov/download/prinseq-lite-0.20.3/*.pl') +
        glob('pathdiscov/scripts/*'),
    package_data = {
        'pathdiscov': ['files/*', 'output_files_templates/*'],
    }
)

def runTask():
    cmd ="paver prepare"
    return subprocess.Popen(cmd, shell=True).communicate()


def main():
    import os
    try:
        import paver.tasks
    except ImportError:
        if os.path.exists("paver-minilib.zip"):
            import sys
            sys.path.insert(0, "paver-minilib.zip")
        else:
            raise ValueError("No paver on the path")
        import paver.tasks
    paver.tasks.main()

if __name__ == '__main__':
    main()
