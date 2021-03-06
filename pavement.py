# -*- coding: utf-8 -*-


from __future__ import print_function

import os
import sys
import time
import subprocess
from paver.easy import *
from os.path import islink,isfile,join,basename,dirname,exists,relpath
from paver.setuputils import setup
try:
    from paver.virtual import bootstrap,virtualenv
except ImportError, e:
    info(
        "VirtualEnv must be installed to enable 'paver bootstrap'. If you need this command, run: pip install virtualenv")

# Import parameters from the setup file.
sys.path.append('.')
from setup import (
    setup_dict, get_project_files, print_success_message,
    print_failure_message, _lint, _test, _test_all,
    CODE_DIRECTORY, DOCS_DIRECTORY, TESTS_DIRECTORY, PYTEST_FLAGS)

from paver.easy import options, task, needs, consume_args
from paver.setuputils import install_distutils_tasks


options(setup=setup_dict,
        bwa=Bunch(
            sdir=path('pathdiscov/download'),
            bindir=path('pathdiscov/bin')
        ),
         prinseq=Bunch(
            sdir=path('pathdiscov/download/prinseq-lite-0.20.3'),
            bindir=path('pathdiscov/bin')
        ),
        minilib=Bunch(
            # extra_files=['doctools','virtual']
        ),
         samtools=Bunch(
              sdir=path('pathdiscov/download'),
              bindir=path('pathdiscov/bin')
          ),
        Ray=Bunch(
            src = path('pathdiscov/download'),
            sfile =path('pathdiscov/download/ray'),
            sfile2 =path('pathdiscov/download/Ray2'),
            olink =path('pathdiscov/bin')
        ),
        CAP3=Bunch(
            sfile = path('pathdiscov/download/CAP3/cap3'),
            olink =path('pathdiscov/bin')

        ),
        FastQC=Bunch(
            url='http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.2.zip',
            downloads=path('pathdiscov/download'),
            installdir=join(sys.prefix,'lib')
        ),

        wkhtmltopdf=Bunch(
            sfile = path('pathdiscov/download/wkhtmltopdf'),
            olink =path('pathdiscov/bin')
        ),

        bowtie2=Bunch(
            sfile =path('pathdiscov/download/bowtie2/bowtie2*'),
            olink =path('pathdiscov/bin')

        ),
        getorf=Bunch(
            src=path('pathdiscov/download'),
            sfile =path('pathdiscov/download/EMBOSS-6.6.0'),
            olink =path('pathdiscov/bin')

        ),
        perl_modules=Bunch(        
            srcdir=path('pathdiscov/Local_Module'),
            dstdir=path(sys.prefix) / 'lib'
        ),
        snap=Bunch(
            sfile =path('pathdiscov/download/snap')
        ),

        #"""
        #local_lib=Bunch(
            #sdir=path('pathdiscov/download/local-lib-2.000011'),
            #bindir=path('pathdiscov/bin')
        #),
        #"""
        settings=Bunch(
            shell_file=path('pathdiscov/files/settings.sh'),
            shell_file_bk=path('pathdiscov/files/settings.sh.base'),
            bash_rc =path('pathdiscov/files/bashrc'),
            config =path('pathdiscov/files/config.yaml'),
            config_bk =path('pathdiscov/files/config.yaml.base'),
            dist_dir =path('.'),
            param_base =path('pathdiscov/files/sample.param.base'),
            param_work =path('pathdiscov/files/sample.param')
        ),
        #virtualenv=Bunch(
            #packages_to_install=['http://bitbucket.org/ianb/pip/get/2cb1db7b2baf.gz#egg=pip', 'urlgrabber', 'jstools', 'virtualenv'],
            #dest_dir='./',
            #install_paver=True,
            #script_name='bootstrap.py',
            #no_site_packages=True,
            #paver_command_line='post_bootstrap'
        #virtualenv=dict(
        #script_name="bootstrap.py",
        #packages_to_install = [
            ## Project dependencies
            #],
        #paver_command_line="init",
        virtualenv=Bunch(
            packages_to_install=[],
            no_site_packages=True
            )
        )


INSTRUCTIONS = """
Run
   $ source pathdiscov/bin/activate
to enter the virtual environment and
   $ deactivate
to exit the environment.
"""


install_distutils_tasks()


# Miscellaneous helper functions

def print_passed():
    # generated on http://patorjk.com/software/taag/#p=display&f=Small&t=PASSED
    print_success_message(r'''  ___  _   ___ ___ ___ ___
 | _ \/_\ / __/ __| __|   \
 |  _/ _ \\__ \__ \ _|| |) |
 |_|/_/ \_\___/___/___|___/
''')


def print_failed():
    # generated on http://patorjk.com/software/taag/#p=display&f=Small&t=FAILED
    print_failure_message(r'''  ___ _   ___ _    ___ ___
 | __/_\ |_ _| |  | __|   \
 | _/ _ \ | || |__| _|| |) |
 |_/_/ \_\___|____|___|___/
''')


class cwd(object):

    """Class used for temporarily changing directories. Can be though of
    as a `pushd /my/dir' then a `popd' at the end.
    """

    def __init__(self, newcwd):
        """:param newcwd: directory to make the cwd
        :type newcwd: :class:`str`
        """
        self.newcwd = newcwd

    def __enter__(self):
        self.oldcwd = os.getcwd()
        os.chdir(self.newcwd)
        return os.getcwd()

    def __exit__(self, type_, value, traceback):
        # This acts like a `finally' clause: it will always be executed.
        os.chdir(self.oldcwd)


# Task-related functions

def _doc_make(*make_args):
    """Run make in sphinx' docs directory.

    :return: exit code
    """
    if sys.platform == 'win32':
        # Windows
        make_cmd = ['make.bat']
    else:
        # Linux, Mac OS X, and others
        make_cmd = ['make']
    make_cmd.extend(make_args)

    # Account for a stupid Python "bug" on Windows:
    # <http://bugs.python.org/issue15533>
    with cwd(DOCS_DIRECTORY):
        retcode = subprocess.call(make_cmd)
    return retcode


# Tasks
#@task
#@virtualenv(dir="test2")
#def venv():
    #sh("paver bootstrap")
@task
def init():
    """Initializing everything so you can start working"""
    info( "virtual environment successfully bootstrapped.")
    info(INSTRUCTIONS)
'''
@task
@needs(['prepareEncoders', 'prepareDetectors', 'prepareTraffic', 'prepareUtilities'])
def prepared(options):
    pass

@task
def prepareEncoders(options):
    pass

@task
def prepareDetectors(options):
    pass

@task
def prepareTraffic(options):
    pass

def prepareUtilities(options):
    pass

@task
def clean(options):
     pass
'''
"""
Try, else:
    wget http://peak.telecommunity.com/dist/virtual-python.py
    python virtual-python.py --no-site-packages --prefix=`pwd`
    wget http://peak.telecommunity.com/dist/ez_setup.py
    ~/bin/python ez_setup.py
    ~/local/bin/easy_install virtualenv
    ~/local/bin/virtualenv --no-site-packages
"""
@task
def bootstrap(options):
    """create virtualenv in ./bootstrap """
    try:
        import virtualenv
    except ImportError, e:
        raise RuntimeError("Virtualenv is needed for bootstrap")
    options.virtualenv.no_site_packages = False
    options.bootstrap.no_site_packages = False
    #options.virtualenv.paver_command_line='prepared'
    call_task('paver.virtual.bootstrap')

@task
def source_shell(options):
    """ source global variabls """
    info("Setting PATH variable")
    currwd = os.getcwd()
    settings= path(currwd) / options.settings.shell_file
    sh('source %s' %(settings) )

@task
def download_install_fastqc(options):
    import zipfile
    from glob import glob
    dlpath = join(options.FastQC.downloads,'fastqc_v*.zip')
    fastqczip = glob(dlpath)
    # No need to redownload
    if not len(fastqczip):
        info("Downloading FastQC from %s" % options.FastQC.url)
        dlcmd = 'cd %s && [ ! -e fastqc*.zip ] && wget %s' % (options.FastQC.downloads,options.FastQC.url)
        sh(dlcmd)
    else:
        info("FastQC Already downloaded")
    fastqczip = glob(dlpath)
    fqcdir = join(options.FastQC.installdir,'FastQC')
    # Check to see if it is extracted already
    if not exists(fqcdir):
        info("Unpacking FastQC")
        zfh = zipfile.ZipFile(fastqczip[-1])
        zfh.extractall(options.FastQC.installdir)
        zfh.close()
    else:
        info("FastQC already unpacked")
    # Make symlink to bin
    src = relpath(join(fqcdir,'fastqc'),join(sys.prefix,'bin'))
    dst = join(sys.prefix,'bin','fastqc')
    if not exists(dst):
        info("Installing fastqc symlink")
        os.symlink(src,dst)
        os.chmod(dst, 0755)
    else:
        info("fastqc symlink already exists")

@task
def download_compile_bwa(options):
    """installs the current package"""
    bwabin=join(sys.prefix,'bin','bwa')
    if not exists(bwabin):
        info("Compiling BWA...")
        currwd = os.getcwd()
        sdir = path(currwd) / options.bwa.sdir
        sh('(cd %s; wget https://github.com/lh3/bwa/archive/0.7.10.tar.gz -O- | tar xzf -; mv bwa-* bwa; cd bwa; make; cd %s)' % (sdir, sdir))

@task
def download_compile_samtools(options):
    """installs the current package"""
    samtoolsbin=join(sys.prefix,'bin','samtools')
    if not exists(samtoolsbin):
        info("Compiling samtools....")
        currwd = os.getcwd()
        sdir = path(currwd) / options.samtools.sdir
        sh('(cd %s; wget https://github.com/samtools/htslib/archive/1.1.tar.gz -O- | tar xzf -; mv htslib-* htslib; wget https://github.com/samtools/samtools/archive/1.1.tar.gz -O- | tar xzf -; mv samtools-* samtools; cd samtools; make; cd %s)' % (sdir, sdir))

@task
def refRay(options):
    """Install  Ray assembler """
    info("Compiling Ray Assembler")
    currwd = os.getcwd()
    src = path(currwd) / options.Ray.src
    sfile = path(currwd) / options.Ray.sfile
    sfile2 = path(currwd) / options.Ray.sfile2
    olink = path(currwd) / options.Ray.olink
    if os.path.isdir("/usr/lib64/openmpi"):
        sh('cd %s; test ! -d ray && tar -xzvf ray.tar.gz; test ! -d RayPlatform && tar -xzvf RayPlatform.tar.gz; which mpicxx && mpicxx=$(which mpicxx) || mpicxx=/usr/lib64/openmpi/bin/mpicxx; export LD_LIBRARY_PATH=/usr/lib64/openmpi:/usr/lib64/openmpi/lib:$LD_LIBRARY_PATH; cd %s;make PREFIX=build2000 MPICXX=$mpicxx;' % (src,sfile))
    elif os.path.isdir("/usr/lib/openmpi"):
        sh('cd %s; test ! -d ray && tar -xzvf ray.tar.gz; test ! -d RayPlatform && tar -xzvf RayPlatform.tar.gz; which mpicxx && mpicxx=$(which mpicxx) || mpicxx=/usr/lib64/openmpi/bin/mpicxx; export LD_LIBRARY_PATH=/usr/lib/openmpi:/usr/lib/openmpi/lib:$LD_LIBRARY_PATH; export PATH=$PATH:/usr/lib/openmpi/bin; cd %s; make PREFIX=build2000 MPICXX=$mpicxx' % (src,sfile))
    else:
        info("Ray is not installed, ... please install `openmpi and openmpi-devel` and try again")
        sys.exit()

@task
def perl_modules(options):
    sh('cp -R {0} {1}'.format(options.srcdir, options.dstdir))

@task
def getorf(options):
    """Install  EMBOSS getorf """
    getorf=join(sys.prefix,'bin','getorf')
    if not exists(getorf):
        info("Compiling EMBOSS...")
        currwd = os.getcwd()
        src = path(currwd) / options.getorf.src
        sfile = path(currwd) / options.getorf.sfile
        sh('(cd %s; ./configure CC="cc" --without-x; make)' %(sfile))

def ensure_line_in_file(filepath, line):
    with open(filepath,'r+') as fh:
        found = False
        # Loop all lines in file
        for file_line in fh:
            # If found flag to not do anything
            if line in file_line:
                found = True
                break
        if not found:
            # Write the line to the end of the file
            fh.write(line)

@task
def installSnap(options):
    """Install  snap aligner """
    snap=join(sys.prefix,'bin','snap')
    if not exists(snap):
        info("Compiling snap aligner...")
        currwd = os.getcwd()
        sfile = path(currwd) / options.snap.sfile
        ddir = dirname(sfile)
        sh('(cd %s; git clone https://github.com/amplab/snap.git; cd snap;git checkout v0.15; make)' %(ddir))

@task
def modifyBashRC():
    "Append the content of setting.sh to .bashrc"
    import fileinput
    import re
    import subprocess
    info("Append path info to your .bashrc ")
    #sfile = os.path.abspath(options.settings.shell_file)
    #sfilebk = os.path.abspath(options.settings.shell_file_bk)
    #sh("cp %s %s" %(sfilebk, sfile))
    #info(sfile)
    #bashrc = os.path.expanduser("~/.bashrc")
    #bashrcbk = os.path.expanduser("~/.bashrc.bak")
    #appdir = os.getcwd()
    #info(sfile)
    #bashrcTemp = os.path.abspath(options.settings.bash_rc)
    #for line in fileinput.input(sfile, inplace=True, backup='.bak'):
            #line = re.sub(r'CWDAPP',  appdir, line.rstrip())
            #line = re.sub(r'GENOMEDIR',  dbdir, line.rstrip())
            #info(line)
    #info(bashrc)
    #if isfile(bashrc):
    #    sourceline = 'source %s' % sfile
    #    ensure_line_in_file(bashrc, sourceline)
    #else:
    #    info ("Creating a new bashrc file...")
    #    #sh("cat %s > %s" %(bashrcTemp, bashrc))
    #    sh("cat %s >> %s" %(bashrcTemp, bashrc))
    #    sh("echo 'source %s' >> %s" %(sfile, bashrc))


@task
def setupConfigFile():
    import fileinput
    import re
    dbdir = os.path.expanduser("~/databases")
    info("Set config file ....")
    conf = os.path.abspath(options.settings.config)
    info(conf)
    conf_bk = os.path.abspath(options.settings.config_bk)
    info(conf_bk)
    sh("cat  %s >  %s" %(conf_bk, conf))
    for line in fileinput.input(conf, inplace=True, backup='.bak'):
        line = re.sub(r'GENOMEDIR', dbdir, line.rstrip())
        info(line)

@task
@needs('install_python_dependencies','install_other_dependencies')
def install_dependencies():
    pass

@task
@needs('download_compile_bwa','download_compile_samtools','refRay','getorf','download_install_fastqc','installSnap','perl_modules')
def install_other_dependencies():
    pass

@task
def install_python_dependencies():
    sh('pip install -r requirements-dev.txt --download-cache pathdiscov/download/.pip_cache')

@task
@needs('install_dependencies')
def prepare():
    """Prepare complete environment
    """
    pass

@task
@needs('prepare','setuptools.command.install')
def install():
    pass

@task
@needs('install')
def develop():
    pass

@task
@needs('prepare', 'doc_html', 'setuptools.command.sdist')
def sdist():
    """Build the HTML docs and the tarball."""
    pass

@task
@consume_args
def unit(args):
    import nose
    nose.run_exit(
        argv=["nosetests"] + args
    )


@task
def test():
    """Run the unit tests."""
    raise SystemExit(_test())


@task
def lint():
    # This refuses to format properly when running `paver help' unless
    # this ugliness is used.
    ('Perform PEP8 style check, run PyFlakes, and run McCabe complexity '
     'metrics on the code.')
    raise SystemExit(_lint())


@task
def test_all():
    """Perform a style check and run all unit tests."""
    retcode = _test_all()
    if retcode == 0:
        print_passed()
    else:
        print_failed()
    raise SystemExit(retcode)


@task
@consume_args
def run(args):
    """Run the package's main script. All arguments are passed to it."""
    # The main script expects to get the called executable's name as
    # argv[0]. However, paver doesn't provide that in args. Even if it did (or
    # we dove into sys.argv), it wouldn't be useful because it would be paver's
    # executable. So we just pass the package name in as the executable name,
    # since it's close enough. This should never be seen by an end user
    # installing through Setuptools anyway.
    from pathdiscov.main import main
    raise SystemExit(main([CODE_DIRECTORY] + args))


@task
def commit():
    """Commit only if all the tests pass."""
    if _test_all() == 0:
        subprocess.check_call(['git', 'commit'])
    else:
        print_failure_message('\nTests failed, not committing.')


@task
def coverage():
    """Run tests and show test coverage report."""
    try:
        import pytest_cov  # NOQA
    except ImportError:
        print_failure_message(
            'Install the pytest coverage plugin to use this task, '
            "i.e., `pip install pytest-cov'.")
        raise SystemExit(1)
    import pytest
    pytest.main(PYTEST_FLAGS + [
        '--cov', CODE_DIRECTORY,
        '--cov-report', 'term-missing',
        TESTS_DIRECTORY])


@task  # NOQA
def doc_watch():
    """Watch for changes in the docs and rebuild HTML docs when changed."""
    try:
        from watchdog.events import FileSystemEventHandler
        from watchdog.observers import Observer
    except ImportError:
        print_failure_message('Install the watchdog package to use this task, '
                              "i.e., `pip install watchdog'.")
        raise SystemExit(1)

    class RebuildDocsEventHandler(FileSystemEventHandler):

        def __init__(self, base_paths):
            self.base_paths = base_paths

        def dispatch(self, event):
            """Dispatches events to the appropriate methods.
            :param event: The event object representing the file system event.
            :type event: :class:`watchdog.events.FileSystemEvent`
            """
            for base_path in self.base_paths:
                if event.src_path.endswith(base_path):
                    super(RebuildDocsEventHandler, self).dispatch(event)
                    # We found one that matches. We're done.
                    return

        def on_modified(self, event):
            print_failure_message('Modification detected. Rebuilding docs.')
            # Strip off the path prefix.
            # import os
            # if event.src_path[len(os.getcwd()) + 1:].startswith(
            #         CODE_DIRECTORY):
            # sphinx-build doesn't always pick up changes on code files,
            # even though they are used to generate the documentation. As
            # a workaround, just clean before building.
            doc_html()
            print_success_message('Docs have been rebuilt.')

    print_success_message(
        'Watching for changes in project files, press Ctrl-C to cancel...')
    handler = RebuildDocsEventHandler(get_project_files())
    observer = Observer()
    observer.schedule(handler, path='.', recursive=True)
    observer.start()
    try:
        while True:
            time.sleep(1)
    except KeyboardInterrupt:
        observer.stop()
        observer.join()


@task
@needs('doc_html')
def doc_open():
    """Build the HTML docs and open them in a web browser."""
    doc_index = os.path.join(DOCS_DIRECTORY, 'build', 'html', 'index.html')
    if sys.platform == 'darwin':
        # Mac OS X
        subprocess.check_call(['open', doc_index])
    elif sys.platform == 'win32':
        # Windows
        subprocess.check_call(['start', doc_index], shell=True)
    elif sys.platform == 'linux2':
        # All freedesktop-compatible desktops
        subprocess.check_call(['xdg-open', doc_index])
    else:
        print_failure_message(
            "Unsupported platform. Please open `{0}' manually.".format(
                doc_index))


@task
def get_tasks():
    """Get all paver-defined tasks."""
    from paver.tasks import environment
    for task in environment.get_tasks():
        print(task.shortname)

@task
@needs('install_python_dependencies')
def doc_man():
    ''' Build man page '''
    retcode = _doc_make('man')

    if retcode:
        raise SystemExit(retcode)

@task
@needs('install_python_dependencies')
def doc_html():
    """Build the HTML docs."""
    retcode = _doc_make('html')

    if retcode:
        raise SystemExit(retcode)


@task
def doc_clean():
    """Clean (delete) the built docs."""
    retcode = _doc_make('clean')

    if retcode:
        raise SystemExit(retcode)
