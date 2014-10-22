# -*- coding: utf-8 -*-


from __future__ import print_function

import os
import sys
import time
import subprocess
from paver.easy import *
from os.path import islink,isfile
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
            sdir=path('usamriidPathDiscov/download/bwa'),
            bindir=path('usamriidPathDiscov/bin')
        ),
        minilib=Bunch(
            # extra_files=['doctools','virtual']
        ),
         samtools=Bunch(
              sdir=path('usamriidPathDiscov/download/samtools'),
              bindir=path('usamriidPathDiscov/bin')
          ),
        Ray=Bunch(
            src = path('usamriidPathDiscov/download'),
            sfile =path('usamriidPathDiscov/download/ray'),
            sfile2 =path('usamriidPathDiscov/download/Ray2'),
            olink =path('usamriidPathDiscov/bin')
        ),
        CAP3=Bunch(
            sfile = path('usamriidPathDiscov/download/CAP3/cap3'),
            olink =path('usamriidPathDiscov/bin')

        ),
        FASTQC=Bunch(
            sfile = path('usamriidPathDiscov/download/FastQC/fastqc'),
            olink =path('usamriidPathDiscov/bin')
        ),

        wkhtmltopdf=Bunch(
            sfile = path('usamriidPathDiscov/download/wkhtmltopdf'),
            olink =path('usamriidPathDiscov/bin')
        ),

        bowtie2=Bunch(
            sfile =path('usamriidPathDiscov/download/bowtie2/bowtie2*'),
            olink =path('usamriidPathDiscov/bin')

        ),
        getorf=Bunch(
            src=path('usamriidPathDiscov/download'),
            sfile =path('usamriidPathDiscov/download/EMBOSS-6.6.0'),
            olink =path('usamriidPathDiscov/bin')

        ),

        #"""
        #local_lib=Bunch(
            #sdir=path('usamriidPathDiscov/download/local-lib-2.000011'),
            #bindir=path('usamriidPathDiscov/bin')
        #),
        #"""
        settings=Bunch(
            shell_file=path('usamriidPathDiscov/files/settings.sh'),
            shell_file_bk=path('usamriidPathDiscov/files/settings.sh.base'),
            bash_rc =path('usamriidPathDiscov/files/bashrc'),
            config =path('config.yaml'),
            config_bk=path('config.yaml.base'),
            dist_dir =path('.'),
            param_base =path('usamriidPathDiscov/files/sample.param.base'),
            param_work =path('usamriidPathDiscov/files/sample.param')
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
   $ source usamriidPathDiscov/bin/activate
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
def install_bwa(options):
    """installs the current package"""
    info("Installing bwa ....\n")
    currwd = os.getcwd()
    sdir = path(currwd) / options.bwa.sdir
    bindir = path(currwd) / options.bwa.bindir
    sh('cd %s; make ; cp bwa %s; cd %s' % (sdir, bindir, sdir))


@task
def install_samtools(options):
    """installs the current package"""
    info("Installing bwa ....\n")
    currwd = os.getcwd()
    sdir = path(currwd) / options.samtools.sdir
    bindir = path(currwd) / options.samtools.bindir
    sh('cd %s; make ; cp samtools %s; cd %s' % (sdir, bindir, sdir))

@task
def refRay(options):
    """Install  Ray assembler """
    info ("Installing  `Ray`  and copy to the bin dir")
    currwd = os.getcwd()
    src = path(currwd) / options.Ray.src
    sfile = path(currwd) / options.Ray.sfile
    sfile2 = path(currwd) / options.Ray.sfile2
    olink = path(currwd) / options.Ray.olink
    if os.path.isdir("/usr/lib64/openmpi/lib"):
        info("install  `Ray` and copy to bin  ....")
        sh(' cd %s; tar -xzvf ray.tar.gz; tar -xzvf RayPlatform.tar.gz; export LD_LIBRARY_PATH=/usr/lib64/openmpi/lib:$LD_LIBRARY_PATH; cd %s;make PREFIX=build2000 MPICXX=/usr/lib64/openmpi/bin/mpicxx; cp Ray %s' % (src,sfile, olink))
        sh('cp %s %s ' %(sfile2, olink))
    elif os.path.isdir("/usr/lib/openmpi/lib"):
        info("install  `Ray` and copy to bin  ....")
        sh(' cd %s; tar -xzvf ray.tar.gz; tar -xzvf RayPlatform.tar.gz; export LD_LIBRARY_PATH=/usr/lib64/openmpi/lib:$LD_LIBRARY_PATH; cd %s;make PREFIX=build2000 MPICXX=/usr/lib64/openmpi/bin/mpicxx; cp Ray %s' % (src,sfile, olink))
        sh('cp %s %s ' %(sfile2, olink))

    else:
        info("Ray is not installed, ... please install `openmpi and openmpi-devel` and try again")
        sys.exit()
@task
def copyWkhtmltopdf(options):
    """copy cap3 to bin filder"""
    info("Copy cap3 binary to bin folder")
    sfile =options.wkhtmltopdf.sfile
    olink = options.wkhtmltopdf.olink
    if isfile(sfile):
        sh('cp %s %s' %(sfile, olink))

@task
def refCAP3(options):
    """copy cap3 to bin filder"""
    info("Copy cap3 binary to bin folder")
    sfile =options.CAP3.sfile
    olink = options.CAP3.olink
    if isfile(sfile):
        sh('cp %s %s' %(sfile, olink))

@task
def getorf(options):
    """Install  EMBOSS getorf """
    info("Installing `getorf` and copy the binary to bin folder")
    currwd = os.getcwd()
    src = path(currwd) / options.getorf.src
    sfile = path(currwd) / options.getorf.sfile
    olink = path(currwd) / options.getorf.olink
    binfile = sfile + "/bin/getorf"
    symlink = olink + "/getorf"
    info("install  `EMBOSS getorf` and copy to bin  ....")
    sh('cd %s; tar -xzvf EMBOSS-6.6.0.tar.gz; cd %s;./configure CC="cc"; ./configure --prefix=%s;make;make install' %(src, sfile, sfile))
    if islink(symlink):
        sh('unlink %s' %(symlink))
    if isfile(binfile):
        sh('ln -s  %s %s/getorf' %(binfile, olink))

@task
def bowtie2(options):
    """copy bowtie to bin filder"""
    info("Copy bowtie binary to bin folder")
    sfile =options.bowtie2.sfile
    olink = options.bowtie2.olink
    sh('cp %s %s' %(sfile, olink))
@task
def refFastQC(options):
    """copy fastqc to bin filder"""
    info("Copy fastqc binary to bin folder")
    sfile = os.path.abspath(options.FASTQC.sfile)
    olink = os.path.abspath(options.FASTQC.olink)
    if islink(olink + "/fastqc"):
        sh ('unlink %s/fastqc' %(olink))
    if isfile(sfile):
        sh('ln -s %s %s/fastqc' %(sfile, olink))
@task
def modifyBashRC():
    "Append the content of setting.sh to .bashrc"
    import fileinput
    import re
    import subprocess
    info("Append path info to your .bashrc ")
    sfile = os.path.abspath(options.settings.shell_file)
    sfilebk = os.path.abspath(options.settings.shell_file_bk)
    sh("cp %s %s" %(sfilebk, sfile))
    info(sfile)
    bashrc = os.path.expanduser("~/.bashrc")
    bashrcbk = os.path.expanduser("~/.bashrc.bak")
    appdir = os.getcwd()
    dbdir = os.path.expanduser("~/databases")
    info(sfile)
    bashrcTemp = os.path.abspath(options.settings.bash_rc)
    for line in fileinput.input(sfile, inplace=True, backup='.bak'):
            line = re.sub(r'CWDAPP',  appdir, line.rstrip())
            line = re.sub(r'GENOMEDIR',  dbdir, line.rstrip())
            info(line)
    info(bashrc)
    if isfile(bashrc):
        sh("cat %s > %s" %(bashrc, bashrcbk))
        cmd = "grep 'files/settings' %s" %(sfile)

        if not (subprocess.call(cmd, shell=True)): # check if the shell return code is "1" or "0"
            pass
        else:
            #sh("cat %s >> %s" %(sfile, bashrc))
            info("Source the path information ....")
            sh("echo 'source %s' >> %s" %(sfile, bashrc))
    else:
        info ("Creating a new bashrc file...")
        #sh("cat %s > %s" %(bashrcTemp, bashrc))
        sh("cat %s >> %s" %(bashrcTemp, bashrc))
        sh("echo 'source %s' >> %s" %(sfile, bashrc))

    info("Set config file ....")
    conf = os.path.abspath(options.settings.config)
    info(conf)
    conf_bk = os.path.abspath(options.settings.config_bk)
    info(conf_bk)
    sh("cat  %s >  %s" %(conf_bk, conf))
    for line in fileinput.input(conf, inplace=True, backup='.bak'):
        line = re.sub(r'GENOMEDIR',  dbdir, line.rstrip())
        info(line)
#@task
#def installPerlLocalLib():
     #"""installs the current local.lib  package"""
     #info("Installing bwa ....\n")
     #currwd = os.getcwd()
     #sdir = path(currwd) / options.local_lib.sdir
     #sh("cd %s; perl Makefile.PL --bootstrap ;make test && make install;mkdir -p ~/perl5/bin" % (sdir))
     #sh("echo [ $SHLVL -eq 1 ] &&  eval '$(perl -I$HOME/perl5/lib/perl5 -Mlocal::lib)' >>~/.bashrc")

@task
def setConfig():
    import glob
    dist_dir = os.path.abspath(options.settings.dist_dir)
    info(dist_dir)
    install_dir =glob.glob(dist_dir +"/usamriidPathDiscov/lib/python*/site-packages/usamriidPathDiscov")
    install_dir =install_dir[0]
    info(install_dir)
    conf = os.path.abspath(options.settings.config)
    sh("cp %s %s" %(conf, install_dir))

@task
def setParam():
    import glob
    dist_dir = os.path.abspath(options.settings.dist_dir)
    info(dist_dir)
    install_dir =glob.glob(dist_dir +"/usamriidPathDiscov/lib/python*/site-packages/usamriidPathDiscov")
    install_dir =install_dir[0]
    info(install_dir)
    sh("mkdir -p %s/files" %(install_dir))
    install_dir = install_dir + "/files"
    info(install_dir)
    param_base = os.path.abspath(options.settings.param_base)
    param_work = os.path.abspath(options.settings.param_work)
    sh("cp %s %s" %(param_base, install_dir))
    sh("cp %s %s" %(param_work, install_dir))

@task
def install_dependencies():
    sh('pip install  -r requirements-dev.txt ')


@task
@needs('install_dependencies', 'source_shell', 'install_bwa', 'install_samtools','refRay','bowtie2','refCAP3' ,'refFastQC','getorf','copyWkhtmltopdf','modifyBashRC','setConfig', 'setParam')
def prepare():
    """Prepare complete environment
    """
    pass

@task
@needs('doc_html', "minilib", 'prepare', 'setuptools.command.sdist')
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
    from usamriidPatDescov.main import main
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
