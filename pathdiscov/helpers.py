import os
import sys
import time
import datetime
import logging
import logging.handlers
from argparse import ArgumentParser
import StringIO
import subprocess
from ruffus.proxy_logger import *
from ruffus import *
from pkg_resources import resource_filename
import yaml
from os.path import *
import shutil
import fileinput
import re
import multiprocessing
#import ast


class Result(object):

    def __init__(self, infiles, outfiles, log=None, stdout=None, stderr=None,
                 desc=None, failed=False, cmds=None):
        """
        A Result object encapsulates the information going to and from a task.
        Each task is responsible for determining the following arguments, based
        on the idiosyncracies of that particular task:

            `infiles`: Required list of input files to the task
            `outfiles`: Required list of output files to the task
            `failed`: Optional (but recommended) boolean indicating whether the
                      task failed or not.  Default is False, implying that the
                      task will always work
            `log`: Optional log file from running the task's commands
            `stdout`: Optional stdout that will be printed in the report if the
                      task failed
            `stderr`: Optional stderr that will be printed in the report if the
                      task failed
            `desc`: Optional description of the task
            `cmds`: Optional string of commands used by the task. Can be very
                    useful for debugging.
        """
        if isinstance(infiles, basestring):
            infiles = [infiles]
        if isinstance(outfiles, basestring):
            outfiles = [outfiles]
        self.infiles = infiles
        self.outfiles = outfiles
        self.log = log
        self.stdout = stdout
        self.stderr = stderr
        self.elapsed = None
        self.failed = failed
        self.desc = desc
        self.cmds = cmds

    def report(self, logger_proxy, logging_mutex):
        """
        Prints a nice report
        """
        with logging_mutex:
            if not self.desc:
                self.desc = ""
            logger_proxy.info(' Task: %s' % self.desc)
            logger_proxy.info('     Time: %s' % datetime.datetime.now())
            if self.elapsed is not None:
                logger_proxy.info('     Elapsed: %s' % nicetime(self.elapsed))
            if self.cmds is not None:
                logger_proxy.debug('     Commands: %s' % str(self.cmds))
            for output_fn in self.outfiles:
                output_fn = os.path.normpath(os.path.relpath(output_fn))
                logger_proxy.info('     Output:   %s' % output_fn)
            if self.log is not None:
                logger_proxy.info('     Log:      %s' % self.log)
            if self.failed:
                logger_proxy.error('=' * 80)
                logger_proxy.error('Error in %s' % self.desc)
                if self.cmds:
                    logger_proxy.error(str(self.cmds))
                if self.stderr:
                    logger_proxy.error('====STDERR====')
                    logger_proxy.error(self.stderr)
                if self.stdout:
                    logger_proxy.error('====STDOUT====')
                    logger_proxy.error(self.stdout)
                if self.log is not None:
                    logger_proxy.error('   Log: %s' % self.log)
                logger_proxy.error('=' * 80)
                sys.exit(1)
            logger_proxy.info('')


def nicetime(seconds):
    """Convert seconds to hours-minutes-seconds"""
    m, s = divmod(seconds, 60)
    h, m = divmod(m, 60)
    elapsed = "%dh%02dm%02.2fs" % (h, m, s)
    return elapsed


def timeit(func):
    """Decorator to time a single run of a task"""
    def wrapper(*arg, **kw):
        t0 = time.time()
        res = func(*arg, **kw)
        t1 = time.time()
        res.elapsed = (t1 - t0)
        if not res.desc:
            res.desc = func.__name__
        return res
    return wrapper


def run(options):
    """
    Run the pipeline according to the provided options (which were probably
    created by get_options())
    """
    if options.just_print:
        pipeline_printout(sys.stdout,
                          options.target_tasks,
                          options.forced_tasks,
                          verbose=options.verbose)

    elif options.flowchart:
        pipeline_printout_graph(open(options.flowchart, "w"),
                                os.path.splitext(options.flowchart)[1][1:],
                                options.target_tasks,
                                options.forced_tasks,
                                size=(11, 8),
                                dpi = 40,
                                user_colour_scheme = {
                                    "colour_scheme_index": 6},
                                no_key_legend=not options.key_legend_in_graph)
    else:
        pipeline_run(options.target_tasks,
                     options.forced_tasks,
                     multiprocess=options.jobs,
                     logger=stderr_logger,
                     touch_files_only=False,
                     verbose=options.verbose)


def get_options():
    """
    Standard set of options to be passed to the command line.  These can be in
    turn passed to run() to actually run the pipeline.
    """
    config = parse_config()
    parser = ArgumentParser()
    '''
    parser.add_argument("-v", "--verbose", dest="verbose",
                        action="count", default=0,
                        help="Print more verbose messages for each "
                        "additional verbose level.")
    parser.add_argument("-L", "--log_file", dest="log_file",
                        metavar="FILE",
                        type=str,
                        help="Name and path of log file")
    parser.add_argument("-t", "--target_tasks", dest="target_tasks",
                        action="append",
                        default=list(),
                        metavar="JOBNAME",
                        type=str,
                        help="Target task(s) of pipeline.")
    parser.add_argument("-n", "--just_print", dest="just_print",
                        action="store_true", default=False,
                        help="Don't actually run any commands; just print "
                             "the pipeline.")
    parser.add_argument("--flowchart", dest="flowchart",
                        metavar="FILE",
                        type=str,
                        help="Don't actually run any commands; just print "
                             "the pipeline as a flowchart.")
    parser.add_argument("--key_legend_in_graph",
                        dest="key_legend_in_graph",
                        action="store_true",
                        default=False,
                        help="Print out legend and key for dependency "
                             "graph.")
    parser.add_argument("--forced_tasks", dest="forced_tasks",
                        action="append",
                        default=list(),
                        metavar="JOBNAME",
                        type=str,
                        help="Pipeline task(s) which will be included "
                             "even if they are up to date.")
    '''
    parser.add_argument(
        "--use-sge",
        dest='sge',
        default=False,
        action="store_true",
        help="Flag to use sge[Default: %(default)s"
    )
    parser.add_argument(
        '--outdir',
        required=True,
        help='output directory'
    )
    parser.add_argument(
        '--R1',
        '-R1',
        required=True,
        help="Path to forward fastq file"
    )
    parser.add_argument(
        '--R2',
        '-R2',
        default=None,
        help ="Path to reverse fastq file"
    )
    parser.add_argument(
        '--param',
        action='store_true',
        help="Generate sample param.txt file and edit after generating " \
                "directory tree"
    )
    parser.add_argument(
        '--noparam',
        action='store_false',
        help="Use modified param.txt file after running with --param"
    )
    parser.add_argument(
        '-c',
        '--cpuNum',
        dest="cpuNum",
        default=config['NODE_NUM'],
        type=int,
        help="Number of CPU to use[Default: %(default)s]"
    )
    parser.add_argument(
        '--blast-unassembled',
        '--blast_unassembled',
        type=int,
        default=config['blast_unassembled'],
        help="Number of reads to use for all iterative_blast_phlo " \
                "steps[Default: %(default)s]" 
    )

    return parser.parse_args()

def runCommand(cmd, printCMD):
    import subprocess
    if printCMD == 'T' or printCMD is True:
        print cmd
    p = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True
    )
    sout, serr = p.communicate()
    return sout,serr

def run_cmd(cmd_str):
    """
    Throw exception if run command fails
    """
    process = subprocess.Popen(cmd_str, stdout = subprocess.PIPE,
                                stderr = subprocess.PIPE, shell = True)
    stdout_str, stderr_str = process.communicate()
    if process.returncode != 0:
        raise Exception("Failed to run '%s'\n%s%sNon-zero exit status %s" %
                            (cmd_str, stdout_str, stderr_str, process.returncode))

def which(program):
    """Check if the excutable exists in the users path
    """
    import os

    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
        else:
            for path in os.environ["PATH"].split(os.pathsep):
                path = path.strip('"')
                exe_file = os.path.join(path, program)
                if is_exe(exe_file):
                    return exe_file

    return None


def make_logger(options, file_name):
    """
    Sets up logging for the pipeline so that messages are synchronized across
    multiple processes
    """

    MESSAGE = 15
    logging.addLevelName(MESSAGE, "MESSAGE")

    def setup_std_logging(logger, log_file, verbose):
        """
        set up logging using programme options
        """
        class debug_filter(logging.Filter):

            """
            Ignore INFO messages
            """

            def filter(self, record):
                return logging.INFO != record.levelno

        class NullHandler(logging.Handler):

            """
            for when there is no logging
            """

            def emit(self, record):
                pass

        # We are interested in all messages
        logger.setLevel(logging.DEBUG)
        has_handler = False

        # log to file if that is specified
        if log_file:
            handler = logging.FileHandler(log_file, delay=False)
            handler.setFormatter(
                logging.Formatter(
                    "%(asctime)s - %(name)s - %(levelname)6s - %(message)s"))
            handler.setLevel(MESSAGE)
            logger.addHandler(handler)
            has_handler = True

        # log to stderr if verbose
        if verbose:
            stderrhandler = logging.StreamHandler(sys.stderr)
            stderrhandler.setFormatter(
                logging.Formatter("    %(message)s")
            )
            stderrhandler.setLevel(logging.DEBUG)
            if log_file:
                stderrhandler.addFilter(debug_filter())
            logger.addHandler(stderrhandler)
            has_handler = True

        # no logging
        if not has_handler:
            logger.addHandler(NullHandler())

    # set up log
    logger_name = os.path.splitext(os.path.basename(file_name))[0]
    logger = logging.getLogger(logger_name)
    setup_std_logging(logger, options.log_file, options.verbose)

    # Allow logging across Ruffus pipeline
    def get_logger(logger_name, args):
        return logger

    (logger_proxy,
     logging_mutex) = make_shared_logger_and_proxy(get_logger,
                                                   logger_name,
                                                   {})

    return logger_proxy, logging_mutex

def parse_config():
    '''
    Parse config file from installation and do anything
    necessary to it such as joining database path to other paths
    and expanding ~ and $
    '''
    config_file = resource_filename(__name__, 'files/config.yaml')
    config = yaml.load(open(config_file).read())

    config['PHRED_OFFSET'] = str(config['PHRED_OFFSET'])
    config['NODE_NUM'] = str(config['NODE_NUM'])
    databases = expanduser(expandvars(config['databases']))
    config['databases'] = databases
    config['host_dna'] = join(databases, config['host_dna'])
    host_rna = config.get('host_rna','')
    if host_rna:
        config['host_rna'] = join(databases, host_rna)
    config['nt_db']= join(databases, config['nt_db'])
    if config['diamond_db']:
        config['diamond_db'] = join(databases, config['diamond_db'])
    config['tax_nodes'] = join(databases, config['tax_nodes'])
    config['tax_names'] = join(databases, config['tax_names'])
    config['blast_unassembled'] = str(config['blast_unassembled'])
    if config['temp_dir']:
        config['temp_dir'] = config['temp_dir']
    if config['nr_db']:
        config['nr_db'] = join(databases, config['nr_db'])
    return config

def setup_shell_environment(config):
    '''
    Set all the shell variables using config values

    Replaces need for user or script to source settings.sh
    '''
    # This will be wherever python setup.py install installs to which
    # is probably pathdiscov/
    installdir = sys.prefix

    # TODO: Investigate if this is actually needed
    os.environ['INNO_PHRED_OFFSET'] = config['PHRED_OFFSET']
    os.environ['INNO_SEQUENCE_PLATFORM'] = config['SEQUENCE_PLATFORM']
    os.environ['INNO_NODE_NUM'] = config['NODE_NUM']
    os.environ['INNO_BOWTIE_HUMAN_GENOME_DB'] = config['host_dna']
    os.environ['INNO_BOWTIE_HUMAN_TRAN_DB'] = config['host_rna']
    os.environ['INNO_BLAST_NT_DB'] = config['nt_db']
    os.environ['INNO_TAX_NODES'] = config['tax_nodes']
    os.environ['INNO_TAX_NAMES'] = config['tax_names']

    os.environ['INNO_SCRIPTS_PATH'] = installdir
    os.environ['PERL5LIB'] = os.path.join(installdir, 'Local_Module')
    os.environ['R_LIBS'] = os.path.join(installdir, 'scripts')
    # Set LD_LIBRARY_PATH
    # TODO: This needs to be better written as it may be 32 bit arch
    if 'LD_LIBRARY_PATH' not in os.environ:
        os.environ['LD_LIBRARY_PATH'] = '/usr/lib64/openmpi/lib'
    else:
        os.environ['LD_LIBRARY_PATH'] += os.pathsep +  '/usr/lib64/openmpi/lib'
    # Set PATH
    # TODO: This needs to be better written as these paths are assumed and
    # really scripts path should not be necessary
    os.environ['PATH'] = installdir + os.pathsep + \
        os.path.join('/usr','lib64','openmpi','bin') + os.pathsep + \
        os.path.join(installdir,'bin') + os.pathsep + \
        os.path.join(installdir,'scripts') + \
        os.pathsep + os.path.join(installdir,'step1') + \
        os.pathsep + os.environ['PATH']

def setup_param(config, paramoutputfile):
    '''
    Setup sample.param file in installation directory

    TODO:
    Not sure why we are doing this every time main is called though
    should only need to be done once during installation

    Uses parsed config.yaml
    '''
    # Copy base sample.param.base to sample.param file
    baseFile = resource_filename(__name__, 'files/sample.param.base')
    shutil.copyfile(baseFile, paramoutputfile)
    # replace some globals in the sample.param file such as db names
    for line in fileinput.input(paramoutputfile, inplace=True):
        line = re.sub(r'SEQPLATFORM',config['SEQUENCE_PLATFORM'], line.rstrip() )
        line = re.sub('NUMINST', config['NODE_NUM'], line.rstrip() )
        line = re.sub(r'HOST_DNA', config['host_dna'], line.rstrip())
        line = re.sub(r'HOST_RNA', config['host_rna'], line.rstrip())
        line = re.sub(r'BLAST_NT', config['nt_db'], line.rstrip())
        line = re.sub(r'TAX_NODES', config['tax_nodes'], line.rstrip())
        line = re.sub(r'TAX_NAMES', config['tax_names'], line.rstrip())
        if config['diamond_db']:
            line = re.sub(r'DIAMOND_NR', config['diamond_db'], line.rstrip())
        if config['temp_dir']:
            line = re.sub(r'TEMPDIR', config['temp_dir'], line.rstrip())
        if config['nr_db']:
            line = re.sub(r'BLASTNR', config['nr_db'], line.rstrip())
        print (line)


def isGzip(input):
    if input.endswith(".gz"):
        return True
    else:
        return False

def symlink(src, dst):
    '''
    Create symlink src -> dst
    Overwrite dst if exists
    Do not create if src does not exist
    '''
    rel = relpath(src, dirname(dst))
    if exists(dst):
        os.unlink(dst)
    if not exists(src):
        return
    os.symlink(rel, dst)

def create_new_project(projpath):
    '''
    Create a new project(just empty dir)
    and/or backup existing project if exists
    '''
    projbk = projpath + '.bk'
    if exists(projpath):
        if exists(projbk):
            shutil.rmtree(projbk)
        os.rename(projpath, projbk)
    os.makedirs(projpath)
