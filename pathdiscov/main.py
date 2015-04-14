#!/usr/bin/env python
import sys
from ruffus import *
import yaml
import tasks
import helpers
import time
import os
import re
import distutils.spawn
import fileinput
from helpers import runCommand,isGzip
from pkg_resources import resource_filename
from os.path import (
    join, expanduser, expandvars,
    splitext, basename, dirname, exists
)
is_64bits = sys.maxsize > 2**32
if not is_64bits:
    print "Please upgrade your operating system to 64 bit, application such as diamond don't run on 32 bit"
    sys.exit(0)
options = helpers.get_options()
#logger_proxy, logging_mutex = helpers.make_logger(options, __file__)

basedir = os.path.relpath('./')
project_dir = options.outdir  # set output dir
pro_dir= os.path.basename(project_dir)
sge = options.sge
R1 = os.path.abspath(options.R1)
R2 = options.R2
if R2:
    R2 = os.path.abspath(options.R2)
blast_unassembled = options.blast_unassembled
# Do all initial setup
config = helpers.parse_config()
helpers.setup_shell_environment(config)

# Setup initial inputs
input_dir = join(project_dir, "input")
logs_dir = join(project_dir, "logs")
results_dir = join(project_dir, "results")
output_dir = join(project_dir, 'output')
quality_analysis_dir = join(results_dir, 'quality_analysis')
sample1 = join(results_dir, "sample1")
sample2 = join(results_dir, "sample2")
paramFile = join(input_dir, "param.txt")

def report(result):
    """Wrapper around Result.report"""
    result.report(logger_proxy, logging_mutex)
    print result


@follows(mkdir(project_dir, input_dir, results_dir, logs_dir))
@originate([paramFile])
def createPram(output_file):
    '''
    Create param.txt inside projectdir/input
    '''
    config['NODE_NUM'] = str(options.cpuNum)
    helpers.setup_param(config, output_file)

def qa_outfile(input):
    ''' Generate quality_analysis output file for input file '''
    if input is not None:
        return join(
            quality_analysis_dir,
            splitext(basename(input).replace('.gz',''))[0] + '_fastqc.zip'
        )
    else:
        return None

@follows(createPram)
@follows(mkdir(quality_analysis_dir))
@files([
    [R1, qa_outfile(R1)],
    [R2, qa_outfile(R2)],
])
def fastQC(input, output):
    '''
    Runs fastqc on both input files(or only R1 if only -R1)
    '''
    if input is not None:
        tasks.createQuality(input, dirname(output))

@follows(fastQC)
@files(
    [R1, R2],
    [join(results_dir,'analysis.log')]
)
def priStage(input, output):
    '''
    Run run_standard.pl with all supplied options
    '''
    result = tasks.priStage(
        input, project_dir, paramFile,
        blast_unassembled, sge, results_dir
    )
    return result

def verify_standard_stages_files(projectpath):
    ''' Hardcoded verification of standard stages '''
    #templates.append(resource_filename(__name__, tfile))
    from verifyproject import STAGES, verify_standard_stages_files
    projectname = basename(projectpath)
    import yaml
    from termcolor import colored
    templatesdir = resource_filename(__name__, 'output_files_templates')
    missingfiles = verify_standard_stages_files(projectpath, templatesdir)
    if missingfiles:
        for path, reason in missingfiles:
            # Don't care about 0 sized files
            if reason == 'Size zero':
                continue
            fname = basename(path)
            if fname == "param.txt":
                print colored("WARNING! :  Unable to create param.txt under the inpute directory", "red")
                sys.exit(1)
            elif  fname == "quality_filter.R1":
                print colored("WARNING! :  Unable to run quality filter step, please check if prinseq is installed and running or the log file under host_map if bowtie2 succeeded", "red")
                sys.exit(1)
            elif fname == "out.bam":
                print colored("WARNING! :  Unable to map the read to the ref genome, please check if bowtie2 is installed or the ref ~/databases exist", "red")
                sys.exit(1)
            elif fname == "out.cap.fa":
                print colored("WARNING! : Unable find out.cap.fa file, please check if cap3 program is running and  ray2_assembly error log", "red")
                sys.exit(1)
            elif fname == "out.ray.fa":
                print colored("WARNING! : Unable to run Ray assembly, please check if Ray2 program is running and possible error in the log file", "red")
                sys.exit(1)
            elif fname == "1.contig.blast":
                print colored("WARNING! : Unable to run iterative_blast_phylo_1, please check the log file under iterative_blast_phylo_1", "red")
                sys.exit(1)
            elif fname == "1.contig.fasta":
                print colored("WARNING! : Unable to run iterative_blast_phylo_1, please check the log file under iterative_blast_phylo_2", "red")
                sys.exit(1)
            elif fname == "contig." + pro_dir + ".top.smallreport.txt":
                print fname
                print colored("WARNING! : Unable to run iterative_blast_phylo , please check the log file under iterative_blast_phylo_1 or iterative_blast_phylo_2", "red")
                sys.exit(1)
    else:
        print colored(" SUCCESS! the tasks completed successfully", "green")

symlink_files = [
    [join(results_dir,'output'), join(project_dir, 'output')],
    [join(results_dir,'iterative_blast_phylo_1','reports'), join(project_dir,'contig_reports')],
    [join(results_dir,'iterative_blast_phylo_2','reports'), join(project_dir,'unassembled_read_reports')],
    [join(results_dir,'step1','R1.count'), join(project_dir, 'R1.count')],
    [join(results_dir,'quality_analysis'), join(project_dir, 'quality_analysis')],
    [join(results_dir,'analysis.log'), join(project_dir, 'analysis.log')]
]
if R2:
    symlink_files += [[join(results_dir,'step1','R2.count'), join(project_dir, 'R2.count')]]

@follows(priStage)
@files(symlink_files)
def symlink(src, dst):
    '''
    Create symlink src -> dst
    Overwrite dst if exists
    Do not create if src does not exist
    '''
    helpers.symlink(src, dst)

def commands(commandlist, index):
    '''
    Just return the correct commands from commandlist based on index given

    :param list commandlist: [(commandpath, commandfunc), ...]
    :param int index: Index of item in each tuple to extract
    '''
    return map(lambda x: x[index], commandlist)

def main():
    from helpers import which
    t0 = time.time()
    print (" Starting time ..... :") + str(t0)

    # Will be all the commands to run
    pipeline_commands = [
        (__name__ + '.createPram', createPram),
    ]

    print "....................." + basedir + "/" + project_dir

    if options.param:
        helpers.create_new_project(project_dir)
    elif not options.noparam:
        pipeline_commands = [
            (__name__ + '.fastQC', fastQC),
            (__name__ + '.priStage', priStage),
            (__name__ + '.symlink', symlink)
        ]
    else:
        helpers.create_new_project(project_dir)
        pipeline_commands += [
            (__name__ + '.fastQC', fastQC),
            (__name__ + '.priStage', priStage),
            (__name__ + '.symlink', symlink),
        ]

    pipeline_printout(sys.stdout, commands(pipeline_commands, 1), verbose=6)
    pipeline_run(commands(pipeline_commands, 0), multiprocess=6)
    pipeline_get_task_names()

    import datetime
    from termcolor import colored
    if not options.param:
        verify_standard_stages_files(project_dir)
    elapsedTime = int((time.time()) - t0)
    elapsedTime = str(datetime.timedelta(seconds=elapsedTime))
    print("Time to complete the task ....." ) + colored (elapsedTime, "red")
