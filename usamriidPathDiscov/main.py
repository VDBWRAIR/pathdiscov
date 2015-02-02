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
from helpers import runCommand
from pkg_resources import resource_filename
from os.path import join, expanduser, expandvars

options = helpers.get_options()
#logger_proxy, logging_mutex = helpers.make_logger(options, __file__)

basedir = os.path.relpath('./')
project_dir = options.outdir  # set output dir
R1 = os.path.abspath(options.R1)
R2 = "none"
if options.R2:
    R2 = os.path.abspath(options.R2)
else:
    R2 = "none"

# Do all initial setup
config = helpers.parse_config()
helpers.setup_shell_environment(config)
helpers.setup_param(config)

# Setup initial inputs
input = project_dir + "/" + "input"
logs = project_dir + "/" + "logs"
results = project_dir + "/" + "results"
sample1 = results + "/" + "sample1"
sample2 = results + "/" + "sample2"
paramFile = input + "/param.txt"

def report(result):
    """Wrapper around Result.report"""
    result.report(logger_proxy, logging_mutex)
    print result

@follows(mkdir(project_dir, input, results, logs))
@originate([paramFile])
def createPram(output_file):
    result = tasks.createParam(output_file)
    return result

#@graphviz(height=1.8, width=2, label="Prepare\nanalysis")

#define file
if options.R2:
    pram1= [
          [R1, join(project_dir, 'input', 'F.fastq')],
          [R2, join(project_dir, 'input', 'R.fastq')],
          ]
else:
    pram1= [
          [R1, join(project_dir, 'input', 'F.fastq')],
          ]


@follows(createPram)
@files(pram1)
def prepare_analysis(input, output):
    """copy the mapfile to analyiss dir

    Arguments:
    - `file_to_copy`: The mapfile from 454 sequencer
    - `file_copy`: a copy of mapfile
    """
    result = tasks.copy_map_file(input, output)
    return result

@follows(mkdir(results + "/quality_analysis"))
@transform(prepare_analysis, formatter("F.fastq", "R.fastq"), results + "/quality_analysis")
def fastQC(input, output):
    result = tasks.createQuality(input, output)
    return result

if R2!="none":
    param2 = [
           [[project_dir + "/input/F.fastq", project_dir + "/input/R.fastq"], results]
         ]
else:
    param2 = [
             [project_dir + "/input/F.fastq", results]
             ]
#@graphviz(height=1.8, width=2, label="Processing\nstep1")
@follows(prepare_analysis)
@files(param2)
#@transform(prepare_analysis, formatter("F.fastq", "R.fastq"), results)
def priStage(input, output):

    result = tasks.priStage(
            input, project_dir, paramFile, config['blast_unassembled'], output)
    return result

def verify_standard_stages_files(projectpath, templatedir):
    ''' Hardcoded verification of standard stages '''
    #templates.append(resource_filename(__name__, tfile))
    from verifyproject import STAGES, verify_project
    projectname = basename(projectpath)
    # Fetch template files for each stage from inside
    # of templatedir
    templates = []
    for stage in STAGES:
        tfile = os.path.join(templatedir,stage+'.lst')
    return verify_project(projectpath, projectname, templates)


def generateSymLink():
    final_out = 'results/output'
    final_out_link = os.path.abspath(project_dir + "/output")
    cmd = "ln -s %s  %s" %(final_out, final_out_link)
    runCommand(cmd, "T")

    final_out = "results/iterative_blast_phylo_1/reports"
    final_out_link = os.path.abspath(project_dir + "/contig_reports")
    cmd = "ln -s %s  %s" %(final_out, final_out_link)
    runCommand(cmd, "T")

    final_out = "results/iterative_blast_phylo_2/reports"
    final_out_link = os.path.abspath(project_dir + "/unassembledread_reports")
    cmd = "ln -s %s  %s" %(final_out, final_out_link)
    runCommand(cmd, "T")

    final_out = "results/step1/R1.count"
    final_out_link = os.path.abspath(project_dir + "/R1.count")
    cmd = "ln -s %s  %s" %(final_out, final_out_link)
    runCommand(cmd, "T")

    final_out = "results/step1/R2.count"
    final_out_link = os.path.abspath(project_dir + "/R2.count")
    cmd = "ln -s %s  %s" %(final_out, final_out_link)
    runCommand(cmd, "T")

    final_out = "results/quality_analysis"
    final_out_link = os.path.abspath(project_dir + "/quality_analysis")
    cmd = "ln -s %s  %s" %(final_out, final_out_link)
    runCommand(cmd, "T")

    final_out = "results/analysis.log"
    final_out_link = os.path.abspath(project_dir + "/analysis.log")
    cmd = "ln -s %s  %s" %(final_out, final_out_link)
    runCommand(cmd, "T")
    # check the process

    cmd = 'verifyproject %s' % project_dir
    import subprocess
    # We can output from this
    subprocess.Popen(cmd, shell=True).wait()
    # print time elapsed to complete the task
    return


def main():

    from helpers import which
    t0 = time.time()
    print (" Starting time ..... :") + str(t0)
    print "print default argument whether to generate default param.txt file ..." +  str(options.param)
    if options.param:
        dir_bak = project_dir + ".bak"
        try:
            if os.path.exists(project_dir):
                tasks.copyDir(project_dir,dir_bak )
                tasks.rmdir(dir)
        except:
            pass
        print "Task names: ", pipeline_get_task_names()

        print "....................." + basedir + "/" + project_dir
        pipeline_printout(sys.stdout, [createPram, prepare_analysis], verbose=6)
        pipeline_run(
            ['usamriidPathDiscov.main.createPram',
            'usamriidPathDiscov.main.prepare_analysis'], multiprocess=6)
        pipeline_get_task_names()  # return task names


    elif options.noparam is False:
        #if options.R2:

        print "....................." + basedir + "/" + project_dir
        pipeline_printout(sys.stdout, [createPram, prepare_analysis,fastQC,priStage], verbose=6)
        pipeline_run(
            ['usamriidPathDiscov.main.fastQC',
            #'usamriidPathDiscov.main.convertToPdf',
            'usamriidPathDiscov.main.priStage'], multiprocess=6)
        pipeline_get_task_names()  # return task names

        generateSymLink()

    else:
        print """ **********************************************************

            Processing Pathogen discovery ....
            *****************************************************************

        """
        dir_bak = project_dir + ".bak"
        try:
            if os.path.exists(project_dir):
                tasks.copyDir(project_dir,dir_bak )
                tasks.rmdir(dir)
        except:
            pass
        print "Task names: ", pipeline_get_task_names()
        print "....................." + basedir + "/" + project_dir
        pipeline_printout(sys.stdout, [createPram, prepare_analysis,fastQC,priStage], verbose=6)
        pipeline_run(
            ['usamriidPathDiscov.main.createPram',
            'usamriidPathDiscov.main.prepare_analysis',
            'usamriidPathDiscov.main.fastQC',
            #'usamriidPathDiscov.main.convertToPdf',
            'usamriidPathDiscov.main.priStage'], multiprocess=6)
        pipeline_get_task_names()
        generateSymLink()

    import datetime
    from termcolor import colored
    elapsedTime = int((time.time()) - t0)
    elapsedTime = str(datetime.timedelta(seconds=elapsedTime))
    print("Time to complete the task ....." ) + colored (elapsedTime, "red")


if __name__ == '__main__':
    main()
