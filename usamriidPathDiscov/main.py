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

print "Set up command line option handling, logger creation, and load config file"
options = helpers.get_options()
#logger_proxy, logging_mutex = helpers.make_logger(options, __file__)
config = yaml.load(open(os.path.join(os.path.dirname(__file__),'config.yaml')).read())
#print yaml.dump(config)
#THIS = os.path.dirname( os.path.abspath(__file__ ) )
#path = os.path.abspath(os.path.join(os.path.dirname(sys.argv[0])))
curdir = os.getcwd()
os.chdir(curdir)
# print mydir
basedir = os.path.relpath('./')
print basedir
project_dir = options.outdir  # set output dir
R1 = options.R1
R2 = options.R2
##################################################
#    Setup databases and few globals             #
#                                                #
##################################################

databases = config['databases']
human_dna = config['human_dna']
h_sapiens_rna = config['h_sapiens_rna']
nt_db = config['nt_db']
tax_nodes = config['tax_nodes']
tax_names = config['tax_names']
blast_unassembled = config['blast_unassembled']
##################################################
#    Seq to process                             #
#                                                #
###############################################
#R1 = config['R1']
#R2 = config['R2']
R1 = os.path.abspath(R1)
R2 = os.path.abspath(R2)
tasks.comment()
print "R1: ", R1
print "R2: ", R2
print "human_dna: ", human_dna
print "h_sapiens_rna: ", h_sapiens_rna
print "tax_nodes: ", tax_nodes
tasks.comment()

#p1 = re.match(r'.+/(.+).fastq', R1)
#pat1 = p1.group(1)
#print pat1

#p2 = re.match(r'.+/(.+).fastq', R2)
#pat2 = p2.group(1)
#print pat2
input = project_dir + "/" + "input"
logs = project_dir + "/" + "logs"
results = project_dir + "/" + "results"
sample1 = results + "/" + "sample1"
sample2 = results + "/" + "sample2"

pram = [
    [R1, project_dir + "/input/F.fastq"],
    [R2, project_dir + "/input/R.fastq"]
]
print "******************************", yaml.dump(pram)

F_fastq = os.path.abspath(project_dir + "/input/F.fastq")
R_fastq = os.path.abspath(project_dir + "/input/R.fastq")

# Copy base sample.param.base to sample.param file
# os.path.abspath("./usamriidPathDiscov/files/sample.param.base")
CWD = os.path.abspath(os.path.join(os.path.dirname(__file__)))
baseFile = os.path.join(CWD,'files','sample.param.base')
# os.path.abspath("./usamriidPathDiscov/files/sample.param")
sampleParam = os.path.join(CWD,'files',"sample.param")
tasks.copy_map_file(baseFile, sampleParam)
# replace some globals in the sample.param file such as db names
for line in fileinput.input(sampleParam, inplace=True, backup='.bak'):
    line = re.sub(r'HUMAN_DNA',  human_dna, line.rstrip())
    line = re.sub(r'H_SAPIENS_RNA',  h_sapiens_rna, line.rstrip())
    line = re.sub(r'BLAST_NT',  nt_db, line.rstrip())
    line = re.sub(r'TAX_NODES',  tax_nodes, line.rstrip())
    line = re.sub(r'TAX_NAMES',  tax_names, line.rstrip())

    print (line)


def report(result):
    """Wrapper around Result.report"""
    result.report(logger_proxy, logging_mutex)
    print result

print "....................." + basedir + "/" + project_dir


pathogenScript = distutils.spawn.find_executable("pathogen.pl")
param3 = [[pathogenScript, input + "/param.txt"]]
paramFile = input + "/param.txt"
print paramFile

print "project_dir ", project_dir
print "paramFile: ", paramFile


@follows(mkdir(project_dir, input, results, logs))
@files(param3)
def createPram(input, output):
    result = tasks.createParam(input, output)
    return result


#@graphviz(height=1.8, width=2, label="Prepare\nanalysis")
@follows(createPram)
@files(pram)
def prepare_analysis(input, output):
    """copy the mapfile to analyiss dir

Arguments:
    - `file_to_copy`: The mapfile from 454 sequencer
    - `file_copy`: a copy of mapfile
    """
    print "copy ", input, " to ", output
    result = tasks.copy_map_file(input, output)
    return result
print tasks.comment()
print F_fastq
print R_fastq
print project_dir
print paramFile
print blast_unassembled
print results
print tasks.comment()


@follows(mkdir(results + "/quality_analysis"))
@transform(prepare_analysis, formatter("F.fastq", "R.fastq"), results + "/quality_analysis")
def fastQC(input, output):
    result = tasks.createQuality(input, output)
    return result

#@transform(fastQC, formatter("F_fastqc.html", "R_fastqc.html"), ".pdf")
#def convertToPdf(input, output):
    #result = tasks.convertHtmlToPDF(input, output)
    #return result
param5 = [
       [[project_dir + "/input/F.fastq", project_dir + "/input/R.fastq"], results]
]


#@graphviz(height=1.8, width=2, label="Processing\nstep1")
@follows(prepare_analysis)
@files(param5)
#@transform(prepare_analysis, formatter("F.fastq", "R.fastq"), results)
def priStage(input, output):
    print "Running step 1..."
    result = tasks.stage1(
        input, project_dir, paramFile, blast_unassembled, output)
    return result


def main():
    t0 = time.time()
    print (" Starting time ..... :") + str(t0)
    dir_bak = project_dir + ".bak"
    try:
        if os.path.exists(project_dir):
            tasks.copyDir(project_dir,dir_bak )
            tasks.rmdir(dir)
    except:
        pass
    print "Task names: ", pipeline_get_task_names()
    #tasks_torun = [createPram, prepare_analysis,fastQC,convertToPdf, priStage]

    #pipeline_printout_graph(
        #'snake_eater.ps', 'ps', tasks_torun, user_colour_scheme={
            #"colour_scheme_index": 6}, no_key_legend=False,
        #pipeline_name="Pathogen Discovery", size=(11, 8), dpi=30,
        #draw_vertically=True, ignore_upstream_of_target=False)
    print "....................." + basedir + "/" + project_dir
    pipeline_printout(sys.stdout, [createPram, prepare_analysis,fastQC,priStage], verbose=6)
    pipeline_run(
        ['usamriidPathDiscov.main.createPram',
         'usamriidPathDiscov.main.prepare_analysis',
         'usamriidPathDiscov.main.fastQC',
         #'usamriidPathDiscov.main.convertToPdf',
         'usamriidPathDiscov.main.priStage'], multiprocess=6)
    pipeline_get_task_names()  # return task names

    #helpers.run(options)

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
    print("End time ....." ) + str((time.time()) - t0)


if __name__ == '__main__':
    main()
