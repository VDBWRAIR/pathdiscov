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

print "Set up command line option handling, logger creation, and load config file"
options = helpers.get_options()
#logger_proxy, logging_mutex = helpers.make_logger(options, __file__)
config_file = resource_filename(__name__, 'files/config.yaml')
config = yaml.load(open(config_file).read())

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
#    Setup environ vars                          #
# All come from old settings.sh                  #
# Effectively replaces the need to source        #
#  settings.sh
##################################################
# This will be wherever python setup.py install installs to which
installdir = sys.prefix

os.environ['INNO_PHRED_OFFSET'] = '33'
os.environ['INNO_SEQUENCE_PLATFORM'] = 'illumina'	# choices are: illumina 454
os.environ['INNO_NODE_NUM'] = '10'
os.environ['INNO_BOWTIE_HUMAN_GENOME_DB'] = human_dna
os.environ['INNO_BOWTIE_HUMAN_TRAN_DB'] = h_sapiens_rna
os.environ['INNO_BLAST_NT_DB'] = nt_db
os.environ['INNO_TAX_NODES'] = tax_nodes
os.environ['INNO_TAX_NAMES'] = tax_names

os.environ['INNO_SCRIPTS_PATH'] = installdir
os.environ['PERL5LIB'] = os.path.join(installdir, 'Local_Module')
os.environ['R_LIBS'] = os.path.join(installdir, 'scripts')
# Set LD_LIBRARY_PATH
if 'LD_LIBRARY_PATH' not in os.environ:
    os.environ['LD_LIBRARY_PATH'] = '/usr/lib64/openmpi/lib'
else:
    os.environ['LD_LIBRARY_PATH'] += os.pathsep +  '/usr/lib64/openmpi/lib'
# Set PATH
os.environ['PATH'] = installdir + os.pathsep + \
    os.path.join(installdir,'bin') + os.pathsep + \
    os.path.join(installdir,'scripts') + \
    os.pathsep + os.path.join(installdir,'step1') + \
    os.pathsep + os.environ['PATH']

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
baseFile = resource_filename(__name__, 'files/sample.param.base')
sampleParam = baseFile.replace('.base','')
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

def main():
    from helpers import which
    print which('pathogen.pl')
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

    from verifyproject import verify_standard_stages_files
    from pprint import pprint
    templatesdir = resource_filename(__name__, 'output_files_templates')
    missingfiles = verify_standard_stages_files(project_dir, templatesdir)
    if missingfiles:
        pprint(missingfiles)

    print("End time ....." ) + str((time.time()) - t0)


if __name__ == '__main__':
    main()
