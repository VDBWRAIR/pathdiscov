import os
import sys
import yaml
import shutil
from helpers import get_options, make_logger, runCommand, which, run_cmd
from pkg_resources import resource_filename

options = get_options()
#logger_proxy, logging_mutex = make_logger(options, __file__)
#config = yaml.load(open(options.config).read())
print 'sys.argv[0] =', sys.argv[0]
config_file = resource_filename(__name__, 'files/config.yaml')
config = yaml.load(open(config_file).read())
#print yaml.dump(config)
#logger_proxy, logging_mutex = helpers.make_logger(options, __file__)
blast_unassembled =config['blast_unassembled']
def work_dir(dirname):
    if not dirname:
        currpath = os.getcwd()
        dir_to_create = currpath + "/" + dirname
        print dir_to_create
        os.mkdir(dir_to_create, 0755)
        return


def rmdir(dirname):
    if os.path.exists(dirname):
        #shutil.rmtree(dirname, ignore_errors=True)
        shutil.rmtree(dirname)
        return


def mkdir(dirname):
    if not os.path.exists(dirname):
        os.mkdir(dirname)
        return


def copy_map_file(file_to_copy, file_copy):
    """copy the mapfile to analyiss dir

 Arguments:
     - `file_to_copy`: The mapfile from 454 sequencer
     - `file_copy`: a copy of mapfile
     """

    shutil.copyfile(file_to_copy, file_copy)
    return
def copyDir(file_to_copy, file_copy):
    """copy the mapfile to analyiss dir

 Arguments:
     - `file_to_copy`: The mapfile from 454 sequencer
     - `file_copy`: a copy of mapfile
     """

    shutil.copytree(file_to_copy, file_copy)
    return


def createParam(input,output):
    cmd = '%s --example  > %s' %(input,output)
    p=runCommand(cmd, "F")
    return p

def createQuality(input,output):
    """Check quality of the two fastq files
    Arguments:
          -`input`: The two fastq files
          -`output`: The output folder for the analysis
    """
    cmds = [
        'fastqc',
        input,
        '-o', output,
    ]
    cmds = '  '.join(cmds)
    cmds += "| tee  -a "  + output + "/analysis_quality.log"
    print cmds
    p = runCommand(cmds, "F")
    return

def convertHtmlToPDF(input,output):
    """Check quality of the two fastq files
    Arguments:
          -`input`: The two fastq files
          -`output`: The output folder for the analysis
    """
    cmds = [
        'wkhtmltopdf',
        input,
        output,
    ]
    cmds = '  '.join(cmds)
    print cmds
    p = runCommand(cmds, "F")
    return

def comment():
    return "*" * 80
#pathogen.pl --sample sample1 --command step1 --paramfile param.txt --outputdir ../results/sample1 --R1 /my/data/R1.fastq.gz --R2 /my/data/R2.fastq.gz > ../logs/out1.o 2> ../logs/out1.e &
#input=[pathDescTest/input/F.fastq, pathDescTest/input/R.fastq, pathDescTest, step1,pathDescTest/input/param.txt, pathDescTest/logs/out.o, pathDescTest/logs/out.e]
def stage1(input, project_dir, paramFile,numreads, output):
    """run step  1 of pathogen discovery
    Arguments:
        -`input`: list of F and R fastq file
        -`output`: output directory
        -`pramfile`: parameter file name
        -`output`: outdir for the output
        -`logdir`: log dir
    """
    #numreads = blast_unassembled * 4
    #cmds = ['pathogen.pl',
            #'--sample', project_dir,
            #'--command', "step1 quality_filter  host_map ray2_assembly",
            #'--paramfile', paramFile,
            #'--outputdir', output,
            #'--R1', F_fastq,
            #'--R2', R_fastq,
            ##'>', output + "/step1/logs/out1.o",
            ##'2>', output + "/step1/logs/out1.e",

            #]
    ffastq,rfastq = input
    cmds = [
        'run_standard_stable4.pl',
        '--sample', project_dir,
        '--paramfile', paramFile,
        '--outputdir', output,
        '--R1', ffastq,
        '--R2', rfastq,
        '--blast_unassembled', str(numreads),
    ]
    cmds = '  '.join(cmds)
    cmds += "| tee  -a "  + output + "/analysis.log"
    print cmds
    p = runCommand(cmds, "F")
    return

