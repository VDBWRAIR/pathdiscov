import os
import sys
import yaml
import shutil
from helpers import get_options, make_logger, runCommand, which, run_cmd
from pkg_resources import resource_filename, resource_stream
import subprocess
import os.path

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
    cmds = ' '.join(cmds)
    cmds += " 2>&1 | tee -a " + output + "/analysis_quality.log"
    runCommand(cmds, True)

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
    cmds = ' '.join(cmds)
    print cmds
    p = runCommand(cmds, False)
    return

def comment():
    return "*" * 80

#pathogen.pl --sample sample1 --command step1 --paramfile param.txt --outputdir ../results/sample1 --R1 /my/data/R1.fastq.gz --R2 /my/data/R2.fastq.gz > ../logs/out1.o 2> ../logs/out1.e &
#input=[pathDescTest/input/F.fastq, pathDescTest/input/R.fastq, pathDescTest, step1,pathDescTest/input/param.txt, pathDescTest/logs/out.o, pathDescTest/logs/out.e]
def priStage(input, project_dir, paramFile,numreads,sge, output):
    """run step  1 of pathogen discovery
    Arguments:
        -`input`: list of F and R fastq file
        -`output`: output directory
        -`pramfile`: parameter file name
        -`output`: outdir for the output
        -`logdir`: log dir
    """
    ffastq, rfastq = input
    cmds = [
        'run_standard_stable4.pl',
        '--sample', os.path.basename(project_dir),
        '--paramfile', paramFile,
        '--outputdir', output,
        '--R1', ffastq,
    ]
    if rfastq is not None:
        cmds += [
            '--R2', rfastq
        ]

    if sge:
        cmds += [
            '--SGE', '1'
        ]

    cmds += ['--blast_unassembled', str(numreads)]
    cmds = ' '.join(cmds)
    cmds += " 2>&1 | tee -a " + output + "/analysis.log"

    runCommand(cmds, True)
