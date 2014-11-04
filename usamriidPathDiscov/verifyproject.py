import string
import os
from os.path import exists, basename

# Default stages
# Each has it's own template under output_files_templates
stages = [
    'quality_analysis',
    'step1',
    'host_map_1',
    'quality_filter',
    'ray2_assembly_1',
    'orf_filter',
    'iterative_blast_phylo_1',
    'iterative_blast_phylo_2',
]

def verify_project(projectpath, projname, templates):
    '''
    Verify project has all files listed in all templates
    templates is a list of template paths
    '''
    # Compiled failed list
    failed = []
    for template in templates:
        failed += verify_files(projectpath, template)
    return failed

def filestemplate_to_listing(projpath, projname, template):
    '''
    Convert string template to list of file/directory paths
    '''
    return string.Template(template).substitute(
        projpath=projpath, projname=projname
    ).splitlines()

def verify_files(projectpath, files_template):
    '''
    Use a files_template to generate a listing of expected files
    $projpath will be replaced with projectpath
    $projname will be replaced with basename(projpath)

    Ensure each file is present and stat(file).st_size > 0
    Ensure each directory is present
    Symbolic links will be checked for existance and that they point to valid files

    return an empty list if all is well otherwise return list of tuples of
        (filedirpath, 'Reason why did not pass')
    '''
    # Get the basename(last thing after the last /)
    projname = basename(projectpath)
    
    # Make sure template file actually exists
    template = None
    try:
        # Read contents into template variable
        with open(files_template) as fh:
            template = fh.read()
    except IOError as e:
        raise ValueError('{0} does not exist'.format(files_template))
    
    # Get the template listing
    expected_listing = filestemplate_to_listing(projectpath, projname, template)

    # List of failed files with reasons
    failed = []
    # Now check all items in expected listing
    for item in expected_listing:
        itemsize = os.stat(item).st_size
        if not exists(item):
            failed.append(
                (item, 'Missing')
            )
        elif itemsize == 0:
            failed.append(
                (item, 'Size zero')
            )
    return failed
