import string
import os
from os.path import exists, basename, normpath
import sys
import pprint

# Default stages
# Each has it's own template under output_files_templates
STAGES = [
    'input',
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
    projpath = normpath(projpath)
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
    projpath = normpath(projectpath)
    # Get the basename(last thing after the last /)
    projname = basename(normpath(projpath))

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
        if not exists(item):
            failed.append(
                (item, 'Missing')
            )
        elif os.stat(item).st_size == 0:
            failed.append(
                (item, 'Size zero')
            )
    return failed

def verify_standard_stages_files(projectpath, templatedir):
    ''' Hardcoded verification of standard stages '''
    #templates.append(resource_filename(__name__, tfile))
    projectname = basename(normpath(projectpath))
    # Fetch template files for each stage from inside
    # of templatedir
    templates = []
    for stage in STAGES:
        tfile = os.path.join(templatedir,stage+'.lst')
        templates.append(tfile)
    return verify_project(projectpath, projectname, templates)

def main():
    import argparse
    from pkg_resources import resource_filename
    parser = argparse.ArgumentParser(
        'Verify all stages have necessary files'
    )

    parser.add_argument(
        dest='projectpath',
        help='Path to project to check'
    )

    templatesdir = resource_filename(__name__, 'output_files_templates')
    parser.add_argument(
        '--templatedir',
        default=templatesdir,
        help='Path to directory that contains templates of files for each stage' \
            '[Default: %(default)s]'
    )

    args = parser.parse_args()

    missingfiles = verify_standard_stages_files(args.projectpath, args.templatedir)
    pprint.pprint(missingfiles)

if __name__ == '__main__':
    main()
