import subprocess
from os.path import *
import os
import sys

def install_diamond(options):
    '''
    Installs diamond

    :param object options.diamond.src: Path to diamond tarball
    :param object options.diamond.install_to: Path to create bin/diamond
    '''
    # Path to diamond bin install path
    bin_path = join(options.diamond.install_to,'bin')
    # Path to diamond executable
    diamond_install_path = join(bin_path, 'diamond')
    sys.stdout.write("Installing diamond into {0}\n".format(options.diamond.install_to))
    # If not exist or not executable
    if os.access(diamond_install_path, os.X_OK):
        return
    cmd = "mkdir -p {1}; tar xzf {0} -C {1}".format(options.diamond.src, bin_path)
    subprocess.check_call(cmd, shell=True)
