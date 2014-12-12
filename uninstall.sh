#!/usr/bin/env bash

# Should perform all necessary steps to uninstall the pipeline

# Clean out compiled binaries so they are recopied
function clean_bin() {
    for f in bwa cap3 fastqc getorf Ray Ray2 samtools bowtie* blast*
    do
        rm -f usamriidPathDiscov/bin/${f}
    done
}

# Full clean build and egg
rm -rf build
rm -rf usamriidPathDiscov.egg-info

# Should we uninstall everything or just the usamriidPathDiscov package
# Takes a long time to reinstall numpy/matplotlib
if [ "$1" == "-full" ]
then
    # Remove all things that installer touches
    rm -rf usamriidPathDiscov/{lib,lib64,include,bin,download}
    # Replace download from git
    git checkout usamriidPathDiscov/download
else
    # Use pip to uninstall(If multiple version, may need to call more than once)
    # Ensure activated
    . usamriidPathDiscov/bin/activate
    while pip uninstall -y usamriidPathDiscov; do sleep 1; done;
    clean_bin
fi

echo "The pipeline should be uninstalled"
