#!/usr/bin/env bash

# Should perform all necessary steps to uninstall the pipeline

# Remove binaries copied from install as well as build and possibly altered files
rm -rf usamriidPathDiscov/{files,build}
# Checkout files from git to replace removed ones
git checkout usamriidPathDiscov/files
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
fi

echo "The pipeline should be uninstalled"
