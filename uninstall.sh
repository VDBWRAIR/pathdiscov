#!/usr/bin/env bash

# Should perform all necessary steps to uninstall the pipeline

# Deactivate virtualenv if needed
deactivate
# Remove data files that may have been altered
rm usamriidPathDiscov/files/*
# Checkout files from git to replace removed ones
git checkout usamriidPathDiscov/files
# Remove virtualenv installation
rm -rf usamriidPathDiscov/{lib,lib64,include,bin}

echo "The pipeline should be uninstalled. Make sure you deactivate prior to reinstalling if needed"
