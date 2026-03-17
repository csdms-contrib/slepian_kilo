#!/bin/bash

# Get staged MATLAB files, exclude externals/
FILES=$(git diff --cached --name-only --diff-filter=ACM | grep '\.m$' | grep -v '^externals/')

if [ -z "$FILES" ]; then
    echo "No MATLAB files staged."
    exit 0
fi

echo "Running MBeautifier on staged MATLAB files..."

# Build MATLAB cell array of files
MATLAB_FILES=$(printf "'%s'," $FILES)

# Run MATLAB in batch mode
matlab -batch "addpath('externals/MBeautifier'); files={${MATLAB_FILES}}; scripts.format_staged(files);"

# Re-add formatted files
git add $FILES

