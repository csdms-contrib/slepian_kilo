#!/bin/bash

# Get staged MATLAB files
FILES=$(git diff --cached --name-only --diff-filter=ACM | grep '\.m$')

# Exit if none
if [ -z "$FILES" ]; then
  exit 0
fi

# Exclude externals/
FILES=$(echo "$FILES" | grep -v '^externals/')

if [ -z "$FILES" ]; then
  exit 0
fi

echo "Formatting MATLAB files..."

# Build MATLAB cell array string
MATLAB_FILES=$(printf "'%s'," $FILES)

matlab -batch "addpath('externals/MBeautifier'); \
files = {${MATLAB_FILES}}; \
for k=1:length(files), MBeautify.formatFile(files{k}); end"

# Re-add formatted files
git add $FILES
