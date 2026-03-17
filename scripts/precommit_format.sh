#!/bin/bash
FILES=$(git diff --cached --name-only --diff-filter=ACM | grep '\.m$' | grep -v '^externals/')

if [ -z "$FILES" ]; then
    echo "No MATLAB files staged."
    exit 0
fi

echo "Running MBeautifier on staged MATLAB files..."

# Convert relative paths to absolute paths
ABS_FILES=()
for f in $FILES; do
    ABS_FILES+=("$(pwd)/$f")
done

# Build MATLAB cell array
MATLAB_FILES=$(printf "'%s'," "${ABS_FILES[@]}")

# Run MATLAB batch
matlab -batch "addpath('externals/MBeautifier'); scripts.format_staged({${MATLAB_FILES}});"

# Re-add files
git add "${FILES[@]}"
