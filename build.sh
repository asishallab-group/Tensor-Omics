#!/bin/bash
# build.sh | Optimized build script for FPM with dynamic alignment
# Build with selected profile and alignment parameter:
# Default fallback alignment for the most likely situation:

source build_utils.sh

init "$@"

mkdir -p build

# trigger clean build on branch switch
if [[ $(command -v git) ]]; then
  current_branch=$(git branch --show-current 2>/dev/null || true)
  filename="build/.${current_branch//\//_.SLASH._}.branch" # replace / by _.SLASH._, as _.SLASH._ is very likely never being part of a branch name
  if [[ ! -f "$filename" ]]; then
    CLEAN_BUILD=1
    rm -f build/.*.branch
    : > "$filename"
  fi
fi

# Clean build directory if it exists
if [[ "$CLEAN_BUILD" ]]; then
  rm -rf build/${COMPILER}_*
fi

# Build with FPM first
utils_fpm build

check_exit_code "Build with fpm failed"

# Remove outdated .so file -> no accidental reuse
rm -f build/libtensor-omics.so

# Retrieve output path for .so from fpm and copy to build directory
while IFS= read -r line; do
    if [[ $line == *libtensor-omics\.so* ]]; then
        # remove leading whitespaces
        tensoromics_so="${line#"${line%%[![:space:]]*}"}"
    fi
done <<< "$(utils_fpm list 2>&1)"

cp "${tensoromics_so}" build 2>/dev/null
check_exit_code "No .so file created"

echo "Build complete with compiler: $COMPILER"
