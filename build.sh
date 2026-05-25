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
  filename=".${COMPILER}.${current_branch}.branch"
  filename=build/${filename//\//_.SLASH._} # replace / by _.SLASH._, as _.SLASH._ is very likely never being part of a branch name
  if [[ ! -f "$filename" ]]; then
    TOX_CLEAN_BUILD=1
    rm -f build/.*.branch
    : > "$filename"
  fi
fi

# Clean build directory if it exists
if [[ "$TOX_CLEAN_BUILD" ]]; then
  rm -rf build/${COMPILER}_*
fi

# clean output directories for safety, so no wrong libs will be linked accidentally
rm -f build/*.so
rm -f external/*.a

# Build with FPM first
# dependencies
cd external/loess_netlib
root=../..
fpm build --compiler "$COMPILER"
find_and_mv_libs "$(fpm build --compiler "$COMPILER" --list 2>&1)" "$root/external"
cd $root
# tox
utils_fpm build

check_exit_code "Build with fpm failed"

# Remove outdated .so file -> no accidental reuse
rm -f build/libtensor-omics.so

# Retrieve output path for .so from fpm and copy to build directory
find_and_mv_libs "$(utils_fpm list 2>&1)" "build"

check_exit_code "No .so file created"

stderr "
${COLOR_GREEN}Build complete with compiler${COLOR_CREAM}: $(echo_compiler $COMPILER)
"
