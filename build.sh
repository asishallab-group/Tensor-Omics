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
    TOX_CLEAN_BUILD=1
    rm -f build/.*.branch
    : > "$filename"
  fi
fi

# Clean build directory if it exists
if [[ "$TOX_CLEAN_BUILD" ]]; then
  rm -rf build/${COMPILER}_*
fi

# Build with FPM first
utils_fpm build

check_exit_code "Build with fpm failed"

# Remove outdated .so file -> no accidental reuse
rm -f build/libtensor-omics.so
utils_fpm list 2>&1
# Retrieve output path for .so from fpm and copy to build directory
while IFS= read -r line; do
    if [[ $line == *\.so* ]]; then
        # remove leading whitespaces
        so="${line#"${line%%[![:space:]]*}"}"
        cp "${so}" build 2>/dev/null
    fi
done <<< "$(utils_fpm list 2>&1)"

check_exit_code "No .so file created"

stderr "
${COLOR_GREEN}Build complete with compiler${COLOR_CREAM}: $(echo_compiler $COMPILER)
"
