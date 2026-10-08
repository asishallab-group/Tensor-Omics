#!/bin/bash
# build.sh | Optimized build script for FPM with dynamic alignment
# Build with selected profile and alignment parameter:
# Default fallback alignment for the most likely situation:

source build_utils.sh

init "$@"

mkdir -p build

# a clean build whenever something changed that fpm cannot see -- see check_build_state
check_build_state

# Clean build directory if it exists
if [[ "$TOX_CLEAN_BUILD" ]]; then
  rm -rf build/${COMPILER}_*
fi

# clean output directories for safety, so no wrong libs will be linked accidentally
rm -f build/*.so
rm -f external/*.a

# Bring the generated sources up to date before anything reads them
generate_code

# Build with FPM first
# dependencies
cd external/loess_netlib
root=../..
fpm build --compiler "$COMPILER"
find_and_mv_libs "$(fpm build --compiler "$COMPILER" --list 2>&1)" "$root/external"
cd $root
# tox
# a library older than its own objects is one an earlier, failed build never linked
remove_stale_libraries "$(utils_fpm list 2>&1)" "$(fpm_package_name)"
utils_fpm build

check_exit_code "Build with fpm failed"

# Retrieve output path for .so from fpm and copy to build directory
find_and_mv_libs "$(utils_fpm list 2>&1)" "build"

stderr "
${COLOR_GREEN}Build complete with compiler${COLOR_CREAM}: $(echo_compiler $COMPILER)
"
