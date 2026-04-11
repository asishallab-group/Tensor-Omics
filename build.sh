#!/bin/bash
# build.sh | Optimized build script for FPM with dynamic alignment
# Build with selected profile and alignment parameter:
# Default fallback alignment for the most likely situation:

source build_utils.sh

COMPILER=$(get_compiler)
FLAGS=$(get_flags)
ALIGN=$(get_alignment)
handle_args "$@"

# trigger clean build on branch switch
if [[ $(which git) ]]; then
  git branch --show-current 2>/dev/null 1> build/.branch.tmp || true
  touch build/.branch
  if [[ $(diff build/.branch.tmp build/.branch) ]]; then
    CLEAN_BUILD=1
  fi
  mv build/.branch.tmp build/.branch
fi

# # Clean build directory if it exists
if [[ "$CLEAN_BUILD" ]]; then
  rm -rf build/${COMPILER}_*
fi

# Build with FPM first
generate_fpm_toml .fpm.toml $COMPILER > fpm.toml
fpm build --compiler $COMPILER --flag "$FLAGS $DIRECTIVES" --flag "-DDEFAULT_ALIGNMENT=$ALIGN" --flag "$MAX_PERF_FLAG"

check_exit_code "Build with fpm failed"

rm fpm.toml

# Copy .mod, .o and .so files from FPM build directories to build
rm -f build/*.o build/*.mod

# Copy .so, .mod, .o files if they exist
find build/"${COMPILER}"_*/ \( -name "*.so" -o -name "*.mod" -o -name "*.o" \) -exec cp {} build/ \;

echo "Build complete with compiler: $COMPILER, alignment: $ALIGN bytes"
