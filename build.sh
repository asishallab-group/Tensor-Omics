#!/bin/bash
# build.sh | Optimized build script for FPM with dynamic alignment
# Build with selected profile and alignment parameter:
# Default fallback alignment for the most likely situation:

source build_utils.sh

init "$@"

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
utils_fpm build

check_exit_code "Build with fpm failed"

if [[ -z "$KEEP_FPM_TOML" ]]; then
  rm fpm.toml
fi

# Copy .mod, .o and .so files from FPM build directories to build
rm -f build/*.o build/*.mod

# Copy latest .so
cp "$(find build/"${COMPILER}"_*/ \( -name "*.so" \) -printf "%T@ %p\n" | sort -n | tail -1 | sed "s/^[^ ]* //" )" build
check_exit_code "No .so file found"

echo "Build complete with compiler: $COMPILER, alignment: $ALIGN bytes"
