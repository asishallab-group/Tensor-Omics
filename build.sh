#!/bin/bash
# build.sh | Optimized build script for FPM with dynamic alignment
# Build with selected profile and alignment parameter:
# Default fallback alignment for the most likely situation:

source build_utils.sh

COMPILER=$(get_compiler)
FLAGS=$(get_flags)
ALIGN=$(get_alignment)

# Detect --max-performance flag
MAX_PERF_FLAG=""
for arg in "$@"; do
  if [[ "$arg" == "--max-performance" ]]; then
    MAX_PERF_FLAG="-DMAX_PERFORMANCE"
  fi
done

# Clean build directory if it exists
if [ -d "build" ]; then
  rm -rf build
fi

# Build with FPM first
fpm build --compiler $COMPILER --flag "$FLAGS" --flag "-DDEFAULT_ALIGNMENT=$ALIGN" --flag "$MAX_PERF_FLAG"

check_exit_code "Build with fpm failed"

# Move .mod, .o and .so files from FPM build directories to root
for compiler_dir in build/${COMPILER}_*; do
  if [ -d "$compiler_dir" ]; then
    echo "Processing FPM directory: $compiler_dir"
    # Move .so files if they exist
    mv "$compiler_dir"/*.so build/ 2>/dev/null || true
    # Move .mod files if they exist
    find "$compiler_dir" -name "*.mod" -exec mv {} build/ \; 2>/dev/null || true
    # Move .o files from subdirectories if they exist
    find "$compiler_dir" -name "*.o" -exec mv {} build/ \; 2>/dev/null || true
    # Remove the compiler-specific directory
    rm -rf "$compiler_dir"
  fi
done

echo "Build complete with compiler: $COMPILER, alignment: $ALIGN bytes"
