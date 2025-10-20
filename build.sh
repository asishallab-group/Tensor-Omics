#!/bin/bash
# build.sh | Optimized build script for FPM with dynamic alignment
# Build with selected profile and alignment parameter:
# Default fallback alignment for the most likely situation:
ALIGN=32
# Detect capabilities in order of descending priority:
if lscpu | grep -q amx; then
ALIGN=128
elif lscpu | grep -q avx512; then
ALIGN=64
elif lscpu | grep -q avx2; then
ALIGN=32
elif lscpu | grep -q sse2; then
ALIGN=16
fi
# Detect compiler and choose appropriate profile:
if [[ "$FC" == "ifx" || "$FC" == "ifort" ]]; then
  FLAGS="-O3 -qopenmp -xHost -align array64byte -qopt-zmm-usage=high -qopt-prefetch=3 -qopt-matmul -fPIC"
  COMPILER="ifx"
else
  FLAGS="-O3 -march=native -mtune=native -fopenmp -funroll-loops -ftree-vectorize -fPIC"
  COMPILER="gfortran"
fi

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
export FC
fpm build --compiler $COMPILER --flag "$FLAGS" --flag "-DDEFAULT_ALIGNMENT=$ALIGN" --flag "$MAX_PERF_FLAG"

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
echo "Flags used: $FLAGS $MAX_PERF_FLAG"

# Verify that we have the necessary files for test_runner.sh
mod_count=$(find build -name "*.mod" | wc -l)
obj_count=$(find build -name "*.o" | wc -l) 
so_count=$(find build -name "*.so" | wc -l)

if [ $mod_count -eq 0 ] || [ $obj_count -eq 0 ]; then
  echo "Warning: Missing .mod or .o files needed for test compilation"
fi