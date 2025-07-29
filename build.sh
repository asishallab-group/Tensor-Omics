#!/bin/bash
# build.sh | Optimized build script for FPM with dynamic alignment
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
  FLAGS="-O3 -march=native -mtune=native -fopenmp -ffast-math -funroll-loops -ftree-vectorize -fassociative-math -fPIC"
  COMPILER="gfortran"
fi

# Detect --max-performance flag
MAX_PERF_FLAG=""
for arg in "$@"; do
  if [[ "$arg" == "--max-performance" ]]; then
    MAX_PERF_FLAG="-DMAX_PERFORMANCE"
  fi
done

# Informative output:
# Build with selected profile and alignment parameter:
export FC
fpm build --compiler $COMPILER --flag "$FLAGS" --flag "-DDEFAULT_ALIGNMENT=$ALIGN" --flag "$MAX_PERF_FLAG"

echo "Build complete with compiler: $COMPILER, alignment: $ALIGN bytes, optimization flags: $FLAGS, max performance: $MAX_PERF_FLAG"

# Find .so file in the build directory
sofile=$(find build -name 'libtensor-omics.so' | head -n 1)

# Create symbolic link in the build directory with relative path
if [[ -n "$sofile" ]]; then
  # Extract the relative path from build/ directory
  relative_path=$(realpath --relative-to=build "$sofile")
  ln -sf "$relative_path" build/libtensor-omics.so
  echo "Created symlink: build/libtensor-omics.so -> $relative_path"
else
  echo "Warning: libtensor-omics.so not found in build directory"
fi