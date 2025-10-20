#!/bin/bash

SOURCE_DIR="src"
TEST_DIR="test"
BUILD_DIR="build"
EXECUTABLE="$BUILD_DIR/run_tests"

# Detect alignment - FIX: Initialize with default value
ALIGN=32
if lscpu | grep -q amx; then
  ALIGN=128
elif lscpu | grep -q avx512; then
  ALIGN=64
elif lscpu | grep -q avx2; then
  ALIGN=32
elif lscpu | grep -q sse2; then
  ALIGN=16
fi

echo "Detected alignment: $ALIGN"

# Detect compiler and flags
if [[ "$FC" == "ifx" || "$FC" == "ifort" ]]; then
  FLAGS="-O3 -qopenmp -xHost -align array64byte -qopt-zmm-usage=high -qopt-prefetch=3 -qopt-matmul -fPIC"
  MODULE_FLAG="-module $BUILD_DIR"
  COMPILER="ifx"
else
  FLAGS="-O3 -march=native -mtune=native -fopenmp -funroll-loops -ftree-vectorize -fPIC"
  MODULE_FLAG="-J$BUILD_DIR"
  COMPILER="gfortran"
fi

MAX_PERF_FLAG=""
TEST_ARGS=()
BUILD_ARGS=()
for arg in "$@"; do
  if [[ "$arg" == "--max-performance" ]]; then
    MAX_PERF_FLAG="-DMAX_PERFORMANCE"
    BUILD_ARGS+=("$arg")
  else
    TEST_ARGS+=("$arg")
  fi
done

echo "Compiling src/"
bash build.sh "${BUILD_ARGS[@]}"

echo "Using compiler: $COMPILER"

mkdir -p $BUILD_DIR

# Clean up any existing run_tests file/directory
if [ -e "$EXECUTABLE" ]; then
  echo "Removing existing $EXECUTABLE..."
  rm -rf "$EXECUTABLE"
fi

echo "Compiling test modules..."
# Then compile test/ modules using .mod files from build/
$COMPILER $FLAGS $MODULE_FLAG -DDEFAULT_ALIGNMENT=$ALIGN $MAX_PERF_FLAG \
  -I$BUILD_DIR -I$SOURCE_DIR -I$TEST_DIR \
  -c $TEST_DIR/*.f90

compilation_result=$?
echo "Test compilation exit code: $compilation_result"


# Move object files to build/
mv *.o $BUILD_DIR/ 2>/dev/null || true
mv *.mod $BUILD_DIR/ 2>/dev/null || true


echo "Linking executable..."
# Finally link everything together
$COMPILER $FLAGS -I$BUILD_DIR \
  $BUILD_DIR/*.o -o $EXECUTABLE

linking_result=$?
echo "Linking exit code: $linking_result"

echo "Running tests..."
# Run the executable with filtered arguments (excluding --max-performance)
$EXECUTABLE "${TEST_ARGS[@]}"

# Clean up test files if they exist
rm -f test_*arr*.bin