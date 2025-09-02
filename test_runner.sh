#!/bin/bash

SOURCE_DIR="src"
TEST_DIR="test"
BUILD_DIR="build"
EXECUTABLE="$BUILD_DIR/run_tests"

source build_utils.sh

COMPILER=$(get_compiler)
FLAGS=$(get_flags)
ALIGN=$(get_alignment)
MODULE_FLAG=$(get_module_flag $BUILD_DIR)

echo "Detected alignment: $ALIGN"
echo "Compiling src/"
bash build.sh
check_exit_code "Build failed"

echo "Using compiler: $COMPILER"

MAX_PERF_FLAG=""
for arg in "$@"; do
  if [[ "$arg" == "--max-performance" ]]; then
    MAX_PERF_FLAG="-DMAX_PERFORMANCE"
  fi
done

mkdir -p $BUILD_DIR

# Clean up any existing run_tests file/directory
if [ -e "$EXECUTABLE" ]; then
  echo "Removing existing $EXECUTABLE..."
  rm -rf "$EXECUTABLE"
fi

echo "Compiling test modules..."
# Then compile test/ modules using .mod files from build/
$COMPILER $FLAGS $MODULE_FLAG -DDEFAULT_ALIGNMENT=$ALIGN $MAX_PERF_FLAG \
  -I$BUILD_DIR \
  -c $TEST_DIR/*.f90

# Move object files to build/
mv *.o $BUILD_DIR/ 2>/dev/null || true
mv *.mod $BUILD_DIR/ 2>/dev/null || true

check_build "*.o" "*.mod"
check_exit_code "Module compilation failed"


echo "Linking executable..."
# Finally link everything together
$COMPILER $FLAGS -I$BUILD_DIR \
  $BUILD_DIR/*.o -o $EXECUTABLE

check_exit_code "Executable compilation failed"

echo "Running tests..."
# Run the executable
$EXECUTABLE "$@"