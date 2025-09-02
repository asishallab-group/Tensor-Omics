#!/bin/bash

SOURCE_DIR="src/parallelization_experiment"
BUILD_DIR="build"
EXECUTABLE="$BUILD_DIR/benchmark"

source build_utils.sh

COMPILER=$(get_compiler)
FLAGS=$(get_flags)
ALIGN=$(get_alignment)
MODULE_FLAG=$(get_module_flag $BUILD_DIR)

echo "Detected alignment: $ALIGN"
echo "Compiling src/"
bash build.sh
check_exit_code "Build failed"

check_build parallelization_experiment.mod
check_exit_code "Missing module"

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

echo "Compiling modules..."
$COMPILER $FLAGS $MODULE_FLAG -DDEFAULT_ALIGNMENT=$ALIGN $MAX_PERF_FLAG \
  -I$BUILD_DIR -c $SOURCE_DIR/*.[fF]90

# Move object files to build/
mv *.o $BUILD_DIR/ 2>/dev/null || true
mv *.mod $BUILD_DIR/ 2>/dev/null || true

check_build "*parallelization_experiment*.o" "*benchmark*.o" "*benchmark*.o"
check_exit_code "Module compilation failed"

echo "Linking executable..."
# Finally link everything together
$COMPILER $FLAGS -I$BUILD_DIR \
  $BUILD_DIR/*parallelization_experiment*.o \
  $BUILD_DIR/*benchmark*.o -o $EXECUTABLE

check_exit_code "Executable compilation failed"

for arg in "$@"; do
  if [[ "$arg" == "norun" ]]; then
    exit
  fi
done

echo "Running benchmark..."
# Run the executable
$EXECUTABLE "$@" > results/benchmark_$COMPILER.csv