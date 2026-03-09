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

handle_args "$@"

echo "Detected alignment: $ALIGN"

if [[ -z "$SKIP_KINDS_TEST" ]]; then
  bash -s -- "$@" <<'EOF'
  function get_directives() {
    echo "-D'OPEN_PAREN=(' -D'CLOSE_PAREN=)' -D'$1(KIND)=KIND(KIND)' -D'$2(KIND)=$1 OPEN_PAREN KIND CLOSE_PAREN' -D'$3(KIND)=$1 OPEN_PAREN 2 CLOSE_PAREN'"
  }

  failed=0
  directives=()
  directives+=("-DTEST_KIND_MISMATCH_C_INT $(get_directives integer int32 c_int) -D'c_int64_t(KIND)=integer OPEN_PAREN KIND CLOSE_PAREN' -D'c_size_t(KIND)=integer OPEN_PAREN KIND CLOSE_PAREN'")
  directives+=("-DTEST_KIND_MISMATCH_C_DOUBLE $(get_directives real real64 c_double)")
  directives+=("-DTEST_KIND_MISMATCH_C_DOUBLE_COMPLEX $(get_directives complex real64 c_double_complex)")
  for d in "${directives[@]}"; do
    test_directive=${d%% *}
    test_directive=${test_directive#-DTEST_KIND_MISMATCH_}
    echo -en "Testing safeguard for mismatch for $test_directive: "
    if [[ $(bash build.sh "$@" "${directives}" 1>kinds.out 2>/dev/null; grep "Divi.*zero" kinds.out) ]]; then
      echo "success"
    else
      echo "failure"
      cat kinds.out
      failed=1
    fi
  done
  rm kinds.out
  exit $failed
EOF
  check_exit_code "Kind Mismatch Test failed"

  echo "Compiling src/"
  bash build.sh --clean-build "$@"
  check_exit_code "Build failed"
else
  bash build.sh "$@"
  check_exit_code "Build failed"
fi

echo "Using compiler: $COMPILER"

mkdir -p $BUILD_DIR

# Clean up any existing run_tests file/directory
if [ -e "$EXECUTABLE" ]; then
  echo "Removing existing $EXECUTABLE..."
  rm -f "$EXECUTABLE"
fi

# Compile mod_test_suite.f90 first to generate mod_test_suite.mod
$COMPILER $FLAGS $MODULE_FLAG -DDEFAULT_ALIGNMENT=$ALIGN $MAX_PERF_FLAG \
  -I$BUILD_DIR \
  -c $TEST_DIR/mod_test_suite.f90
check_exit_code "mod_test_suite compilation failed"

echo "Compiling test modules..."
# Then compile test/ modules using .mod files from build/
$COMPILER $FLAGS $MODULE_FLAG -DDEFAULT_ALIGNMENT=$ALIGN $MAX_PERF_FLAG \
  -I$BUILD_DIR \
  -c $TEST_DIR/*.[fF]90
check_exit_code "Test compilation failed"

# Move object files to build/
mv *.o $BUILD_DIR/

check_build "*.o" "*.mod"
check_exit_code "Module compilation failed"


echo "Linking executable..."
# Finally link everything together
$COMPILER $FLAGS -I$BUILD_DIR \
  $BUILD_DIR/*.o -o $EXECUTABLE

check_exit_code "Executable compilation failed"

rm -f test_*.bin
rm -f test_*.zip
rm -f manifest.txt

echo "Running tests..."
# Run the executable
$EXECUTABLE $ARGS

check_exit_code "Tests failed"

if [[ -z "$KEEP_BIN" && -z "$KEEP_FILES" ]]; then
  rm -f test_*.bin
fi
if [[ -z "$KEEP_ZIP" && -z "$KEEP_FILES" ]]; then
  rm -f test_*.zip
fi
if [[ -z "$KEEP_TXT" && -z "$KEEP_FILES" ]]; then
  rm -f test_*.txt
  rm -f manifest.txt
fi
