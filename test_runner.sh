#!/bin/bash

source build_utils.sh

init "$@"

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
  bash build.sh --clean-build --keep-fpm-toml "$@"
  check_exit_code "Build failed"
else
  bash build.sh --keep-fpm-toml "$@"
  check_exit_code "Build failed"
fi

echo "Using compiler: $COMPILER"

rm -f test_*.bin
rm -f test_*.zip
rm -f manifest.txt

echo "Running tests..."

# By default same behavior as compiling manually, as fpm has some struggles sometimes with correct linking, e.g. some routine changes in src, but the tests use the old implementation.
if [[ -z $REUSE_MOD_FILES ]]; then
  rm build/$COMPILER_*/**/*mod_test*
fi

# Run the executable
utils_fpm test ${TEST_TARGET:-run_tests}

check_exit_code "Tests failed"

rm fpm.toml

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
