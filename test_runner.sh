#!/bin/bash

source build_utils.sh

init "$@"

generate_code

if [[ -z "$TOX_SKIP_KINDS_TEST" ]]; then
  bash -s -- "$@" <<'EOF'
  source build_utils.sh
  echo 

  function get_directives() {
    echo '"--directive='OPEN_PAREN=\('" "--directive='CLOSE_PAREN=\)'" "--directive='$1\(KIND\)=KIND\(KIND\)'" "--directive='$2\(KIND\)=$1 OPEN_PAREN KIND CLOSE_PAREN'" "--directive='$3\(KIND\)=$1 OPEN_PAREN 2 CLOSE_PAREN'"'
  }
  failed=0
  c_int=("--directive=TEST_KIND_MISMATCH_C_INT" $(get_directives integer int32 c_int) "--directive='c_int64_t(KIND)=integer OPEN_PAREN KIND CLOSE_PAREN'" "--directive='c_size_t(KIND)=integer OPEN_PAREN KIND CLOSE_PAREN'")
  c_double=("--directive=TEST_KIND_MISMATCH_C_DOUBLE" $(get_directives real real64 c_double))
  c_double_complex=("--directive=TEST_KIND_MISMATCH_C_DOUBLE_COMPLEX" "$(get_directives complex real64 c_double_complex)")
  # c_char needs no get_directives: the other three override a kind that modules elsewhere
  # declare with, so those declarations have to be macro-rewritten too or they fail first and
  # for the wrong reason. Overriding c_char breaks nothing before safeguard, which is compiled
  # first by design, so its own guard is what fails.
  # These five need no get_directives. The three above override a kind that modules elsewhere
  # declare with, so those declarations must be macro-rewritten too or they fail first and for
  # the wrong reason. f42_safeguard depends on nothing and so compiles before any of them,
  # which is what lets its own guard be the failure in every case below.
  c_char=("--directive=TEST_KIND_MISMATCH_C_CHAR")
  c_bool=("--directive=TEST_KIND_MISMATCH_C_BOOL")
  c_size_t=("--directive=TEST_KIND_MISMATCH_C_SIZE_T")
  c_int64_t=("--directive=TEST_KIND_MISMATCH_C_INT64_T")
  c_signed_char=("--directive=TEST_KIND_MISMATCH_C_SIGNED_CHAR")
  for d in c_int c_double c_double_complex c_char c_bool c_size_t c_int64_t c_signed_char; do
    msg_prefix="Testing safeguard for mismatch for $COLOR_COPPER$d"
    declare -n directives="$d"
    # these builds are meant to fail in the preprocessor, so regenerating for each of them
    # would only cost time -- the build below does it once for the run
    if [[ $(bash build.sh "$@" --skip-code-generation "${directives[@]}" 1>kinds.out 2>/dev/null ; grep "Divi.*zero" kinds.out) ]]; then
      stderr "$msg_prefix$COLOR_CREAM: ${COLOR_GREEN}success"
    else
      stderr "$msg_prefix$COLOR_CREAM: ${COLOR_RED}failure"
      cat kinds.out >&2
      failed=1
    fi
  done
  rm kinds.out
  exit $failed
EOF
  check_exit_code "Kind Mismatch Test failed"

  stderr "Compiling src/"
  bash build.sh --clean-build "$@" --skip-code-generation --compiler="$COMPILER"
  check_exit_code "Build failed"
else
  bash build.sh "$@" --compiler="$COMPILER" --skip-code-generation
  check_exit_code "Build failed"
fi

rm -f *.test.*
rm -f manifest.txt

stderr "Running tests..."

# By default same behavior as compiling manually, as fpm has some struggles sometimes with correct linking, e.g. some routine changes in src, but the tests use the old implementation.
if [[ -z $TOX_REUSE_MOD_FILES ]]; then
  rm -f build/$COMPILER_*/**/*mod_test*
fi

# Run the executable
utils_fpm test ${TOX_TEST_TARGET:-run_tests}

check_exit_code "Tests failed"

if [[ -z "$TOX_KEEP_FILES" ]]; then
  for type in zip txt bin; do
    keep=TOX_KEEP_${type^^}
    if [[ -z "${!keep}" ]]; then
      rm -f *.test.$type
    fi
  done

  if [[ -z "$TOX_KEEP_TXT" ]]; then
    rm -f manifest.txt
  fi
fi

stderr "
${COLOR_GREEN}All tests passed with compiler${COLOR_CREAM}: $(echo_compiler $COMPILER)
"
