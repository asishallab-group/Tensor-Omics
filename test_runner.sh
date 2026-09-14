#!/bin/bash

source build_utils.sh

init "$@"

generate_code

if [[ -z "$TOX_SKIP_KINDS_TEST" ]]; then
  bash -s -- "$@" <<'EOF'
  source build_utils.sh
  echo 

  failed=0
  # One directive per kind, and nothing else. Each TEST_KIND_MISMATCH_* redefines its kind
  # inside f42_safeguard only, and the safeguard depends on nothing, so it compiles before
  # any module that declares with a C kind -- which is what makes its own guard, and not
  # some unrelated declaration, the error this looks for. (Macro-rewrites of int32, real64
  # and c_int used to ride along here. They never reached the compiler: the quoting split
  # them into stray positional arguments, and the guards fired without them anyway.)
  for kind in c_int c_double c_double_complex c_char c_bool c_size_t c_int64_t c_signed_char; do
    directive="--directive=TEST_KIND_MISMATCH_${kind^^}"
    msg_prefix="Testing safeguard for mismatch for $COLOR_COPPER$kind"
    # these builds are meant to fail in the preprocessor, so regenerating for each of them
    # would only cost time -- the build below does it once for the run
    bash build.sh "$@" --skip-code-generation "$directive" 1>kinds.out 2>/dev/null
    if file_matches 'Divi.*zero' kinds.out; then
      stderr "$msg_prefix$COLOR_CREAM: ${COLOR_GREEN}success"
    else
      stderr "$msg_prefix$COLOR_CREAM: ${COLOR_RED}failure"
      printf '%s\n' "$(<kinds.out)" >&2
      failed=1
    fi
  done
  rm kinds.out
  exit $failed
EOF
  check_exit_code "Kind Mismatch Test failed"

  stderr "Compiling src/"
  bash build.sh "$@" --skip-code-generation --compiler="$COMPILER"
  check_exit_code "Build failed"
else
  bash build.sh "$@" --compiler="$COMPILER" --skip-code-generation
  check_exit_code "Build failed"
fi

rm -f *.test.*
rm -f manifest.txt

stderr "Running tests..."

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
