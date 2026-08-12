#!/bin/bash
# Every suite in one run: the Fortran tests through test_runner.sh, then the Python and the R
# bindings against the library that build just produced.
#
# Exits non-zero if ANY of the three fails. It did not used to: `./test_runner.sh` ran with its
# status ignored, each binding loop swallowed its own with `&&`/`||`, and the script's exit code
# was whatever the trailing `rm` returned. A Fortran suite that failed to *compile* therefore
# reported a clean run -- hit for real on 2026-08-07, with `<ERROR> stopping due to failed
# compilation` sitting in the log and every Python and R suite printed green beneath it.

source build_utils.sh

failed=0

./test_runner.sh "$@" || failed=1

echo
for i in python/test/*.py; do
  if python3 "$i" >/dev/null 2>/dev/null; then
    cecho "$i: ${COLOR_GREEN}success"
  else
    cecho "$i: ${COLOR_RED}failed"
    failed=1
  fi
done

echo
for i in r/test/*.R; do
  if Rscript "$i" >/dev/null 2>/dev/null; then
    cecho "$i: ${COLOR_GREEN}success"
  else
    cecho "$i: ${COLOR_RED}failed"
    failed=1
  fi
done

# the suites leave these behind; -f so an empty run is not itself a failure
rm -f *.zip *.bin

exit $failed
