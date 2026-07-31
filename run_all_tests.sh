#!/bin/bash
./test_runner.sh "$@"

source build_utils.sh
echo
for i in  python/test/*.py; do python3 $i >/dev/null 2>/dev/null && cecho "$i: ${COLOR_GREEN}success" || cecho "$i: ${COLOR_RED}failed"; done
echo
for i in  r/test/*.R; do Rscript $i >/dev/null 2>/dev/null && cecho "$i: ${COLOR_GREEN}success" || cecho "$i: ${COLOR_RED}failed"; done

rm *.zip *.bin