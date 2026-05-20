#!/bin/bash
./test_runner.sh "$@"
echo
for i in  python/test/*.py; do python3 $i >/dev/null 2>/dev/null && echo -en success: || echo -en failed:; echo " $i"; done
echo
for i in  rcpp/test/*.R; do Rscript $i >/dev/null 2>/dev/null && echo -en success: || echo -en failed:; echo " $i"; done