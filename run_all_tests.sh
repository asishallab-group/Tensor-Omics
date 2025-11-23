#!/bin/bash
./test_runner.sh "$@"
for i in  python/test/*.py;do echo; echo $i; python3 $i; done
for i in  r/test/*.R;do echo; echo $i; Rscript $i; done