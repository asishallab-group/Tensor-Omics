#!/bin/bash
# Two compilation units, no -flto / -ipo, so the copies cannot be elided.
set -e
cd "$(dirname "$0")"
build_and_run() {
    local cc="$1"; shift
    local label="$1"; shift
    rm -f ./*.mod ./*.o
    echo "############ $cc $label"
    $cc "$@" -c bench_kernels.f90 -o /tmp/bk.o
    $cc "$@" /tmp/bk.o bench_logical.f90 -o /tmp/bench.x
    /tmp/bench.x
    echo
}
build_and_run gfortran "-O0 (default build)"          -O0
build_and_run gfortran "-O3 (--max-performance)"      -O3 -march=native -mtune=native -funroll-loops -ftree-vectorize
build_and_run ifx      "-O0 -heap-arrays"             -O0 -heap-arrays
build_and_run ifx      "-O3 -xHost -heap-arrays"      -O3 -xHost -heap-arrays
rm -f ./*.mod ./*.o
