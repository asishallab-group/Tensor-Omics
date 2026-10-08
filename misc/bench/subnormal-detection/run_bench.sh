#!/bin/bash
# The helpers and the kernels are separate compilation units built into a shared library with
# -fPIC, as libtensor_omics.so is, and the driver links against it: no -flto / -ipo, so a call to
# the helper costs what it costs in the library. A configuration whose compiler is missing or
# fails is reported and skipped; the rest still run. nvfortran is included when it is on PATH or
# named by NVFORTRAN=/path/to/nvfortran. Compiler output goes to a log printed only on failure.
cd "$(dirname "$0")" || exit 1
here=$PWD
build=/tmp/subnormal-detection-bench.$$
trap 'rm -rf "$build" "$build.log"' EXIT
build_and_run() {
    local fc="$1"; shift
    local label="$1"; shift
    echo "############ $fc $label"
    if ! command -v "$fc" >/dev/null; then
        echo "not found, skipped"; echo; return
    fi
    rm -rf "$build"; mkdir -p "$build"
    if ! (cd "$build" &&
          $fc "$@" -fPIC -c "$here/bench_helpers.f90" &&
          $fc "$@" -fPIC -c "$here/bench_kernels.f90" &&
          $fc "$@" -shared -o libbench.so bench_helpers.o bench_kernels.o &&
          $fc "$@" -c "$here/bench_driver.f90" &&
          $fc "$@" -o bench.x bench_driver.o -L. -lbench -Wl,-rpath,"$build") >"$build.log" 2>&1; then
        echo "build failed:"; cat "$build.log"; rm -f "$build.log"; echo; return
    fi
    rm -f "$build.log"
    "$build/bench.x" || echo "run failed"
    echo
}
nv="${NVFORTRAN:-nvfortran}"
build_and_run gfortran "-O0 (default build)"      -O0
build_and_run gfortran "-O3 (--max-performance)"  -O3 -march=native -mtune=native -fopenmp -funroll-loops -ftree-vectorize
build_and_run ifx      "-O0, precise (default build on #202/#220)" -O0 -assume protect_parens -fp-model precise
build_and_run ifx      "-O3 -xHost, main's flags (--max-performance)" -O3 -xHost -align array64byte -qopt-zmm-usage=high -qopt-prefetch=3 -qopt-matmul
build_and_run ifx      "-O3 -xHost, precise (--max-performance on #202/#220)" -O3 -xHost -align array64byte -qopt-zmm-usage=high -qopt-prefetch=3 -qopt-matmul -assume protect_parens -fp-model precise
build_and_run "$nv"    "-O0 (default build)"      -O0
build_and_run "$nv"    "-O3, without -stdpar"     -O3
build_and_run "$nv"    "-O3 -Mconcur -stdpar=multicore (--max-performance)" -O3 -Mconcur -fopenmp -stdpar=multicore
rm -f "$build.log"
