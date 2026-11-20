#!/bin/bash
# Why f42_bit_masks is a core of procedures on plain words and not a derived type.
#
#   ./run_bench.sh           probes, then the benchmark
#   ./run_bench.sh probes    only the probes: what compiles, what leaks, what crashes
#   ./run_bench.sh bench     only the timings
#
# The probes also run under nvfortran when it is on PATH. If it needs a different environment
# to link, pass it in, e.g. NVFORTRAN="env PATH=/usr/bin:/opt/nvidia/.../compilers/bin nvfortran"
# (the value is split on spaces).
set -e
cd "$(dirname "$0")"

build_dir="${TMPDIR:-/tmp}/bit_masks_bench.$$"
mkdir -p "$build_dir"
trap 'rm -rf "$build_dir"' EXIT

NVFORTRAN="${NVFORTRAN:-nvfortran}"
GFORTRAN_O3="-O3 -march=native -mtune=native -funroll-loops -ftree-vectorize"
IFX_O3="-O3 -xHost -align array64byte -qopt-zmm-usage=high"
IFX_PROJECT="-assume protect_parens -fp-model precise"   # fpm.toml's default.ifx.flags

# compile <compiler> <flags...> -- <sources...> : objects and modules go to $build_dir
compile() {
    local compiler="$1"; shift
    local flags=()
    while [[ $1 != -- ]]; do flags+=("$1"); shift; done
    shift
    local module_dir=(-J "$build_dir")
    [[ $compiler != gfortran ]] && module_dir=(-module "$build_dir")
    $compiler "${flags[@]}" "${module_dir[@]}" -I "$build_dir" "$@"
}

first_error() {
    local line
    while IFS= read -r line; do
        if [[ $line == *rror* ]]; then echo "$line"; return; fi
    done
}

# run_probe <program> : runs it under the default 8 MB stack, without core dumps, output
# indented, and reports a nonzero exit status
run_probe() {
    local status=0 line
    # the outer 2>/dev/null only silences bash's own "Segmentation fault" job report
    { ( ulimit -s 8192 -c 0; exec "$1" ) >"$build_dir/out" 2>&1; } 2>/dev/null || status=$?
    while IFS= read -r line; do echo "    $line"; done <"$build_dir/out"
    if (( status != 0 )); then echo "    ABORTED: exit status $status"; fi
}

run_probes() {
    local compilers=(gfortran ifx)
    $NVFORTRAN --version >/dev/null 2>&1 && compilers+=("$NVFORTRAN")
    local compiler variant
    for compiler in "${compilers[@]}"; do
        echo "############ probes: ${compiler##* }"
        for variant in VIEW_OF_INOUT VIEW_OF_INPUT; do
            echo "probe_pure_view -D$variant"
            if compile $compiler -D$variant -- probes/probe_pure_view.F90 -o "$build_dir/probe" \
                    >"$build_dir/log" 2>&1; then
                run_probe "$build_dir/probe"
            else
                echo "    REJECTED: $(first_error <"$build_dir/log")"
            fi
        done
        for variant in INTRINSIC_ASSIGNMENT DEFINED_ASSIGNMENT; do
            echo "probe_owning_mask -D$variant"
            compile $compiler -D$variant -- probes/probe_owning_mask.F90 -o "$build_dir/probe" \
                >"$build_dir/log" 2>&1 || { cat "$build_dir/log"; continue; }
            run_probe "$build_dir/probe"
        done
        echo "probe_assignment_stack"
        { compile $compiler -O2 -- -c mask_candidates.f90 -o "$build_dir/mask_candidates.o" &&
          compile $compiler -O2 -- "$build_dir/mask_candidates.o" probes/probe_assignment_stack.f90 \
              -o "$build_dir/probe"; } >"$build_dir/log" 2>&1 || { cat "$build_dir/log"; exit 1; }
        run_probe "$build_dir/probe"
        echo
    done
}

# bench <compiler> <label> <flags...> : three units, no -flto / -ipo (../README.md)
bench() {
    local compiler="$1" label="$2"; shift 2
    echo "############ $compiler $label"
    compile $compiler "$@" -- -c mask_candidates.f90 -o "$build_dir/mask_candidates.o"
    compile $compiler "$@" -- -c bench_kernels.f90 -o "$build_dir/bench_kernels.o"
    compile $compiler "$@" -- "$build_dir/mask_candidates.o" "$build_dir/bench_kernels.o" \
        bench_bit_masks.f90 -o "$build_dir/bench"
    # pack:assignment at 1e7 bits overflows ifx's default stack (see probe_assignment_stack)
    ( ulimit -s unlimited; "$build_dir/bench" )
    echo
}

run_bench() {
    bench gfortran "-O0 (default build)"     -O0
    bench gfortran "-O3 (--max-performance)" $GFORTRAN_O3
    bench ifx      "-O0 (default build)"     -O0 $IFX_PROJECT
    bench ifx      "-O3 (--max-performance)" $IFX_O3 $IFX_PROJECT
}

case "${1:-all}" in
    probes) run_probes ;;
    bench)  run_bench ;;
    all)    run_probes; run_bench ;;
    *)      echo "usage: $0 [probes|bench]" >&2; exit 2 ;;
esac
