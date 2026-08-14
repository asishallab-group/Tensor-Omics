# What the `c_bool` marshalling copy costs

Run 2026-08-14, WSL2, gfortran 16.1 and ifx 2026.1, best of three trials per measurement.
Reproduce with `./run_bench.sh`.

Every generated C wrapper taking or returning a logical array declares an automatic array of
default `logical` and copies elementwise, because `c_bool` is one byte and a default `logical`
is four:

```fortran
logical, dimension(n_vectors) :: vectors_selection_mask_f
vectors_selection_mask_f = vectors_selection_mask
```

31 such array temporaries across 16 wrappers (plus 13 scalar ones, which cost nothing).
Declaring the implementation's dummy `logical(c_bool)` would remove the temporary and the copy.

## Numbers, ns per element

`in`/`out` are the two marshalling directions, `auto` is the same copy into an automatic array
(what a wrapper really pays), `mk` is computing a mask, `rd` is reading one with `count`.

### gfortran `-O0` — the default build

| n | in | out | auto | mk:logical | mk:c_bool | rd:logical |
|---:|---:|---:|---:|---:|---:|---:|
| 10 000 | 1.375 | 1.287 | 1.379 | 3.473 | 3.382 | 0.548 |
| 100 000 | 1.363 | 1.284 | 1.383 | 3.606 | 3.461 | 0.541 |
| 1 000 000 | 1.465 | 1.358 | 1.440 | 3.989 | 3.713 | 0.569 |
| 10 000 000 | 3.382 | 3.273 | 4.438 | 2.847 | 2.561 | 2.411 |

### gfortran `-O3 -march=native -funroll-loops -ftree-vectorize` — `--max-performance`

| n | in | out | auto | mk:logical | mk:c_bool | rd:logical |
|---:|---:|---:|---:|---:|---:|---:|
| 10 000 | 0.611 | 0.633 | 0.615 | 0.213 | 0.323 | 0.039 |
| 100 000 | 0.653 | 0.654 | 0.650 | 0.300 | 0.335 | 0.058 |
| 1 000 000 | 0.711 | 0.692 | 0.725 | 0.673 | 0.523 | 0.071 |
| 10 000 000 | 0.988 | 0.826 | 1.940 | 1.054 | 0.751 | 0.212 |

### ifx `-O0 -heap-arrays`

| n | in | out | auto | mk:logical | mk:c_bool | rd:logical |
|---:|---:|---:|---:|---:|---:|---:|
| 10 000 | 3.736 | 3.831 | 3.720 | 6.322 | 6.339 | 1.928 |
| 100 000 | 3.528 | 3.495 | 3.485 | 6.482 | 6.504 | 1.673 |
| 1 000 000 | 3.639 | 3.606 | 3.584 | 6.914 | 6.844 | 1.687 |
| 10 000 000 | 6.076 | 6.067 | 7.046 | 5.126 | 5.102 | 3.955 |

### ifx `-O3 -xHost -heap-arrays`

| n | in | out | auto | mk:logical | mk:c_bool | rd:logical |
|---:|---:|---:|---:|---:|---:|---:|
| 10 000 | 0.140 | 0.084 | 0.084 | 0.132 | 0.144 | 0.071 |
| 100 000 | 0.169 | 0.112 | 0.117 | 0.223 | 0.191 | 0.073 |
| 1 000 000 | 0.180 | 0.125 | 0.179 | 0.512 | 0.366 | 0.080 |
| 10 000 000 | 0.494 | 0.381 | 1.222 | 0.797 | 0.591 | 0.226 |

## What it means

**1. The copy costs more than a pass over the mask, and far less than a call.** It runs at
2–16× the cost of one `count` over the same array, so a wrapper whose kernel reads the mask a
handful of times does spend a visible slice of its time marshalling. But in absolute terms it
is small: a 20 000-gene mask at `--max-performance` is **13 µs**, once per call. Against any
kernel doing real work per gene that is noise.

**2. `c_bool` inside the kernels is not the cost that was feared.** Computing a mask into
`c_bool` rather than default `logical` is *cheaper* at large n on every configuration — 29%
cheaper at 10M under gfortran `-O3`, because it moves a quarter of the bytes — and dearer only
at small n under the optimising builds (+52% at 10 000 on gfortran `-O3`), where the kind
conversion is not amortised by memory traffic. Crossover is somewhere between 100 000 and
1 000 000 elements. The second-order objection does not hold.

**3. The real finding is not performance, it is the stack.** `logical :: dst(n)` with `n` a
runtime extent is an automatic array. At n = 10 000 000 that is 40 MB, and **ifx segfaults on
it at default settings** — an 8 MB stack. gfortran survives, because it moves large automatic
arrays to the heap; ifx does not unless told, and `-heap-arrays` is not in the ifx profile in
`fpm.toml`. The allocation is also what makes `auto` diverge from `in` at 10M: 1.94 vs 0.99
ns/element under gfortran `-O3`, 1.22 vs 0.49 under ifx.

This is reachable from the published API: `serialize_logical`/`deserialize_logical` take
`arr(n_elements)` of the caller's choosing, and `is_outlier(n_genes)` is an output copy of the
same shape. Demonstrated for the exact construct the wrappers use; not yet reproduced end to
end through a real ifx-built binding.

## Outcome

**The stack finding was fixed directly, not with a flag.** Adding `-heap-arrays` to the ifx
profile was rejected: it fixes the crash per compiler, by a flag anyone can drop, and leaves
the failure unchecked. The generated locals became `allocatable` through `M_ALLOCATE` instead,
which is heap on every compiler and returns `ERR_ALLOC_FAIL` rather than aborting. Landed on
`131-codegen-new` as *Allocate the C layer's converted locals instead of declaring them
automatic*.

**`logical(c_bool)` was adopted, on memory rather than speed.** The decision (FES) rests on
removing the copy buffer entirely and on three bytes saved per logical throughout, with the
benchmark asked only whether computation would get significantly *slower*. It does not — and
at the sizes where the memory argument bites, `c_bool` is faster for the same reason it is
smaller: 29% at 10M elements under gfortran `-O3`, 26% under ifx. The two arguments are one
argument. The only measured cost is below ~100 000 elements, where computing a mask into
`c_bool` runs ~50% dearer per element — 1.1 µs on a 10 000-element mask, which is not a
decision input.

The "all-or-nothing sweep" objection raised against it does not survive either: a
`logical(c_bool)` actual against a default-`logical` dummy is a compile error, which is
ordinary kind discipline exactly as for `int32` and `real64`, and the errors are the sweep's
to-do list rather than a risk.

**Still worth profiling if masks ever look hot.**
`tox_data_integration_per_family` has six 2-D mask temporaries of `n_neighbors x n_points`. At
50 x 1000 that is 300 000 element-copies per call, ~210 µs, and if the wrapper is entered per
family across thousands of families it is the only site measured here where the copy could
reach seconds.

## Notes on the benchmark itself

Two compilation units, and it must not be built with `-flto` / `-ipo`. With everything in one
program ifx sees that a copied mask is only ever counted, rewrites `count(dst)` as `count(src)`
and reports a copy costing a tenth of the truth.

Best of three trials, discarding non-positive results. Partly the usual noise floor — the
single-shot `-O0` figures came out 3x too high — and partly necessary, because **ifx's
`system_clock` is wall-clock based** (its count is the Unix epoch in microseconds, where
gfortran's is monotonic uptime), so an NTP step lands as a backwards jump mid-measurement.
