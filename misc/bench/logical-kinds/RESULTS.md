# What `logical(c_bool)` costs, and what the marshalling copy cost before it

Run 2026-08-14 (marshalling) and re-run in full 2026-08-24 (read side), WSL2, gfortran 16.1
and ifx 2026.1. Reproduce with `./run_bench.sh`.

The benchmark was written for a question that is now settled: every generated C wrapper taking
or returning a logical array used to declare an automatic array of default `logical` and copy
elementwise, because `c_bool` is one byte and a default `logical` is four:

```fortran
logical, dimension(n_vectors) :: vectors_selection_mask_f
vectors_selection_mask_f = vectors_selection_mask
```

31 such array temporaries across 16 wrappers. `logical(c_bool)` was adopted throughout, the
copy is gone, and with it the meaning of the `in`/`out`/`auto` columns — **they measure a
copy the framework no longer performs.** They are kept because they are what the 2026-08-14
decision rested on, and because they still price the boundary for anyone reintroducing one.

The live question is what the *implementations* pay for holding their masks in `c_bool`, and
that is a same-column comparison between the two kinds. Writing a mask (`mk:`) was measured
from the start. **Reading one was not** — there was only `rd:logical`, with nothing to compare
it against, which made it the one number in the table that answered nothing. The two read
columns below close that gap, and they matter more than `mk:` does: a mask is written once
and read repeatedly.

## Numbers, ns per element

`in`/`out` are the two marshalling directions, `auto` the same copy into an automatic array
(what a wrapper really paid), `mk` computing a mask, `rd` reading a whole one with `count`,
`br` reading it one element at a time to guard work — `if (mask(i)) ...`, which is what the
implementations actually do (`tox_tissue_versatility_impl.F90:73`,
`tox_data_integration_jsd_impl.F90:182`, and every other mask site in `src/tox`).

### gfortran `-O0` — the default build

| n | in | out | auto | mk:logical | mk:c_bool | rd:logical | rd:c_bool | br:logical | br:c_bool |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 10 000 | 1.223 | 1.096 | 1.241 | 3.254 | 3.091 | 0.491 | 0.487 | 1.890 | 1.893 |
| 100 000 | 1.221 | 1.103 | 1.227 | 3.280 | 3.196 | 0.489 | 0.485 | 2.191 | 2.185 |
| 1 000 000 | 1.286 | 1.164 | 1.298 | 3.643 | 3.444 | 0.505 | 0.483 | 2.436 | 2.416 |
| 10 000 000 | 3.036 | 2.909 | 4.029 | 2.611 | 2.392 | 2.243 | 2.108 | 3.623 | 3.658 |

### gfortran `-O3 -march=native -funroll-loops -ftree-vectorize` — `--max-performance`

| n | in | out | auto | mk:logical | mk:c_bool | rd:logical | rd:c_bool | br:logical | br:c_bool |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 10 000 | 0.560 | 0.581 | 0.560 | 0.204 | 0.296 | 0.036 | 0.068 | 0.624 | 0.620 |
| 100 000 | 0.587 | 0.587 | 0.585 | 0.281 | 0.311 | 0.058 | 0.068 | 0.962 | 0.962 |
| 1 000 000 | 0.614 | 0.629 | 0.610 | 0.597 | 0.491 | 0.064 | 0.071 | 1.178 | 1.133 |
| 10 000 000 | 0.953 | 0.774 | 1.828 | 0.947 | 0.673 | 0.193 | 0.094 | 1.797 | 1.827 |

### ifx `-O0 -heap-arrays`

| n | in | out | auto | mk:logical | mk:c_bool | rd:logical | rd:c_bool | br:logical | br:c_bool |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 10 000 | 3.373 | 3.444 | 3.363 | 5.624 | 5.701 | 1.886 | 1.853 | 1.810 | 1.829 |
| 100 000 | 3.133 | 3.114 | 3.105 | 5.672 | 5.829 | 1.517 | 1.525 | 2.195 | 2.189 |
| 1 000 000 | 3.189 | 3.161 | 3.200 | 6.006 | 6.090 | 1.552 | 1.518 | 2.432 | 2.400 |
| 10 000 000 | 5.360 | 5.265 | 6.363 | 4.584 | 4.561 | 3.734 | 3.474 | 3.906 | 3.884 |

### ifx `-O3 -xHost -heap-arrays`

| n | in | out | auto | mk:logical | mk:c_bool | rd:logical | rd:c_bool | br:logical | br:c_bool |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 10 000 | 0.127 | 0.076 | 0.074 | 0.117 | 0.127 | 0.062 | 0.046 | 0.252 | 0.199 |
| 100 000 | 0.150 | 0.097 | 0.107 | 0.194 | 0.160 | 0.065 | 0.045 | 0.272 | 0.228 |
| 1 000 000 | 0.156 | 0.110 | 0.156 | 0.466 | 0.321 | 0.069 | 0.048 | 0.572 | 0.415 |
| 10 000 000 | 0.439 | 0.330 | 1.144 | 0.758 | 0.536 | 0.224 | 0.074 | 0.684 | 0.568 |

A second full pass agreed cell-for-cell: median deviation 1.2%, 136 of 144 cells within 5%,
worst 15% (ifx `-O3`, `in`, n = 1 000 000).

## What the read side says

**Reading a whole mask (`rd`, a `count` reduction): `c_bool` wins wherever the array does not
fit in cache, and the one place it loses is free.** At 10 M elements it is 2.1× faster under
gfortran `-O3` (0.094 vs 0.193) and 3.0× under ifx `-O3` (0.074 vs 0.224) — a quarter of the
bytes, read at the same bandwidth. Under ifx `-O3` it is ahead at *every* size (1.35–1.44× at
10 k–1 M). The single adverse cell is gfortran `-O3` at small n, where `count` over default
`logical` vectorises to 0.036 ns/element and over `c_bool` to 0.068: 1.9× relative, and
**0.3 µs absolute on a 20 000-gene mask.** Unoptimised builds show no difference either way.

**Reading it element by element (`br`, the pattern the implementations use): gfortran is
indifferent, ifx prefers `c_bool`.** Under both gfortran profiles the two kinds are within 4%
at every size — the loop is branch-bound, not bandwidth-bound, and the load width does not
show. Under ifx `-O3`, `c_bool` is 19–38% faster (0.415 vs 0.572 at 1 M). There is no size and
no compiler at which the branch loop is measurably slower in `c_bool`.

So the kind change is a win on the read side and never a loss on it. That is the opposite of
the concern it was raised against, and it holds for the read pattern the code actually uses.

## What it means

**1. The marshalling copy cost more than a pass over the mask, and far less than a call.** It
ran at 2–16× the cost of one `count` over the same array, so a wrapper whose implementation
read the mask a handful of times did spend a visible slice of its time marshalling. In
absolute terms it was small — a 20 000-gene mask at `--max-performance` is 12 µs, once per
call — but it is now zero.

**2. `c_bool` inside the implementations is not the cost that was feared, on either side.**
Writing a mask into `c_bool` is *cheaper* at large n on every configuration (29% at 10 M under
gfortran `-O3`) and dearer only at small n under the optimising builds (+45% at 10 000 on
gfortran `-O3`), where the kind conversion is not amortised by memory traffic. Reading follows
the same shape but more favourably, and under ifx has no small-n penalty at all. Crossover for
both is between 100 000 and 1 000 000 elements. The second-order objection does not hold.

**3. The real finding was never performance, it was the stack.** `logical :: dst(n)` with `n` a
runtime extent is an automatic array. At n = 10 000 000 that is 40 MB, and **ifx segfaults on
it at default settings** — an 8 MB stack. gfortran survives, because it moves large automatic
arrays to the heap; ifx does not unless told, and `-heap-arrays` is not in the ifx profile in
`fpm.toml`. The allocation is also what makes `auto` diverge from `in` at 10 M: 1.83 vs 0.95
ns/element under gfortran `-O3`, 1.14 vs 0.44 under ifx.

This was reachable from the published API: `serialize_logical`/`deserialize_logical` take
`arr(n_elements)` of the caller's choosing, and `is_outlier(n_genes)` is an output copy of the
same shape. Demonstrated for the exact construct the wrappers used; not reproduced end to end
through a real ifx-built binding.

## Outcome

**The stack finding was fixed directly, not with a flag.** Adding `-heap-arrays` to the ifx
profile was rejected: it fixes the crash per compiler, by a flag anyone can drop, and leaves
the failure unchecked. The generated locals became `allocatable` through `M_ALLOCATE` instead,
which is heap on every compiler and returns `ERR_ALLOC_FAIL` rather than aborting. Landed on
`131-codegen-new` as *Allocate the C layer's converted locals instead of declaring them
automatic*.

**`logical(c_bool)` was adopted, on memory rather than speed.** The decision (FES) rested on
removing the copy buffer entirely and on three bytes saved per logical throughout, with the
benchmark asked only whether computation would get significantly *slower*. It does not — and
at the sizes where the memory argument bites, `c_bool` is faster for the same reason it is
smaller. The two arguments are one argument.

**The read-side measurement, added 2026-08-24, is the one that actually tests that claim**, and
it strengthens it: up to 3× faster on a whole-array read at 10 M, 19–38% faster on ifx's
element-wise read, and never slower on the element-wise read on any compiler at any size. The
only measured cost of `c_bool` anywhere in the table is below ~100 000 elements under
gfortran `-O3` — sub-microsecond per call on a 10 000-element mask, which is not a decision
input.

The "all-or-nothing sweep" objection does not survive either: a `logical(c_bool)` actual
against a default-`logical` dummy is a compile error, which is ordinary kind discipline exactly
as for `int32` and `real64`, and the errors are the sweep's to-do list rather than a risk.

**Still worth profiling if masks ever look hot.**
`tox_data_integration_per_family` has six 2-D masks of `n_neighbors x n_points`. At 50 x 1000
that is 300 000 elements per call; read element-wise at `--max-performance` that is ~340 µs,
and across thousands of families it is the only site measured here that could reach seconds.

## Notes on the benchmark itself

Two compilation units, best of three trials, and never `-flto` / `-ipo`. Both rules were
learned here and both produced confidently wrong numbers first; they apply to every benchmark
in this tree, so the reasoning lives in [`../README.md`](../README.md) rather than here.

**The 2026-08-24 re-run also replaced `system_clock` with `cpu_time`, and the earlier ifx
tables should be read with that in mind.** ifx's `system_clock` is wall-clock based — its
count is the Unix epoch in microseconds, where gfortran's is monotonic uptime — so an NTP step
under WSL2 lands as a backwards jump mid-measurement. That does not merely add noise: it makes
a trial look *faster*, and best-of-N then keeps precisely the corrupted trial. It showed as
single ifx cells an order of magnitude below their neighbours, in a different cell on every
run (`br:logical` 0.200 between neighbours of 2.39 and 4.27). `cpu_time` is monotonic on both
compilers, the benchmark is single-threaded so its CPU time is its wall time, and every
measurement runs for 0.1 s or more — far above the clock's resolution. With it, two full
passes agree to a 1.2% median.
