# How f42 should recognise a subnormal, and what the test costs in the length kernel

Run twice in full on 2026-10-06, WSL2 on an Intel Core i5-10500H (12 threads), gfortran 16.1.0,
ifx 2026.1.1 and nvfortran 26.9. Reproduce with `./run_bench.sh`; nvfortran is included when it is
on `PATH` or named by `NVFORTRAN=/path/to/nvfortran ./run_bench.sh`.

This benchmark was written for a rule f42 considered and then dropped (see Decision below): reading
a subnormal value as zero wherever a length or a direction is derived from it, so that a normal
build would give the same answer as a host that flushes subnormals (FTZ/DAZ). Such a rule needs a
reliable test for "is this subnormal", and the obvious one is a floating-point comparison:

```fortran
abs(val) < tiny(val)
```

and ifx's default floating-point model folds it. Under `-O3 -xHost` without `-fp-model precise`,
which is `main`'s `--max-performance` profile, it is `.false.` for a real subnormal when the process
runs in gradual-underflow mode, the mode Python and R leave the CPU in. The same build starts its
own programs with FTZ and DAZ on, so its Fortran tests read subnormals as zero anyway and cannot see
the fold; only a call from Python or R can. ifx `-O0` and the `-fp-model precise` builds start with
both off and do not fold the comparison.

Three spellings read the representation instead of comparing values. Checked separately, with the
values built from their bit patterns and gradual underflow forced, each classifies the largest and
the smallest subnormal and their negatives as subnormal, and `tiny`, `±0`, `1` and NaN as not, and
keeps `-0.0` as `-0.0`, under every profile the project builds and under ifx `-fp-model fast=2` and
nvfortran `-fast` as well.

| spelling | test |
|---|---|
| exponent | `exponent(val) < minexponent(val)` |
| class | `ieee_class(val) == ieee_positive_denormal .or. ieee_class(val) == ieee_negative_denormal` |
| bits | `0 < magnitude < transfer(tiny(1.0_real64), 0_int64)`, with `magnitude = iand(transfer(val, 0_int64), huge(0_int64))` |

The bits test rests on positive IEEE doubles ordering like their bit patterns. Being integer
arithmetic on the representation, it also classifies a subnormal correctly when denormals-are-zero
is on, as gfortran's `-ffast-math` startup sets it, where exponent and class miss it. Under DAZ with
ifx (`main`'s `-O3` flags, `-fp-model fast=2`) class still classifies correctly and exponent does
not; nvfortran, which starts every program with FTZ and DAZ on, even at `-O0`, gets all three right.
Missing it costs nothing in practice: with DAZ on, the arithmetic reads the value as zero whatever
the test says.

## What is measured

`scaled_length`'s two loops (`f42_vector_impl`), the largest magnitude and then the scaled sum of
squares, with the test applied to every component read in both. The vectors hold ordinary values in
[0.5, 2), so no variant changes a single length and all seven compute the same result. Seven
kernels:

- `plain`: today's loops, no test;
- `… call`: each spelling as a `pure` function in another module, which is what an fpm build gives a
  caller, since nothing is inlined across files without `-flto`/`-ipo`;
- `… inline`: each spelling written into the loop itself, which separates the cost of the test from
  the cost of the call.

The helpers and the kernels are separate compilation units built into a shared library with
`-fPIC`, as `libtensor_omics.so` is, and the driver links against it. The pass count is calibrated
per vector length so that a trial of `plain` takes at least 0.15 s, and is then held fixed for all
seven kernels. The input is perturbed every pass, and the summed lengths are compared across
kernels.

**The clock is the process CPU time from `clock_gettime(CLOCK_PROCESS_CPUTIME_ID)`, not
`cpu_time`.** nvfortran's `cpu_time` is `gettimeofday` underneath, a wall clock, and this machine
steps its wall clock back by about a second every half minute: the trap `../README.md` describes for
ifx's `system_clock`. A first run timed with `cpu_time` produced two nvfortran cells 2.2 and 3.9 times
lower than in the tables below, which did not reproduce. On gfortran and ifx the two clocks measure the same thing.

**Each figure below is the best of ten trials**: best of five in each of two full runs. Across the
147 cells outside nvfortran's threaded profile the two runs agree to a median of 4.3%, 132 cells
within 10% and the worst 21%. Other work ran on the machine meanwhile; the load average reached 13
during nvfortran's threaded profile, which runs on all 12 threads itself.

## Numbers, ns per element

### gfortran `-O0` — the default build

| dims | plain | exponent call | class call | bits call | exponent inline | class inline | bits inline |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 8 | 5.203 | 15.281 | 30.408 | 14.082 | 11.058 | 21.586 | 16.852 |
| 64 | 4.642 | 14.102 | 25.513 | 12.079 | 9.948 | 20.978 | 15.686 |
| 1024 | 4.521 | 13.291 | 25.738 | 11.897 | 9.479 | 21.280 | 15.177 |

### gfortran `-O3 -march=native` — `--max-performance`

| dims | plain | exponent call | class call | bits call | exponent inline | class inline | bits inline |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 8 | 1.354 | 9.038 | 13.124 | 4.565 | 7.016 | 9.959 | 1.656 |
| 64 | 0.780 | 8.970 | 12.535 | 4.676 | 6.429 | 9.395 | 1.269 |
| 1024 | 1.189 | 8.598 | 12.631 | 4.810 | 6.632 | 9.499 | 1.307 |

### ifx `-O0 -fp-model precise` — the default build on #202's and #220's branches

| dims | plain | exponent call | class call | bits call | exponent inline | class inline | bits inline |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 8 | 5.303 | 24.351 | 28.963 | 11.370 | 20.429 | 25.477 | 9.577 |
| 64 | 5.244 | 21.772 | 26.039 | 11.152 | 17.320 | 22.995 | 9.054 |
| 1024 | 5.208 | 21.168 | 26.882 | 10.200 | 16.720 | 23.067 | 8.341 |

### ifx `-O3 -xHost`, `main`'s flags — `--max-performance`

| dims | plain | exponent call | class call | bits call | exponent inline | class inline | bits inline |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 8 | 1.362 | 20.932 | 23.511 | 5.300 | 15.891 | 20.182 | 1.733 |
| 64 | 0.476 | 18.827 | 21.616 | 4.786 | 14.006 | 18.474 | 0.681 |
| 1024 | 0.329 | 18.662 | 21.455 | 4.636 | 13.713 | 18.101 | 0.572 |

### ifx `-O3 -xHost -assume protect_parens -fp-model precise` — `--max-performance` on #202's and #220's branches

| dims | plain | exponent call | class call | bits call | exponent inline | class inline | bits inline |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 8 | 1.553 | 18.538 | 23.486 | 4.676 | 15.238 | 21.119 | 2.217 |
| 64 | 1.551 | 17.159 | 21.969 | 4.641 | 14.257 | 19.193 | 1.885 |
| 1024 | 1.942 | 16.998 | 21.959 | 4.588 | 13.878 | 18.850 | 2.002 |

### nvfortran `-O0` — the default build

| dims | plain | exponent call | class call | bits call | exponent inline | class inline | bits inline |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 8 | 5.640 | 12.144 | 23.489 | 22.875 | 10.624 | 22.195 | 21.768 |
| 64 | 5.121 | 11.164 | 20.391 | 20.541 | 9.783 | 18.960 | 19.008 |
| 1024 | 4.925 | 11.045 | 20.123 | 19.841 | 9.513 | 18.521 | 18.363 |

### nvfortran `-O3`, without `-stdpar`

| dims | plain | exponent call | class call | bits call | exponent inline | class inline | bits inline |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 8 | 1.444 | 7.030 | 22.621 | 16.041 | 5.194 | 18.497 | 14.777 |
| 64 | 0.546 | 6.466 | 21.100 | 15.195 | 5.271 | 16.884 | 13.340 |
| 1024 | 0.464 | 6.633 | 21.313 | 16.121 | 5.178 | 17.371 | 14.002 |

### nvfortran `-O3 -Mconcur -stdpar=multicore` — `--max-performance`

Not tabulated: this profile threads every `do concurrent`, including the reduction inside
`scaled_length`, so each call of an 8-element kernel pays for starting and joining threads. The
process CPU time, summed over the threads, came to 40 000–270 000 ns per element at 8 dims,
6 000–23 000 at 64 and 190–3 400 at 1024, varying up to fourfold between the two runs and with no
consistent order among the kernels. The threading cost dominates every column, so the profile says
nothing about the spellings and a good deal about itself. nvfortran also warns that `-Mconcur` is
deprecated.

## What they mean

**The class test is the most expensive spelling in every profile**, except at nvfortran `-O0`, where
it ties with bits. Every compiler implements it as library calls: gfortran calls
`__ieee_arithmetic_MOD_ieee_class_type_eq` twice per element, ifx `for_ieee_class_k8` and
`for_ieee_class_eq`.

**The exponent test is cheap only on nvfortran.** gfortran lowers `exponent()` to a call of libm's
`frexp`, ifx to `for_exponent8_v`, nvfortran to `pgf90_expondx`, which is the fastest of the three:
6.5–7.0 ns per element called at `-O3`, against 8.6–9.0 on gfortran and 17–21 on ifx.

**The bits test, written inline, is close to free on gfortran and ifx at `-O3`**: +10–63% over
`plain` on gfortran, +27–74% with `main`'s ifx flags and +3–43% with `-fp-model precise`, where the
sum runs sequentially anyway and the integer work hides in its latency. The spread is wide because
`plain` itself is under 2 ns; the absolute cost is at most 0.7 ns per element. On nvfortran the bits
test is the slow one, 13–16 ns per element at `-O3` whether called or inline, because nvfortran
turns `transfer` into a call of `pghpf_transfer_i8` even at `-O3`.

**At `-O3`, calling the test from another module costs 2–5 ns per element more than writing it
inline** on gfortran and ifx: 2.5–4.1 ns for bits, 2.0–5.0 for exponent, 2.4–3.4 for class. For
bits, the inline loop vectorizes and a call in its body stops that. Exponent and class vectorize in
neither form, being library calls already, so for them the difference is the call itself.

**No spelling is cheapest on all three compilers.** Called at `-O3`, bits is 1.8–2.0× cheaper than
exponent on gfortran and 3.7–4.0× on ifx, and exponent is 2.3–2.4× cheaper than bits on nvfortran
(1.8–1.9× at `-O0`). At `-O0` gfortran is indifferent (1.1–1.2×) and ifx still favours bits by
2.0–2.1×.

**In absolute terms the cost is small at the sizes the library sees.** The length kernel reads
genes × axes components per pass, a few hundred thousand. For 300 000 that is about 1.5 ms per
pass for bits called, and about 6.5 ms for class called on ifx.

**The loop construct does not matter.** `scaled_length` sums with `do concurrent ... reduce`, and so
do these kernels. Written as a plain `do` instead, the three `--max-performance` profiles of
gfortran and ifx give the same figures within the run-to-run noise and in the same order, and ifx's
checksum differences stay: neither compiler uses the `reduce` clause to vectorize the sum; the
floating-point model decides.

## Decision

f42 does not detect subnormals at all (decided 2026-10-07). Its algorithms are to give a valid
result whether a subnormal reads as itself or as zero instead: a length or a direction is computed
in power-of-two scaled coordinates, so that no intermediate falls into the subnormal range, and no
test or division depends on a subnormal being non-zero. An input made only of subnormals can then
give different results in a build that flushes them — a correct direction where they exist, the
zero vector where they do not — and nothing tries to make those agree. Such inputs do not occur in
expression data, and a subnormal added to any value above about 2e-292 is below half an ulp of it,
so whether it reads as zero changes no result beyond rounding.

The numbers above are what that decision was weighed against: called from the module that would
hold it, the cheapest reliable test costs 4.6–7.0 ns per component read at `-O3`, against 0.3–1.9
for the loop without it, and only the bits test written inline comes close to free; and it would
have had to be applied at every read in every kernel to give consistent answers. Had the rule been kept, the plan was the bits test, and the
exponent test where `__NVCOMPILER` is defined.

## Things the numbers show on the way

- **Bit-identity depends on `-fp-model precise` on ifx.** With `main`'s flags the summed lengths
  differ between kernels in the last bits at 1024 dims, in both runs. Switching off vectorization or
  fused multiply-add contraction, alone or together, still leaves differences, because ifx's
  default model also reassociates scalar sums; only `-fp-model precise` makes every kernel give the
  same sum. nvfortran's
  threaded profile differs in the last bits too, from its threaded partial sums.
- **`-fp-model precise` makes today's loop 3.3× slower on ifx at 64 dims and 5.9× at 1024**
  (0.476 → 1.551, 0.329 → 1.942 ns per element), because it forbids reassociating the sum, which is
  what vectorizing it needs.
