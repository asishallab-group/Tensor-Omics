# Benchmarks

One directory per question. Each is self-contained and answers something a design decision
turned on, so that the decision can be re-checked later against a new compiler or a new
machine rather than re-argued.

| directory | question |
|---|---|
| [`logical-kinds/`](logical-kinds/) | What does the `c_bool` marshalling copy at the C boundary cost, and does holding a mask in `c_bool` rather than the default kind cost anything to write or to read? |
| [`r-binding-backend/`](r-binding-backend/) | Should the R binding go through Rcpp, cpp11, pure C `.Call`, or `.C`? |

## Layout of a benchmark

```
<question>/
    RESULTS.md          the numbers, when they were taken, on what, and what they meant
    run_bench.sh        one entry point that builds and runs every configuration,
      (or bench.R)      named for whatever language drives it
    bench_*.f90         the sources
    ...
```

`RESULTS.md` is the point of the directory. A benchmark whose numbers live only in a terminal
that has since scrolled away has to be run again before it can be cited, and by then nobody
remembers which flags it used.

Nothing here is built by `fpm` — `misc/` is outside its source directories — so a benchmark can
be a plain program with its own `main`, and adding one cannot break the library build.

## Two things that apply to every benchmark here

Both were learned the hard way in `logical-kinds/`, and both produced confidently wrong numbers
before they were noticed.

**Put the timed operation in its own compilation unit, and never build with `-flto` or `-ipo`.**
With everything in one program, ifx saw that a copied mask was only ever counted, rewrote
`count(dst)` as `count(src)`, and reported a copy costing a tenth of the truth — the tell was
that a copy *with* an allocation came out faster than the same copy without one. A separately
compiled callee is also the more honest model, since that is what a real caller faces.

**Time with `cpu_time`, not `system_clock`, and take the best of several trials.** The trials
are the usual noise floor — the single-shot `-O0` figures came out three times too high. The
clock is a correctness matter: **ifx's `system_clock` is wall-clock based**, its count being
the Unix epoch in microseconds where gfortran's is monotonic uptime, so an NTP step lands as a
jump backwards in the middle of a measurement. On WSL2 that happens often enough to see in a
single run.

That combination is worse than either part alone, and it is the trap: a backwards jump makes a
trial look *faster*, so best-of-N does not reject the corrupted trial — it selects it.
Discarding non-positive elapsed times catches only the jumps large enough to invert the
interval; a smaller step survives as a plausible number. It showed in `logical-kinds/` as
single ifx cells an order of magnitude below both their neighbours, landing in a different cell
on every run. `cpu_time` is monotonic on both compilers and is the same quantity as wall time
for a single-threaded benchmark whose measurements each run for 0.1 s or more.

Two smaller ones worth repeating:

- **Consume the result and vary the input.** Accumulate into a checksum the program prints at
  the end, and perturb one element per repetition, or the loop is dead code.
- **Measure both profiles.** `build.sh` defaults to `-O0` and only gives `-O3` under
  `--max-performance`, and conclusions differ between them: the copy is 4–6 ns per element at
  `-O0` and well under 1 ns at `-O3`.
