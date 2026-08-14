# Benchmarks

One directory per question. Each is self-contained and answers something a design decision
turned on, so that the decision can be re-checked later against a new compiler or a new
machine rather than re-argued.

| directory | question |
|---|---|
| [`logical-kinds/`](logical-kinds/) | What does the `c_bool` marshalling copy at the C boundary cost, and does computing a mask in `c_bool` cost more than in the default kind? |

## Layout of a benchmark

```
<question>/
    RESULTS.md          the numbers, when they were taken, on what, and what they meant
    run_bench.sh        builds and runs every configuration; the whole entry point
    bench_*.f90         the sources
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

**Take the best of several trials and discard non-positive results.** Partly the usual noise
floor: the single-shot `-O0` figures came out three times too high. But it is also necessary
for correctness on this toolchain, because **ifx's `system_clock` is wall-clock based** — its
count is the Unix epoch in microseconds, where gfortran's is monotonic uptime — so an NTP step
lands as a jump backwards in the middle of a measurement. On WSL2 that happens often enough to
see in a single run.

Two smaller ones worth repeating:

- **Consume the result and vary the input.** Accumulate into a checksum the program prints at
  the end, and perturb one element per repetition, or the loop is dead code.
- **Measure both profiles.** `build.sh` defaults to `-O0` and only gives `-O3` under
  `--max-performance`, and conclusions differ between them: the copy is 4–6 ns per element at
  `-O0` and well under 1 ns at `-O3`.
