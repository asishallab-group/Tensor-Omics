# Why `f42_bit_masks` is a core of procedures on plain words, not a derived type

Run 2026-09-15 on an i5-10500H (AVX2, no AVX-512) under WSL2, with gfortran 16.1, ifx 2026.1
and nvfortran 25.7. Reproduce with `./run_bench.sh` (`probes` or `bench` runs one half).

## The question

Issue #167 asks for a bit-mask module that replaces the hand-rolled int32 bit masks in
`tox_paralog_analysis_impl` and serves Shatter's per-seed masks (#205). Two designs were
weighed:

- **A derived type** `bit_mask` that behaves like a logical array: `.and.`, `.or.`, `.not.`,
  assignment to and from `logical(c_bool)` arrays, `all`/`any`/`count`. The first sketch held
  the words in a pointer component, so a mask could also be a view of words a caller owns.
  Lazy evaluation was considered too: `a .and. b` building a node that computes a bit only when
  it is asked for, instead of materialising the result.
- **A core** of pure procedures on plain int32 words plus `n_bits`: set, test, fill, in-place
  `and`/`or`/`and_not`/`xor`/`not`, copy, count, fused `and_count`, set-bit iteration, pack and
  unpack. A column of an `n_bits x n_masks` mask is simply `words(:, i_mask)`.

The core was chosen. The probes show that every form of the type breaks a rule or a compiler,
and the timings show that the type's convenient forms are the slow ones.

`mask_candidates.f90` holds both designs as the module could contain them: the core, and the
type in three forms (`owned_mask` with an allocatable component, an eager `.and.` and a defined
assignment from flags; `view_mask` with a pointer component; `mask_expression`, lazy).

## The probes: what compiles, what leaks, what crashes

| probe | gfortran | ifx | nvfortran |
|---|---|---|---|
| `probe_pure_view -DVIEW_OF_INOUT`: a pure procedure points a view at words it may modify | compiles | compiles | compiles |
| `probe_pure_view -DVIEW_OF_INPUT`: the same for `intent(in)` words | **rejected** | **rejected** (#7145) | compiles (non-standard) |
| `probe_owning_mask -DINTRINSIC_ASSIGNMENT`: pointer component, owns flag, final | **SIGSEGV** | **SIGSEGV** | runs |
| `probe_owning_mask -DDEFINED_ASSIGNMENT`: the same, `=` copies the words | clean | clean | **leaks 102 of 105 buffers** |
| `probe_assignment_stack`: `mask = flags` at 10^7 bits, 8 MB stack | runs | **SIGSEGV** | runs |

**A pure procedure cannot make a view of its input.** Pointing at an `intent(in)` dummy is
forbidden in pure code (F2018 C1594; gfortran: "Bad target in pointer assignment in PURE
procedure"). TensorOmics implementations are pure, and most masks reach them as `intent(in)`
arguments. So a view type could only wrap masks the procedure writes, and a lazy node, which has
to point at its operands, cannot be built at all. nvfortran accepts the construct, so a design
checked only there would break on the other two.

**No owning pointer type is correct on all three compilers.** Operators have to return new
masks, so the type must own some of its words and free them in a final procedure.
- With the default `=`, which copies the pointer, the finalizer of the function result frees
  the words the left-hand side now points at: gfortran and ifx crash on the first
  `a = new_mask(...)`.
- With a defined `=` that copies the words, gfortran and ifx are clean, but nvfortran 25.7
  never finalizes a function result and leaks every operator result.

**A defined assignment from a logical array can crash ifx.** The standard passes the
right-hand side of a defined assignment as if it were in parentheses, and ifx copies it onto the
stack first: 10 MB at 10^7 bits, beyond the default 8 MB. `call core_pack(...)` on the same data
runs. `-Warray-temporaries` does not flag the copy.

An allocatable component avoids all three, and is what the timings use for `owned_mask`. But
every operator result and every whole-mask assignment is then a heap allocation, and the
implementation layer forbids allocation (`helper/codegen/README.md`): its check looks for the
`allocatable` attribute on locals and dummies, so a local `type(owned_mask)` would pass it
while allocating.

## The timings

Each operation is timed in every form the two designs offer, with `logical(c_bool)` arrays,
today's masks, as the reference (`bench_kernels.f90`):

| form | what it is |
|---|---|
| `flags` | `logical(c_bool)` arrays |
| `core` | plain words, in-place or fused procedures: the chosen design |
| `eager` | `owned_mask`'s `.and.`, which allocates its result |
| `view` | `view_mask`, the words behind a pointer component |
| `lazy_bits`, `lazy_words` | `mask_expression`, evaluated per bit or per word on demand |
| `owned`, `assignment` | `owned_mask`'s intrinsic `=` copy and its defined `mask = flags` |

The operations are the ones Shatter and the paralog code perform: `count(a .and. b)` (Shatter's
merge test), `c = a .and. b` into an existing mask, the per-iteration backup copy, and packing a
`logical(c_bool)` array at the API edge. The last line is Shatter's whole merge step: `count`
of the conjunction of every pair of 100 masks of 10^5 bits. Density 0.5 (merge: 0.3).

### gfortran -O0 — the default build

| ns per bit | 10^3 | 10^5 | 10^7 |
|---|---:|---:|---:|
| `and_count:flags` | 1.3936 | 3.4096 | 3.2247 |
| `and_count:core` | 0.1576 | 0.1224 | 0.1118 |
| `and_count:eager` | 0.2519 | 0.1825 | 0.1707 |
| `and_count:lazy_bits` | 12.6388 | 20.2630 | 19.6603 |
| `and_count:lazy_words` | 0.6028 | 0.4188 | 0.4080 |
| `and:flags` | 1.2126 | 1.1810 | 1.2184 |
| `and:core` | 0.1018 | 0.0690 | 0.0698 |
| `and:eager` | 0.1339 | 0.0724 | 0.0721 |
| `and:view` | 0.1454 | 0.1063 | 0.1095 |
| `copy:flags` | 0.6670 | 0.6212 | 0.6444 |
| `copy:core` | 0.0933 | 0.0626 | 0.0638 |
| `copy:owned` | 0.0378 | 0.0023 | 0.0038 |
| `pack:core` | 2.1991 | 7.0557 | 6.9270 |
| `pack:assignment` | 2.1909 | 6.8856 | 6.9781 |

Merge, ms: `merge:flags` 902.0, `merge:core` 53.5

### gfortran -O3 — `--max-performance`

| ns per bit | 10^3 | 10^5 | 10^7 |
|---|---:|---:|---:|
| `and_count:flags` | 0.3929 | 3.1162 | 3.1967 |
| `and_count:core` | 0.0252 | 0.0112 | 0.0116 |
| `and_count:eager` | 0.0373 | 0.0120 | 0.0148 |
| `and_count:lazy_bits` | 5.6850 | 10.2873 | 10.2802 |
| `and_count:lazy_words` | 0.2752 | 0.1535 | 0.1569 |
| `and:flags` | 0.7614 | 3.4555 | 3.4514 |
| `and:core` | 0.0149 | 0.0035 | 0.0066 |
| `and:eager` | 0.0286 | 0.0040 | 0.0066 |
| `and:view` | 0.0282 | 0.0134 | 0.0138 |
| `copy:flags` | 0.0197 | 0.0185 | 0.0608 |
| `copy:core` | 0.0136 | 0.0014 | 0.0038 |
| `copy:owned` | 0.0243 | 0.0021 | 0.0039 |
| `pack:core` | 0.5460 | 3.5677 | 3.5183 |
| `pack:assignment` | 0.5618 | 3.4960 | 3.5335 |

Merge, ms: `merge:flags` 1182.9, `merge:core` 5.5

### ifx -O0 — the default build

| ns per bit | 10^3 | 10^5 | 10^7 |
|---|---:|---:|---:|
| `and_count:flags` | 2.0626 | 4.4259 | 4.4249 |
| `and_count:core` | 0.1139 | 0.0768 | 0.0769 |
| `and_count:eager` | 0.3062 | 0.1386 | 0.2062 |
| `and_count:lazy_bits` | 16.6381 | 20.5605 | 20.0346 |
| `and_count:lazy_words` | 1.0804 | 0.3912 | 0.3805 |
| `and:flags` | 1.8919 | 1.8379 | 1.8895 |
| `and:core` | 0.0949 | 0.0610 | 0.0632 |
| `and:eager` | 0.1813 | 0.0636 | 0.0668 |
| `and:view` | 0.1126 | 0.0727 | 0.0759 |
| `copy:flags` | 1.6887 | 1.6679 | 1.7084 |
| `copy:core` | 0.0717 | 0.0428 | 0.0434 |
| `copy:owned` | 0.0292 | 0.0015 | 0.0039 |
| `pack:core` | 2.4009 | 7.7948 | 7.6162 |
| `pack:assignment` | 4.0817 | 9.2297 | 9.5087 |

Merge, ms: `merge:flags` 1372.0, `merge:core` 36.5

### ifx -O3 — `--max-performance`

| ns per bit | 10^3 | 10^5 | 10^7 |
|---|---:|---:|---:|
| `and_count:flags` | 0.1401 | 0.1389 | 0.1741 |
| `and_count:core` | 0.0205 | 0.0071 | 0.0075 |
| `and_count:eager` | 0.1179 | 0.0141 | 0.0757 |
| `and_count:lazy_bits` | 7.3713 | 8.9533 | 9.2588 |
| `and_count:lazy_words` | 0.8041 | 0.1462 | 0.1359 |
| `and:flags` | 0.0270 | 0.0360 | 0.1788 |
| `and:core` | 0.0103 | 0.0036 | 0.0066 |
| `and:eager` | 0.0822 | 0.0070 | 0.0103 |
| `and:view` | 0.0134 | 0.0036 | 0.0068 |
| `copy:flags` | 0.0167 | 0.0206 | 0.0711 |
| `copy:core` | 0.0097 | 0.0013 | 0.0039 |
| `copy:owned` | 0.0108 | 0.0013 | 0.0039 |
| `pack:core` | 0.6536 | 3.9084 | 4.0723 |
| `pack:assignment` | 0.7086 | 4.0442 | 4.2587 |

Merge, ms: `merge:flags` 76.6, `merge:core` 3.7

A second full pass agreed cell for cell: median deviation 1.8%, 161 of 176 cells within 5%,
all within 10% (worst: `pack:assignment` at 10^3 bits, gfortran `-O0`, 9%).

### What the timings say

Ratios below are at `--max-performance`, gfortran / ifx.

**Lazy evaluation is the slowest way to do anything.** Evaluated per bit, `count(a .and. b)`
costs about 10 ns per bit: 900× / 1300× the core's fused `and_count` at 10^5 bits, and slower
even than `logical(c_bool)`. Evaluated per word, it is still 14× / 21× the core. What it saves
is one result mask, n_bits/8 bytes, which a fused procedure saves too, at full speed.

**The eager operator pays an allocation per call, and that dominates small masks.** At 10^3
bits `c = a .and. b` costs 1.9× / 8.0× the in-place core, and `count(a .and. b)` 1.5× / 5.8× the
fused one. Shatter calls these once per seed and growth iteration, so a fixed cost per call
adds up. At 10^7 bits ifx's eager `and_count` is 10×
the core again: a fresh 1.25 MB result on every call.

**The pointer component costs gfortran its vectorization.** `and:view` is 1.9–3.8× `and:core`
under gfortran; ifx is unaffected. The loops are identical except that `view_and_into` reaches
the words through pointer components, and gfortran's `-fopt-info-vec` reports
`core_and_into`'s loop as vectorized (32-byte vectors) and says nothing about `view_and_into`'s.

**Not every form of the type is slow.** The intrinsic `=` copy of `owned_mask` (`copy:owned`) is
as fast as the core's copy at `-O3`, and faster at `-O0`, where it is a `memcpy` against an
unoptimized loop. `pack:assignment` runs the same packing loop as `pack:core` and costs the same
at `-O3`. Neither decides anything: the copy is only possible with the allocatable component the
implementation layer forbids, and the assignment's problem is the crash, not its speed.

**Why bit masks at all.** Against `logical(c_bool)`, the core's `and_count` is 280× / 20× faster
at 10^5 bits, and Shatter's merge 215× / 21× (1183 ms against 5.5 ms under gfortran). gfortran
compiles `count(a .and. b)` and `c = a .and. b` over `logical(c_bool)` as a branch per element
(`-fopt-info-vec`: "not vectorized: unsupported control flow in loop"), which at density 0.5
mispredicts half the time; ifx avoids the branch (`and:flags` 0.036 ns per bit at 10^5). A mask also takes an eighth of the memory:
Shatter's 10^6 x 1000 result is 125 MB instead of 1 GB.

## Outcome

**`f42_bit_masks_impl` is a core of pure procedures on plain int32 words, with no mask type**
(FES, 2026-09-15). The type fails on rules before it fails on speed:

- a pure implementation cannot make a view of an `intent(in)` mask, so a view type and every
  lazy node are out;
- an owning pointer type is wrong on every compiler in one of its two forms;
- an allocatable component allocates on every operator and whole-mask assignment, which the
  implementation layer forbids and its allocation check cannot see;
- a defined assignment from a logical array can crash ifx.

And where the type's forms run at all, the operator allocation, the pointer component and lazy
evaluation are each measurably slower than the in-place and fused procedures that the core
offers.

For the core itself, the numbers set two rules:

- Whole-mask operations are in place (`and_into`, `copy_into`) or fused (`and_count`); nothing
  materialises a result nobody keeps.
- The packing loop here, bit by bit with a branch, costs 3.5–4 ns per bit. The module's pack
  should build each word branch-free.

## Notes on the benchmark itself

- Three compilation units, never `-flto`/`-ipo`, `cpu_time`, and every measurement repeated
  until it runs for at least 0.1 s, best of 3 (see [`../README.md`](../README.md)).
- Every form of an operation is checked against the `logical(c_bool)` result before it is timed.
- One bit of every left operand is flipped per repetition, so nothing can be hoisted.
- `run_bench.sh` runs the benchmark with `ulimit -s unlimited`, because `pack:assignment` at
  10^7 bits would otherwise crash ifx (the probe above).
- nvfortran is not timed: it is not a CI compiler, and the probes already rule the type out
  there. On this machine it links only with `/usr/bin` first on its PATH, which is why the
  probes were run with `NVFORTRAN="env PATH=/usr/bin:/opt/nvidia/hpc_sdk/Linux_x86_64/2025/compilers/bin nvfortran"`.
