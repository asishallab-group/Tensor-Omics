# Designing the kernel layer

Why the generator produces the `foo` and `foo_alloc` wrappers from a hand-written kernel, how
the kernel is annotated, and what was rejected on the way.

Companion documents: [`c-layer.md`](c-layer.md) for the Fortran/C wrappers,
[`language-layers.md`](language-layers.md) for Python and R, [`../README.md`](../README.md)
for the generator itself.

---

## The rule the layer follows

**The author writes only the kernel — the implementation — and annotates it. The generator
writes everything the kernel's signature and docs already imply: the validating wrapper and
the allocating wrapper.**

Everything below follows from that, plus one constraint: **the kernel is the single source of
truth.** A fact stated on the kernel is never restated by hand anywhere else; if the generator
cannot derive something, the kernel gains an annotation for it rather than the wrapper gaining
hand-written code.

---

## The three layers, and who writes them

The framework splits a procedure into three:

- `foo_kernel` — the implementation. **Hand-written.** Lives in `src/kernel/**/tox_*_kernel.F90`
  (module `tox_X_kernel`, procedure `X_kernel`). Not exported; takes no `ierr`.
- `foo` — validates its inputs, then calls the kernel. **Generated.**
- `foo_alloc` — allocates the kernel's work arrays, builds and sorts its permutations, calls
  its "recommend" sizing routines, then calls the kernel. **Generated.**

Both generated wrappers land in `src/tox/tox_X.F90` (`module tox_X`, `use tox_X_kernel`) and
carry `M_EXPORT_C`, so the ordinary binding pipeline wraps them to C/Python/R with no special
casing. The clean name — `tox_X`, `loess_fit`, without the `_kernel` suffix — *is* the
recommended public API.

### Why `_kernel`, and why a separate tree

The suffix was `_helper` and is now `_kernel`. `_helper` was overloaded: genuine helpers (the
serde flat-implementation shims, small private subroutines) also end in `_helper`. `_kernel`
names exactly one thing — the implementation that the generator turns into an API — and it is
unmissable at every call site, which matters because calling a raw kernel means calling
something with no input validation.

The suffix *and* the `src/kernel/` location together are the generation trigger. There is no
separate marker macro: a procedure is a generatable kernel iff it is a `_kernel` procedure in
`src/kernel/`. This keeps the trigger in the one place a reader already looks — the name.

**Rejected — a `M_KERNEL` marker macro.** It would restate what the suffix already says. The
suffix cannot be forgotten (it is part of every call); a marker can.

### Why f42 is untouched

`src/f42/` is generic infrastructure, not the Tensor Omics kernel. Whether an f42 procedure is
exported is a per-case judgement, so f42 keeps explicit `M_EXPORT_C` exactly as today. The
kernel-generation feature is tox-scoped; it reads only `src/kernel/`.

---

## `src/tox/` is generated-only

Because the wrappers land in `src/tox/`, the simplest and least error-prone arrangement is for
`src/tox/` to contain *nothing but* generated files. Then the generator treats it exactly like
`src/bindings/c`: Ford excludes the whole directory (`exclude_dir`), and `clean` removes the
whole directory (`rmtree`). No per-file markers, no glob that might spare or delete the wrong
thing.

**Rejected — a mixed `src/tox/` (generated wrappers beside a few hand-written modules).** It
would force per-file Ford exclusion and a marker-based clean that must delete only generated
files and never the hand-written ones — permanent complexity and a sharp edge in `clean`, in
exchange for avoiding a one-time reorg.

Reaching generated-only required relocating the modules in `src/tox/` that are neither kernels
nor wrappers:

| Module(s) | New home | Why it was misplaced |
|---|---|---|
| `tox_errors` | `src/f42/` | no exports; `use`d by f42 itself (72 sites) — it always sat *below* f42 in the dependency stack, so it was never tox application code |
| `tox_conversions` | `src/f42/` | no exports; C-interop glue for the binding layer (27 sites) |
| `tox_data_archive`, `tox_data_tools`, `tox_data_read_write` | `src/io/` | hand-written, exported, but file-I/O and external-library (`libzip`, netlib) code — not numeric kernels. Their bindings still auto-generate from `M_EXPORT_C` wherever they live |
| `tox_data_accessors` | `test/` (or delete) | unused anywhere in `src/` |
| recommend routines (`*_required_workspace`, `calc_*_size`, `calc_neighborhood_size`) | `src/kernel/` (public in the kernel module) | called by the generated allocs and by the `_expert` bindings |
| mode/enum params (`MODE_*`, `METHOD_*`, `*_PATTERN`) | `src/kernel/` (kernel module) | referenced by kernel signatures and `DM_REQUIRED_IF_MODE` |

The `tox_data_*` modules are hand-written for now; the intended end state is that they become
adapters over generic CSV and zip modules — a separate day's work.

---

## What the author annotates on the kernel

The generator already reads, and this feature reuses unchanged: explicit `intent`; kinded
types; per-argument `!!` docs; `tmp_<name>` work arrays; `<arg>_shape`, `n_selected_<arg>`,
`mode`/`method`; `DM_DEFAULT`, `DM_REQUIRED_IF_MODE`, `DM_OPTIONAL_OUTPUT`, `DM_RESULT_SIZE_IS`;
and `DM_OUTPUT_FROM(size_arg, recommend_proc, tox_X_kernel, AUTO)`.

### Validation ranges

Bounds are opt-in; **finiteness is not.** A real value that is NaN or infinite is a bug unless
the kernel says otherwise, so finiteness is the framework's default contract: every real
argument is checked, and an argument opts *out* per failure mode.

```
DM_MIN(EXPR)         The minimum valid value is `EXPR`.
DM_MAX(EXPR)         The maximum valid value is `EXPR`.
DM_SENTINEL(EXPR)    The value `EXPR` is additionally accepted.
DM_ALLOW_NAN         NaN is permitted for this value.
DM_ALLOW_INFINITE    Infinite values are permitted for this value.
```

These compose onto the existing `tox_errors` validators — the generator never writes an ad-hoc
check:

- Extent scalars → `validate_dimension_size` (automatic, from the argument's extent role; no
  macro needed).
- Every real data argument → `validate_in_range_real` / `validate_all_in_range_real`, which
  reject NaN/Inf by default. `DM_MIN`/`DM_MAX`/`DM_SENTINEL` add `min=`/`max=`/`sentinel=`
  (exclusive bounds via `above(0.0_real64)` / `below(…)`, written verbatim). `DM_ALLOW_NAN` /
  `DM_ALLOW_INFINITE` pass the corresponding new optional logical.
- Integer scalars/arrays → validated only if they carry a bound (integers cannot be NaN/Inf).
- `arg_pos` is the 1-based position of the argument in the *emitted* procedure's dummy list, so
  `foo` and `foo_alloc` number independently.

This requires one small extension to `tox_errors`: `validate_in_range_real` and
`validate_all_in_range_real` gain optional `allow_nan` / `allow_infinite` logicals, defaulting
`.false.` — i.e. today's behaviour, so nothing hand-written changes.

**Rejected — a `DM_FINITE` opt-in.** It puts the burden the wrong way round: the common case
(a value that must be finite) would need an annotation, and the rare, dangerous case (NaN
tolerated) would be silent. Finite-by-default with explicit opt-out states the tolerance where
the tolerance actually lives.

### Sorting and work arrays (`foo_alloc`)

`foo_alloc` exists only when the kernel needs work arrays, a pre-sorted permutation, or a
recommend-sized buffer; a kernel with none of these generates just `foo`. It drops from its
signature every argument the caller should not have to supply — `tmp_*` work arrays,
permutations, and recommend-computed sizes — allocates them locally, prepares them, and calls
the kernel.

Permutations follow a convention: an argument named `<base>_perm` is, in `foo_alloc`,
allocated, seeded with `init_perm`, and heapsorted against `<base>`. An explicit override
macro covers what the convention cannot express (a bare `perm`, or quicksort with its stack
arrays).

`foo_alloc` calls the **kernel directly**, not `foo`. It has just built the permutation and
work arrays itself, so re-running `foo`'s validation — in particular the O(n)
`validate_all_in_range` on a permutation it knows is `[1..n]` — would be wasted work.

**Rejected — `foo_alloc` → `foo`.** It is simpler for the generator (validation lives in one
place) and is what the tox families do today, but it pays an O(n) re-validation of data the
alloc layer just constructed correctly. The f42 families call the kernel directly for exactly
this reason; the generator follows them.

### Recommend / work-array sizing

A work-array size documented with `DM_OUTPUT_FROM(size_arg, recommend_proc, tox_X_kernel, AUTO)`
already tells the binding languages to call `recommend_proc` themselves. The same annotation
now tells `foo_alloc` to call it in Fortran: size the buffer, `M_ALLOCATE` it, call the kernel.
One annotation, two consumers — the `_expert` binding still exposes the size argument and calls
the routine; `foo_alloc` internalises it.

**Rejected — a new `RECOMMENDED` mode alongside `AUTO`.** `AUTO` already carries everything
needed (the producer and how to supply its inputs, resolved by the existing output-from pass).
A second mode would be a second name for the same fact.

### Mode-dependent kernels (opt-in split)

The house style is one procedure per mode — `detect_dosage_effect`,
`detect_subfunctionalization` — not a single procedure taking a runtime `mode`. A kernel opts
into this by giving its `mode`/`method` argument a mode table **with a procedure-name column**:
the generator then emits one `{plain, _alloc}` pair per mode value, names each from the table,
drops the `mode` dummy and fixes it, and treats each `DM_REQUIRED_IF_MODE(mode, …, Vᵢ)`
argument as a mandatory dummy in `Vᵢ`'s wrapper and absent from the others. Mode-independent
optionals stay optional in every wrapper.

Without that name column, nothing changes: a single runtime-mode procedure, with
`DM_REQUIRED_IF_MODE` keeping its "optional, required-in-mode" meaning. The column *is* the
opt-in — the per-mode names could not be mechanically derived anyway (`detect_dosage_effect` is
not spelled by `DOSAGE_PATTERN`), so the place that supplies the names is the natural place to
signal the split.

**Rejected — splitting every mode argument by default.** Some mode arguments are genuinely
runtime choices; forcing them into separate procedures would be wrong. Opt-in keeps both shapes
available. **Rejected — a separate split marker macro:** it would duplicate the signal the name
column already carries.

---

## How the generator builds it

Only the frontend reads Fortran; everything after it works on the in-memory IR, which is
hand-constructible by design. So the wrappers are **synthesised as IR and injected into the
project** before analysis — the C/Python/R pipeline then treats them like any procedure it
read, and a new Fortran emitter renders them to `src/tox/`.

1. **Directives.** The new `DM_*` macros are wired end-to-end (`src/macros.h`,
   `ir/directives.py`, `frontend/macros.py`) like the existing ones; the mode-table reader
   gains the optional procedure-name column; `tox_errors` gains the two `allow_*` logicals.
2. **Synthesis** (`synthesize.py`, run between parse and `analyse_project`). For each kernel —
   or once per mode value when it opts into mode-split — clone the kernel's arguments (carrying
   type, dimension, intent, optionality, docs, directives) into `foo` (+ synthesised `ierr`)
   and `foo_alloc` (minus the dropped arguments), give both `Meta(category = the export
   category)` and the frontend's conventions, and put both in one synthesised `module tox_X`.
   Placing the pair in one module makes the existing `_alloc`↔`_expert` naming produce
   `foo_alloc → foo_c` and `foo → foo_expert_c` with no special handling. A side table records,
   per wrapper, the kernel, the dropped arguments, and each permutation's base and sort — the
   recommend calls are read post-analysis off the resolved output-from plan.
3. **Emission** (`emit/fortran_wrapper.py`, modelled on `emit/fortran_c.py`). Renders
   `module tox_X` with its `use … only:` imports (kernel, recommend routines, the `tox_errors`
   validators actually used, and `f42_utils` only when permutations are present), the
   declaration-ordering hoist that keeps an extent declared before the array it sizes, the
   validation block, and — for `foo_alloc` — the recommend calls, `M_ALLOCATE`s, `init_perm`/
   sort, and the kernel call. These are ordinary library sources: no `NO_C_BINDING` guard.
4. **Orchestration.** A new `fortran` target and `_fortran_files` builder in `generate.py`;
   `kernel_src_dir` and `tox_out_dir` in `config.py`; `src/tox` added to Ford's `exclude_dir`
   so the generated wrappers are never read back; `clean` `rmtree`s `tox_out_dir`. `--check`
   diffs the emitted `.F90` against disk for free, like every other target.

**Rejected — a two-pass round-trip** (emit the wrappers, then re-run Ford over `src/tox` to
read them for binding generation). On the first run, or after any kernel change, the on-disk
wrappers are stale, so the bindings would be generated from the old API. Synthesising in memory
keeps the kernel the single source of truth and needs Ford only once.

**Rejected — a `tox_X_variants` module name.** The generated module takes the clean `tox_X`
name precisely because it is the API callers should use; a decorated name would push the
recommended entry point behind an odd import.

---

## Rollout

The generator feature is independent of converting existing families, so it lands first and is
proven on pilots before the bulk migration:

1. Build the feature (macros, synthesis, emitter, target, wiring) with unit tests and an
   end-to-end test that compiles and runs an emitted wrapper.
2. Relocate the six blocker modules (green build at each step), reaching generated-only
   `src/tox`.
3. Pilot A — a validation-only family (e.g. `tox_shift_vectors`).
4. Pilot B — an alloc + permutation + recommend family (e.g. `tox_get_outliers`).
5. Pilot C — a mode-split family (`tox_paralog_analysis`).
6. Full rollout, ending with `tox_data_integration`, whose parent-interface + submodule
   structure dissolves: the kernels become plain `src/kernel` procedures, the wrappers
   regenerate into one `src/tox/tox_data_integration.F90`.

Verification at each stage: the generator test suite plus `--check` idempotency; `./build.sh`
producing the `_c`/`_call` symbols for the wrappers in one library; and the full Python and R
suites after the relocations and pilots.
