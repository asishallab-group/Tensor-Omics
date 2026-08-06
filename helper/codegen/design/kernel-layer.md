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

Both generated wrappers land in `src/generated/tox/tox_X.F90` (`module tox_X`, `use tox_X_kernel`) and
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

### A family too big for one file: re-export, not submodules

A family whose kernels do not fit comfortably in one file splits into several kernel modules
gathered by a **parent that holds no procedures of its own and only `use`s its children** —
the shape `f42_serde_arrays_deserialize` already uses:

```fortran
module tox_data_integration_kernel
    use tox_data_integration_preprocessing_kernel
    use tox_data_integration_jsd_kernel
    ...
end module tox_data_integration_kernel
```

The generator mirrors it: a procedure-less kernel module that uses other kernel modules
generates the matching parent over the *wrappers* (`use tox_data_integration_jsd`, …), so
`use tox_data_integration` still reaches the whole family and the split stays an
implementation detail of the kernel tree. Only children that actually generate something are
re-exported — a kernel module of nothing but constants or recommend routines has no generated
counterpart to `use` — and a parent may gather another parent.

**Rejected — `submodule`s (what `tox_data_integration` used to be).** The parent had to
restate every procedure's full signature and documentation inside an `interface` block, and
each submodule then repeated the declarations a second time: two copies of every signature to
keep in step, for no benefit the plain split does not give. It also cost the generator a
special case, since Ford reports the implementing submodule as having no routines and the
interface declaration is the only place the signature can be read from. Nothing in the
framework needs submodules; a re-exporting parent is the same encapsulation with one copy of
each signature.

### Why f42 is untouched

`src/f42/` is generic infrastructure, not the Tensor Omics kernel. Whether an f42 procedure is
exported is a per-case judgement, so f42 keeps explicit `M_EXPORT_C` exactly as today. The
kernel-generation feature is tox-scoped; it reads only `src/kernel/`.

---

## The tree says who writes what

```
src/
  macros.h            included by every source, by this path
  f42/                infrastructure, library-agnostic
    utils/              f42_utils re-exports f42_math, _sort, _random, _vector, _stats
    serde/              likewise, per element type
  kernel/             the kernels -- the API's source of truth, hand-written
  data/               the hand-written data-set API (tox_data_*)
  generated/          NOTHING here is hand-written
    tox/                the wrappers
    bindings/c/         the Fortran C wrappers
    bindings/r/         the R `.Call` shims
```

A family gets a directory once it is more than one module, and its parent re-exports the
rest — the shape `f42/serde/` established, that `f42/utils/` now follows, and that the
generator mirrors for a split kernel family (above). A directory holding one family's single
module would name nothing the module name does not.

The rule a reader needs is one line: **edit anything outside `src/generated/`.** That is also
what `.gitattributes` marks, what a review can collapse, and what an editor can be told to
ignore — one path prefix rather than a list that has to be kept in step with the generator.

Reaching that required two things. The generated wrappers had to *become* the whole of their
directory, which meant relocating what used to sit in `src/tox` beside them:

| Module(s) | New home | Why it was misplaced |
|---|---|---|
| `tox_errors` | `src/f42/` | no exports; `use`d by f42 itself (72 sites) — it always sat *below* f42 in the dependency stack, so it was never tox application code |
| `tox_conversions` | `src/f42/` | no exports; C-interop glue for the binding layer (27 sites) |
| the `tox_data_*` family | `src/data/` | hand-written, exported, but not numeric kernels: the data-set API — archive, readers, validation, accessors. Their bindings still auto-generate from `M_EXPORT_C` wherever they live |
| recommend routines (`*_required_workspace`, `calc_*_size`, `calc_neighborhood_size`) | `src/kernel/` (public in the kernel module) | called by the generated allocs and by the `_expert` bindings |
| mode/enum params (`MODE_*`, `METHOD_*`, `*_PATTERN`) | `src/kernel/` (kernel module) | referenced by kernel signatures and `DM_REQUIRED_IF_MODE` |

The `tox_data_*` modules are hand-written for now; the intended end state is that they become
adapters over generic CSV and zip modules — a separate day's work.

Then the three generated trees moved under one root. Fortran locates modules by name and fpm
scans `src/` recursively, so no `use` statement changes — the move is a `git mv` plus the
`Paths` fields in `config.py`.

**Rejected — a mixed generated/hand-written directory.** Sorting a directory's files by who
owns them is a per-file rule, and every tool that has to respect it (Ford's exclusion, `clean`,
review, search) needs its own copy. A path prefix is one rule that all of them already
understand.

**Rejected — grouping the hand-written trees under a `tox/` parent** (`src/tox/kernel`,
`src/tox/io`), leaving `src/f42`, `src/tox` and `src/generated` at the top. Tidier on paper,
but it renames the tree the author edits most, and `src/tox` would come to mean hand-written
code after years of meaning the opposite. The generated/hand-written split is the one that
earns its churn; this one does not.

**Rejected — deleting `src/generated` wholesale in `clean`, now that it is generated-only.**
`clean` removes exactly the files the generator is about to write, derived from the kernel
tree, and Ford excludes that same list. One list, two consumers, no way for them to disagree
about what is generated. `rmtree` would be a second, independent answer to the same question.

---

## What the author annotates on the kernel

The generator already reads, and this feature reuses unchanged: explicit `intent`; kinded
types; per-argument `!!` docs; `tmp_<name>` work arrays; `<arg>_shape`, `n_selected_<arg>`,
`mode`/`method`; `DM_DEFAULT`, `DM_REQUIRED_IF_MODE`, `DM_RESULT_SIZE_IS`;
and `DM_OUTPUT_FROM(size_arg, recommend_proc, tox_X_kernel, AUTO)`.

### The argument position an `ierr` carries

An `ierr` packs the position of the argument it blames (`arg_pos * 10000 + code`). A kernel may
report a genuine runtime error, and when it does the generated wrapper **clears that position**
before returning — `call clear_err_arg_pos(ierr)` after every kernel call, and after every
prologue call. The code survives; the position does not.

The position is the *kernel's* numbering, and the wrapper's dummy list is a different list: it
drops work arrays, permutations and recommend-sized values, appends `ierr` last, and a
mode-split wrapper drops the mode argument too. Worse, a position propagates unchanged through
every call that does not rewrite it, so what arrives is often not even the kernel's own
numbering but that of some private helper three frames down — `compute_family_scaling` returns
`calc_percentile`'s positions 1–3, which in *its* signature are `n_genes`, `n_families` and
`distances`. Reporting any of that to a caller is worse than saying nothing, and "not argument
related" is what position 0 means.

This is safe because `ierr` is provably OK at the moment of the call: validation ends with
`if (is_err(ierr)) return`, and so do the prologue and recommend calls. Nothing the clear
discards can be a position the wrapper itself set. **A change that lets validation fall through
without returning would break that, and silently.**

The cost is real but small: six kernel-side positions are correct today and become 0. Three are
already unreachable from Python and R, because the C wrapper converts the mode *string* and
rejects an unknown one before Fortran is entered; only a direct Fortran caller ever saw them.

**Rejected — remapping the kernel's position onto the wrapper's by name.** This is what the
hand-written per-mode wrappers did, in runs of a dozen `map_err_arg_pos` calls apiece. It
assumes the packed position refers to the kernel's own dummy, which is false wherever the
position propagated in from a helper — the majority of sites.

**Rejected — keeping the position when the leading arguments coincide.** `detect_dosage_effect`'s
first four arguments *are* the kernel's first four, so a helper-borrowed 3 or 4 would survive
looking perfectly plausible as `n_genes` or `n_dims`. There is no prefix a generator can trust.

**Rejected — just deleting `arg_pos=` from the kernels.** It cannot reach a position that
propagates in from f42, which is exactly the `compute_family_scaling` case above. The kernels
keep their annotations: they are correct for anyone calling a kernel directly, and the test
suite does.

Two of the six are recovered rather than lost: a mode argument the wrapper still takes is now
checked against the values its own table names (see below), which puts the position back — the
wrapper's own this time, not the kernel's.

`map_err_arg_pos` stays for hand-written code where the two dummy lists really are known to
correspond. The recommend-routine call site would need the same clear, but no recommend routine
reachable from generated code takes an `ierr` today, so nothing is emitted there.

### Mode membership

A mode argument the wrapper still takes — one whose table has no per-mode procedure column, so
it was not split — is checked against exactly the values that table names:

```fortran
if (baseline_mode /= MODE_RAW .and. baseline_mode /= MODE_MEAN .and. baseline_mode /= MODE_MIN) &
    call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=4_int32)
```

The generator parsed the table, so it knows the accepted set exactly. An optional mode is
guarded by `present()`: absent means "use the documented default", which is one of the accepted
values by construction, and reading an absent optional to discover that is not a Fortran program.

A chain of `/=` rather than `all(x /= [...])`, so the check costs a few integer comparisons and
no temporary array.

**Rejected — annotating the mode argument with `DM_MIN`/`DM_MAX`.** The accepted set need not be
contiguous, and a hand-written range duplicates knowledge the generator already holds and goes
stale the moment a mode is added, silently admitting the new value everywhere it was not updated.

**Rejected — leaving it to the kernel.** The kernel's own `case default` is what used to report
this, and its `arg_pos` is in the kernel's numbering — which the wrapper now clears. The check
has to be the wrapper's for the position to mean anything. A direct Fortran caller is the one
this protects: the C layer already rejects an unknown mode *string* before Fortran is entered.

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

**Rejected — `DM_VALIDATE`, a hook for a kernel's own validation helper.** The proposal was an
argument-level directive letting the generated wrapper call a hand-written validation routine
inside the kernel, passing *the wrapper's* position for that argument — the one number the
kernel cannot know, and the reason `clear_err_arg_pos` has to zero what a kernel returns. It had
exactly two customers, both the same check: `clock_hand_angle_between_vectors` needed its three
axis indices to be mutually distinct, a relation between an argument and itself that no
per-argument bound can state.

That check is gone, because the argument that needed it was the wrong argument. The signed angle
now takes an orientation *vector* rather than three axis indices, so the rule is a bound like any
other. The general lesson, and the reason this stays rejected rather than deferred: **a procedure
that needs a validation nothing in the contract can express is a procedure whose signature is
wrong.** The generator's vocabulary is the API review. If a genuine new kind of constraint ever
turns up, it earns a *convention* — a directive with a name, a meaning, and every kernel that
qualifies using it — not a per-procedure escape hatch that hides one bad signature and stops
anyone from noticing the next.

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

**Rejected — chaining producers, so one recommend routine feeds another.** The motivating case
is a workspace whose exact size depends on a count the kernel only discovers at run time:
`normalize_by_std_dev` fits LOESS through the `n_valid` genes that carry variance, not all
`n_genes`, so the exact workspace needs an `n_valid` producer whose output becomes the
`max_neighborhood_size` input of `tox_loess_required_workspace`. Supporting that means the
generator resolves a dependency graph over producers and emits them in topological order,
allocating each producer's own inputs before its call and the buffers that depend on it after —
and it still cannot help a producer whose input is itself a `tmp_` buffer that does not exist
yet.

The kernels do not need it. Every such size has a cheap upper bound that is already an argument,
and the recommend formulas are monotone in it, so the buffer is sized for the bound and sliced at
the call: `compute_family_scaling_kernel` sizes for `n_families` and passes
`loess_x(1:n_valid)`, and `normalize_by_std_dev_kernel` does the same with `n_genes`. What is
over-allocated is the fraction of points the fit drops; what is bought is that `foo_alloc` stays
a flat sequence of calls and allocations. Should a case appear where the bound is genuinely far
from the truth, the topological version is the answer — not before.

### The prologue is the allocating tier's sugar, and has no scope

`DM_PROLOGUE(PROCEDURE, MODULE)` names a routine the *allocating* wrapper runs after the work
arrays are prepared and before the kernel; it may write the outputs and report `handled`, and
the kernel is then skipped. It is what `foo_alloc` derives that `foo` lets a caller pass in --
the same relation the `<base>_perm` convention already has, where `foo_alloc` seeds and
heapsorts and a caller wanting another sort reaches for `foo`. A prologue is that convention's
general form.

**Rejected -- a `SCOPE` of EXPERT / ALLOC / BOTH.** The macro carried one, and two of its three
values were mistakes:

- `BOTH` is inlinable. Both wrappers call the kernel, so anything that must run in both, before
  the kernel body, *is* the first thing in the kernel body -- `handled` becomes a local and the
  early return a plain statement. The prologue has no access the kernel lacks, so the directive
  bought nothing. `loess_degenerate_fit` was the only user, and it moved into both LOESS kernels;
  the third caller (`normalize_by_std_dev_inplace_helper`) had been performing the same check by
  hand, which is what a precondition-every-caller-must-remember looks like when it should have
  been the procedure's own contract.
- `EXPERT` contradicts what the expert tier is. FES: *expert is the entry point for full control
  over what reaches the kernel -- a specific threshold, a specific initialised permutation --
  while `_alloc` derives the threshold from a percentile and sorts with heapsort.* A prologue
  running in the expert tier would override exactly the control that tier exists to give.

So the scope is not derived, it is fixed, and the placement with it: always in `foo_alloc`,
always below the work arrays. A kernel with no work arrays generates no `foo_alloc`, so a
prologue on one would never run -- an error, not a silence.

**Nothing validated a prologue before this**: the emitter was the directive's only consumer, and
an emitter has no author's line to point at. Every way of getting one wrong therefore passed --
a name that resolved to nothing produced a wrapper with no prologue; a dummy that matched nothing
was dropped from the keyword call; a prologue without `handled` left the wrapper branching on an
undefined logical. All are errors now, in `ir/validate.py`.

**Rejected -- a rename table for prologue dummies.** A `DM_OUTPUT_FROM` producer needs one
because it is published and its parameter names cannot move. A prologue is internal to the kernel
module, so a mismatch is fixed by renaming its dummy, and the diagnostic says so.

### No allocation in a kernel module

Nothing in a kernel module allocates: every buffer is a `tmp_` argument, so the generated
`_alloc` owns the memory and an expert caller can hand in buffers it already holds. The rule
covers the module's helpers as well as its kernels — a kernel that allocates nothing itself but
calls a helper that does gives its expert caller nothing.

It is enforced on the **declaration**: a local declared `allocatable` is refused. The generator
never reads a body, and an `M_ALLOCATE` needs an allocatable to allocate into, so the declaration
is a complete proxy for this tree. (A `pointer` local that is `allocate`d would slip through. No
kernel does that, and the gap is cheaper than teaching the frontend to read statements.) Ford
hides a procedure's own variables by default, so the frontend sets `proc_internals`.

A `pointer` local is explicitly fine. Aliasing a buffer the kernel was handed allocates nothing,
and the `target` attribute it needs on a dummy never reaches the wrapper — the emitter carries
intent, type and dimension, not attributes. `normalization_pipeline_kernel` relies on both: it
points its LOESS scratch at the still-unwritten columns of its own output buffer.

Two naming rules fall out of the same place. A kernel may not carry `M_EXPORT_C` — its wrapper is
the entry point, and exporting it beside the wrapper publishes an unvalidated twin under a name
(`foo_kernel`) a binding caller cannot tell apart from `foo`. And a kernel may not be named
`<x>_alloc_kernel`: `_alloc` is the generator's suffix, and such a kernel generates a `foo_alloc`
that allocates nothing, with no expert tier behind it. Support routines in a kernel module — the
recommend routines above — are neither, and stay exported.

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
read, and a new Fortran emitter renders them to `src/generated/tox/`.

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
   `kernel_src_dir` and `tox_out_dir` in `config.py`; `src/generated/tox` added to Ford's `exclude_dir`
   so the generated wrappers are never read back; `clean` `rmtree`s `tox_out_dir`. `--check`
   diffs the emitted `.F90` against disk for free, like every other target.

**Rejected — a two-pass round-trip** (emit the wrappers, then re-run Ford over `src/generated/tox` to
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
   `src/generated/tox`.
3. Pilot A — a validation-only family (e.g. `tox_shift_vectors`).
4. Pilot B — an alloc + permutation + recommend family (e.g. `tox_get_outliers`).
5. Pilot C — a mode-split family (`tox_paralog_analysis`).
6. Full rollout, ending with `tox_data_integration`, whose parent-interface + submodule
   structure dissolves: the kernels become plain `src/kernel` procedures, the wrappers
   regenerate into one `src/generated/tox/tox_data_integration.F90`.

Verification at each stage: the generator test suite plus `--check` idempotency; `./build.sh`
producing the `_c`/`_call` symbols for the wrappers in one library; and the full Python and R
suites after the relocations and pilots.
