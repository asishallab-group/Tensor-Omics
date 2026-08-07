# Designing the implementation layer

Why the generator produces the `foo` and `foo_expert` wrappers from a hand-written
implementation, how the implementation is annotated, and what was rejected on the way.

Companion documents: [`c-layer.md`](c-layer.md) for the Fortran/C wrappers,
[`language-layers.md`](language-layers.md) for Python and R, [`../README.md`](../README.md)
for the generator itself.

---

## The rule the layer follows

**The author writes only the implementation and annotates it. The generator writes everything
the implementation's signature and docs already imply: the entry point a caller reaches for
first, and — where that entry point has anything to prepare — the expert one beside it.**

Everything below follows from that, plus one constraint: **the implementation is the single
source of truth.** A fact stated on the implementation is never restated by hand anywhere
else; if the generator cannot derive something, the implementation gains an annotation for it
rather than the wrapper gaining hand-written code.

---

## The tiers, and who writes them

The framework splits a procedure into an implementation and the wrappers it implies:

- `foo_impl` — the implementation. **Hand-written.** Lives in a module `tox_X_impl`, in a file
  of that same name. Not exported; takes an `ierr` only where it can genuinely fail.
- `foo` — validates its inputs, calls the "recommend" sizing routines, allocates the work
  arrays, seeds and heapsorts the permutations, runs the prologue, then calls the
  implementation. **Generated.**
- `foo_expert` — validates its inputs and calls the implementation with what the caller
  supplies. Allocates nothing, prepares nothing. **Generated only where the two signatures
  would differ.**

**The plain name always goes to the entry point a caller should reach for first.** Where the
implementation has work arrays, permutations or `DM_OUTPUT_FROM(…, AUTO)` values for `foo` to
take over — or a prologue that asks for an argument of its own — the two signatures differ,
both wrappers are generated, and the validate-and-call one becomes `foo_expert`. Where neither
happens the two would be the same procedure twice, so a *single* wrapper is generated and it
is called `foo`; there is no `foo_expert`, because there is no second tier for the name to
distinguish it from. That is 39 of the 87 procedures the tox tree generates today: an
implementation does not imply a pair, it implies at most one.

The pair is a **sugar/control** split rather than a fast/slow one: `foo` derives what
`foo_expert` lets a caller pass in -- the heapsorted permutation, a prologue's threshold or
short-circuit, the recommend-computed workspace sizes. That is why a prologue belongs to `foo`
alone (below), and why the expert tier is only published to Python and R where `foo` does one
of those things: elsewhere those languages allocate the work arrays for both tiers, so
`foo_expert` would be the same call under a name promising control it cannot give. Fortran and
C always get both, because there the expert tier really does hand the buffers over. Two of the
24 expert tiers generated today are published to Python and R, and both are the ones that sort
a permutation.

The wrappers land in the mirror of the implementation's own path — `src/tox/tox_X_impl.F90`
generates `src/generated/tox/tox_X.F90` (`module tox_X`, `use tox_X_impl`) — and carry the
export category `M_EXPORT_C` stands for, so the ordinary binding pipeline wraps them to
C/Python/R with no special casing. The category is set on the synthesised IR rather than
written into the emitted source: nothing ever parses the generated tree again, so a marker in
it would be decoration. The clean name — `tox_X`, `loess_fit`, without the `_impl` suffix —
*is* the recommended public API.

### Why the plain name belongs to the tier that prepares

It did not always. The generated pair used to be `foo` for the validating wrapper and
`foo_alloc` for the allocating one, on the argument that a suffix should say what a wrapper
*does* and that the bare name should therefore go to the tier that assumes least. That
argument is superseded, and it was reasoning about the wrong thing. A plain name is not a
description, it is a recommendation: it is what a programmer types first, what an example
shows, and what a reader takes for the ordinary way in. The ordinary way in is the tier that
prepares what the implementation needs, so that tier takes the plain name and the one for a
caller who wants to prepare it themselves says so — `foo_expert`.

The rename also ends a discrepancy the framework had carried for as long as it has had
bindings. Python and R have always published `foo` and `foo_expert`;
`abi/c_abi.py:stripped_name` produced exactly that from a Fortran `foo`/`foo_alloc` pair, so
the two halves of the same procedure went by different names depending on which language you
read. Now Fortran says what Python and R say, and one name means one thing in all four.

The consequence worth recording is what that does to `stripped_name`: for generated code it is
now the **identity**. Synthesis names the tiers `foo` and `foo_expert` outright, and the C
symbols come out `foo_c` and `foo_expert_c` by adding a suffix rather than by translating a
name. The translating branch survives only for f42's hand-written pairs (below) and retires
itself when those are converted. Because the translation was already producing these names,
the exported C symbols and the published Python and R functions came through the rename
unchanged, byte for byte; the whole of it is internal to Fortran and to the generator.

**Rejected — naming the lone wrapper `foo_expert` unconditionally.** It is the tidier rule
(one implementation, one suffix, no case analysis) and it would generate, compile and bind
perfectly well. It would also rename 39 of the 87 entry points and every C symbol under them,
in a change that is supposed to publish nothing new — and rename them *towards* a name
promising control over a tier that does not exist. `synthesize._wrappers_for` therefore asks
whether there is anything to take over, and only then adds the suffix.

**Rejected — diagnostics for the retired names.** A check that recognised the old suffixes and
pointed at their replacements is what a released API would owe its callers. Nobody outside
this branch ever wrote them: the names were introduced and retired between two of its commits,
so the diagnostic's whole audience is the branch's own history. The rename is a rename.

### Why `_impl`, and why the module name is the whole trigger

The suffix was `_helper` before it was `_impl`, and `_helper` was overloaded: genuine helpers
(the serde flat-implementation shims, small private subroutines) end in it too. `_impl` names
exactly one thing — the implementation the generator turns into an API — and it is unmissable
at every call site, which matters because calling one directly means calling something with no
input validation.

The **module name is the whole trigger**: a module `tox_X_impl` generates `tox_X`, and the
procedures in it whose names carry the same suffix are what it generates from. It used to be
the suffix *and* the one directory implementations were required to sit in. The directory half
is gone; where a file sits now decides only where its wrappers are *written*, by a mechanical
mirror (`synthesize.generated_path_for`):

```
src/<rest>/<module>_impl.F90   ->   src/generated/<rest>/<module>.F90
```

so `src/tox/tox_loess_impl.F90` generates `src/generated/tox/tox_loess.F90`, and
`src/tox/data_integration/tox_data_integration_jsd_impl.F90` generates
`src/generated/tox/data_integration/tox_data_integration_jsd.F90`. Sub-directories are
mirrored rather than flattened, which they used to be, so a family that earned a directory in
the source tree keeps it in the generated one.

The fixed directory was dropped so that a layer other than tox can write an implementation
without the generator learning a second case. Under the old rule an f42 implementation was
either impossible or a special case in the placement, the Ford exclusion and the cleaner
alike; under the mirror it is an ordinary file whose wrappers appear next to their own layer.
That the tox layer's implementations all happen to live in `src/tox` is now an observation
about the tree rather than a rule in the generator.

Two consequences follow, and both are enforced:

- **`_impl` is a reserved module suffix across the whole of `src/`.** Any module so named
  acquires generated wrappers and is held to the implementation-module rules — allocates
  nothing, not exported, its prologue checked (`validate.validate_impl_module`). That is the
  price of dropping the directory the rules used to be scoped to, and it is the cheaper half
  of the trade: a rule that follows a name is one a reader can apply without knowing where
  they are in the tree.
- **An implementation module lives in a file named after it**
  (`validate._check_module_is_named_for_its_file`). Two independent derivations meet here: the
  wrapper's *contents* come from the module name, which synthesis reads, while its *path*
  comes from the file name, which the cleaner scans without parsing anything. Let the two
  disagree and the generator emits `module tox_a` into a file called `tox_b.F90` — legal
  Fortran, no error anywhere, and a tree in which nothing can be found by the name it has.

**Rejected — a marker macro on the implementation.** It would restate what the suffix already
says. The suffix cannot be forgotten, because it is part of every call to the thing; a marker
can.

---

## What an implementation may reach

An implementation may `use` another `_impl` module, or one of a whitelist held in
`Conventions.impl_import_whitelist`: the intrinsic modules, `tox_errors`, `tox_conversions`,
`f42_config`, `f42_safeguard`, `f42_utils`. Anything else is refused
(`validate._check_impl_imports`).

The rule exists because **it is what makes the allocation rule mean anything beyond one
file.** `_check_impl_allocates` reads declarations — it sees a procedure that allocates for
itself, and nothing else. An implementation that allocates nothing but calls a helper in some
other module that does would pass it, and the guarantee an expert caller relies on (that
handing in buffers really avoids the allocation) would quietly not hold. Bounding the imports
is the only place that can be seen: another `_impl` module is held to the same rule, and the
whitelist is a curated rest.

It also fixes the **direction** of the dependency, which nothing else did. An implementation
could `use` a *generated* wrapper — a layering inversion, since the wrapper is generated from
implementations, and within one family a module cycle that surfaces as a build error naming
neither the cause nor the rule it broke.

**A procedure-level `use` counts.** Fortran allows one inside a procedure as well as in the
module header, and a rule that reads only the header would be one indentation level away from
being bypassed. `Procedure.uses` exists for this; the frontend fills it from the same Ford
attribute it already read for modules.

**It is a curated boundary, not a proof.** `f42_utils` re-exports `f42_stats`, whose
hand-written `_alloc` halves allocate; the other four children (`f42_sort`, `f42_math`,
`f42_random`, `f42_vector`) allocate nothing at all. Whitelisting only those four was
considered and rejected: implementations really do call `calc_percentile`,
`calc_percentile_helper` and `compute_scaled_distance_quantile`, all of which live in
`f42_stats`. So the entry stands, and it stops being an exception when f42 is converted to
`_impl` — at which point its allocating halves are generated rather than written.

**Rejected — having the generator verify each whitelisted module is allocation-free.** It
would turn the list from an assertion into a claim the build tests, which is the better shape.
It would also fail on `f42_stats` from the first run, so it would ship as a standing warning
about work that is already scheduled. Worth revisiting after the f42 conversion, when it would
pass and would then keep the list honest.

**An `allocatable` dummy is an allocation too**, and is refused with the locals. It looks like
the caller's memory — and the caller does receive it — but whoever fills it is the one who
allocated it, and the expert tier's promise is broken just the same. A `tmp_` argument
expresses the same intent without an allocatable in the signature, which the C layer could not
carry across the ABI in any case.

**Rejected — a "layer" in the output path** (`src/<layer>/**` → `src/generated/<layer>/`,
naming `tox` and `f42` as the layers). It is the same mirror for every file that exists today
and reads as more intentional. But mirroring the *whole* relative path is total: it has an
answer for a file directly under `src/`, for a family three directories deep, and for a layer
nobody has invented yet, and it needs no list to be kept in step with the tree. A rule with an
edge case is a rule someone eventually lands on.

### A family too big for one file: re-export, not submodules

A family whose implementations do not fit comfortably in one file splits into several
implementation modules gathered by a **parent that holds no procedures of its own and only
`use`s its children** — the shape `f42_serde_arrays_deserialize` already uses:

```fortran
module tox_data_integration_impl
    use tox_data_integration_preprocessing_impl
    use tox_data_integration_jsd_impl
    ...
end module tox_data_integration_impl
```

The generator mirrors it: a procedure-less implementation module that uses other
implementation modules generates the matching parent over the *wrappers*
(`use tox_data_integration_jsd`, …), so `use tox_data_integration` still reaches the whole
family and the split stays an implementation detail. Only children that actually generate
something are re-exported — an implementation module of nothing but constants or recommend
routines has no generated counterpart to `use` — and a parent may gather another parent. The
parent's own file mirrors like any other, so the generated family lands in one directory
alongside the parent that gathers it.

**Rejected — `submodule`s (what `tox_data_integration` used to be).** The parent had to
restate every procedure's full signature and documentation inside an `interface` block, and
each submodule then repeated the declarations a second time: two copies of every signature to
keep in step, for no benefit the plain split does not give. It also cost the generator a
special case, since Ford reports the implementing submodule as having no routines and the
interface declaration is the only place the signature can be read from. Nothing in the
framework needs submodules; a re-exporting parent is the same encapsulation with one copy of
each signature.

### f42: the mechanism is available, the conversion is deferred

Nothing scopes the generator to tox any more. `src/f42/utils/f42_stats_impl.F90` would
generate `src/generated/f42/utils/f42_stats.F90` with no special case anywhere — that is the point
of the mirror. f42 has nevertheless **not** been converted here. It still writes its two tiers
by hand as `compute_edf` / `compute_edf_alloc`, and `stripped_name` still translates that
hand-written shape into the published `compute_edf_expert` / `compute_edf`.

So `_alloc` remains a live convention, for **hand-written pairs only**. The generator never
emits it, and an implementation may not be named for a wrapper either:
`_check_impl_is_not_named_for_a_wrapper` refuses `foo_alloc_impl` and `foo_expert_impl` alike.
The first would generate a `foo_alloc` whose published Python name is `foo`, colliding with a
real `foo` in the same family; the second would put two procedures called `foo_expert` in one
module, and the emitter would strip the suffix and call `foo_impl` from the wrong one — wrong
code that compiles, because `foo_impl` exists.

Converting f42's four hand-written pairs is deliberately a **later, separate job**, and the
cost is what makes it one. Its pairs sit in *mixed* modules: `f42_stats` also holds
`loess_smooth_2d`, `compute_scaled_distance_quantile` and other ordinary exports, and
`f42_kd_tree` the same. A generated module is whole-file — everything in `f42_stats_impl`
would generate into `src/generated/f42/utils/f42_stats.F90` — so each conversion is first a split of
the module into the part that becomes an implementation and the part that stays a hand-written
export. That is a change to f42's public structure and belongs in its own change, with its own
review.

**Rejected — deleting `alloc_suffix` with the rest of the old naming.** The generator no
longer emits `_alloc`, so the field looks like dead vocabulary. It is not: it is how those
four hand-written pairs are *recognised*, and with it gone `compute_edf_alloc` would publish
as `compute_edf_alloc` and `compute_edf` as `compute_edf` — a silent rename of four published
Python and R functions in a change whose whole claim is that it publishes nothing new. The
field stays, documented as the hand-written convention it now is, and goes when f42 goes.

---

## The tree says who writes what

```
src/
  macros.h            included by every source, by this path
  f42/                infrastructure, library-agnostic
    utils/              f42_utils re-exports f42_math, _sort, _random, _vector, _stats
    serde/              likewise, per element type
  tox/                the tox implementations -- the API's source of truth, hand-written
    data_integration/   a family too big for one file
  data/               the hand-written data-set API (tox_data_*)
  generated/          NOTHING here is hand-written
    tox/                the wrappers, mirroring src/tox
    bindings/c/         the Fortran C wrappers
    bindings/r/         the R `.Call` shims
```

A family gets a directory once it is more than one module, and its parent re-exports the
rest — the shape `f42/serde/` established, that `f42/utils/` now follows, and that the
generator mirrors on both sides for a split family (above). A directory holding one family's
single module would name nothing the module name does not.

The rule a reader needs is one line: **edit anything outside `src/generated/`.** That is also
what `.gitattributes` marks, what a review can collapse, and what an editor can be told to
ignore — one path prefix rather than a list that has to be kept in step with the generator.

Reaching that required two things. The generated wrappers had to *become* the whole of their
directory, which meant relocating what used to sit in `src/tox` beside them:

| Module(s) | New home | Why it was misplaced |
|---|---|---|
| `tox_errors` | `src/f42/` | no exports; `use`d by f42 itself (72 sites) — it always sat *below* f42 in the dependency stack, so it was never tox application code |
| `tox_conversions` | `src/f42/` | no exports; C-interop glue for the binding layer (27 sites) |
| the `tox_data_*` family | `src/data/` | hand-written, exported, but not numeric implementations: the data-set API — archive, readers, validation, accessors. Their bindings still auto-generate from `M_EXPORT_C` wherever they live |
| recommend routines (`*_required_workspace`, `calc_*_size`, `calc_neighborhood_size`) | `src/tox/` (public in the implementation module) | called by the generated allocating wrappers and by the `_expert` bindings |
| mode/enum params (`MODE_*`, `METHOD_*`, `*_PATTERN`) | `src/tox/` (the implementation module) | referenced by implementation signatures and `DM_REQUIRED_IF_MODE` |

The `tox_data_*` modules are hand-written for now; the intended end state is that they become
adapters over generic CSV and zip modules — a separate day's work.

Then the three generated trees moved under one root, and the implementations moved into the
name that emptied — `src/tox`, which now means the hand-written tox tree and nothing else.
Fortran locates modules by name and fpm scans `src/` recursively, so neither move touched a
`use` statement: each is a `git mv` plus a field in `config.py`, where `generated_dir` is now
the single place the whole rule above is spelled. (What did rewrite the `use` lines was the
suffix on the implementation modules, which is a rename and not a relocation.)

**Rejected — a mixed generated/hand-written directory.** Sorting a directory's files by who
owns them is a per-file rule, and every tool that has to respect it (Ford's exclusion, `clean`,
review, search) needs its own copy. A path prefix is one rule that all of them already
understand.

**Rejected — grouping the hand-written trees under a `tox/` parent** (`src/tox/impl`,
`src/tox/io`), leaving `src/f42`, `src/tox` and `src/generated` at the top. Tidier on paper,
but it buys a level of nesting for two directories and pushes the tree the author edits most
one step further down. What the relocation did instead was reuse the freed `src/tox` for the
implementations themselves — the objection that `src/tox` would come to mean hand-written code
"after years of meaning the opposite" expired the moment the wrappers moved to
`src/generated/tox`, and the name was then the obvious one for what is left. `src/data` stayed
a sibling.

**Rejected — deleting `src/generated` wholesale in `clean`, now that it is generated-only.**
`clean` runs per target, and the generated tree holds three targets' output; the Fortran
target removes exactly the files it is about to write, one per implementation module, derived
from the source tree by the same mirror that places them. An `rmtree` of the root would make
`--target fortran` delete the C and R bindings it was not asked about.

**Ford excludes the generated tree by directory, not by file name.** It used to name the files
— the same derived list `clean` uses — which is one fact with one definition and looks
strictly better. It is not, for two reasons the layer-free mirror created. Ford matches a bare
exclude as `**/<name>`, so once an f42 implementation generates `f42_stats.F90`, excluding
that name would also drop the hand-written `src/f42/utils/f42_stats.F90` from the parse, and
every binding generated from it would vanish with no error at all. And a list derived from the
implementations that *still exist* cannot name the stale wrapper of one that was deleted,
which is exactly the file that must not be read back. One directory covers both, and covers
every target at once.

---

## What the author annotates on the implementation

The generator already reads, and this feature reuses unchanged: explicit `intent`; kinded
types; per-argument `!!` docs; `tmp_<name>` work arrays; `<arg>_shape`, `n_selected_<arg>`,
`mode`/`method`; `DM_DEFAULT`, `DM_REQUIRED_IF_MODE`, `DM_RESULT_SIZE_IS`;
and `DM_OUTPUT_FROM(size_arg, recommend_proc, tox_X_impl, AUTO)`.

### The argument position an `ierr` carries

An `ierr` packs the position of the argument it blames (`arg_pos * 10000 + code`). An
implementation may report a genuine runtime error, and when it does the generated wrapper
**clears that position** before returning — `call clear_err_arg_pos(ierr)` after every
implementation call that takes an `ierr`, and after every prologue call. The code survives; the
position does not.

The position is the *implementation's* numbering, and the wrapper's dummy list is a different
list: it drops work arrays, permutations and recommend-sized values, appends `ierr` last, and
a mode-split wrapper drops the mode argument too. Worse, a position propagates unchanged
through every call that does not rewrite it, so what arrives is often not even the
implementation's own numbering but that of some private helper three frames down —
`compute_family_scaling_impl` returns `calc_percentile`'s positions 1–3, which in *its*
signature are `n_genes`, `n_families` and `distances`. Reporting any of that to a caller is
worse than saying nothing, and "not argument related" is what position 0 means.

This is safe because `ierr` is provably OK at the moment of the call: validation ends with
`if (is_err(ierr)) return`, and so do the prologue and recommend calls. Nothing the clear
discards can be a position the wrapper itself set. **A change that lets validation fall through
without returning would break that, and silently.**

The cost is real but small: six positions that are correct in the implementations today become
0. Three are already unreachable from Python and R, because the C wrapper converts the mode
*string* and rejects an unknown one before Fortran is entered; only a direct Fortran caller
ever saw them.

**Rejected — remapping the implementation's position onto the wrapper's by name.** This is what
the hand-written per-mode wrappers did, in runs of a dozen `map_err_arg_pos` calls apiece. It
assumes the packed position refers to the callee's own dummy, which is false wherever the
position propagated in from a helper — the majority of sites.

**Rejected — keeping the position when the leading arguments coincide.** `detect_dosage_effect`'s
first four arguments *are* the implementation's first four, so a helper-borrowed 3 or 4 would
survive looking perfectly plausible as `n_genes` or `n_dims`. There is no prefix a generator
can trust.

**Rejected — just deleting `arg_pos=` from the implementations.** It cannot reach a position
that propagates in from f42, which is exactly the `compute_family_scaling_impl` case above.
The implementations keep their annotations: they are correct for anyone calling one directly,
and the test suite does.

Two of the six are recovered rather than lost: a mode argument the wrapper still takes is now
checked against the values its own table names (see below), which puts the position back — the
wrapper's own this time, not the implementation's.

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

**Rejected — leaving it to the implementation.** Its own `case default` is what used to report
this, and its `arg_pos` is in the implementation's numbering — which the wrapper now clears.
The check has to be the wrapper's for the position to mean anything. A direct Fortran caller is
the one this protects: the C layer already rejects an unknown mode *string* before Fortran is
entered.

### Validation ranges

Bounds are opt-in; **finiteness is not.** A real value that is NaN or infinite is a bug unless
the implementation says otherwise, so finiteness is the framework's default contract: every
real argument is checked, and an argument opts *out* per failure mode.

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
  `foo` and `foo_expert` number independently — they have to, since `foo` is the shorter list.

This requires one small extension to `tox_errors`: `validate_in_range_real` and
`validate_all_in_range_real` gain optional `allow_nan` / `allow_infinite` logicals, defaulting
`.false.` — i.e. today's behaviour, so nothing hand-written changes.

**Rejected — a `DM_FINITE` opt-in.** It puts the burden the wrong way round: the common case
(a value that must be finite) would need an annotation, and the rare, dangerous case (NaN
tolerated) would be silent. Finite-by-default with explicit opt-out states the tolerance where
the tolerance actually lives.

**Rejected — `DM_VALIDATE`, a hook for an implementation's own validation helper.** The
proposal was an argument-level directive letting the generated wrapper call a hand-written
validation routine inside the implementation module, passing *the wrapper's* position for that
argument — the one number the implementation cannot know, and the reason `clear_err_arg_pos`
has to zero what it returns. It had exactly two customers, both the same check:
`clock_hand_angle_between_vectors` needed its three axis indices to be mutually distinct, a
relation between an argument and itself that no per-argument bound can state.

That check is gone, because the argument that needed it was the wrong argument. The signed angle
now takes an orientation *vector* rather than three axis indices, so the rule is a bound like any
other. The general lesson, and the reason this stays rejected rather than deferred: **a procedure
that needs a validation nothing in the contract can express is a procedure whose signature is
wrong.** The generator's vocabulary is the API review. If a genuine new kind of constraint ever
turns up, it earns a *convention* — a directive with a name, a meaning, and every implementation
that qualifies using it — not a per-procedure escape hatch that hides one bad signature and stops
anyone from noticing the next.

### Sorting and work arrays (the allocating tier)

`foo_expert` exists only when the implementation needs work arrays, a pre-sorted permutation,
or a recommend-sized buffer — or when its prologue asks for an argument of its own; an
implementation with none of these generates just `foo`, and that `foo` is the plain
validate-and-call body. Where the pair does appear, `foo` drops from its signature every
argument the caller should not have to supply — `tmp_*` work arrays, permutations, and
recommend-computed sizes — allocates them locally, prepares them, and calls the
implementation. `foo_expert` keeps the full list, which is also where the emitter reads the
implementation's whole picture from when it writes `foo`'s body.

Permutations follow a convention: an argument named `<base>_perm` is, in `foo`, allocated,
seeded with `init_perm`, and heapsorted against `<base>` — and only while `<base>` is an
argument too, since a name that merely ends in `_perm` and orders something the implementation
never receives is the caller's own data. There is no override macro: an ordering the default
sort cannot give is expressed by a prologue that declares the permutation `intent(out)`
(below), and an argument the convention does not recognise — a bare `perm` — stays the
caller's own.

`foo` calls the **implementation directly**, not `foo_expert`. It has just built the
permutation and work arrays itself, so re-running the expert tier's validation — in particular
the O(n) `validate_all_in_range` on a permutation it knows is `[1..n]` — would be wasted work.

**Rejected — `foo` calling `foo_expert`.** It is simpler for the generator (validation lives in
one place) and is what the tox families did before they were generated, but it pays an O(n)
re-validation of data the allocating tier just constructed correctly. The f42 families call
their shared helper directly for exactly this reason; the generator follows them.

### Recommend / work-array sizing

A work-array size documented with `DM_OUTPUT_FROM(size_arg, recommend_proc, tox_X_impl, AUTO)`
already tells the binding languages to call `recommend_proc` themselves. The same annotation
now tells `foo` to call it in Fortran: size the buffer, `M_ALLOCATE` it, call the
implementation. One annotation, two consumers — the `foo_expert` binding still exposes the
size argument and calls the routine; `foo` internalises it.

**Rejected — a new `RECOMMENDED` mode alongside `AUTO`.** `AUTO` already carries everything
needed (the producer and how to supply its inputs, resolved by the existing output-from pass).
A second mode would be a second name for the same fact.

**Rejected — chaining producers, so one recommend routine feeds another.** The motivating case
is a workspace whose exact size depends on a count the implementation only discovers at run
time: `normalize_by_std_dev` fits LOESS through the `n_valid` genes that carry variance, not
all `n_genes`, so the exact workspace needs an `n_valid` producer whose output becomes the
`max_neighborhood_size` input of `tox_loess_required_workspace`. Supporting that means the
generator resolves a dependency graph over producers and emits them in topological order,
allocating each producer's own inputs before its call and the buffers that depend on it after —
and it still cannot help a producer whose input is itself a `tmp_` buffer that does not exist
yet.

The implementations do not need it. Every such size has a cheap upper bound that is already an
argument, and the recommend formulas are monotone in it, so the buffer is sized for the bound
and sliced at the call: `compute_family_scaling_impl` sizes for `n_families` and passes
`loess_x(1:n_valid)`, and `normalize_by_std_dev_impl` does the same with `n_genes`. What is
over-allocated is the fraction of points the fit drops; what is bought is that `foo` stays a
flat sequence of calls and allocations. Should a case appear where the bound is genuinely far
from the truth, the topological version is the answer — not before.

### The prologue is the allocating tier's sugar, and has no scope

`DM_PROLOGUE(PROCEDURE, MODULE)` names a routine `foo` runs after the work arrays are prepared
and before the implementation; it may write the outputs and report `handled`, and the
implementation is then skipped. It is what `foo` derives that `foo_expert` lets a caller pass
in -- the same relation the `<base>_perm` convention already has, where `foo` seeds and
heapsorts and a caller wanting another sort reaches for `foo_expert`. A prologue is that
convention's general form.

**Rejected -- a `SCOPE` of EXPERT / ALLOC / BOTH.** The macro carried one, and two of its three
values were mistakes:

- `BOTH` is inlinable. Both wrappers call the implementation, so anything that must run in
  both, before its body, *is* the first thing in that body -- `handled` becomes a local and the
  early return a plain statement. The prologue has no access the implementation lacks, so the
  directive bought nothing. `loess_degenerate_fit` was the only user, and it moved into both
  LOESS implementations; the third caller (`normalize_by_std_dev_inplace_helper`) had been
  performing the same check by hand, which is what a precondition-every-caller-must-remember
  looks like when it should have been the procedure's own contract.
- `EXPERT` contradicts what the expert tier is. FES: *the expert tier is the entry point for
  full control over what reaches the implementation -- a specific threshold, a specific
  initialised permutation -- while the plain one derives the threshold from a percentile and
  sorts with heapsort.* A prologue running in the expert tier would override exactly the
  control that tier exists to give.

So the scope is not derived, it is fixed, and the placement with it: always in `foo`, always
below the work arrays. An implementation with nothing to take over generates only `foo`, and
that `foo` is the validating body, which runs no prologue -- so a prologue on such an
implementation would never be called anywhere. That is an error
(`validate._check_prologue_runs_somewhere`), not a silence, and its note says where the work
belongs instead: at the top of the implementation, where every caller gets it.

**Nothing validated a prologue before this**: the emitter was the directive's only consumer, and
an emitter has no author's line to point at. Every way of getting one wrong therefore passed --
a name that resolved to nothing produced a wrapper with no prologue; a dummy that matched nothing
was dropped from the keyword call; a prologue without `handled` left the wrapper branching on an
undefined logical. All are errors now, in `ir/validate.py`.

**Rejected -- letting the presence of a prologue switch off the `<base>_perm` sorting.** More
explicit on the face of it, but at the wrong granularity: a prologue is per-implementation and
a permutation is per-argument, so an implementation with two of them would lose both defaults,
and adding a prologue to derive a threshold would silently change how an unrelated argument is
prepared. What went in instead is per-argument and derived from the two signatures: **a
permutation the prologue declares `intent(out)` is the prologue's**, and `foo` then allocates
it and stops. `intent(inout)` keeps the default and hands the prologue an order to refine. This
is what makes a non-default ordering expressible without the argument becoming `tmp_`, which
would take it out of the expert tier's signature -- the one place a caller is supposed to be
able to supply their own. `synthesize.sorted_permutations` is the single definition; the
emitter and `validate` both read it, so the sort that is emitted and the sort the read-first
check reasons about cannot diverge.

**What the prologue writes and the implementation reads is a wrapper local.** An argument that
is `intent(out)` on the prologue and `intent(in)` on the implementation involves nobody outside
the wrapper, so it joins the `tmp_`/perm/recommend drop set: gone from `foo`'s signature,
declared as a local there (allocated from the implementation's extents when it is an array),
and still a dummy of `foo_expert`, where the caller supplies it themselves. The two intents are
the whole annotation.

**Rejected -- requiring a `tmp_` name for it.** That was the first answer, and it forces the
implementation to declare `intent(inout)` for a value it only reads, because `_check_temporary`
rejects a `tmp_` argument that is `intent(in)`. A naming convention that makes an intent lie is
worse than no convention: the intents already carry the fact, and deriving from them keeps
`tmp_` meaning one thing (scratch) instead of two.

The one case that stays an error is an argument that would be dropped but also gives an extent of
something the caller still passes or receives -- the same brake `_sizes_a_kept_argument` puts on a
recommend-sized value, for the same reason. It has to stay a dummy, and Fortran will not let the
wrapper hand an `intent(in)` dummy to something that writes it.

`taken_over_arguments` therefore takes the prologue, and both the signature (built in
`synthesize`) and the body (written in the emitter) pass it, so the two cannot disagree about what
the caller passes.

**A prologue dummy the implementation does not have becomes an argument of `foo`.** What a
prologue derives *from* is the allocating tier's own vocabulary -- a threshold's `percentile` --
and
the implementation, which takes the threshold, has no use for it. So it joins that wrapper's
signature, after the implementation's own arguments and before `ierr`, and is validated there
like any other argument. `foo_expert` is untouched: it takes the derived value directly. This
is also why the two tiers now appear whenever their signatures would differ, rather than only
when something is taken over -- a prologue with an argument of its own is reason enough.

The cost is that a misspelling has nowhere left to be caught, since it reads as a new argument.
So a name one edit from an implementation argument is refused with "looks like a misspelling of
...": `n_gene` beside `n_genes` would otherwise leave the prologue and the implementation
working from different numbers, silently. One edit only -- two names differing by more than
that are two names.

**Rejected -- a rename table for prologue dummies.** A `DM_OUTPUT_FROM` producer needs one
because it is published and its parameter names cannot move. A prologue is internal to the
implementation module, so a mismatch is fixed by renaming its dummy, and the diagnostic says so.

### Validation is compiled out as a whole, or not at all

Both wrappers' checks sit behind `#ifndef NO_INPUT_VALIDATION`, so a build can drop them. Two
things decide the shape of that:

**`call set_ok(ierr)` is outside the guard.** It is not a check -- it is what leaves `ierr`
defined on the path where nothing goes wrong. Inside the guard, a validation-free build would
hand every caller an undefined `ierr` and an implementation's own runtime errors (§ the `ierr`
an implementation may keep) would be unreadable.

**Rejected -- an optional `validate` argument, or a per-procedure opt-out.** Whether to check
inputs is not a property of one call: a caller either trusts its data or does not, and mixing
the two within a build gives the worst of both -- the cost of the checks and none of the
confidence. A whole-build directive also costs nothing at run time, where an argument would
branch on every call. It follows `NO_C_BINDING` and `NO_R_BINDING`, which are the same kind of
decision (see [`c-layer.md`](c-layer.md)).

The C layer's null checks are *not* under this directive. They prevent a segfault rather than
reject a bad value, and a caller who has established that their inputs are good has established
nothing about a binding language passing a null pointer.

### No allocation in an implementation module

Nothing in an implementation module allocates: every buffer is a `tmp_` argument, so the
generated `foo` owns the memory and an expert caller can hand in buffers it already holds. The
rule covers the module's helpers as well as its implementations — one that allocates nothing
itself but calls a helper that does gives its expert caller nothing.

It is enforced on the **declaration**: a local declared `allocatable` is refused. The generator
never reads a body, and an `M_ALLOCATE` needs an allocatable to allocate into, so the declaration
is a complete proxy for this tree. (A `pointer` local that is `allocate`d would slip through. No
implementation does that, and the gap is cheaper than teaching the frontend to read statements.)
Ford hides a procedure's own variables by default, so the frontend sets `proc_internals`.

A `pointer` local is explicitly fine. Aliasing a buffer the implementation was handed allocates
nothing, and the `target` attribute it needs on a dummy never reaches the wrapper — the emitter
carries intent, type and dimension, not attributes. `normalization_pipeline_impl` relies on
both: it points its LOESS scratch at the still-unwritten columns of its own output buffer.

Two naming rules fall out of the same place. An implementation may not carry `M_EXPORT_C` — its
wrapper is the entry point, and exporting it beside the wrapper publishes an unvalidated twin
under a name (`foo_impl`) a binding caller cannot tell apart from `foo`. And it may not be named
for one of the wrappers it generates, `foo_expert_impl` or `foo_alloc_impl` alike (above).
Support routines in an implementation module — the recommend routines — are neither, and stay
exported.

### Mode-dependent implementations (opt-in split)

The house style is one procedure per mode — `detect_dosage_effect`,
`detect_subfunctionalization` — not a single procedure taking a runtime `mode`. An
implementation opts into this by giving its `mode`/`method` argument a mode table **with a
procedure-name column**: the generator then emits one entry point per mode value — with its
expert tier where there is one — names each from the table, drops the `mode` dummy and fixes
it, and treats each `DM_REQUIRED_IF_MODE(mode, …, Vᵢ)` argument as a mandatory dummy in `Vᵢ`'s
wrapper and absent from the others. Mode-independent optionals stay optional in every wrapper.

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
read, and a Fortran emitter renders them into `src/generated/`.

1. **Directives.** The `DM_*` macros are wired end-to-end (`src/macros.h`,
   `ir/directives.py`, `frontend/macros.py`) like the existing ones; the mode-table reader
   gains the optional procedure-name column; `tox_errors` gains the two `allow_*` logicals.
2. **Synthesis** (`synthesize.py`, run between parse and `analyse_project`). For each
   implementation — or once per mode value when it opts into mode-split — clone its arguments
   (carrying type, dimension, intent, optionality, docs, directives) into the validating
   wrapper (+ synthesised `ierr`) and, where anything is taken over or the prologue asks for an
   argument of its own, into the allocating one (minus the dropped arguments, plus the
   prologue's own); give both `Meta(category = the export category)` and the frontend's
   conventions, and put both in one synthesised `module tox_X`. The names are decided here and
   are final: the allocating wrapper is `foo` and its twin `foo_expert`, or, with nothing to
   take over, the lone wrapper is `foo` — so `abi.c_abi.stripped_name` has nothing left to
   translate and the C symbols follow by suffix alone. A side table records, per wrapper, only
   what its own signature cannot give: the implementation it calls, and the mode a per-mode
   wrapper fixes. The dropped arguments and the permutations to sort are not recorded — the
   emitter recomputes them through the same `taken_over_arguments` and `sorted_permutations`
   that shaped the signature, so the two cannot disagree — and the recommend calls are read
   post-analysis off the resolved output-from plan.
3. **Emission** (`emit/fortran_wrapper.py`, modelled on `emit/fortran_c.py`). Renders
   `module tox_X` with its `use … only:` imports (the implementation, the recommend routines,
   the `tox_errors` validators actually used, and `f42_sort` only when permutations are
   present), the declaration-ordering hoist that keeps an extent declared before the array it
   sizes, the validation block, and — for the allocating tier — the recommend calls,
   `M_ALLOCATE`s, `init_perm`/sort, the prologue and the implementation call. The plain entry
   point leads both the `public ::` block and the subroutine order, so a reader meets the tiers
   in the order they should be reached for. These are ordinary library sources: no
   `NO_C_BINDING` guard.
4. **Orchestration.** A `fortran` target and `_fortran_files` builder in `generate.py`, placing
   each module by the mirror `synthesize.generated_path_for` defines; `generated_dir` in
   `config.py` as the single path the rule is spelled in; `src/generated` added to Ford's
   `exclude_dir` so the generator never reads its own output back; `clean` removing the
   Fortran target's own files. `--check` diffs the emitted `.F90` against disk for free, like
   every other target.

**Rejected — a two-pass round-trip** (emit the wrappers, then re-run Ford over `src/generated`
to read them for binding generation). On the first run, or after any change to an
implementation, the on-disk wrappers are stale, so the bindings would be generated from the old
API. Synthesising in memory keeps the implementation the single source of truth and needs Ford
only once.

**Rejected — a `tox_X_variants` module name.** The generated module takes the clean `tox_X`
name precisely because it is the API callers should use; a decorated name would push the
recommended entry point behind an odd import.

---

## Rollout

The generator feature is independent of converting existing families, so it landed first and
was proven on pilots before the bulk migration:

1. Build the feature (macros, synthesis, emitter, target, wiring) with unit tests and an
   end-to-end test that compiles and runs an emitted wrapper.
2. Relocate the six blocker modules (green build at each step), reaching a generated-only
   `src/generated/tox`.
3. Pilot A — a validation-only family (e.g. `tox_shift_vectors`).
4. Pilot B — an allocating family with a permutation and a recommend routine (e.g.
   `tox_get_outliers`).
5. Pilot C — a mode-split family (`tox_paralog_analysis`).
6. Full rollout, ending with `tox_data_integration`, whose parent-interface + submodule
   structure dissolves: the submodules become plain implementation modules under
   `src/tox/data_integration/` gathered by a re-exporting parent, and the wrappers regenerate
   into the mirrored `src/generated/tox/data_integration/`.

Verification at each stage: the generator test suite plus `--check` idempotency; `./build.sh`
producing the `_c`/`_call` symbols for the wrappers in one library; and the full Python and R
suites after the relocations and pilots.
