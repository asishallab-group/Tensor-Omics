# Writing code the generator can wrap

TensorOmics generates its public API. There are two ways in, and the guide is split along them:

- **The implementation path** (Part I) — for a numeric procedure of the pipeline. You write one
  annotated implementation; the generator writes the entry point, its expert sibling where there
  is one, and the C, Python and R bindings, with the documentation, input validation and error
  handling derived from your implementation rather than restated four times by hand.
- **The export path** (Part II) — for IO, infrastructure, and everything that is not a numeric
  procedure of the pipeline. You write the whole procedure, validation included, and mark it
  `M_EXPORT_C`; the generator writes the three bindings.

Either way this guide is the author's side of the contract: **what to write so the generator can
do its half**, case by case, with a real example for each. Every snippet here is taken from the
current `src/` tree.

Related reading, once this is not enough:

| Document | What it answers |
|---|---|
| [`helper/codegen/README.md`](helper/codegen/README.md) | what the generator produces, how to run it, how it is built and tested |
| [`helper/codegen/design/impl-layer.md`](helper/codegen/design/impl-layer.md) | *why* the implementation layer works this way, and every alternative that was rejected |
| [`helper/codegen/design/c-layer.md`](helper/codegen/design/c-layer.md) | the Fortran↔C wrapper decisions |
| [`helper/codegen/design/language-layers.md`](helper/codegen/design/language-layers.md) | the Python and R decisions |
| [`src/macros.h`](src/macros.h) | the macros themselves, each with its contract in a comment |

---

## Contents

- [1. The one rule](#1-the-one-rule)
- [2. Where code lives, and which path you are on](#2-where-code-lives-and-which-path-you-are-on)
  - [The names, and who each is for](#the-names-and-who-each-is-for)

**Part I — the implementation path** (a numeric procedure of the TOX pipeline)

- [3. An implementation end to end](#3-an-implementation-end-to-end)
- [4. What every implementation must have](#4-what-every-implementation-must-have)
- [5. Case by case](#5-case-by-case)
  - [5.1 A value with a valid range](#51-a-value-with-a-valid-range)
  - [5.2 A real that may be NaN or infinite](#52-a-real-that-may-be-nan-or-infinite)
  - [5.3 Sizes, masks and counts the caller never passes](#53-sizes-masks-and-counts-the-caller-never-passes)
  - [5.4 A distance matrix](#54-a-distance-matrix)
  - [5.5 An optional argument and its default](#55-an-optional-argument-and-its-default)
  - [5.6 An output filled only partially](#56-an-output-filled-only-partially)
  - [5.7 Work arrays, and how the expert tier appears](#57-work-arrays-and-how-the-expert-tier-appears)
  - [5.8 A permutation](#58-a-permutation)
  - [5.9 A workspace sized by a recommend routine](#59-a-workspace-sized-by-a-recommend-routine)
  - [5.10 A mode argument](#510-a-mode-argument)
  - [5.11 One procedure per mode](#511-one-procedure-per-mode)
  - [5.12 An argument only one mode needs](#512-an-argument-only-one-mode-needs)
  - [5.13 Work the plain wrapper does beyond allocating](#513-work-the-plain-wrapper-does-beyond-allocating)
  - [5.14 A genuine runtime error](#514-a-genuine-runtime-error)
  - [5.15 A family too big for one file](#515-a-family-too-big-for-one-file)

**Part II — the export path** (IO, infrastructure, and everything else)

- [6. Exporting a hand-written procedure](#6-exporting-a-hand-written-procedure)
  - [6.1 When this is the right path](#61-when-this-is-the-right-path)
  - [6.2 What you write: `M_EXPORT_C`](#62-what-you-write-m_export_c)
  - [6.3 What you now do yourself: validate](#63-what-you-now-do-yourself-validate)
  - [6.4 What works exactly as in Part I](#64-what-works-exactly-as-in-part-i)
  - [6.5 An `_alloc` pair by hand — gone](#65-an-_alloc-pair-by-hand--gone)
  - [6.6 An array of any rank through one signature](#66-an-array-of-any-rank-through-one-signature)
  - [6.7 Never export an implementation](#67-never-export-an-implementation)

**Both paths**

- [7. Documentation that survives four languages](#7-documentation-that-survives-four-languages)
- [8. What the generator refuses](#8-what-the-generator-refuses)
- [9. The workflow](#9-the-workflow)
- [10. Before you commit](#10-before-you-commit)

---

## 1. The one rule

> **Write the implementation. Annotate it. Write nothing else.**

Everything the generator emits follows from your implementation's signature and its
documentation. So a fact belongs in exactly one place — the implementation — and if the
generator cannot derive something, the implementation gains an annotation, never the wrapper
hand-written code.

The corollary is the part worth internalising, because it will decide arguments you have one
day:

> **If the generator cannot express your procedure's contract, the procedure is wrong — not the generator.**

The annotation vocabulary is the API review. A signature that needs a rule none of the macros
can state is a signature that is asking its callers for the wrong thing. Redesign it; that is
almost always cheap, because only the tests depend on it. If a genuinely new *kind* of
constraint ever turns up, it earns a new convention — a macro with a name, a meaning, and every
implementation that qualifies using it — never a per-procedure escape hatch.

---

## 2. Where code lives, and which path you are on

```
src/
  macros.h            included by every source, by this path
  f42/                infrastructure, library-agnostic
    utils/              f42_utils_impl re-exports f42_{math,sort,random,vector,stats}_impl
    serde/              likewise, per element type
  tox/                the tox implementations -- the API's source of truth, hand-written
    data_integration/   a family split over several files (§5.15)
  data/               the hand-written data-set API (tox_data_*), incl. the zip archive
  generated/          NOTHING here is hand-written
    tox/                the wrappers, mirroring src/tox/ sub-directories and all
    f42/                likewise for src/f42/
    bindings/c/         the Fortran C wrappers
    bindings/r/         the R `.Call` shims
```

The rule inside `src/` is one line: **edit anything outside `src/generated/`.** That prefix is
what the generator deletes and rewrites on every run — a change you make there survives exactly
until the next build.

Three more trees are generated in full and live outside `src/` only because that is where their
languages expect them:

```
python/tensor_omics/    the Python package
r/tensor_omics/         the R package
snippets/               the VS Code snippets (except the hand-written toxdev_snippets.json)
```

All four are marked `linguist-generated` in [`.gitattributes`](.gitattributes), so a review
collapses them and sees the implementation that changed rather than the fan-out it produced. The
test suites under `python/test/` and `r/test/` are hand-written — those you do edit.

Two names are load-bearing, and together they are the generation trigger:

- the module is named **`<something>_impl`**
- the procedure inside it is named **`<name>_impl`**

That is all. There is no marker macro to forget, and **no fixed directory either**: the module
name is the whole trigger, so an implementation generates wherever under `src/` it is written.
The suffix is part of every call site, which matters, because calling an implementation directly
means calling something with no input validation. A procedure in an implementation module that
does *not* end in `_impl` is not one — the recommend routines of §5.9 and private helpers live
there untouched.

Because the suffix alone triggers, **`_impl` is reserved across the whole of `src/`**: any module
so named acquires generated wrappers and is held to §4's rules, whether it sits under `src/tox/`
or anywhere else. The file must be named for the module, too (`tox_shift_vectors_impl` in
`tox_shift_vectors_impl.F90`) — see §4 for why.

Where the wrappers are written is a mechanical mirror of where the implementation sits, with no
knowledge of layers in it:

```
src/<rest>/<module>_impl.F90   ->   src/generated/<rest>/<module>.F90
```

so `src/tox/tox_loess_impl.F90` generates `src/generated/tox/tox_loess.F90`, and
`src/tox/data_integration/tox_data_integration_jsd_impl.F90` generates
`src/generated/tox/data_integration/tox_data_integration_jsd.F90`. The generated module drops the
suffix and keeps everything else: `compute_shift_vector_field_impl` in `tox_shift_vectors_impl`
becomes `compute_shift_vector_field` in `tox_shift_vectors`. **That clean name is the public API**
— in Fortran, C, Python and R alike.

### The names, and who each is for

One implementation puts two or three names into the tree. They are not variants of one function;
each is addressed at a different caller.

| Name | Who it is for | What it does |
|---|---|---|
| `foo_impl` | **nobody** — it is never exported | the implementation. No validation, no `ierr`. Reached only through a wrapper |
| `foo` | almost everyone | validates, calls the recommend routines, allocates the work arrays, seeds and heapsorts the permutations, runs the prologue, then calls `foo_impl` |
| `foo_expert` | a caller who wants control over what reaches the implementation | validates and calls `foo_impl` with what you supply. Allocates nothing, prepares nothing |

A wrapper is **`pure` exactly when everything it calls is** — the implementation, the
`DM_OUTPUT_FROM` producers and the prologue. `tox_errors`' validators and `f42_sort`'s
`init_perm` / `sort_array_heapsort` already are, and `allocate(..., stat=)` is permitted in a
pure procedure, so the allocating tier is no less eligible than the expert one. You declare
nothing: write `pure` on the implementation and the wrapper follows. That is what lets a
caller reach a validated entry point from `do concurrent` instead of reaching past it.

**The plain name always goes to the entry point a caller should reach for first**, and it is the
same name in Fortran, C, Python and R: nothing is renamed on the way out, and the Fortran source
and the published API can be grepped for with one string.

**`foo_expert` exists only where there is something for `foo` to take over** — work arrays, a
`<base>_perm` permutation, a `DM_OUTPUT_FROM(..., AUTO)` value, an argument the prologue asks for
of its own. Where there is nothing, one wrapper is generated and it is called `foo`, because
there is no second tier for a suffix to distinguish it from. That is 39 of the 87 generated
procedures today, so do not read — or write — as though every implementation produced a pair.

The pair, where there is one, is a **sugar/control** split, not a fast/slow one. `foo` derives
what the expert tier lets you pass: the heapsorted permutation of §5.8, the threshold or
short-circuit a prologue computes (§5.13), the workspace sizes of §5.9. Want a different sort
order, your own reused buffers, or the computation without a degenerate-input policy? That is
what the expert tier is.

**The expert tier is not published to Python and R unless it offers something.** Where `foo` only
validates and allocates, those languages allocate the work arrays for *both* tiers anyway, so
`foo_expert` would be the same call under a name claiming otherwise — the generator emits only
`foo` there. Fortran and C always get both: there the expert tier really does hand the buffers
over. Today two generated procedures keep a Python/R expert tier — `pool_means` and
`determine_shared_residual_range` — each because it seeds and sorts a permutation you may supply
yourself. Where both are published, **each docstring says what the other does** -- which
permutation the plain one seeds and sorts, which prologue it runs -- so a reader does not have to
work out why there are two.

### Which path is yours

Not everything is an implementation, and a third of the public API is not — 39 hand-written
exports against 87 generated wrappers today. File and archive IO, the f42 trees and statistics,
the whole serde family, and the sizing routines the implementations themselves call are all
**hand-written and exported**: the generator wraps them to C, Python and R, but does not write
their Fortran.

| Your procedure | Path | Where | Marker |
|---|---|---|---|
| a numeric procedure of the TOX pipeline | **Part I — implementation** | `src/tox/` | the `_impl` name |
| the data-set API: archive, CSV/TSV, validation, accessors | **Part II — export** | `src/data/` | `M_EXPORT_C` |
| library-agnostic infrastructure | **Part II — export** | `src/f42/` | `M_EXPORT_C` |
| a sizing / utility routine an implementation or a caller needs | **Part II — export** | the implementation module, `public` | `M_EXPORT_C` |

`src/tox/` in the first row is convention, not mechanism: the trigger is the name, so an `_impl`
module under `src/f42/` would generate exactly the same way, into `src/generated/f42/`. Nothing
there has been converted yet (§6.5 says why), but the mechanism no longer knows the difference.

The two paths share everything about the *bindings* — naming conventions, documentation, type
rules, the refusals in §8. They differ in exactly one thing: **who writes the validation.** On
the implementation path the generator does. On the export path you do (§6.3).

---

# Part I — the implementation path

---

## 3. An implementation end to end

The whole of `src/tox/tox_shift_vectors_impl.F90`, minus the body:

```fortran
#include <src/macros.h>

!> Implementation of the shift vector field for all genes.
!|
!| Hand-written implementation only. The generator turns this into the entry point
!| [[tox_shift_vectors(module):compute_shift_vector_field]] in module `tox_shift_vectors`.
module tox_shift_vectors_impl
    use, intrinsic :: iso_fortran_env, only: real64, int32
    M_IMPLICIT_NONE
    private
    public :: compute_shift_vector_field_impl
contains

    !> summary: Compute the shift vector field for all genes.
    !| AUTHOR_ALEXANDER_SCHWARZPAUL
    !| Computes the shift vectors by subtracting the corresponding family centroid from the expression vector.
    pure subroutine compute_shift_vector_field_impl(n_tissues, n_genes, n_families, expression_vectors, &
                                                    family_centroids, gene_to_fam, shift_vectors)
        integer(int32), intent(in) :: n_tissues
            !! Expression vector dimension
        integer(int32), intent(in) :: n_genes
            !! Total number of genes
        integer(int32), intent(in) :: n_families
            !! Total number of families
        real(real64), intent(in) :: expression_vectors(n_tissues, n_genes)
            !! Gene expression matrix
        real(real64), intent(in) :: family_centroids(n_tissues, n_families)
            !! Family centroid matrix
        integer(int32), intent(in) :: gene_to_fam(n_genes)
            !! M_GENE_TO_FAM_DOC(expression_vectors)
            !! DM_MIN(1_int32)
            !! DM_MAX(n_families)
            !! DM_SENTINEL(M_GENE_TO_FAM_SENTINEL)
        real(real64), intent(out) :: shift_vectors(n_tissues, 2, n_genes)
            !! Output, real matrix array. For each gene it holds two vectors: ...

        ! ... the implementation, and nothing else
    end subroutine compute_shift_vector_field_impl
end module tox_shift_vectors_impl
```

Note what is *absent*: no `ierr`, no validation, no `set_ok`, no `M_EXPORT_C`. Note the three
annotation lines on `gene_to_fam` — the only thing in the file that is not implementation.

**What the generator writes** into `src/generated/tox/tox_shift_vectors.F90`:

```fortran
    subroutine compute_shift_vector_field(&
            n_tissues, n_genes, n_families, expression_vectors, family_centroids, gene_to_fam, shift_vectors, ierr)
        ...                                              ! your declarations and your docs
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_tissues, ierr, arg_pos=1_int32)
        call validate_dimension_size(n_genes, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_families, ierr, arg_pos=3_int32)
        call validate_all_in_range_real(expression_vectors, n_tissues * n_genes, ierr, arg_pos=4_int32)
        call validate_all_in_range_real(family_centroids, n_tissues * n_families, ierr, arg_pos=5_int32)
        call validate_all_in_range_int(gene_to_fam, n_genes, ierr, arg_pos=6_int32, &
                                       min=1_int32, max=n_families, sentinel=0_int32)
        if (is_err(ierr)) return
#endif

        call compute_shift_vector_field_impl(...)
    end subroutine compute_shift_vector_field
```

There is nothing to take over here — no work array, no permutation, no recommend-sized buffer —
so this one wrapper is the whole of it, and it takes the plain name. Five of those six checks
you never asked for: the extents are validated because they *are* extents, and both real matrices
are checked for NaN and infinity because that is the framework's default contract (§5.2). The
`#ifndef` is what `--directive=NO_INPUT_VALIDATION` compiles out (§9).

**What a caller gets**, from the same source, without another line written:

```python
compute_shift_vector_field(expression_vectors, family_centroids, gene_to_fam)   # extents derived
```
```r
compute_shift_vector_field(expression_vectors, family_centroids, gene_to_fam)
```

with your `summary:` as the docstring, your argument docs as the parameter docs, and an
`ierr` that raises a typed error naming the offending argument.

---

## 4. What every implementation must have

| Requirement | Why |
|---|---|
| `#include <src/macros.h>` at the top of the file | a `DM_` that never expands is an error, as is a misspelt `M_`/`CM_`/`DM_` in any doc comment |
| `M_IMPLICIT_NONE` in the module | `implicit none (type, external)` — a typo'd call is a compile error, not a link-time surprise |
| `private` + an explicit `public ::` list | the implementation is API by way of its wrapper, not by accident |
| **The module's own name as the file name** (`tox_x_impl` in `tox_x_impl.F90`) | two independent things read that name and must agree: generation triggers on the *module* name, while the cleaner and the Ford exclusion scan *file* names without parsing. Let them diverge and the generator writes a wrapper nothing cleans and nothing excludes — so the next run parses its own output and defines the module twice |
| Explicit `intent` on every dummy | decides constness in C and R, and input/output/in-out everywhere |
| Kinded numeric types (`real(real64)`, `integer(int32)`) | a default kind has no defensible C mapping |
| `!> summary: ...` on the procedure | becomes the docstring in Python and R |
| **A `!>` block on the module, written as API prose** | it is carried verbatim onto the generated module, so it is the Python module docstring and the Ford page for the published API. Say what the family is for and name its entry points; do *not* describe the implementation, and do not explain what the generator will do with it (§7) |
| An author tag (`AUTHOR_*` from [`authors.h`](authors.h)) | attribution, rendered into the Ford docs |
| A `!!` doc on **every** argument | inherited by the C wrapper and by both language layers |
| **No `ierr` for validation** | validation is the wrapper's job; see §5.14 for the one case where an implementation keeps an `ierr` |
| **No `M_EXPORT_C` on the implementation** | the generated wrapper is what the bindings call; exporting the implementation beside it publishes an unvalidated twin under a name a caller cannot tell apart. Support routines in the same module — the recommend routines of §5.9 — *are* exported: they have no wrapper |
| **No `_alloc` or `_expert` in the implementation's name** | both are wrapper suffixes, and neither is yours to choose. `foo_expert_impl` would generate a second procedure called `foo_expert`; `foo_alloc_impl` would generate `foo_alloc`, which every author here reads as the allocating tier while being an ordinary second procedure beside the generated `foo` that is one (§6.5). Name the implementation for what it computes; whether a second tier appears at all is decided by its `tmp_` arguments (§5.7) |
| **No allocation, anywhere in the module** | every buffer is a `tmp_` argument, so the generated `foo` owns the memory and an expert caller can hand in buffers it already has. The rule covers the module's helpers too: an implementation that allocates nothing itself but calls a helper that does is no better off. Enforced on the declaration — a local *or a dummy* declared `allocatable` is refused. A `pointer` local is fine: aliasing a buffer you were handed allocates nothing |
| **Only implementations and infrastructure may be `use`d** | another `_impl` module, or one of `impl_import_whitelist` — the intrinsic modules, `tox_errors`, `tox_conversions`, `f42_config`, `f42_safeguard`. That bound is what makes the rule above hold *across* modules: the check reads declarations, so only the import list can see a helper elsewhere that allocates. It also fixes the direction — an implementation cannot reach a generated wrapper, which would invert the layering and, within one family, be a module cycle. A `use` inside a procedure counts the same as one in the module header |

Extents are recognised on sight: `n_tissues` in `expression_vectors(n_tissues, n_genes)` is an
extent, is validated as one, and is never asked of a Python or R caller — they pass the array,
and the binding reads its shape.

---

## 5. Case by case

Each case is: **when** it applies, **what you write**, **what you get**.

### 5.1 A value with a valid range

**When** an argument is only meaningful inside bounds — an index, a count, a probability, an
angle.

**Write** `DM_MIN` / `DM_MAX` in the argument's doc. The expression is Fortran source and may
name other arguments or module constants. Wrap it in `above(...)` / `below(...)` for an
exclusive bound. `DM_SENTINEL` adds one value that is accepted regardless — a marker that is not
a datum.

```fortran
integer(int32), intent(in) :: gene_to_fam(n_genes)
    !! M_GENE_TO_FAM_DOC(expression_vectors)
    !! DM_MIN(1_int32)
    !! DM_MAX(n_families)
    !! DM_SENTINEL(M_GENE_TO_FAM_SENTINEL)

real(real64), intent(in), optional :: gain_gamma
    !! positive magnitude gain for dosage effect
    !! DM_MIN(above(0.0_real64))
```

**You get** a `validate_all_in_range_int` / `validate_in_range_real` call in the wrapper,
carrying exactly the keywords you documented, blaming this argument's position. And the bound
appears in the Python and R documentation, because it was written as prose in the first place.

> An **integer** with no bound is not validated at all — integers cannot be NaN, and there is
> nothing else to check. If an integer has a meaningful range, say so; nobody else will.

### 5.2 A real that may be NaN or infinite

**When** an argument legitimately carries non-finite values — a masked-out mean, a distance to
something that is not there.

**Write** the opt-out. Finiteness is the *default*: every real argument is checked unless it
says otherwise, and each failure mode opts out separately.

```fortran
real(real64), intent(in) :: distances(n_genes)
    !! Array of Euclidean distances for each gene
    !! DM_ALLOW_NAN
    !! DM_ALLOW_INFINITE
```

**You get** `allow_nan=` / `allow_infinite=` passed to the validator — and a tolerance stated
where the tolerance actually lives, which is the point. The burden is deliberately this way
round: the common case (must be finite) is free, and the rare, dangerous case is visible.

> Counter-example worth knowing: netlib LOESS cannot fit a non-finite sample. It dies inside the
> decomposition with `svddc failed in l2fit` rather than reporting. So loess arguments must
> **not** carry these macros — the default rejection is what keeps that from happening.

### 5.3 Sizes, masks and counts the caller never passes

**When** an argument is an extent, or the count of `.true.` values in a mask.

**Write** the naming convention and nothing else:

| Name | Meaning |
|---|---|
| an extent used in a declaration (`vec(n_dims)`) | validated with `validate_dimension_size`, derived by the binding from the array |
| `<arg>_mask`, `<arg>_selection_mask` | a selection mask |
| `n_selected_<arg>` | the count of `.true.` in that mask |

```fortran
logical, dimension(n_axes), intent(in) :: axes_selection_mask
    !! Logical array (n_axes), .TRUE. for axes to include in calculation
integer(int32), intent(in) :: n_selected_axes
    !! Number of selected axes (count of .TRUE. in axes_selection_mask)
    !! DM_MIN(1_int32)
```

**You get** a cross-check Fortran itself cannot make, plus a floor if you asked for one:

```fortran
call validate_in_range_int(n_selected_axes, ierr, arg_pos=7_int32, min=1_int32)
if (count(axes_selection_mask, kind=int32) /= n_selected_axes) &
    call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=7_int32)
```

and Python/R callers who pass the mask alone — `n_selected_axes` is computed from it at the call.

### 5.4 A distance matrix

**When** a square real matrix holds pairwise distances.

**Write** the name: the argument is called `distances` or `distance_matrix`, or ends in
`_distances` / `_distance_matrix`.

```fortran
real(real64), dimension(n_points, n_points), intent(inout) :: distances
```

**You get** `validate_distance_matrix` instead of the plain finiteness sweep — symmetry,
non-negativity and a zero diagonal. (A *vector* named `distances` is not a distance matrix and
gets the ordinary treatment.)

### 5.5 An optional argument and its default

**When** an argument may be omitted.

**Write** `intent(in), optional` plus `DM_DEFAULT`, whose value must be a **constant
expression** — the generator evaluates it at generation time so the binding languages can pass
it. Use `M_DEFAULT_VAL` in the body to resolve it into a local.

```fortran
real(real64), intent(in), optional :: span
    !! Span parameter for LOESS smoothing
    !! DM_DEFAULT(CM_FAMILY_SPAN_DEFAULT)
```
```fortran
real(real64) :: actual_span
M_DEFAULT_VAL(span, actual_span, CM_FAMILY_SPAN_DEFAULT)
```

**You get** an optional that Python and R expose *with the real default in the signature and in
the docs*, evaluated: `CM_FAMILY_SPAN_DEFAULT` reaches a Python caller as `0.7`.

> **Never pass an optional dummy straight into a mandatory dummy of something you call.** That
> is a null dereference when it is absent, and it is exactly how `loess_fit_robust_impl`
> segfaulted. Resolve it with `M_DEFAULT_VAL` first and pass the local.

An optional **without** a default is nullable: absent in Fortran, `None`/`NULL` from the
language layers. An optional **with** one is always passed, which keeps the wrapper flat.

### 5.6 An output filled only partially

**When** a result array is sized for the worst case but only its leading elements are filled —
a selection whose size is not known until the work is done.

**Write** `DM_RESULT_SIZE_IS(<argument>)` on the result, naming the `intent(out)` scalar integer
your implementation sets to the number of elements it actually filled:

```fortran
real(real64), dimension(n_genes), intent(out) :: results
    !! the results.
    !! DM_RESULT_SIZE_IS(n_results)
integer(int32), intent(out) :: n_results
    !! how many leading elements of `results` were filled
```

**You get** a Python and R caller who never sees the padding, and never sees the count either —
the array they get back is already that long:

```python
results[..., :n_results.value]          # Python
```
```r
utils::head(.result$results, .result$n_results)   # R
```

The counting argument is dropped from the return, because the returned array already carries the
answer; returning both would invite the caller to slice a second time.

This is a **binding-level** trim: the generated Fortran wrapper passes both arguments through
untouched, so a Fortran caller still receives the full buffer plus the count, which is what a
Fortran caller wants. Size the array for the worst case as usual — the language layers hide it.

> The example is `fx_masked` in [`helper/codegen/tests/fixtures/src/fx_edges.F90`](helper/codegen/tests/fixtures/src/fx_edges.F90).
> No implementation in `src/` uses this yet — f42's hand-written `compute_edf` does, on the
> export path (§6.4) — so the fixture is the worked case for this path.

### 5.7 Work arrays, and how the expert tier appears

**When** your implementation needs scratch space it does not want to allocate itself (because a
hot loop should not allocate, and because a Fortran caller may want to reuse a buffer).

**Write** the argument with a `tmp_` prefix and an `intent(out)` or `intent(inout)`:

```fortran
integer(int32), intent(out) :: tmp_stack_left(n_genes)
    !! Stack array for left indices during sorting
```

**You get** a *second* wrapper, and the two split the name. `foo` drops the work arrays from its
signature, declares them locally as allocatables, `M_ALLOCATE`s them and calls the
implementation:

```fortran
integer(int32), dimension(:), allocatable :: tmp_stack_left
...
M_ALLOCATE(tmp_stack_left(n_genes))
```

`foo_expert` keeps them — that is the expert entry point, for a caller who manages buffers. The
language layers need no translation for this: `foo` is `foo` and `foo_expert` is `foo_expert`
everywhere, subject only to §2's rule that Python and R drop an expert tier with nothing to
offer. An implementation with no work arrays, permutations or recommend-sized buffers — and no
prologue asking for arguments of its own (§5.13) — generates only `foo`.

> A `tmp_` argument that is `intent(in)` is an error. A work array is an output or an in-out;
> nothing else makes sense.

This is the *only* way an implementation module gets scratch space: an `allocatable` local
anywhere in it is refused (§4). A buffer whose size is not an expression over the other arguments
is still a `tmp_` argument — `DM_OUTPUT_FROM(..., AUTO)` (§5.9) names the routine that sizes it.
Where even that routine cannot be called ahead of time because the size depends on something only
the implementation discovers, size the buffer at the **upper bound** and slice it:
`normalize_by_std_dev_impl` takes its LOESS workspace for all `n_genes` and hands the fit
`tmp_loess_x(1:n_valid)`, because dropping the zero-variance genes only ever makes the fit
smaller.

### 5.8 A permutation

**When** the implementation needs an index permutation sorted against one of its arrays.

**Write** the argument as `<base>_perm`, where `<base>` is **an argument too** — the wrapper has
to have the array in hand to sort against it:

```fortran
real(real64), intent(in) :: pooled_means(pool_size)
    !! Pooled means
integer(int32), intent(in) :: pooled_means_perm(pool_size)
    !! Sorting permutation for `pooled_means`
```

**You get**, in `foo` only:

```fortran
M_ALLOCATE(pooled_means_perm(pool_size))
call init_perm(pooled_means_perm)
call sort_array_heapsort(pooled_means, pooled_means_perm)
```

`foo_expert` still takes it, because an expert caller may already have it sorted. A `<base>` that
is not an argument means the wrapper cannot build the permutation, so the argument stays an
ordinary one the caller supplies — and a `tmp_`-prefixed permutation is the implementation's own
scratch, allocated and left alone (`tmp_perm` in `compute_family_scaling_impl` is seeded inside
the implementation).

**For a different ordering** — descending, stable, by magnitude — let the prologue (§5.13) build
it: a permutation the prologue declares `intent(out)` is the prologue's, so `foo` allocates it and
stops. The argument keeps its `<base>_perm` name, so the expert tier still accepts an order from
the caller; only what `foo` builds by default changes. Declare it `intent(inout)` instead and the
wrapper still seeds and sorts first, handing the prologue an order to refine.

`foo` calls the **implementation directly**, not `foo_expert` — it just built that permutation
itself, so re-running an O(n) validation over `[1..n]` would be wasted work. A bare `perm` with no
base name is not a permutation by this convention; it stays an ordinary argument.

### 5.9 A workspace sized by a recommend routine

**When** a buffer's size can only be computed by calling something — netlib's LOESS workspaces
are the standing example.

**Write** `DM_OUTPUT_FROM(<size_arg>, <producer>, <its module>, AUTO)` on the size argument.
Producer inputs are matched to your arguments **by name**; where the two spell a quantity
differently, or where the producer wants a constant, a table says so:

```fortran
integer(int32), intent(in) :: int_workspace_size
    !! Length of integer workspace.
    !! DM_OUTPUT_FROM(int_workspace_size, tox_loess_required_workspace, tox_loess_impl, AUTO)
    !!
    !! | Producer input        | Supplied by |
    !! |-----------------------|-------------|
    !! | n_dim                 | 1_int32     |
    !! | max_neighborhood_size | n_families  |
    !! | save_factorization    | .false.     |
integer(int32), intent(out) :: tmp_int_workspace(int_workspace_size)
    !! Integer workspace array
```

**You get** the producer called for you in `foo`, before the allocation that needs it:

```fortran
call tox_loess_required_workspace(n_dim = 1_int32, max_neighborhood_size = n_families, &
                                  int_workspace_size = int_workspace_size, &
                                  real_workspace_size = real_workspace_size, &
                                  save_factorization = .false.)
M_ALLOCATE(tmp_int_workspace(int_workspace_size))
```

— and, in the `foo_expert` binding, the size argument still exposed with its documentation
pointing at the routine that computes it. One annotation, two consumers.

The producer must be exported (`M_EXPORT_C`), so that a wrapper exists to call; recommend
routines therefore live in the implementation module, public, tagged — they are not
implementations themselves and the `_impl` rules of §4 leave them alone. Use `JUST_INFO` instead
of `AUTO` when the caller genuinely has to make the call themselves — then the docs say where to
get the value, and nothing is called for them.

### 5.10 A mode argument

**When** one procedure does a related thing several ways, chosen at run time.

**Write** an integer argument named `mode`, `method`, or ending in `_mode` / `_method`, compared
against `MODE_*` / `METHOD_*` parameters, and document the accepted values in a table. The
argument's own name decides the prefix, and the table's first column header must agree with it:

```fortran
integer(int32), intent(in), optional :: mode
    !! Mode for LOESS fitting
    !! DM_DEFAULT(CM_FAMILY_MODE_DEFAULT)
    !!
    !! | Mode | Value |
    !! |------|-------|
    !! | Plain LOESS fitting  | [[tox_loess_impl(module):MODE_PLAIN(variable)]]  |
    !! | Robust LOESS fitting | [[tox_loess_impl(module):MODE_ROBUST(variable)]] |
```

**You get** two things. Python and R callers pass a **string** — the parameter name without its
prefix, lower-cased, so `MODE_PLAIN` is `"plain"` — which the C wrapper maps back to the
integer, rejecting an unknown one before Fortran is entered. And the wrapper checks membership
against exactly the values the table names:

```fortran
if (present(mode)) then
    if (mode /= MODE_PLAIN .and. mode /= MODE_ROBUST) &
        call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=11_int32)
end if
```

Do **not** bound a mode with `DM_MIN`/`DM_MAX`: the accepted set need not be contiguous, and a
hand-written range goes stale the moment a mode is added — silently admitting the new value
everywhere it was not updated. The table is the source of truth, and the generator already reads
it.

### 5.11 One procedure per mode

**When** the modes are really different procedures wearing one signature — the house style is
`detect_dosage_effect` and `detect_subfunctionalization`, not `detect(pattern_mode)`.

**Write** the same mode table with a **third column naming the procedure per value**. That
column *is* the opt-in; nothing else changes.

```fortran
integer(int32), intent(in) :: pattern_mode
    !! used pattern for detection
    !!
    !! |         Mode         |                            Value                             |         Procedure           |
    !! |----------------------|--------------------------------------------------------------|-----------------------------|
    !! |    Dosage Effect     | [[tox_paralog_analysis_impl(module):MODE_DOSAGE_PATTERN(variable)]]  | detect_dosage_effect        |
    !! | Subfunctionalization | [[tox_paralog_analysis_impl(module):MODE_SUBFUNC_PATTERN(variable)]] | detect_subfunctionalization |
```

**You get** one entry point per mode value — named from the column, with the `mode` dummy dropped
and fixed internally — plus its `_expert` sibling wherever that mode's wrapper has something to
take over (§5.7): `detect_dosage_effect`, `detect_dosage_effect_expert`,
`detect_subfunctionalization`, `detect_subfunctionalization_expert`, and (nothing to take over, so
no second tier) `filter_paralogs_by_pattern_dosage_effect` and
`filter_paralogs_by_pattern_subfunctionalization`.

Each validates its own arguments at its own positions, so no position remapping is needed
anywhere — which is exactly what the hand-written per-mode wrappers used to spend a dozen
`map_err_arg_pos` calls apiece on.

Without the column, nothing changes: one procedure, mode chosen at run time. Both shapes stay
available, because some mode arguments genuinely are runtime choices.

### 5.12 An argument only one mode needs

**When** an argument is meaningless outside one mode.

**Write** `DM_REQUIRED_IF_MODE(<mode arg>, <module>, <mode parameter>)` on an optional argument:

```fortran
real(real64), intent(in), optional :: gain_gamma
    !! positive magnitude gain for dosage effect
    !! DM_REQUIRED_IF_MODE(pattern_mode, tox_paralog_analysis_impl, MODE_DOSAGE_PATTERN)
    !! DM_DEFAULT(0.1_real64)
    !! DM_MIN(above(0.0_real64))
```

**You get** behaviour that depends on whether the mode splits (§5.11):

- **Split**: the argument is a **mandatory** dummy in that mode's wrapper and **absent** from
  every other mode's. Mode-independent optionals stay optional in all of them. A `DM_DEFAULT`
  alongside is fine here — it applies within the mode that has the argument.
- **Not split**: it stays optional, and the documentation states in which mode it is required.
  A `DM_DEFAULT` alongside is an **error** here: the binding always passes a defaulted argument
  on, so "required in that mode" would say nothing.

### 5.13 Work the plain wrapper does beyond allocating

**When** `foo` should *derive* something the expert tier lets a caller pass in — a threshold taken
from a percentile of the data, a permutation built and sorted a particular way — or should decide
the input is too degenerate to compute on and answer it directly.

This is the same idea as the `<base>_perm` convention of §5.8, which is a prologue hard-coded for
one case: `foo` seeds and heapsorts, and a caller who wants a different sort reaches for
`foo_expert`. A prologue is the general form.

**Write** `DM_PROLOGUE(<procedure>, <module>)` in the implementation's own doc block. The prologue
takes the implementation's arguments by name — the `tmp_` work arrays included — plus `handled`
and `ierr`.

**A dummy the implementation does not have becomes an argument of `foo`.** What the prologue
derives *from* is usually no business of the implementation's: a threshold comes from a
`percentile`, and the implementation takes the threshold. So `percentile` joins the plain
wrapper's signature — that one alone, since only it runs the prologue — and is validated there
like any other argument, its range macros included. `foo_expert` takes the derived value directly
and never sees it:

```
foo         (values, n, result, percentile, ierr)               ← you supply the percentile
foo_expert  (values, n, threshold, tmp_scratch, result, ierr)   ← you supply the threshold
```

The prologue's own arguments come after the implementation's and before `ierr`.

A name that is *nearly* one the implementation has is refused as a misspelling — `n_gene` beside
`n_genes` would otherwise become a new argument silently, leaving the prologue and the
implementation working from different numbers.

```fortran
!> summary: Flag directional outliers
!| AUTHOR_A_N_AUTHOR
!| DM_PROLOGUE(prepare_ranking, tox_outliers_impl)
```
```fortran
pure subroutine prepare_ranking(z_scores, n_genes, tmp_valid_perm, threshold, handled, ierr)
    logical, intent(out) :: handled
        !! `.true.` when the data was degenerate and the outputs already hold the answer
```

**What the prologue writes and the implementation only reads never reaches a caller.** An argument
that is `intent(out)` on the prologue and `intent(in)` on the implementation is dropped from
`foo`'s signature and declared as a local there — a scalar as a plain one, an array allocated from
the implementation's own extents. Both intents already say this, so no `tmp_` name is needed, and
the implementation keeps the honest `intent(in)` for something it only reads. `foo_expert` still
takes it, because `foo_expert` has no prologue: the expert tier is where a caller supplies the
value themselves.

The other two pairings mean different things and are left alone. An implementation `intent(out)`
makes the two *alternative producers* of one output — the prologue's value is what the caller gets
on the `handled` path, the implementation overwrites it otherwise. An implementation
`intent(inout)` means the caller supplies a value the prologue then refines.

**You get** the prologue called in `foo` — after the work arrays are prepared, before the
implementation — and the implementation skipped when it says so:

```fortran
M_ALLOCATE(tmp_valid_perm(n_genes))
call prepare_ranking(z_scores = z_scores, n_genes = n_genes, &
                     tmp_valid_perm = tmp_valid_perm, threshold = threshold, &
                     handled = handled, ierr = ierr)
if (is_err(ierr)) return
if (handled) return
```

**There is no scope**, and that is the point. `foo_expert` is the tier that gives a caller full
control over what reaches the implementation; a prologue running there would override exactly the
control that tier exists to provide. And work that *both* tiers need is work every caller of the
implementation needs — both call it — so it belongs at the top of the implementation, where it is
three lines and no directive. That is where LOESS's degenerate-input check lives.

**What is refused**, all of it silent before:

- a `DM_PROLOGUE` naming a procedure that does not exist
- a prologue with no `handled`, or one that is not a scalar `logical, intent(out)` — the wrapper
  returns early on it regardless, so without it that branch reads an undefined value
- a dummy that is one edit from an implementation argument — a misspelling, not a new argument
- a dummy that is the mode argument, or one scoped to a mode, of an implementation that splits per
  mode (§5.11) — the wrappers for the other modes do not have it
- a prologue on an implementation with nothing to take over *and* no arguments of its own, which
  generates a single wrapper with no prologue in it, so the call is emitted nowhere
- a value the prologue writes and the implementation reads that *also* sizes something the caller
  still passes or receives — it cannot become a local then, and Fortran will not hand an
  `intent(in)` dummy to something that writes it
- a value the prologue produces that the allocations, the permutation sorts or the recommend
  calls *above* it read — the name resolves either way, so it would compile and compute rubbish.
  (A permutation the prologue itself builds is not sorted above it, so rewriting that one's
  `<base>` is fine.)

There is deliberately **no rename table**: unlike a `DM_OUTPUT_FROM` producer, a prologue is
internal to the implementation module, so the fix is to rename its dummy.

Reach for this only when the work is genuinely the plain wrapper's. A prologue that merely
*validates* is a smell — that is what the range macros are for, and a pre-validation pass that
duplicates them is a bug rather than a feature.

### 5.14 A genuine runtime error

**When** the implementation can fail for a reason no input check could have foreseen: an external
library reporting failure, a division by zero, a configuration of otherwise-valid inputs that
has no answer.

**Write** an `ierr` on the implementation after all, set it with `set_err_once`, and say so in the
documentation:

```fortran
logical :: undefined_sign
call set_ok(ierr)
call clock_hand_angle_between_vectors_helper(v1, v2, n_dims, orientation_reference, &
                                             signed_angle, undefined_sign)
! not a bad argument on its own -- the reference only fails to orient *this* rotation,
! which no check on any single argument could have foreseen
if (undefined_sign) call set_err_once(ierr, ERR_INVALID_INPUT)
```

**You get** the code propagated to the caller — and the **argument position cleared**
(`clear_err_arg_pos`) on the way out. That is deliberate: a position numbers the
*implementation's* dummy list, the wrapper's list is a different one, and a position often
propagates in unchanged from a private helper three frames down. A wrong argument name is worse
than none, and position 0 means "not argument related".

Keep `arg_pos=` in your implementation anyway. It is correct for a direct Fortran caller, and the
Fortran test suite is one.

### 5.15 A family too big for one file

**When** one family's implementations no longer fit comfortably in a single file.

**Write** several implementation modules under `src/tox/<family>/`, plus a **parent that holds no
procedures of its own and only `use`s its children**:

```fortran
module tox_data_integration_impl
    use tox_data_integration_preprocessing_impl
    use tox_data_integration_jsd_impl
    use tox_data_integration_per_family_impl
    use tox_data_integration_stats_impl
end module tox_data_integration_impl
```

**You get** the same shape mirrored over the wrappers — sub-directory and all, into
`src/generated/tox/<family>/` — so `use tox_data_integration` still reaches the whole family and
the split stays an implementation detail:

```fortran
module tox_data_integration
    use tox_data_integration_jsd
    use tox_data_integration_per_family
    use tox_data_integration_preprocessing
    use tox_data_integration_stats
end module tox_data_integration
```

Only children that actually generate something are re-exported, and a parent may gather another
parent. **Do not use `submodule`s** — they force every signature to be written twice (once in
the parent's `interface`, once in the implementation) for encapsulation this already provides.

---

# Part II — the export path

---

## 6. Exporting a hand-written procedure

The generator wraps any procedure marked `M_EXPORT_C` to C, Python and R, wherever it lives. It
just does not write its Fortran — so the validating wrapper the implementation path gives you for
free is, here, your own procedure's job.

### 6.1 When this is the right path

Because the procedure is not a numeric procedure of the pipeline:

- **`src/data/`** — the `tox_data_*` family: the zip archive (`libzip`), the CSV/TSV and
  OrthoFinder readers, the data-set validation and the accessors. One family, one directory,
  named for what the modules are rather than for the one mechanism only `tox_data_archive`
  performs.
- **`src/f42/`** — library-agnostic infrastructure: the serde family, whose procedures open files.
  Whether an f42 procedure is exported at all is a per-case judgement, which is exactly why the
  marker stays explicit here. It is not a licence: `f42_stats`, `f42_binary_search_tree` and
  `f42_kd_tree` were all hand-written exports until someone asked what a wrapper around them
  would do, and the answer was "validate what they were validating by hand, and allocate and sort
  what their callers were allocating and sorting" — so all three are implementations now.
- **inside an implementation module** — a recommend/sizing routine
  (`tox_loess_required_workspace`, `calc_neighborhood_size`) or a utility a caller genuinely
  needs (`mask_chunk_count`). These *must* be exported: `DM_OUTPUT_FROM(..., AUTO)` needs a
  wrapper to call (§5.9). The test is the same one: a wrapper around `calc_neighborhood_size`
  would validate nothing and prepare nothing, so there is nothing to generate.

If your procedure is a numeric procedure of the pipeline, take Part I instead. "It was easier to
export it directly" is how an unvalidated API gets shipped.

One name is not yours on this path either: **`_impl` is reserved across the whole of `src/`.** A
module named for it acquires generated wrappers wherever it sits, and every procedure in it is
held to §4 — it may allocate nothing, it may `use` only implementations and the whitelisted
infrastructure, its `_impl` procedures may not be exported, and it must live in a file named
after it. Nothing hand-written and exported may carry the suffix.

### 6.2 What you write: `M_EXPORT_C`

The marker goes in the procedure's Ford pre-comment, and that is the whole of it:

```fortran
!> M_EXPORT_C
!| summary: Compute the Empirical Distribution Function (EDF) from pre-sorted permutation
!| AUTHOR_JITU_DABA
!| Returns the sorted unique values and their cumulative frequencies in [0,1].
!| Assumes `values` is already sorted by `values[perm]`. Caller controls sorting algorithm.
pure subroutine compute_edf(values, n_values, perm, unique_values, cdf_values, n_unique, ierr)
```

`M_EXPORT_C` expands to a Ford `category` tag, and the generator reads the category from that
same macro — so the marker and what the generator recognises cannot drift.

**Everything in §4 still applies**: explicit intents, kinded types, a `summary:`, an author, a
`!!` on every argument, `M_IMPLICIT_NONE`, `#include <src/macros.h>`. Those are binding
requirements, not implementation requirements. Untagged procedures are held to none of it.

### 6.3 What you now do yourself: validate

This is the one real difference, and the one that bites.

> **The range macros generate checks only in generated wrappers.** On a hand-written procedure,
> `DM_MIN`, `DM_MAX`, `DM_SENTINEL`, `DM_ALLOW_NAN` and `DM_ALLOW_INFINITE` are **documentation
> only** — they render into the Ford, Python and R docs and validate nothing. There is likewise
> **no default finiteness contract**: a real argument here is checked only if you check it.

So the procedure declares its own `ierr`, opens with `set_ok`, and calls the `tox_errors`
validators itself — the same ones the generator would have used, with `arg_pos=` counted in its
own dummy list:

```fortran
subroutine serialize_int_helper(arr, n_elements, arr_shape, filename, ierr)
    ...
    call set_ok(ierr)
    call validate_in_range_int(n_elements, ierr, min=0_int32, arg_pos=1_int32)
    if (is_err(ierr)) return
```

Because you own the numbering, `arg_pos` here is *correct* and stays — nothing clears it (contrast
§5.14). Use `set_err_once` so the first failure is the one reported, and end validation with
`if (is_err(ierr)) return` before touching the data.

An exported procedure with no `ierr` gets one synthesised, because the binding languages always
need a channel to raise on. A procedure that can fail and does not declare one simply cannot
report — give it one.

### 6.4 What works exactly as in Part I

Everything that is a *binding* rule rather than a wrapper rule, which is most of the contract:

| Convention | Effect here |
|---|---|
| extents (`vec(n_dims)`) | derived by the binding from the array; never asked of a Python or R caller |
| `<arg>_mask` / `n_selected_<arg>` | the count is computed from the mask at the call |
| `tmp_<name>` | a work array: allocated by the binding language, never returned |
| `<arg>_shape` | a serialized array (§6.6) |
| `DM_DEFAULT` | the binding passes the evaluated default — `tox_data_tools` uses `DM_DEFAULT(char(9))` for a tab separator |
| `DM_OUTPUT_FROM(..., AUTO)` | the *language layer* calls the producer for the caller — `tox_data_archive` sizes eleven arguments this way from `get_tox_data_dims` |
| `DM_OUTPUT_FROM(..., JUST_INFO)` | documentation pointing at where the value comes from |
| `DM_RESULT_SIZE_IS` | as in Part I (§5.6): Python and R trim the result, Fortran callers still get the full buffer |
| a mode table | Python and R still pass the mode as a *string*, and the C wrapper still rejects an unknown one |

What is *not* available: the generated validation block, the plain wrapper's automatic allocation
and sorting, the prologue hook, and the mode split. Those are things the generator writes into a
wrapper, and here there is no wrapper for it to write into.

### 6.5 An `_alloc` pair by hand — gone

**Do not write one, and do not expect the generator to know what it is.** Write an `_impl` and
let §2 give you `foo` / `foo_expert`.

The shape was two hand-written procedures in one module, `<p>_alloc` and `<p>`, of which the
allocating one is the one callers want — how this framework spelled the two tiers before any
wrapper was generated. `abi/c_abi.py:stripped_name` used to translate it on the way out, so the
pair published as `foo` / `foo_expert` and matched what the implementation path produces.
`f42_kd_tree` was the last source writing one; it converted on 2026-08-11, and the translation
was retired on 2026-08-12. **Every procedure is now published under its own Fortran name.**

So a `foo_alloc` written today is published as `foo_alloc`. The reserved-name rule still refuses
`foo_alloc_impl`, because that name beside a generated `foo` reads as the allocating tier while
being an ordinary second procedure (see `design/impl-layer.md`).

Nothing about the shape was ever the obstacle to converting one, either. An implementation module
may hold ordinary exported procedures alongside its implementations, so no module ever had to be
split first.

### 6.6 An array of any rank through one signature

**When** a procedure must accept an array of any rank — serialization is the standing case.

**Write** a flat array plus a `<arg>_shape` companion that carries its extents:

```fortran
integer(int32), dimension(n_elements), intent(in) :: arr
    !! Array to be serialized
integer(int32), dimension(:), intent(in) :: arr_shape
    !! Extents of `arr`, one per dimension
```

`<arg>_shape` must be `intent(in)`, rank-1 integer, and **not optional** — the C wrapper reads it
to size `arr` before it may take `c_loc` of anything.

**You get** a Python caller who passes an array of any rank and never writes the shape out:

```python
serialize_int_helper(np.arange(12).reshape(3, 4), path)   # shape derived, flattened Fortran-order
```

Characters work the same way, with the length carried as a leading extent.

### 6.7 Never export an implementation

`M_EXPORT_C` on an `_impl` procedure does **not** replace wrapper generation — it would add to it.
The implementation still gets its wrappers, *and* the raw implementation is exported alongside, so
Python and R end up offering an unvalidated twin of the same call, one `_impl` suffix apart.

That is never what you want, and it never happens: the generator refuses the tag outright on an
`_impl` procedure of an `_impl` module (§8). If an implementation's functionality should be
reachable, it is reachable through its wrapper; if some *other* procedure in an implementation
module should be exported (a recommend routine, a utility), that one is not an `_impl` and §6.1
already covers it.

> Its bindings land in a module named for the file they came from, so the exports of
> `tox_loess_impl` are published as the Python module `tox_loess_impl` beside the wrappers'
> `tox_loess`. That is the one place the suffix is visible to a binding caller — a support routine
> is not part of the tiered API and is not renamed into it.

---

# Both paths

---

## 7. Documentation that survives four languages

Your `!!` and `!|` text is not just Ford docs. It is the C wrapper's comment, the Python
docstring and the R `.Rd` help page. So:

- **The generator never rewrites your prose.** It renders *its own* macro output per language,
  and resolves markup — nothing else. If a sentence is wrong in Python, it was wrong in the
  implementation; fix it there and all four fix at once.
- **The module's `!>` block is published too**, verbatim, as the generated module's own
  documentation — the Python module docstring and the Ford page a caller lands on. So write it
  for a caller: what the family is for, what its entry points are, what a result means. Not
  "Implementations for X", and not a description of what the generator will make of it — the
  generated file says where it came from, and this guide says how. The `_impl` modules were
  rewritten along these lines on 2026-08-12; `tox_get_outliers_impl` is a short worked example.
- **Literals are rendered per language.** `0.7_real64` reaches Python and R as `0.7`; `.true.`
  becomes `True` / `TRUE`.
- **Ford links resolve to what the reader can call.** `[[tox_loess_impl(module):MODE_PLAIN(variable)]]`
  becomes the mode *string* in a Python or R doc, and a link to an implementation becomes a link
  to the binding that wraps it — never to an implementation the reader cannot reach.
- **A link that names nothing is an error.** Name the module that *defines* the thing —
  `[[f42_math_impl(module):is_close(interface)]]`, not the parent that re-exports it — and note
  that a private routine counts as naming nothing, because Ford documents only what a module
  makes public. Linking to a *generated* module is right where the reader is meant to call it:
  `[[f42_binary_search_tree(module):build_bst_index]]` becomes the Python and R binding.
  The two type vocabularies are **not interchangeable**: a component takes `module`, `type`,
  `subroutine`, `function`, `proc`, `interface`, `absinterface`, `submodule`, `program`,
  `block`, `file`, `namelist` (any of them `ext`-prefixed for another project); an item takes
  `variable`, `bound`, `common`, `constructor`, `final`, `modproc` and the procedure/interface
  kinds. Give one the other's type and it is not a link at all.
  To *show* the syntax rather than follow it, wrap it in backticks — Ford 7 leaves inline code
  verbatim, and so does the generator.
- **Write plainly about arguments the reader actually passes.** Terms like "pre-allocated" or
  "workspace" describe the Fortran expert signature; a Python caller who never sees that
  argument only finds them confusing.
- **A `mode` mentioned in prose is misleading in a split wrapper** (§5.11), where that argument
  does not exist. Link to the mode-specific procedure instead.

Errors follow the same principle. An error names the argument the *caller* passed, and points at
the derived one only as a route:

```
invalid input arguments (argument 'vec1', via 'n_elements')
```

You get that for free, from the extent/mask/shape roles — nothing to annotate.

One thing worth knowing when you write a Python test or example against a procedure you
authored: **an array a Python binding returns is read-only.** A result is a value, so
modifying one in place raises, and `.copy()` is the way to a modifiable array. An
`intent(inout)` argument is untouched by this — it is the caller's array and is still modified
in place. R needs no equivalent, being copy-on-modify.

---

## 8. What the generator refuses

Errors. Nothing is written until they are fixed. Each diagnostic points at the line in *your*
source and carries a note saying what to write instead.

| Refused | Why |
|---|---|
| a deferred-length character (`len=:`) | implies allocatable/pointer; not interoperable |
| a character array of rank > 1 | `tox_conversions` converts up to a vector |
| an `allocatable` or `pointer` dummy | no C-interoperable form |
| a derived type | not mapped |
| a missing `intent` | constness and direction are undecidable |
| a kind with no C mapping | an error, never a guess — a wrong guess compiles and lies |
| an optional shape or extent argument | the wrapper must read it before it may take `c_loc` of what it sizes |
| a `tmp_` argument that is `intent(in)` | a work array is an output or an in-out |
| an **optional output** (`intent(out), optional`, and not a `tmp_`) | no binding can honour it — Python would return a dict whose keys vary per call, R a list of varying length. Express it as an optional input *flag* plus a `tmp_` work array, the way `loess_fit_plain_impl` does with `compute_influence`: the work stays skippable and the return type stays fixed |
| `DM_DEFAULT` **and** `DM_REQUIRED_IF_MODE` on a mode that does *not* split | with a runtime mode the argument is always passed on, so "required in that mode" says nothing. On a split mode (§5.11) the pair is meaningful and accepted |
| a mode argument that is not a scalar integer | modes are compared against `MODE_*` parameters |
| a mode table with no values, or a mode string matching no parameter | the table is what the C layer maps the caller's string through |
| a `DM_OUTPUT_FROM` producer input that is neither name-matched nor in the table | the generator will not guess what to pass |
| a misspelt `M_`/`CM_`/`DM_` in any doc comment | an unexpanded macro is a silently wrong document |
| an `M_EXPORT_C` on an **implementation** | its wrappers are the entry points; the export publishes an unvalidated twin beside them (§4) |
| an implementation module in a file not named after it | generation reads the module name while the cleaner and the Ford exclusion read file names — the two would name different files, and the generator would parse its own output next run (§4) |
| an implementation named `<something>_alloc_impl` or `<something>_expert_impl` | both suffixes belong to the wrappers: the first would generate a `foo_alloc` that reads as the allocating tier beside the generated `foo` that is one (§6.5), the second a second procedure called `foo_expert` (§4) |
| an `allocatable` local **or dummy** anywhere in an **implementation module** | the generated `foo` owns the memory, not the implementation (§4, §5.7) |
| a `use` in an **implementation module** that names neither another `_impl` module nor a whitelisted one | the bound on what it may reach is what makes the allocation rule hold across modules, and what keeps it below the wrappers rather than beside them (§4) |
| a `DM_PROLOGUE` naming a procedure that does not exist | the wrapper would be generated with no prologue at all (§5.13) |
| a prologue with no `handled`, or one that is not a scalar `logical, intent(out)` | the wrapper returns early on it regardless, so the branch would read an undefined value (§5.13) |
| a prologue dummy one edit from an implementation argument | a misspelling would otherwise become a new argument, and the two would be different values (§5.13) |
| a prologue dummy that some mode's wrapper does not have | the prologue runs in all of them (§5.13) |
| a prologue on an implementation with nothing to take over and no arguments of its own | only one wrapper is generated, and it has no prologue for this to run in (§5.13) |
| a prologue producing something the setup above it reads | the name resolves either way, so it would compile and compute rubbish (§5.13) |
| a Ford `[[...]]` link naming a module, or a name in it, that does not exist | it silently stops being a link in all four languages, and nothing downstream can tell (§7) |

Warnings never stop generation; errors write nothing. A few things are warnings rather than
errors, and the house bar is **zero warnings**, so treat them as errors anyway:

- a missing `summary:` or `author` meta tag, and an argument with no documentation — the API
  would be generated undocumented
- a value table on an argument that is *not* named as a mode (`baseline` rather than
  `baseline_mode`) — the table is then documentation only, and the value crosses the binding as
  a bare integer instead of its name


---

## 9. The workflow

```sh
./build.sh                          # regenerates, then compiles -- this is the normal loop
./build.sh --skip-code-generation   # compile the tree as it stands

python helper/generate_code.py            # regenerate only
python helper/generate_code.py --check    # write nothing; non-zero if the tree is stale
python helper/generate_code.py --target python
```

**Compiling the checks out.** `./build.sh --directive=NO_INPUT_VALIDATION` builds every
generated wrapper without its input validation — for a caller who has already established that
the inputs are good, typically an inner loop over data it produced itself. It gives up every
diagnostic in this guide, so it is a whole-build decision rather than one to take per call site.
What survives it: `call set_ok(ierr)`, so `ierr` is still defined, and every runtime error an
implementation raises itself (§5.14) — those are not input checks. The C layer's null checks stay
too; they guard against a segfault rather than a bad value.

`build.sh` runs the generator before fpm, so a source change and its generated layers cannot
drift apart in a build. It needs Python and [`ford`](https://forddocs.readthedocs.io); without
them it warns and compiles what is committed. CI checks the same thing with `--check`.

Then the suites, all of which must pass:

```sh
./run_all_tests.sh --skip-kinds-test --reuse-mod-files   # Fortran, then Python, then R
python -m pytest helper/codegen/tests -q                 # the generator's own suite
```

`run_all_tests.sh` exits non-zero if anything fails, and stops at a failed Fortran build rather
than reporting a library the Python and R suites could not have loaded. It prints the Python and
R files as pass/fail only, so run a failing one alone for the actual message:

```sh
./test_runner.sh --skip-kinds-test --reuse-mod-files     # Fortran, with output
python python/test/mod_test_<x>.py
Rscript r/test/mod_test_<x>.R
```

The Fortran suite is **not** part of `build.sh`. It has been left red for a whole conversion
before now; run it.

Two things that will bite when you convert an existing procedure to an implementation:

- **Generated wrappers take `ierr` last**, after the optionals. Every positional `ierr` in an
  old test lands on a mask, a quantile or a status instead. Name it: `ierr=ierr`.
- **`ierr` packs an argument position**, so `assert_equal_int(ierr, ERR_X)` no longer matches.
  Use `assert_err(ierr, code, msg, arg_pos)` from [`test/asserts.F90`](test/asserts.F90), which
  reports code and position separately.

---

## 10. Before you commit

**Both paths**

- [ ] `#include <src/macros.h>`, `M_IMPLICIT_NONE`, `private` + an explicit `public ::`.
- [ ] Every dummy has an `intent`, a kind, and a `!!` doc; the procedure has a `summary:` and an
      author.
- [ ] No optional dummy forwarded into a mandatory one (`M_DEFAULT_VAL` first).
- [ ] `python helper/generate_code.py` — **0 warnings** — then `--check` is clean.
- [ ] `./build.sh` green; Fortran, Python, R and codegen suites green.
- [ ] Nothing generated was edited by hand (`src/generated/`, `python/tensor_omics/`,
      `r/tensor_omics/`, `snippets/`); the regenerated diff is committed with the source change.
- [ ] If you renamed anything published, `python helper/generate_code.py --target snippets` too:
      `snippets/` is **not** in the default targets, so neither a plain run nor `--check`
      touches it, and it drifts silently.

**Part I — an implementation**

- [ ] Module `tox_*_impl` in a file of the same name, procedure `*_impl`, and **not**
      `M_EXPORT_C`.
- [ ] The module's `!>` block reads as documentation of the published API, not of the
      implementation — it is carried onto the generated module verbatim (§7).
- [ ] No `!!` on an ordinary comment inside a procedure *body*: `!!` is Ford's docmark, silent
      while the procedure is unexported and stray documentation lines in the wrapper once it is
      not. `grep -n '^\s*!! ' <file>` before converting one.
- [ ] Every bounded value carries `DM_MIN` / `DM_MAX` / `DM_SENTINEL`. Every real that may be
      non-finite carries the matching `DM_ALLOW_*` — and no other real does.
- [ ] No validation in the implementation; no `ierr` unless it reports a genuine runtime failure.
- [ ] The generated wrapper was read once, and validates what you expected it to — and where two
      were generated, the plain one prepares what you meant it to.

**Part II — a hand-written export**

- [ ] `!> M_EXPORT_C`, and the procedure is genuinely not an implementation (§6.1).
- [ ] It declares `ierr`, opens with `set_ok`, **validates every argument itself**, and ends
      validation with `if (is_err(ierr)) return`.
- [ ] `arg_pos=` on each validator matches this procedure's own dummy list.
- [ ] Any `DM_MIN`/`DM_MAX`/`DM_ALLOW_*` it documents is backed by a check it actually makes —
      here they are prose, not code (§6.3).
