# Writing code the generator can wrap

TensorOmics generates its public API. There are two ways in, and the guide is split along them:

- **The kernel path** (Part I) — for a numeric procedure of the pipeline. You write one annotated
  kernel; the generator writes the validating wrapper, the allocating wrapper, and the C, Python
  and R bindings, with the documentation, input validation and error handling derived from your
  kernel rather than restated four times by hand.
- **The export path** (Part II) — for IO, infrastructure, and everything that is not a kernel. You
  write the whole procedure, validation included, and mark it `M_EXPORT_C`; the generator writes
  the three bindings.

Either way this guide is the author's side of the contract: **what to write so the generator can
do its half**, case by case, with a real example for each. Every snippet here is taken from the
current `src/` tree.

Related reading, once this is not enough:

| Document | What it answers |
|---|---|
| [`helper/codegen/README.md`](helper/codegen/README.md) | what the generator produces, how to run it, how it is built and tested |
| [`helper/codegen/design/kernel-layer.md`](helper/codegen/design/kernel-layer.md) | *why* the kernel layer works this way, and every alternative that was rejected |
| [`helper/codegen/design/c-layer.md`](helper/codegen/design/c-layer.md) | the Fortran↔C wrapper decisions |
| [`helper/codegen/design/language-layers.md`](helper/codegen/design/language-layers.md) | the Python and R decisions |
| [`src/macros.h`](src/macros.h) | the macros themselves, each with its contract in a comment |

---

## Contents

- [1. The one rule](#1-the-one-rule)
- [2. Where code lives, and which path you are on](#2-where-code-lives-and-which-path-you-are-on)

**Part I — the kernel path** (a numeric procedure of the TOX pipeline)

- [3. A kernel end to end](#3-a-kernel-end-to-end)
- [4. What every kernel must have](#4-what-every-kernel-must-have)
- [5. Case by case](#5-case-by-case)
  - [5.1 A value with a valid range](#51-a-value-with-a-valid-range)
  - [5.2 A real that may be NaN or infinite](#52-a-real-that-may-be-nan-or-infinite)
  - [5.3 Sizes, masks and counts the caller never passes](#53-sizes-masks-and-counts-the-caller-never-passes)
  - [5.4 A distance matrix](#54-a-distance-matrix)
  - [5.5 An optional argument and its default](#55-an-optional-argument-and-its-default)
  - [5.6 An output filled only partially](#56-an-output-filled-only-partially)
  - [5.7 Work arrays, and how `_alloc` appears](#57-work-arrays-and-how-_alloc-appears)
  - [5.8 A permutation](#58-a-permutation)
  - [5.9 A workspace sized by a recommend routine](#59-a-workspace-sized-by-a-recommend-routine)
  - [5.10 A mode argument](#510-a-mode-argument)
  - [5.11 One procedure per mode](#511-one-procedure-per-mode)
  - [5.12 An argument only one mode needs](#512-an-argument-only-one-mode-needs)
  - [5.13 Work that must happen before the kernel](#513-work-that-must-happen-before-the-kernel)
  - [5.14 A genuine runtime error](#514-a-genuine-runtime-error)
  - [5.15 A family too big for one file](#515-a-family-too-big-for-one-file)

**Part II — the export path** (IO, infrastructure, and everything that is not a kernel)

- [6. Exporting a hand-written procedure](#6-exporting-a-hand-written-procedure)
  - [6.1 When this is the right path](#61-when-this-is-the-right-path)
  - [6.2 What you write: `M_EXPORT_C`](#62-what-you-write-m_export_c)
  - [6.3 What you now do yourself: validate](#63-what-you-now-do-yourself-validate)
  - [6.4 What works exactly as in Part I](#64-what-works-exactly-as-in-part-i)
  - [6.5 An `_alloc` pair by hand](#65-an-_alloc-pair-by-hand)
  - [6.6 An array of any rank through one signature](#66-an-array-of-any-rank-through-one-signature)
  - [6.7 Never export a `_kernel` procedure](#67-never-export-a-_kernel-procedure)

**Both paths**

- [7. Documentation that survives four languages](#7-documentation-that-survives-four-languages)
- [8. What the generator refuses](#8-what-the-generator-refuses)
- [9. The workflow](#9-the-workflow)
- [10. Before you commit](#10-before-you-commit)

---

## 1. The one rule

> **Write the implementation. Annotate it. Write nothing else.**

Everything the generator emits follows from your kernel's signature and its documentation. So
a fact belongs in exactly one place — the kernel — and if the generator cannot derive
something, the kernel gains an annotation, never the wrapper hand-written code.

The corollary is the part worth internalising, because it will decide arguments you have one
day:

> **If the generator cannot express your procedure's contract, the procedure is wrong — not the generator.**

The annotation vocabulary is the API review. A signature that needs a rule none of the macros
can state is a signature that is asking its callers for the wrong thing. Redesign it; that is
almost always cheap, because only the tests depend on it. If a genuinely new *kind* of
constraint ever turns up, it earns a new convention — a macro with a name, a meaning, and every
kernel that qualifies using it — never a per-procedure escape hatch.

---

## 2. Where code lives, and which path you are on

```
src/
  macros.h            included by every source, by this path
  f42/                infrastructure, library-agnostic
    utils/              f42_utils re-exports f42_math, _sort, _random, _vector, _stats
    serde/              likewise, per element type
  kernel/             the kernels -- the API's source of truth, hand-written
  data/               the hand-written data-set API (tox_data_*), incl. the zip archive
  generated/          NOTHING here is hand-written
    tox/                the wrappers
    bindings/c/         the Fortran C wrappers
    bindings/r/         the R `.Call` shims
```

The rule inside `src/` is one line: **edit anything outside `src/generated/`.** That prefix is
what the generator deletes and rewrites on every run — a change you make there survives exactly
until the next build.

Two more trees are generated in full and live outside `src/` only because that is where their
languages expect them:

```
python/tensor_omics/    the Python package
r/tensor_omics/         the R package
snippets/               the VS Code snippets (except the hand-written toxdev_snippets.json)
```

All four are marked `linguist-generated` in [`.gitattributes`](.gitattributes), so a review
collapses them and sees the kernel that changed rather than the fan-out it produced. The test
suites under `python/test/` and `r/test/` are hand-written — those you do edit.

Two names are load-bearing, and together they are the generation trigger:

- the file is under **`src/kernel/`**, in a module named `tox_<family>_kernel`
- the procedure is named **`<name>_kernel`**

That is all. There is no marker macro to forget — the suffix is part of every call site, which
matters, because calling a raw kernel means calling something with no input validation.

The generated wrapper takes the clean name: `compute_shift_vector_field_kernel` in
`tox_shift_vectors_kernel` becomes `compute_shift_vector_field` in `tox_shift_vectors`. **That
clean name is the public API** — in Fortran, C, Python and R alike.

### Which path is yours

Not everything is a kernel, and a third of the public API is not — 44 hand-written exports
against 85 generated wrappers today. File and archive IO, the f42 trees and statistics, the whole
serde family, and the sizing routines the kernels themselves call are all **hand-written and
exported**: the generator wraps them to C, Python and R, but does not write their Fortran.

| Your procedure | Path | Where | Marker |
|---|---|---|---|
| a numeric procedure of the TOX pipeline | **Part I — kernel** | `src/kernel/` | the `_kernel` name |
| the data-set API: archive, CSV/TSV, validation, accessors | **Part II — export** | `src/data/` | `M_EXPORT_C` |
| library-agnostic infrastructure | **Part II — export** | `src/f42/` | `M_EXPORT_C` |
| a sizing / utility routine a kernel or a caller needs | **Part II — export** | the kernel module, `public` | `M_EXPORT_C` |

The two paths share everything about the *bindings* — naming conventions, documentation, type
rules, the refusals in §8. They differ in exactly one thing: **who writes the validation.** On
the kernel path the generator does. On the export path you do (§6.3).

---

# Part I — the kernel path

---

## 3. A kernel end to end

The whole of `src/kernel/tox_shift_vectors_kernel.F90`, minus the body:

```fortran
#include <src/macros.h>

!> Kernel for computing the shift vector field for all genes.
!|
!| Hand-written implementation only. The generator turns this into the validating wrapper
!| [[tox_shift_vectors(module):compute_shift_vector_field]] in module `tox_shift_vectors`.
module tox_shift_vectors_kernel
    use, intrinsic :: iso_fortran_env, only: real64, int32
    M_IMPLICIT_NONE
    private
    public :: compute_shift_vector_field_kernel
contains

    !> summary: Compute the shift vector field for all genes.
    !| AUTHOR_ALEXANDER_SCHWARZPAUL
    !| Computes the shift vectors by subtracting the corresponding family centroid from the expression vector.
    pure subroutine compute_shift_vector_field_kernel(n_tissues, n_genes, n_families, expression_vectors, &
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
    end subroutine compute_shift_vector_field_kernel
end module tox_shift_vectors_kernel
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
        call validate_dimension_size(n_tissues, ierr, arg_pos=1_int32)
        call validate_dimension_size(n_genes, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_families, ierr, arg_pos=3_int32)
        call validate_all_in_range_real(expression_vectors, n_tissues * n_genes, ierr, arg_pos=4_int32)
        call validate_all_in_range_real(family_centroids, n_tissues * n_families, ierr, arg_pos=5_int32)
        call validate_all_in_range_int(gene_to_fam, n_genes, ierr, arg_pos=6_int32, &
                                       min=1_int32, max=n_families, sentinel=0_int32)
        if (is_err(ierr)) return

        call compute_shift_vector_field_kernel(...)
    end subroutine compute_shift_vector_field
```

Three of those six checks you never asked for: the extents are validated because they *are*
extents, and both real matrices are checked for NaN and infinity because that is the framework's
default contract (§5.2).

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

## 4. What every kernel must have

| Requirement | Why |
|---|---|
| `#include <src/macros.h>` at the top of the file | a `DM_` that never expands is an error, as is a misspelt `M_`/`CM_`/`DM_` in any doc comment |
| `M_IMPLICIT_NONE` in the module | `implicit none (type, external)` — a typo'd call is a compile error, not a link-time surprise |
| `private` + an explicit `public ::` list | the kernel is API by way of its wrapper, not by accident |
| Explicit `intent` on every dummy | decides constness in C and R, and input/output/in-out everywhere |
| Kinded numeric types (`real(real64)`, `integer(int32)`) | a default kind has no defensible C mapping |
| `!> summary: ...` on the procedure | becomes the docstring in Python and R |
| An author tag (`AUTHOR_*` from [`authors.h`](authors.h)) | attribution, rendered into the Ford docs |
| A `!!` doc on **every** argument | inherited by the C wrapper and by both language layers |
| **No `ierr` for validation** | validation is the wrapper's job; see §5.14 for the one case where a kernel keeps an `ierr` |
| **No `M_EXPORT_C` on the kernel** | the generated wrapper is what the bindings call; exporting the kernel beside it publishes an unvalidated twin under a name a caller cannot tell apart. Support routines in the same module — the recommend routines of §5.9 — *are* exported: they have no wrapper |
| **No `_alloc` in the kernel's name** | `_alloc` is the generator's suffix. Name the kernel for what it computes; whether an allocating wrapper appears is decided by its `tmp_` arguments (§5.7) |
| **No allocation, anywhere in the module** | every buffer is a `tmp_` argument, so the generated `_alloc` owns the memory and an expert caller can hand in buffers it already has. The rule covers the module's helpers too: a kernel that allocates nothing itself but calls a helper that does is no better off. Enforced on the declaration — a local declared `allocatable` is refused. A `pointer` local is fine: aliasing a buffer you were handed allocates nothing |

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
> is a null dereference when it is absent, and it is exactly how `loess_fit_robust_kernel`
> segfaulted. Resolve it with `M_DEFAULT_VAL` first and pass the local.

An optional **without** a default is nullable: absent in Fortran, `None`/`NULL` from the
language layers. An optional **with** one is always passed, which keeps the wrapper flat.

### 5.6 An output filled only partially

**When** a result array is sized for the worst case but only its leading elements are filled —
a selection whose size is not known until the work is done.

**Write** `DM_RESULT_SIZE_IS(<argument>)` on the result, naming the `intent(out)` scalar integer
your kernel sets to the number of elements it actually filled:

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
> No kernel in `src/` uses this yet, so the fixture is the worked case.

### 5.7 Work arrays, and how `_alloc` appears

**When** your kernel needs scratch space it does not want to allocate itself (because a hot loop
should not allocate, and because a Fortran caller may want to reuse a buffer).

**Write** the argument with a `tmp_` prefix and an `intent(out)` or `intent(inout)`:

```fortran
integer(int32), intent(out) :: tmp_stack_left(n_genes)
    !! Stack array for left indices during sorting
```

**You get** a *second* wrapper. `foo` keeps the work arrays — that is the expert entry point,
for a caller who manages buffers. `foo_alloc` drops them from its signature, declares them
locally as allocatables, `M_ALLOCATE`s them and calls the kernel:

```fortran
integer(int32), dimension(:), allocatable :: tmp_stack_left
...
M_ALLOCATE(tmp_stack_left(n_genes))
```

In the language layers this pairing is automatic: `foo_alloc` becomes the plain name a Python or
R caller sees, and `foo` becomes `foo_expert`. A kernel with no work arrays, permutations or
recommend-sized buffers generates only `foo`.

> A `tmp_` argument that is `intent(in)` is an error. A work array is an output or an in-out;
> nothing else makes sense.

This is the *only* way a kernel module gets scratch space: an `allocatable` local anywhere in it
is refused (§4). A buffer whose size is not an expression over the other arguments is still a
`tmp_` argument — `DM_OUTPUT_FROM(..., AUTO)` (§5.9) names the routine that sizes it. Where even
that routine cannot be called ahead of time because the size depends on something only the kernel
discovers, size the buffer at the **upper bound** and slice it: `normalize_by_std_dev_kernel`
takes its LOESS workspace for all `n_genes` and hands the fit `tmp_loess_x(1:n_valid)`, because
dropping the zero-variance genes only ever makes the fit smaller.

### 5.8 A permutation

**When** the kernel needs an index permutation sorted against one of its arrays.

**Write** the argument as `<base>_perm`, where `<base>` is the array it orders.

```fortran
real(real64), intent(in) :: distances(n_genes)
integer(int32), intent(out) :: tmp_perm(n_genes)
    !! Permutation array for sorting gene distances
```

**You get**, in `foo_alloc` only: the array allocated, seeded with `init_perm`, and heapsorted
against `<base>`. `foo` still takes it, because an expert caller may already have it sorted.

`foo_alloc` calls the **kernel directly**, not `foo` — it just built that permutation itself, so
re-running an O(n) validation over `[1..n]` would be wasted work. A bare `perm` with no base
name is not a permutation by this convention; it stays an ordinary argument.

### 5.9 A workspace sized by a recommend routine

**When** a buffer's size can only be computed by calling something — netlib's LOESS workspaces
are the standing example.

**Write** `DM_OUTPUT_FROM(<size_arg>, <producer>, <its module>, AUTO)` on the size argument.
Producer inputs are matched to your arguments **by name**; where the two spell a quantity
differently, or where the producer wants a constant, a table says so:

```fortran
integer(int32), intent(in) :: int_workspace_size
    !! Length of integer workspace.
    !! DM_OUTPUT_FROM(int_workspace_size, tox_loess_required_workspace, tox_loess_kernel, AUTO)
    !!
    !! | Producer input        | Supplied by |
    !! |-----------------------|-------------|
    !! | n_dim                 | 1_int32     |
    !! | max_neighborhood_size | n_families  |
    !! | save_factorization    | .false.     |
integer(int32), intent(out) :: tmp_int_workspace(int_workspace_size)
    !! Integer workspace array
```

**You get** the producer called for you in `foo_alloc`, before the allocation that needs it:

```fortran
call tox_loess_required_workspace(n_dim = 1_int32, max_neighborhood_size = n_families, &
                                  int_workspace_size = int_workspace_size, &
                                  real_workspace_size = real_workspace_size, &
                                  save_factorization = .false.)
M_ALLOCATE(tmp_int_workspace(int_workspace_size))
```

— and, in the `_expert` binding, the size argument still exposed with its documentation pointing
at the routine that computes it. One annotation, two consumers.

The producer must be exported (`M_EXPORT_C`), so that a wrapper exists to call; recommend
routines therefore live in the kernel module, public, tagged. Use `JUST_INFO` instead of `AUTO`
when the caller genuinely has to make the call themselves — then the docs say where to get the
value, and nothing is called for them.

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
    !! | Plain LOESS fitting  | [[tox_loess_kernel(module):MODE_PLAIN(variable)]]  |
    !! | Robust LOESS fitting | [[tox_loess_kernel(module):MODE_ROBUST(variable)]] |
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
    !! |    Dosage Effect     | [[tox_paralog_analysis_kernel(module):MODE_DOSAGE_PATTERN(variable)]]  | detect_dosage_effect        |
    !! | Subfunctionalization | [[tox_paralog_analysis_kernel(module):MODE_SUBFUNC_PATTERN(variable)]] | detect_subfunctionalization |
```

**You get** one wrapper per mode value — named from the column, with the `mode` dummy dropped
and fixed internally — plus its `_alloc` where the kernel needs one (§5.7):
`detect_dosage_effect`, `detect_dosage_effect_alloc`, `detect_subfunctionalization`,
`detect_subfunctionalization_alloc`, and (no work arrays, so no `_alloc`)
`filter_paralogs_by_pattern_dosage_effect` and `filter_paralogs_by_pattern_subfunctionalization`.

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
    !! DM_REQUIRED_IF_MODE(pattern_mode, tox_paralog_analysis_kernel, MODE_DOSAGE_PATTERN)
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

### 5.13 Work that must happen before the kernel

**When** something has to run *before* the kernel and may make the kernel unnecessary — a
degenerate input that already determines the answer, or preparation a caller should not have to
do.

**Write** `DM_PROLOGUE(<procedure>, <module>, <scope>)` in the procedure's own doc block, with
scope `EXPERT`, `ALLOC` or `BOTH`. The prologue takes the **kernel's** arguments by name, plus
`handled` and `ierr`:

```fortran
!> summary: Perform plain LOESS fitting
!| AUTHOR_FRANZ_ERIC_SILL
!| DM_PROLOGUE(loess_degenerate_fit, tox_loess_kernel, BOTH)
```
```fortran
pure subroutine loess_degenerate_fit(n, x, y, degree, fitted_values, handled, ierr)
    logical, intent(out) :: handled
        !! `.true.` when the data was degenerate and `fitted_values` already holds the answer
```

**You get** the prologue called first in the wrapper, and the kernel skipped when it says so:

```fortran
call loess_degenerate_fit(n = n, x = x, y = y, degree = degree, &
                          fitted_values = fitted_values, handled = handled, ierr = ierr)
if (is_err(ierr)) return
if (handled) return
```

**Where it runs.** As early as it can. In the allocating wrapper that means *above* the work
arrays, so a prologue that refuses a degenerate input spares the whole setup — half of why a
prologue exists. A prologue that takes one of those work arrays cannot run there (there would be
nothing to hand it), so it runs below them instead. You do not choose this; what the prologue
takes decides it. The one consequence: a late prologue may not produce anything the allocations
or the recommend calls above it read, and the generator refuses it if you try.

```fortran
!| DM_PROLOGUE(prepare_ranking, tox_outliers_kernel, ALLOC)
...
real(real64), intent(out) :: tmp_valid_perm(n_genes)
    !! the prologue fills this, so the prologue runs after it is allocated
```

**What is refused**, all of it silent before: a `DM_PROLOGUE` naming a procedure that does not
exist; a prologue with no `handled`, or one that is not a scalar `logical, intent(out)` (the
wrapper returns early on it regardless, so without it the branch reads an undefined value); a
dummy naming nothing the kernel has. There is deliberately **no rename table** — unlike a
`DM_OUTPUT_FROM` producer, a prologue is internal to the kernel module, so the fix is to rename
its dummy.

Reach for this only when the work is genuinely not the kernel's. A prologue that merely
*validates* is a smell — that is what the range macros are for, and a pre-validation pass that
duplicates them is a bug rather than a feature.

### 5.14 A genuine runtime error

**When** the kernel can fail for a reason no input check could have foreseen: an external
library reporting failure, a division by zero, a configuration of otherwise-valid inputs that
has no answer.

**Write** an `ierr` on the kernel after all, set it with `set_err_once`, and say so in the
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
(`clear_err_arg_pos`) on the way out. That is deliberate: a position numbers the *kernel's*
dummy list, the wrapper's list is a different one, and a position often propagates in unchanged
from a private helper three frames down. A wrong argument name is worse than none, and position
0 means "not argument related".

Keep `arg_pos=` in your kernel anyway. It is correct for a direct Fortran caller, and the
Fortran test suite is one.

### 5.15 A family too big for one file

**When** one family's kernels no longer fit comfortably in a single file.

**Write** several kernel modules under `src/kernel/<family>/`, plus a **parent that holds no
procedures of its own and only `use`s its children**:

```fortran
module tox_data_integration_kernel
    use tox_data_integration_preprocessing_kernel
    use tox_data_integration_jsd_kernel
    use tox_data_integration_per_family_kernel
    use tox_data_integration_stats_kernel
end module tox_data_integration_kernel
```

**You get** the same shape mirrored over the wrappers, so `use tox_data_integration` still
reaches the whole family and the split stays an implementation detail:

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
just does not write its Fortran — so the validating wrapper the kernel path gives you for free
is, here, your own procedure's job.

### 6.1 When this is the right path

Because the procedure is not a numeric kernel of the pipeline:

- **`src/data/`** — the `tox_data_*` family: the zip archive (`libzip`), the CSV/TSV and
  OrthoFinder readers, the data-set validation and the accessors. One family, one directory,
  named for what the modules are rather than for the one mechanism only `tox_data_archive`
  performs.
- **`src/f42/`** — library-agnostic infrastructure: `f42_kd_tree`, `f42_binary_search_tree`,
  `f42_stats`, the serde family. Whether an f42 procedure is exported at all is a per-case
  judgement, which is exactly why the marker stays explicit here.
- **inside a kernel module** — a recommend/sizing routine (`tox_loess_required_workspace`,
  `calc_neighborhood_size`) or a utility a caller genuinely needs (`mask_chunk_count`). These
  *must* be exported: `DM_OUTPUT_FROM(..., AUTO)` needs a wrapper to call (§5.9).

If your procedure is a numeric kernel, take Part I instead. "It was easier to export it directly"
is how an unvalidated API gets shipped.

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
requirements, not kernel requirements. Untagged procedures are held to none of it.

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
| `<p>_alloc` / `<p>` | the same pairing, written by hand (§6.5) |

What is *not* available: the generated validation block, `foo_alloc`'s automatic allocation and
sorting, the prologue hook, and the mode split. Those are things the generator writes into a
wrapper, and here there is no wrapper for it to write into.

### 6.5 An `_alloc` pair by hand

The `_alloc` ↔ `_expert` convention is a binding rule, so it works for hand-written procedures
too: write both procedures in one module, name them `<p>_alloc` and `<p>`, and export what you
want callers to have.

```fortran
!> M_EXPORT_C
!| summary: Build a k-d tree index using a stack-based, non-recursive approach
pure subroutine build_kd_index_alloc(points, n_dimensions, n_points, kd_indices, dimension_order, ierr)

! not exported: the expert form, for a Fortran caller that owns the buffers
pure subroutine build_kd_index(points, n_dimensions, n_points, kd_indices, dimension_order, &
                               tmp_workspace, tmp_value_buffer, tmp_permutation, tmp_recursion_stack, ierr)
```

Export **both** and the C symbols are `build_kd_index_c` (from `_alloc`) and
`build_kd_index_expert_c` (from the plain one), exactly as on the kernel path. Export only the
`_alloc`, as `f42_kd_tree` does, and the buffer-managing form stays a Fortran-only entry point —
a deliberate choice, not an oversight: there is no sensible way for a Python caller to own those
buffers.

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

### 6.7 Never export a `_kernel` procedure

`M_EXPORT_C` on a `_kernel` procedure does **not** replace wrapper generation — it adds to it.
The kernel still gets its validating wrapper, *and* the raw kernel is exported alongside, so
Python and R end up offering an unvalidated twin of the same call, one `_kernel` suffix apart.

That is never what you want. If a kernel's functionality should be reachable, it is reachable
through its wrapper; if some *other* procedure in a kernel module should be exported (a recommend
routine, a utility), that one is not a `_kernel` and §6.1 already covers it.

> Present state: `tox_trajectory_contribution_analysis_kernel` still carries five such tags, left
> over from its conversion. They generate `compute_all_contributions_kernel` and four siblings in
> the Python and R packages, beside the wrappers that validate. Nothing calls them.

---

# Both paths

---

## 7. Documentation that survives four languages

Your `!!` and `!|` text is not just Ford docs. It is the C wrapper's comment, the Python
docstring and the R `.Rd` help page. So:

- **The generator never rewrites your prose.** It renders *its own* macro output per language,
  and resolves markup — nothing else. If a sentence is wrong in Python, it was wrong in the
  kernel; fix it there and all four fix at once.
- **Literals are rendered per language.** `0.7_real64` reaches Python and R as `0.7`; `.true.`
  becomes `True` / `TRUE`.
- **Ford links resolve to what the reader can call.** `[[tox_loess_kernel(module):MODE_PLAIN(variable)]]`
  becomes the mode *string* in a Python or R doc, and a link to a kernel becomes a link to the
  binding that wraps it — never to a kernel the reader cannot reach.
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
| an **optional output** (`intent(out), optional`, and not a `tmp_`) | no binding can honour it — Python would return a dict whose keys vary per call, R a list of varying length. Express it as an optional input *flag* plus a `tmp_` work array, the way `loess_fit_plain` does with `compute_influence`: the work stays skippable and the return type stays fixed |
| `DM_DEFAULT` **and** `DM_REQUIRED_IF_MODE` on a mode that does *not* split | with a runtime mode the argument is always passed on, so "required in that mode" says nothing. On a split mode (§5.11) the pair is meaningful and accepted |
| a mode argument that is not a scalar integer | modes are compared against `MODE_*` parameters |
| a mode table with no values, or a mode string matching no parameter | the table is what the C layer maps the caller's string through |
| a `DM_OUTPUT_FROM` producer input that is neither name-matched nor in the table | the generator will not guess what to pass |
| a misspelt `M_`/`CM_`/`DM_` in any doc comment | an unexpanded macro is a silently wrong document |
| an `M_EXPORT_C` on a **kernel** | its wrapper is the entry point; the export publishes an unvalidated twin beside it (§4) |
| a kernel named `<something>_alloc_kernel` | `_alloc` is the generator's suffix, and such a kernel generates an allocating wrapper that allocates nothing (§4) |
| an `allocatable` local anywhere in a **kernel module** | the generated `_alloc` owns the memory, not the kernel (§4, §5.7) |
| a `DM_PROLOGUE` naming a procedure that does not exist | the wrapper would be generated with no prologue at all (§5.13) |
| a prologue with no `handled`, or one that is not a scalar `logical, intent(out)` | the wrapper returns early on it regardless, so the branch would read an undefined value (§5.13) |
| a prologue dummy naming nothing the kernel has | it would be dropped from the generated call, and Fortran would reject code you did not write (§5.13) |
| a late prologue producing something the allocations above it read | the name resolves either way, so it would compile and compute rubbish (§5.13) |

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

`build.sh` runs the generator before fpm, so a source change and its generated layers cannot
drift apart in a build. It needs Python and [`ford`](https://forddocs.readthedocs.io); without
them it warns and compiles what is committed. CI checks the same thing with `--check`.

Then the suites, all of which must pass:

```sh
./run_all_tests.sh --skip-kinds-test --reuse-mod-files   # Fortran, then Python, then R
python -m pytest helper/codegen/tests -q                 # the generator's own suite
```

`run_all_tests.sh` reports the Python and R files pass/fail only. When one fails, run it alone
for the actual message:

```sh
./test_runner.sh --skip-kinds-test --reuse-mod-files     # Fortran, with output
python python/test/mod_test_<x>.py
Rscript r/test/mod_test_<x>.R
```

The Fortran suite is **not** part of `build.sh`. It has been left red for a whole conversion
before now; run it.

Two things that will bite when you convert an existing procedure to a kernel:

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

**Part I — a kernel**

- [ ] In `src/kernel/`, module `tox_*_kernel`, procedure `*_kernel`, and **not** `M_EXPORT_C`.
- [ ] Every bounded value carries `DM_MIN` / `DM_MAX` / `DM_SENTINEL`. Every real that may be
      non-finite carries the matching `DM_ALLOW_*` — and no other real does.
- [ ] No validation in the kernel; no `ierr` unless it reports a genuine runtime failure.
- [ ] The generated wrapper was read once, and validates what you expected it to.

**Part II — a hand-written export**

- [ ] `!> M_EXPORT_C`, and the procedure is genuinely not a kernel (§6.1).
- [ ] It declares `ierr`, opens with `set_ok`, **validates every argument itself**, and ends
      validation with `if (is_err(ierr)) return`.
- [ ] `arg_pos=` on each validator matches this procedure's own dummy list.
- [ ] Any `DM_MIN`/`DM_MAX`/`DM_ALLOW_*` it documents is backed by a check it actually makes —
      here they are prose, not code (§6.3).
