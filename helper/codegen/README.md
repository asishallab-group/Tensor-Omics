# TensorOmics binding generator

Generates the C, Python and R bindings to the Fortran library from the Fortran sources
themselves. A procedure marked for export becomes a C-callable wrapper, a Python function,
and an R function — with the documentation, error handling and argument validation derived
from the Fortran, not restated by hand.

This document is the whole of it: what it produces, how to run it, how it is built, and —
the part a library author needs — **how to write Fortran the generator can wrap**. The
design rationale (why each choice, and what was rejected) lives in [`design/`](design/);
this points there rather than repeating it.

---

## Contents

- [What it produces](#what-it-produces)
- [Running it](#running-it)
- [Configuration: what is read, and where output goes](#configuration-what-is-read-and-where-output-goes)
- [How it is built](#how-it-is-built)
- [Writing generator-compliant Fortran](#writing-generator-compliant-fortran)
  - [Marking a procedure for export](#marking-a-procedure-for-export)
  - [What every exported procedure needs](#what-every-exported-procedure-needs)
  - [Naming conventions](#naming-conventions)
  - [Mode arguments](#mode-arguments)
  - [Documentation macros (`DM_`)](#documentation-macros-dm_)
  - [Serialized arrays](#serialized-arrays)
  - [What is rejected, and why](#what-is-rejected-and-why)
- [Edge cases handled](#edge-cases-handled)
- [How it is tested](#how-it-is-tested)
- [Extending it](#extending-it)
- [Open items](#open-items)

---

## What it produces

From one exported Fortran procedure, three things:

| Target | Output | What it is |
|---|---|---|
| C | `src/bindings/c/<module>_c.F90` | a `bind(C)` wrapper: a plain-pointer ABI, null validation, type conversion |
| Python | `python/tensor_omics/<module>.py` | a `ctypes` function with a numpydoc docstring |
| R | `src/bindings/r/<module>.c` + `r/tensor_omics/<module>.R` | a C `.Call` shim (marshalling, bundled into the `.so`) under an R function (validation, docs) |
| Snippets | `snippets/<Language>_<root>_snippets.json` | VS Code call/setup snippets, split by language and module root |

Plus, once per project: an error module for each language, generated from `tox_errors`, and
the loader/marshalling scaffolding.

The `snippets` target (`emit/vscode_snippets.py`) emits VS Code snippets split into six
files -- `{Fortran,Python,R}_{f42,tox}_snippets.json` under `snippets/` -- the language in
the file name (so no per-snippet `scope`) and the root keeping the two namespaces apart.
These are **regenerated artifacts, git-ignored** (like the Python/R packages); the
hand-written `snippets/toxdev_snippets.json` sits alongside and is tracked. Per exported procedure: a native Fortran `call` (plus a variant that guards
`ierr`) and a wrapper call for Python and R, arguments rendered as keyword tabstops -- a
`mode`/`method` argument becomes a *choice* of its accepted values. Per module: a Fortran
`use ..., only:` and a Python import (both a choice of that module's procedures). Plus
generic aids: an R loader `source`, an error-handling wrapper per binding language, and
an error-code picker. Every prefix starts with the module's root -- `f42:` or `tox:`, taken
from the module name up to its first underscore.

The generated bindings are consistent by construction: all three read the *same* model of
what a procedure looks like from C (the `abi` layer), so they cannot drift.

---

## Running it

```sh
python helper/generate_code.py            # generate everything into the tree
python helper/generate_code.py --check    # report problems, write nothing, non-zero if stale
python helper/generate_code.py --target python   # one target only
python helper/generate_code.py --help
```

Run from the repository root. `--check` is what a CI guard wants: it exits non-zero if the
committed bindings no longer match the sources.

Requirements: Python 3.11+ (`contextlib.chdir`), and `FORD`, `pcpp`, `numpy` (see
`requirements-dev.txt`). The generator reads the Ford settings from `fpm.toml`.

Diagnostics point at the offending line of the *original* source, with the entity chain and
a note on what to do:

```
error: shape argument 'data_shape' is optional
  --> src/tox_x.F90:90
  argument 'data_shape' in procedure 'x_alloc' in module 'tox_x'
  note: the C wrapper reads it before it may take c_loc of the arrays it describes, so it
        has to be there; see the null validation order in the generator README
```

Warnings never stop generation; errors write nothing.

---

## Configuration: what is read, and where output goes

All paths are relative to the repository root (`--root`, default the current directory) and
live in [`config.py`](config.py) as `Paths`. The defaults:

| Setting | Default | Role |
|---|---|---|
| `src_dir` | `src` | the sources to read (`--src`) |
| `macros_header` | `src/macros.h` | the macro definitions, incl. the `DM_` doc macros |
| `c_binding_dir` | `src/bindings/c` | **output**: the Fortran C wrappers |
| `r_binding_dir` | `src/bindings/r` | **output**: the R C `.Call` shims (fpm bundles these into `libtensor-omics.so`) |
| `python_out_dir` | `python/tensor_omics` | **output**: the Python package |
| `r_out_dir` | `r/tensor_omics` | **output**: the R wrappers + loader |
| `snippets_dir` | `snippets` | **output**: the VS Code snippets (six files, git-ignored) |

So a default run writes:

```
src/bindings/c/<module>_c.F90        # Fortran C wrappers, one per module
src/bindings/r/                      # R C .Call shims -- fpm compiles them into the .so
    tox_marshal.h  init.c  <module>.c
python/tensor_omics/
    __init__.py  library.py  error_handling.py
    <module>.py                       # one per module
r/tensor_omics/
    tox_validate.R   error_handling.R   <module>.R
```

Two things worth knowing:

- **Each target directory is cleaned before writing** (unless `--no-clean`), so a procedure
  that stops being exported leaves no stale wrapper. For R only the generated `src/` and
  `R/` subdirectories are cleaned — a hand-written `DESCRIPTION` or `NAMESPACE` next to them
  is left alone.
- **`--library`** (default `build/libtensor-omics.so`) is not an output path; it is where the
  *generated Python loader* will look for the compiled shared library at runtime. Override it
  for an installed or relocated build, or set `TENSOR_OMICS_LIBRARY` at run time.

The conventions the generator recognises in the sources (prefixes, suffixes, the
`category` tag) are also in `config.py`, as `Conventions` — one place, no naming literals
scattered through the code. They are the source-language contract, documented under
[Writing generator-compliant Fortran](#writing-generator-compliant-fortran).

---

## How it is built

A strict one-way pipeline. Nothing downstream imports Ford; nothing upstream knows a target
language exists.

```
frontend  →  ir  →  abi  →  emit
```

| Stage | Package | Responsibility |
|---|---|---|
| frontend | `frontend/` | the only code that imports Ford; turns its parse tree into the IR |
| ir | `ir/` | a language-neutral model: types, docs, directives, entities, roles, errors, validation |
| abi | `abi/` | one decision layer for "what does this procedure look like from C" |
| emit | `emit/` | one emitter per target, each rendering the ABI model |

The keystone is that the **IR is constructible without Ford**. A test builds a `Procedure`
directly, so every rule is unit-tested without parsing anything — which is why the suite is
large and fast. The Ford frontend is just one producer of IR; the test fixtures are another.

`generate.py` wires the stages; `cli.py` is the command-line shell over it.

See [`design/c-layer.md`](design/c-layer.md) and
[`design/language-layers.md`](design/language-layers.md) for why each layer decides what it
does.

---

## Writing generator-compliant Fortran

This is the dev guide: what to write so a procedure is wrapped correctly. The generator
reads conventions from names and from documentation macros; nothing here requires editing
the generator.

### Marking a procedure for export

The `M_EXPORT_C` macro, in the procedure's Ford pre-comment (needs
`#include <src/macros.h>`):

```fortran
!> M_EXPORT_C
!| summary: Normalizes a vector to unit length in-place
!| author: A Developer
pure subroutine normalize_unit_length(vector, n_dims, ierr)
```

`M_EXPORT_C` expands to a Ford `category` meta tag, and the generator reads the category
value from the *same* macro — so the marker and what the generator recognises cannot drift.
Change the macro in one place and both follow.

Only tagged procedures are wrapped. Everything else is held to none of the rules below.

### What every exported procedure needs

- **Explicit `intent`** on every dummy argument. It decides constness in C and R, and
  whether an argument is an input, an output, or both.
- **An error argument `ierr`**, `integer, intent(out)`. If a procedure has none, the wrapper
  synthesises one — the binding languages always raise on error, so there must be a
  channel. Give it one if it can fail.
- **Kinded numeric types** (`real(real64)`, `integer(int32)`). A default kind has no
  defensible C mapping.
- **Documentation** — a `summary:` meta tag (becomes the docstring) and a `!!` comment on
  each argument (inherited by all three targets).

### Naming conventions

| Pattern | Meaning |
|---|---|
| `<p>_alloc` | the allocating variant; its wrapper is `<p>_c`. Its non-alloc twin `<p>` becomes `<p>_expert_c` |
| `tmp_<name>` | a work array: allocated by the binding language, never returned, never asked for |
| `<arg>_shape` | carries the shape of a flat `<arg>` passed separately (rank-independent serialization) |
| `n_selected_<arg>` | the count of an `<arg>_mask` / `<arg>_selection_mask`; computed from the mask |
| `mode`, `method`, `*_mode`, `*_method` | a mode argument (see below) |

Extents (`n_dims` in `vector(n_dims)`) are recognised automatically and never asked of the
caller — the binding language takes them from the array.

### Mode arguments

An integer argument compared against `MODE_*` / `METHOD_*` parameters. The binding
languages pass a **string**; the C wrapper maps it back. Document the accepted values in a
markdown table — the argument's own name decides the prefix (`mode` → `MODE_`, `method` →
`METHOD_`), and the header must agree:

```fortran
integer(int32), intent(in) :: link_method
    !! how to link clusters
    !!
    !! | Method | Value |
    !! |--------|-------|
    !! | minimises variance | [[tox_clustering(module):METHOD_WARD(variable)]] |
    !! | nearest neighbour   | [[tox_clustering(module):METHOD_SINGLE(variable)]] |
```

The string is the parameter name without its prefix, lower-cased: `METHOD_WARD` → `"ward"`.

### Documentation macros (`DM_`)

Written inside `!!` / `!|` comments, these carry what a signature cannot. They expand to
prose in the rendered Ford docs, so they read naturally; the generator recognises the
expansion. **Requires `#include <src/macros.h>` at the top of the file** — a `DM_` that
never expands is an error (as is a misspelt `M_`/`CM_`/`DM_` anywhere in a doc comment).

| Macro | On | Meaning |
|---|---|---|
| `DM_DEFAULT(VALUE)` | an optional | the value used when omitted. A *constant expression* — evaluated at generation time |
| `DM_REQUIRED_IF_MODE(MODE_ARG, MODULE, MODE_PARAM)` | an optional | required only in one mode; nullable otherwise |
| `DM_OPTIONAL_OUTPUT` | an `intent(out)` | the caller may decline it |
| `DM_RESULT_SIZE_IS(ARG)` | a result array | `ARG` holds how many leading elements are filled; the rest is trimmed |
| `DM_OUTPUT_FROM(OUT, PROC, MODULE, AUTO)` | an input | obtained by calling `PROC`; the caller never supplies it |
| `DM_OUTPUT_FROM(OUT, PROC, MODULE, JUST_INFO)` | an input | the caller supplies it; the doc says where to get it |

An optional **with** a `DM_DEFAULT` is required in C — the binding languages know the
default and pass it, which keeps the wrapper flat. Only an optional with no default is
nullable.

`DM_OUTPUT_FROM(..., AUTO)` matches the producer's inputs to the consumer's arguments by
name. `mask_chunk_count(n_genes, count)` is called wherever the consumer also has `n_genes`.
The producer must be exported, so that there is a wrapper to call; it may live in any
module. Python imports it inside the calling function rather than at the top of the file,
because two modules may size each other's outputs and a module-level import would then be
circular. R needs no import: every wrapper is sourced into one environment.

Where the producer and consumer spell the same quantity differently, the consumer argument
maps them in a table:

```fortran
integer(int32), intent(in) :: n_work
    !! size of the work array.
    !! DM_OUTPUT_FROM(n_work, fx_work_size, fx_edges, AUTO)
    !!
    !! | Producer input | Supplied by |
    !! |----------------|-------------|
    !! | n_values       | n_samples   |
```

Name-matching is tried first, so only the differences need a row. A producer input that is
neither name-matched nor in the table is an error naming it.

### Serialized arrays

A flat array whose shape travels separately, so any rank can be passed through one
signature:

```fortran
real(real64), dimension(:), intent(in) :: data
    !! the values, flat
integer(int32), dimension(:), intent(in) :: data_shape
    !! the extents of `data`, one per dimension
```

`data_shape` must be `intent(in)`, a rank-1 integer, and **not optional** (the wrapper reads
it to size `data`). `data` may be assumed-shape `(:)` or explicit `(n)`; the wrapper makes it
assumed-size and slices it to `product(data_shape)`.

Python then accepts an array of *any* rank and flattens it in Fortran order at the call, so
the shape argument never has to be written out by hand:

```python
serialize_int_helper(np.arange(12).reshape(3, 4), path)   # shape derived
deserialize_int_helper(shape, path)                       # count = product(shape)
```

Characters work the same way: the shape is read off the caller's array before the encode
rebuilds the buffer, and the strings are encoded in Fortran order.

### What is rejected, and why

The generator refuses what it cannot wrap correctly (only for exported procedures):

- a **deferred-length** character (`len=:`) — implies allocatable/pointer, not interoperable
- a character array of **rank > 1** — `tox_conversions` converts up to a vector
- an **allocatable** or **pointer** dummy — no C interoperable form
- a **derived type** — not mapped yet
- a **missing intent**
- a **`tmp_` argument that is `intent(in)`** — a work array is an output or in-out
- a **kind with no C mapping** — an error, never a guess (a wrong guess compiles and lies)
- an **optional shape or extent argument** — the wrapper must read it before `c_loc`

---

## Edge cases handled

Beyond the conventions above, the cases that took care to get right. Several exist *because*
the compiler or a run rejected the first attempt.

- **Null validation order** (`ierr`, then scalars, then arrays) so `c_loc` is never taken of
  a zero-size target; an empty array passes through to the callee's own validation. See
  [`design/c-layer.md`](design/c-layer.md).
- **Nullable optionals** cross as `bind(C)` `OPTIONAL` (interoperable since F2018), a plain
  pointer that is absent when C passes null — no descriptor, no branch.
- **`intent(inout)`** is modified in place in Python but copied-and-returned in R, because R
  is copy-on-modify. The R side clones before Fortran touches the data.
- **Characters** carry their length as a leading C extent. Output buffers are zero-filled
  (Fortran fills them only partially; the nulls terminate the strings), numeric outputs are
  not (wasteful and dishonest).
- **`c_bool`** logicals cross as real booleans; a dummy already `logical(c_bool)` needs no
  copy.
- **Shape cross-checks** the binding language makes but Fortran cannot: two arguments
  sharing an extent must agree, or a wrong answer / segfault results.
- **Argument-named errors**: `ierr` encodes the offending argument's position, so an error
  names it (`invalid input arguments (argument 'n_dims')`).
- **`M_IMPLICIT_NONE`** (`implicit none (type, external)`) in generated modules, so a typo
  in a called helper is a compile error, not a link-time surprise.
- **NA** (R) is rejected only where the check is free — integers, logicals, characters —
  never doubles, whose `NA` is a NaN Fortran already catches.

The full rationale for each is in [`design/language-layers.md`](design/language-layers.md)
and [`design/c-layer.md`](design/c-layer.md), including the alternatives that were rejected.

---

## How it is tested

Three layers, fastest first:

1. **Unit tests** on the IR, built from hand-made fixtures with no Ford — most of the suite.
2. **Frontend tests** parsing small fixture modules through the real Ford, plus a parse of
   the entire real `src/` that must come out with no diagnostics.
3. **End-to-end tests** that generate from the fixtures, **compile** the output
   (`gfortran -std=f2018`), build a shared library, and **call it** — from generated Python
   and from generated R. This is the only thing that proves the output is not just
   plausible but correct: `fx_sum_matrix` really sums a matrix.

The fixture modules in [`tests/fixtures/`](tests/fixtures/) are the generator's
specification — a deliberate, complete set of the constructs it supports, held to a
no-diagnostics bar. See their [README](tests/fixtures/README.md).

Run the suite:

```sh
python -m pytest helper/codegen/tests -q
```

The end-to-end tests skip cleanly without `gfortran` / `R`.

---

## Extending it

- **A new target language** is a new module in `emit/` consuming the `abi` model. It never
  touches the frontend or IR. The three existing emitters are the template.
- **A new convention or `DM_` macro**: add the macro to `src/macros.h`, its recognition to
  `ir/directives.py` (or a naming rule to `config.py`), its meaning to `ir/roles.py`, and a
  validation rule if one applies. Nothing about a target language changes.
- **A new type mapping**: extend the tables in `abi/c_abi.py` and the emitters. An unmapped
  kind is already an error, so a gap is loud.

Every convention lives in `config.py`; there are no naming literals scattered through the
code.

---

## Open items

- **`DM_OUTPUT_FROM(AUTO)` with a producer input that is a *constant*** rather than another
  argument. `tox_loess_required_workspace(n_dim, max_neighborhood_size, save_factorization)`
  is the real case: a consumer can now supply `max_neighborhood_size` by renaming, but
  `n_dim` is always `1` and `save_factorization` always `.false.`, and neither is an
  argument of the consumer to point at. Errors clearly until implemented.
- **ifx**: the F2018 features used (`OPTIONAL` in `bind(C)`, `implicit none (type,
  external)`) are verified with gfortran only. ifx is expected to agree; worth a check.
- **Compile check in CI**: the end-to-end tests need a compiler. They are marked to skip
  without one; wiring them into CI is a configuration step.
