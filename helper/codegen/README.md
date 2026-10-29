# TensorOmics binding generator

Generates the C, Python and R bindings to the Fortran library from the Fortran sources
themselves. A procedure marked for export becomes a C-callable wrapper, a Python function,
and an R function — with the documentation, error handling and argument validation derived
from the Fortran, not restated by hand.

This document covers the generator itself: what it produces, how to run it, how it is built
and tested, and what it refuses. **What to write in the Fortran** is
[`codegen_guide.md`](../../codegen_guide.md) at the repository root; the design rationale (why
each choice, and what was rejected) is in [`design/`](design/). This points at both rather than
repeating them.

---

## Contents

- [What it produces](#what-it-produces)
- [Running it](#running-it)
- [Configuration: what is read, and where output goes](#configuration-what-is-read-and-where-output-goes)
- [How it is built](#how-it-is-built)
- [The source contract](#the-source-contract)
  - [What the generator makes of it](#what-the-generator-makes-of-it)
  - [What is rejected, and why](#what-is-rejected-and-why)
  - [Documentation that names something](#documentation-that-names-something)
- [Edge cases handled](#edge-cases-handled)
- [How it is tested](#how-it-is-tested)
- [Extending it](#extending-it)
- [Open items](#open-items)

---

## What it produces

From one exported Fortran procedure, three things:

| Target | Output | What it is |
|---|---|---|
| C | `src/generated/bindings/c/<module>_c.F90` | a `bind(C)` wrapper: a plain-pointer ABI, null validation, type conversion |
| Python | `python/tensor_omics/<module>.py` | a `ctypes` function with a numpydoc docstring |
| R | `src/generated/bindings/r/<module>.c` + `r/tensor_omics/<module>.R` | a C `.Call` shim (marshalling, bundled into the `.so`) under an R function (validation, docs) |
| Snippets | `snippets/<Language>_<root>_snippets.json` | VS Code call/setup snippets, split by language and module root |

Plus, once per project: an error module for each language, generated from `tox_errors`, and
the loader/marshalling scaffolding.

Most of the exported procedures are themselves generated. A hand-written implementation
module produces the Fortran wrappers that then feed the table above — `src/generated/tox/`
is that output, and everything downstream of the synthesis treats those wrappers as ordinary
sources. See [How it is built](#how-it-is-built).

The `snippets` target (`emit/vscode_snippets.py`) emits VS Code snippets split into six
files -- `{Fortran,Python,R}_{f42,tox}_snippets.json` under `snippets/` -- the language in
the file name (so no per-snippet `scope`) and the root keeping the two namespaces apart.
These are **regenerated artifacts** -- never edited by hand, and committed to the tree like
the Python/R packages; the hand-written `snippets/toxdev_snippets.json` sits alongside them. Per exported procedure: a native Fortran `call` (plus a variant that guards
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
python helper/generate_code.py            # the wrappers and the bindings, into the tree
python helper/generate_code.py --check    # report problems, write nothing, non-zero if stale
python helper/generate_code.py --target python   # one target only
python helper/generate_code.py --help
```

Run from the repository root. `--check` is what a CI guard wants: it exits non-zero if the
committed bindings no longer match the sources -- a file that would change, or one the clean
step would remove and not write back, such as the wrapper of a deleted implementation
(`generate.files_to_remove`, which reads the same per-target list `_clean` deletes from).
Python's `__pycache__` is not counted: it is there after any test run and is no output.

Requirements: Python 3.11+ (`contextlib.chdir`), and `FORD`, `pcpp`, `numpy` (see
`requirements-dev.txt`). The generator reads the Ford settings from `fpm.toml`.

### Build switches the output honours

The generator always emits all of it; these decide what the *build* keeps. Pass them through
`./build.sh --directive=NAME`, which forwards to the preprocessor.

| Directive | Drops | Notes |
|---|---|---|
| `NO_C_BINDING` | the `bind(C)` wrappers, and the R shims with them | they call the `bind(C)` symbols. It does not disarm the kind guards: `f42_safeguard` compiles unconditionally, and the hand-written sources still name C kinds without it |
| `NO_R_BINDING` | the R `.Call` shims alone | keeps the C ABI for Python and direct C use; the build auto-disables this one with a warning when R is not installed |
| `NO_INPUT_VALIDATION` | the generated wrappers' input checks | see below |

`NO_INPUT_VALIDATION` is for a caller who has already established that the inputs are good —
an inner loop over data it produced itself. It gives up every diagnostic the framework offers,
so it is a whole-build decision rather than one to take per call site. What survives it:

- **`call set_ok(ierr)`**, which is not a check. It is what leaves `ierr` defined on the path
  where nothing goes wrong, and an implementation's own runtime errors are still reported
  through it.
- **the C layer's null checks**, which guard against a segfault rather than a bad value.
- **Python's and R's own validation**, which is a runtime layer and cannot be preprocessed out.
  A caller who wants the checks gone calls through C or Fortran.

A wrapper with nothing to check emits no guard at all, rather than an empty `#ifndef`/`#endif`.
The directives are written at column 0 whatever the surrounding indent — gfortran's preprocessor
rejects an indented one outright ("Invalid character in name"), and a subroutine body is
rendered on its own writer before being nested into its module's, so `render/writer.py` keeps
them there through that nesting.

Diagnostics point at the offending line of the *original* source, with the entity chain and
a note on what to do. A generated wrapper inherits the implementation's location, so a rule
broken in code the author never wrote still names the line they did — the chain then reads
the generated module while the file is the implementation's:

```
error: shape argument 'data_shape' is optional
  --> src/tox/tox_x_impl.F90:90
  argument 'data_shape' in procedure 'x_expert' in module 'tox_x'
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
| `generated_dir` | `src/generated` | **output**: everything written back into the source tree, the wrapper modules included |
| `c_binding_dir` | `src/generated/bindings/c` | **output**: the Fortran C wrappers |
| `r_binding_dir` | `src/generated/bindings/r` | **output**: the R C `.Call` shims (fpm bundles these into `libtensor_omics.so`) |
| `python_out_dir` | `python/tensor_omics` | **output**: the Python package |
| `r_out_dir` | `r/tensor_omics` | **output**: the R wrappers + loader |
| `snippets_dir` | `snippets` | **output**: the VS Code snippets (six files, only on `--target snippets`) |

So a default run writes:

```
src/generated/<rest>/<module>.F90              # the wrapper modules, mirroring src/<rest>/
src/generated/bindings/c/<module>_c.F90        # Fortran C wrappers, one per module
src/generated/bindings/r/                      # R C .Call shims -- fpm compiles them into the .so
    tox_marshal.h  init.c  <module>.c
python/tensor_omics/
    __init__.py  library.py  error_handling.py   # __init__: every procedure, code and exception class
    _contract.py                      # the input contract every wrapper checks by
    <module>.py                       # one per module
r/tensor_omics/
    tox_validate.R   error_handling.R   <module>.R
```

Three things worth knowing:

- **There is no setting for where a wrapper module goes**, because it is a mirror rather than
  a destination: an implementation at `src/<rest>/<module>_impl.F90` generates
  `src/generated/<rest>/<module>.F90`, sub-directories and all. Nothing in the rule knows
  about `tox` or `f42`, which is what lets any layer acquire implementations without the
  generator learning a second case — and what makes `_impl` a suffix reserved across the whole
  of `src/`, since a module so named acquires wrappers wherever it sits.
  `synthesize.generated_path_for` is the single place that decides it, and the emitter goes
  through it. The cleaner and the Ford frontend need no path at all: the cleaner removes every
  wrapper module under `generated_dir` that it finds on disk, and Ford excludes `generated_dir`
  as a whole — rather than three literals kept in step by hand, which is how one of them comes
  to write, delete and exclude three different paths.
- **Each target directory is cleaned before writing** (unless `--no-clean`), so a procedure
  that stops being exported leaves no stale wrapper. In the R *package* directory only the
  generated `*.R` are removed — a hand-written `DESCRIPTION` or `NAMESPACE` next to them is
  left alone. The wrapper modules are the one target cleaned file by file rather than by
  directory, because they share `src/generated` with the binding trees: every `.F90` under
  `generated_dir` outside them, read off the disk (`synthesize.generated_wrapper_paths`).
  The list used to come from the implementations that exist, so the wrapper of a *deleted*
  one survived — and fpm compiles all of `src/generated`, so that orphan, still `use`-ing the
  module that went with it, broke every build after, regenerated or not (#221). Reading the
  disk is safe because nothing hand-written lives under `generated_dir`; a missing source tree
  (a mistyped `--src`) owns nothing. `--check` reports such an orphan as a file that would be
  removed, so a committed one fails the CI guard rather than only a build that skips
  regenerating. The snippets target is not cleaned at all: `snippets/`
  holds the hand-written `toxdev_snippets.json` beside the generated files, each of which is
  rewritten whole.
- **`--library`** (default `build/libtensor_omics.so`) is not an output path; it is where the
  *generated Python loader* will look for the compiled shared library at runtime. Override it
  for an installed or relocated build, or set `TENSOR_OMICS_LIBRARY` at run time.

The conventions the generator recognises in the sources (prefixes, suffixes, the
`category` tag) are also in `config.py`, as `Conventions` — one place, no naming literals
scattered through the code. They are the source-language contract, which
[`codegen_guide.md`](../../codegen_guide.md) documents from the author's side.

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

`synthesize.py` is a third, and sits between the frontend and the semantic pass: it reads every
`_impl` module, builds the wrappers each implementation implies, and injects them into the
project before roles are analysed. **The module name is the whole trigger** — a module called
`<x>_impl`, wherever under `src/` it lives — so everything downstream, validation and all three
emitters included, sees the wrappers as ordinary procedures read from source and needs to know
nothing about implementations. Where the file sits decides only where they are written.

The frontend excludes that output from the next parse *by directory* rather than by file name,
which is what it used to be. Ford matches a bare exclusion as `**/<name>`, so once a layer
generates an `f42_stats.F90`, excluding that name would also drop the hand-written
`src/f42/utils/f42_stats.F90` and every binding built from it, with no error at all. One
directory also covers every target at once, a wrapper left behind by a deleted implementation
included — which a list derived from the implementations that still exist would not name.

`generate.py` wires the stages; `cli.py` is the command-line shell over it.

See [`design/impl-layer.md`](design/impl-layer.md) for why the wrappers are synthesised at all,
and [`design/c-layer.md`](design/c-layer.md) and
[`design/language-layers.md`](design/language-layers.md) for why each layer decides what it
does.

---

## The source contract

**[`codegen_guide.md`](../../codegen_guide.md), at the repository root, is the contract of
record** — what to write so a procedure is wrapped correctly, case by case, with a worked
example for each and a real snippet from the current tree. It covers both ways in: the
implementation path, where you write one annotated `_impl` procedure and the generator writes
the wrappers, and the export path, where you write the whole procedure and mark it
`M_EXPORT_C`.

| If you want | Go to |
|---|---|
| the whole contract, task by task | the guide |
| valid ranges, finiteness, masks, distance matrices, optionals | guide Sections 5.1–5.6 |
| work arrays, permutations, recommend sizing | guide Sections 5.7–5.9 |
| modes, per-mode procedures | guide Sections 5.10–5.12 |
| prologues, runtime errors, split families | guide Sections 5.13–5.15 |
| exporting a hand-written procedure | guide Section 6 |
| every `DM_` macro, with its contract in a comment | [`src/macros.h`](../../src/macros.h) |

None of that is repeated here. Two copies of one contract drift, and these two had; what stays
in this document is what is about the *generator* rather than about the sources — what it makes
of the contract, and what it refuses.

### What the generator makes of it

The parts a source author does not need, and a generator maintainer does:

- **`M_EXPORT_C` expands to a Ford `category` meta tag**, and the generator reads the category
  value from that same macro (`frontend.export_category`) rather than hardcoding the string. The
  marker and what the generator looks for therefore cannot disagree — change the macro and both
  the sources and the generator follow. Only tagged procedures are wrapped; everything else is
  held to none of the rules.
- **The plain name goes to the entry point a caller should reach for first**, which is why the
  suffix is never the author's to choose (`synthesize._wrappers_for`). Where the implementation
  has something to take over — a `tmp_` work array, a `<base>_perm` permutation, a
  `DM_OUTPUT_FROM(..., AUTO)` value — or where its prologue asks for an argument of its own,
  two wrappers are generated: `foo`, which prepares all of that and calls the implementation,
  and `foo_expert`, which validates and calls it with exactly what you supply, preparing
  nothing. Where there is nothing to take over there is **one** wrapper and it is called `foo`,
  there being no second tier for an `_expert` name to distinguish it from — nearly half the
  wrappers in the tree today. Naming the lone one `foo_expert` unconditionally would generate,
  compile and bind perfectly well, and quietly rename every entry point that has no work arrays.
- **The generated module carries the implementation module's own documentation**, verbatim,
  and its `Meta` with it (`synthesize`). It is the published API -- what Python imports, what
  the R help pages are built from, what a Fortran caller `use`s -- so an author's prose about
  the family is what a reader of any of those gets, and `Doc.summary` is their own first line.
  It was a two-line "do not edit" banner until 2026-08-12, which left every published module
  described by the same sentence; the argument for it (an implementation module's doc describes
  the implementation) was answered by rewriting eighteen of those docs as API prose instead.
  The "generated from" note is added by each **emitter**, not carried in the IR: it is about
  one file, and carried in the doc, Python printed it beside its own trailer.
- **A name is published as it stands.** Synthesis names the tiers `foo` and `foo_expert`
  outright, so they cross as `foo_c` and `foo_expert_c` with nothing to decide, and a
  hand-written export is published as written. `abi/c_abi.stripped_name` is therefore the
  identity — it was not until 2026-08-12, when it translated a hand-written `foo` / `foo_alloc`
  pair into the published `foo` / `foo_expert`. f42 wrote the last of those, and the rule went
  with them: renaming a procedure behind its author's back needs a live customer.
- **A dummy is referred to by its declared spelling.** Fortran is case-insensitive, so
  `integer :: N` sizes `x(n)`; Python and the R shim are not, and saw two variables, one never
  defined. The frontend rewrites every reference to a dummy -- in extents, character lengths,
  `DM_RESULT_SIZE_IS`, `DM_REQUIRED_IF_MODE` and the `DM_MIN`/`DM_MAX`/`DM_SENTINEL`
  expressions -- to the declared spelling once (`ford_frontend._use_declared_spellings`), so
  every emitter can compare names exactly. The declared names are never changed: they are the
  published API, mixed case included (`neighbor_mask_S1`, `omics_field_RAP_projection`).
- **A default is the value a Fortran caller gets** (`ir/constants.py`). The bindings pass the
  evaluated `DM_DEFAULT`, so it is evaluated with Fortran's arithmetic, not Python's: every
  integer carries its kind (`constants.Integer`; a literal without one is default, kind 4) and
  overflowing that kind is an error, as it is for the compiler; integer `/` and `mod` truncate
  toward zero without going through a float (`9007199254740993_int64 / 1` is exact); an integer
  to a negative power is an integer (`2**-1` is 0); `nint`/`anint` round a half away from zero,
  and not as `floor(x + 0.5)`, whose addition rounds `0.49999999999999994` up; a negative real
  to a real power, a complex result and an infinite real are errors. Module parameters keep the
  kind they are declared with (`Project.constant_values`), and a kinded logical
  (`.true._c_bool`, as the guide asks for) evaluates like a plain one. A sign binds below `**`
  and `*`, as in Fortran: `-2.0_real64**2` is -4. The one rule not followed is a real's kind:
  every real literal is read as real64, so an unkinded `0.95`, which the compiler stores as the
  real32 0.949999988, is passed as 0.95 -- write `_real64` and the two agree.
- **A mode crosses as a string.** The binding languages pass the parameter name without its
  prefix, lower-cased (`METHOD_WARD` → `"ward"`), and the C wrapper maps it back to the integer,
  rejecting an unknown one before Fortran is entered. The Fortran wrapper separately checks
  membership against exactly the values the mode table names.
- **`DM_OUTPUT_FROM(..., AUTO)` is called by the language layer**, so the producer must be
  exported — there has to be a wrapper to call. Python imports it *inside* the calling function
  rather than at the top of the file, because two modules may size each other's outputs and a
  module-level import would then be circular. R needs no import: every wrapper is sourced into
  one environment.
- **A refusal of the binding's own is an input error without a code** (#221). Python raises
  `ToxArgumentError` (`emit.errors_python.ARGUMENT_EXCEPTION`), a `ToxInputError` whose `code`
  is None and which is a `TypeError` and a `ValueError` too, for every shape, type and value it
  refuses before the call; R's `tox_type_error`, `tox_shape_error` and `tox_na_error` inherit
  `tox_input_error` (`.tox_argument_error`), `code` NA. They used to be a bare `ValueError` /
  `TypeError` and a plain `tox_error`, so the same wrong shape needed a different handler when
  the binding caught it than when the Fortran did. Each docstring's `Raises` names it.
- **Every input is checked by one input contract, the same in Python and R** (#221,
  `design/language-layers.md`, "what an argument accepts"). The input's own kind first -- an
  allowlist on numpy's `dtype.kind`, R's `is.numeric` / `is.logical` / `is.character`, never
  `typeof` -- then a conversion accepted only if converting back gives exactly the same values:
  an integer takes whole floats in range, a logical takes numbers that are 0 or 1, a real takes
  any integer or float that float64 holds exactly (`2**53 + 1` is refused, a 64-bit integer is
  compared only beyond 2^53 and never cast back from past its range, a Python int too large for
  numpy by its own arithmetic, the integers of a list numpy would make floats of as they were
  given), and nothing parses a string, counts a datetime or reads a factor's
  codes.
  A string is a `str` (or an `os.PathLike`) or one R string, valid Unicode without a NUL, at most
  its declared length; an array of strings is a sequence of them (an R factor by its labels); a
  mode matches whole after its trailing blanks, never cut to the buffer; a rank is the declared one, nothing is flattened
  (an R 1-D array and table are vectors, a numeric data.frame a matrix). Python's is one generated
  module every wrapper calls, `tensor_omics/_contract.py` (`emit.python_contract`), free for an
  input of the Fortran type. R checks kind, rank and a single value in `tox_validate.R`; every
  element of an integer, logical or character array is checked by the C shim as it converts it
  (`tox_int_in`, `tox_bool_in`, `tox_char_in` in `tox_marshal.h`), which raises the R wrapper's own
  condition from C (`tox_refuse`, from `emit.errors_r`'s classes). An R string crosses by what R
  marks it as, the same in every locale (`tox_utf8`): marked UTF-8 as it stands, marked latin1
  translated to UTF-8, unmarked as it stands if its bytes are UTF-8; one marked "bytes", one that
  is not well-formed and a latin1 byte with no character are refused by name, so the escape text
  R makes up for what it cannot convert (`<e9>`) never reaches Fortran. A returned string is
  marked UTF-8. Trailing blanks are Fortran semantics,
  accepted and dropped, and every generated string parameter says so (`emit.doc_contract`).
- **A single complex is two doubles wide in Python** (#221 WP8). `ctypes` has no complex type,
  and the box was the real part's `c_double`: Fortran wrote 16 bytes into 8 for an output and
  read 16 from it for an input. The box is `(c_double * 2)`, passed as the `c_double*` the argtype
  names (`python_ctypes._is_complex_box`): allocated for an output; for an input and an
  `intent(inout)`, filled from the value `_contract.complex_number` accepted (a real number by the
  real rule, a complex one whose parts float64 holds); an output and an inout one come back as
  `complex(z[0], z[1])`. An array is numpy's `complex128`, and C and R had it right already.
- **Every input axis is compared before the call, by one rule** (#221 A-shape,
  `emit.shape_checks`): an axis that names an extent is compared with that extent's value, except
  the axis the value was read from, and an axis the declaration fixes at a constant (`fields(n, 2,
  m)`, `points(N_AXES, n)`) with the constant. The value is read off the first input that must be
  passed (`shape_checks.extent_provider`, which Python's derivation uses too, so the axis left out
  is the one read), counted from a mask, computed by a producer, or passed by the caller; an
  optional array is compared when given. Before, only an extent two inputs shared was compared, and
  optional arrays never: `build_residual_histograms` read past a 1x1 `neighbor_mask`, an input
  sized `n_selected_<x>` was compared with nothing, and `omics_field_RAP_projection` read a field
  without its target. A message names the array, the axis, where the value came from and the
  declared shape: `'w' has 1 along axis 0, but 'x' implies n == 3; 'w' is declared with shape (n)`.
- **An input axis the declaration computes is given to C and compared there** (#221 A-shape). Only
  Fortran can evaluate `expr(sum(reps_per_tissue), n_genes)` or `velocity(max(0_int32,
  n_timepoints - 1), ...)`, so C is given the axis's actual extent as an argument right after the
  array (`n_expr_elements_dim_1`, `CArgument.checked_against`), which Python and the R shim read
  off the array like an assumed-shape extent. The C wrapper declares the array with it and
  null-checks with it -- the declaration `sum(reps_per_tissue)` dereferenced `reps_per_tissue`
  before its own null check, so a null one crashed the caller -- and, after every null check,
  compares it with the declaration: `ERR_SIZE_MISMATCH` at the array's position. Not under
  `NO_INPUT_VALIDATION`, which drops checks of values; this one keeps a short buffer from being
  read past, as the null checks keep a null one from being dereferenced. An array whose shape
  travels in a shape argument and whose declaration counts it too (`arr(n_strings)` beside
  `arr_shape`) gets the same comparison of `product(arr_shape)`, no argument added: a C caller's
  inconsistent count had the callee write past the view. Which arrays get which is the rule's to
  decide, not a list: on main at #221 it came to three axes in three symbols (`calc_tiss_avg_c`,
  `compute_acceleration_from_velocity[_trajectory]_c`) and the ten serde helpers' counts. The
  binding passes the declared shape of every array the C wrapper compares to `check_err_code`
  (`error_arguments.shapes_of`), so the error says what the array must be: `Array size mismatch
  (argument 'expr'); 'expr' is declared with shape (sum(reps_per_tissue), n_genes)`. Docs keep showing the declared extents (`CArgument.declared_extents`), a
  module constant by its name.
- **A function is published as a subroutine, its result right before `ierr`** (#221 B1). Both
  Fortran tiers of a function implementation are subroutines taking the result as an
  `intent(out)` argument of the result's own name, before the implementation's own `ierr` where
  it declares one (`synthesize.published_arguments`), and the wrapper assigns the call to it:
  `total = total_impl(n = n, x = x)`. From synthesis on it is an ordinary argument -- validated,
  bound, documented and returned like any output -- and the implementation stays a function, so
  other implementations keep calling it as one. An exported function's C wrapper takes its
  result at the same place (`abi.c_abi`), which it used to append after an own `ierr`
  (`hw_ownierr_c(n, x, ierr, total)`). An `elemental` function is wrapped as an elemental
  subroutine is: the wrapper takes scalars.
- **`foo` ends in `ierr` even with a prologue.** A prologue's own arguments are `foo`'s
  vocabulary, so they go after the implementation's -- and before the implementation's own
  `ierr` where it declares one (`synthesize._wrappers_for`). They used to be appended after it,
  publishing `..., ierr, quantile_level` -- the shape `detect_outliers` once had.
- **The drop set of `foo`** -- what it prepares rather than asks for -- is
  `synthesize.taken_over_arguments`: `tmp_` work arrays, `<base>_perm` permutations,
  `DM_OUTPUT_FROM` values, and anything the prologue writes that the implementation only
  reads. The signature is shaped from it in `synthesize` and the body written from it in the
  emitter, so both call it with the same prologue; an argument that also sizes something the
  caller still sees is held back, or nothing could size what comes back. `foo_expert` keeps
  every one of them: that difference in signature is the whole reason for two entry points.
- **Every `use ..., only:` names what the body uses, and nothing more.** An unused `only:`
  name is not a compile error in any Fortran and no compiler warns about one, so a fixed
  import list goes stale in silence -- which is what happened to the C bindings' `tox_errors`
  list, twice: once as recorded here, and again when the mode conversion stopped being able to
  fail and left `is_err` behind. Each conditional import has a predicate beside the code that
  emits it (`_conversion_helpers`, `_emitted_literal_kinds`), with tests that a module needing
  none asks for none -- and those tests are what caught the second one. When auditing by hand,
  expand `src/macros.h` first: most apparent unused imports are named inside `M_CHECK_NON_NULL`.
- **A published pair explains itself** (`emit/doc_tiers.py`). Where both `foo` and
  `foo_expert` reach a language, each docstring says what the other does -- which permutation
  the allocating half seeds and sorts, which prologue it runs -- because the two are otherwise
  indistinguishable to a reader. The parts are stored with identifiers in backticks and
  rendered per language: Python keeps them, R turns them into `\code{}`. Where `foo` only
  validates and allocates, `foo_expert` does not reach a language at all
  (`SynthesisResult.expert_only_in_fortran`, applied once in `generate`): the binding allocates
  the work arrays and computes the `DM_OUTPUT_FROM` sizes for both tiers, so the two would be
  the same call under a name promising control over it. Fortran and C keep both, because there
  the expert tier really does hand the buffers over.
- **The package's top level is the whole Python API** (#221 C). `tensor_omics/__init__.py`
  re-exports every procedure and everything the error module defines: each `ERR_*` and `STAT_*`
  code, `ARG_POS_FACTOR`, every exception class (`ToxError`, one per group, `ToxArgumentError`)
  and `check_err_code`, all listed in `__all__` -- where R has had them all along, as globals of
  `error_handling.R`. Both the error module and the re-export are emitted from one list,
  `emit/errors_python.public_names`, with the error table `generate` already reads, so a code
  added to `tox_errors` is exported with no other edit, and `__init__` cannot name something the
  module lacks. It used to export `ToxError` and `check_err_code` only, by a hard-coded pair.
- **A serialized array is made assumed-size and sliced** to `product(<arg>_shape)` in the wrapper,
  which is why the shape argument may not be optional — the wrapper reads it before it may take
  `c_loc` of what it sizes. Python then accepts an array of any rank and flattens it in Fortran
  order at the call, so the shape is never written out by hand. Characters work the same way: the
  shape is read off the caller's array before the encode rebuilds the buffer.

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
- an **optional output** — no binding can honour it; use an optional input flag plus a `tmp_`
  work array
- an axis with a **lower bound other than 1** (`x(0:n)`, `x(0:)`) — every layer sizes an axis
  from its extent alone, so the validator's count, the C wrapper's null check and the bindings'
  shapes would all paste `0:n` where a number belongs. A written-out lower bound of 1 is
  dropped by the frontend (`x(1:n)` is `x(n)`)
- a **`DM_DEFAULT` that does not fit its argument's kind** (`3000000000_int64` on an
  `integer(int32)`), or that is not a constant Fortran would accept (an overflowing integer, a
  negative real raised to a real power) — the bindings would pass a value no Fortran caller can
- a **kind the bindings cannot carry** — anything but `integer(int32)`, `real(real64)`,
  `complex(real64)`, `logical(c_bool)` and a default `logical`. A numeric kind is compared by
  value (`validate.SUPPORTED_KIND_VALUES`, so `c_int` is `int32` and `integer(8)` is int64); a
  logical kind by name, since `logical(1)` -- c_bool's value -- would still be converted into a
  default logical and fail to compile. One error per cause: the C ABI stays quiet about a kind
  this refuses (`validate.bindings_cannot_carry`). R stores an integer as int32 and a
  number as real64, so the shim handed Fortran the wrong storage for every other kind (an `int64`
  output was written eight bytes per element into a four-byte `INTSXP`, a heap overflow), and
  the generated validators exist for int32 and real64 only. `KIND_MAP` still knows the C
  counterparts of the others; it is the bindings, not C, that cannot take them
- an **argument named like a word a binding language reserves** — Python's hard keywords and
  builtins (a frozen list, not read from the running interpreter), R's reserved words
  (case-sensitively: `NA` is refused, `na` is a name), and C's keywords as R's `gcc -std=gnu2x`
  takes them (`asm` included), and the functions the R shim's header defines (`tox_int_in`,
  `tox_char_in`, ... -- `names.R_SHIM_HELPERS`, which a shim parameter would hide); or named like
  a procedure the wrapper names (its own name, its C symbol, a `DM_OUTPUT_FROM` producer). Not refused: Fortran keywords (Fortran has none
  reserved, and main publishes `value`) and R function names (an argument cannot shadow a call in
  R, and main publishes `mode`, `vector`, `array`, `expression`, `factor`, `points`, `weights`)
- a **procedure named like a Python keyword or builtin, or an R reserved word** — it is a Python
  and an R function of that name. Its C symbols (`foo_c`, `foo_call`) carry a suffix, so C's
  keywords do not apply
- a **procedure or module named like something else at the Python package's top level** -- an
  export of the error module (`check_err_code`, `ERR_OK`, `ToxInputError`, `ARG_POS_FACTOR`), one
  of the package's own submodules (`library`, `error_handling`, `_contract`), or a module, the
  procedure's own included (`emit/python_ctypes.check_package_names`, run with the Python target).
  The package binds each name once, so one would hide the other: a procedure `check_err_code`
  would replace the decoder, one called `library` the submodule `tensor_omics.library`, and a
  module `error_handling` would be written over the error module's file. Compared ignoring case,
  as Fortran compares names, so a second spelling cannot reach the Fortran and R layers either.
  None of it happens on main or the open branches (#221 C6): procedures are snake_case, codes
  `ERR_`/`STAT_`, classes `Tox*`
- an **argument named like something its generated Fortran wrapper uses** — an intrinsic it
  calls (`count(mask, kind=int32)` beside an `n_selected_` count, `size(x, 1)` for an
  assumed-shape array), a name it imports (`set_err_once`, a validator, a kind), a kind it writes
  (`3_int32`, `.true._c_bool`), a local it declares (`handled` beside a prologue, `<x>_f` for a
  converted `x`, `<x>_address` for an array of strings `x`, `<x>_value` for a materialised
  default), or a variable of the module around it (`empty_view_anchor`, where an array of
  strings is viewed). Read off the wrappers as emitted, after
  their macros are expanded (`emit/name_clashes.py`), so it is refused exactly where it collides:
  `count` stays a good argument of `mask_chunk_count`, which counts nothing. The R shim and Python
  have no such clash to refuse -- everything they bind beside the arguments starts with `_`
  (`_v_k`, `_nprot`, `_np`, `_error`), which no Fortran name can. All of these names are built in
  [`names.py`](names.py), and `tests/test_names.py` reads the generated Python, R shim, R wrapper
  and both Fortran wrappers back to hold the emitters to it
- an **`ierr` that is not the last argument**, of a hand-written export or, through the
  synthesised tiers, of an implementation. Every published signature ends in the error code; the C
  ABI checks its own order too (`c_abi._Builder._check_error_is_last`), a function included, whose
  result goes before it
- a **function without a `result(...)` clause**, an implementation or an export -- its result is
  published under its own name, which is then the function's, and the wrapper that calls the
  function names it too: gfortran says "Symbol ... conflicts", ifx error #6401. A result **named
  `ierr`** (two dummies of one name), one **with the `tmp_` prefix** (a work array, never
  returned: the function's only output, thrown away) and one of **assumed length** (`len=*`, whose
  length comes from a declaration no binding writes) are refused as well
  (`validate.function_result_problems`); synthesis leaves such a result out of the wrappers, so
  nothing builds the signature the error is about. Otherwise a result is held to every rule an
  `intent(out)` argument is -- `allocatable`, `pointer`, `len=:`, a kind the bindings cannot
  carry, a reserved word -- and an undocumented one warns like any argument
- the **same exported name twice**, in two modules or in one (a mode table's third column can do
  it), compared ignoring case: it is the C symbol, the R `.Call` entry and the name Python and R
  publish, so the link fails with a duplicate symbol, or the second silently shadows the first
- an **`M_EXPORT_C` on a private procedure** (reported by the frontend, which is the only place
  that sees one: Ford leaves it out of the module's routines, so it used to publish nothing,
  silently)
- an **optional `intent(inout)`** — Python kept no result for it and R wrote into the caller's
  own vector, which breaks copy-on-modify
- an **`intent(inout)` character or logical array** — both cross through a buffer the binding
  builds, and nothing hands the modification back (the R shim did not even compile)
- a **`DM_DEFAULT` on an array** — a default is one value, so an omitted array arrived with one
  element and failed its shape check
- an **extent that is not `intent(in)`, or that only optional arrays carry** — the bindings read
  an extent off the arrays they are passed and drop it from the signature, so a written one is
  neither theirs to read nor to return, and one only an omitted array carries is read off nothing
- an **`n_selected_<x>` count and its mask that are not both `intent(in)`, or an optional mask
  with a count** — the count is derived from the mask before the call, in every layer
- a **`DM_RESULT_SIZE_IS` naming an `intent(in)` count** — the bindings read the count back after
  the call; a size the caller chooses sizes the result directly
- a **module constant sizing an argument that the generator cannot evaluate to an integer** — Python and R
  allocate and pad with its value, so the generator has to know it
- an **input axis computed from arguments that the C wrapper cannot evaluate** (`c_abi._can_check`):
  one reading an optional (absent, it has no value), an `intent(out)` dummy (the procedure's to
  set), a dummy the C wrapper converts (a string, a mode, a default-kind logical: under its name
  there it is not the value the declaration means), or an array whose own extent is computed (the
  two checks would have to be ordered; nothing on main needs it). C is given such an axis's extent
  and the C wrapper compares it with the declaration, see above
- an **extent the bindings compute outside the grammar they compute in** -- the extents of an
  output or a work array, the length of a string, and an input's constant axis, which Python and R
  allocate, pad or compare by before the call (`abi.extents.check_binding_extents`, on what
  reaches the languages). The
  grammar is the part of Fortran that means the same in Python and C: names, module constants,
  integer literals, `+ - * /`, parentheses, `max`, `min`, `size`. Two regex translators used to
  pass anything else through unchanged (#221 B11), and `/` floors in Python where Fortran
  truncates, so `(n - 10)/3` sized a buffer one element short. Refused too: such an extent reading
  an output or an optional argument, which has no value yet, or a name that is neither an
  argument nor a module constant. The note points at `DM_OUTPUT_FROM`

An implementation module is checked separately, because its procedures are read rather than
exported and so are reached by none of the above. These rules follow the `_impl` suffix
wherever it is written — under `src/tox`, under `src/f42`, anywhere — which is the price of
dropping the fixed directory they used to be scoped to:

- a module **not in a file named after it** — the *module* name is what synthesis triggers on
  and what the generated module is called; the *file* name is what its path is mirrored from.
  Let the two disagree and `module tox_x`
  is written into some unrelated file — or into none, since a module named `<x>_impl` in a
  file that is not an `_impl` file has no mirrored path at all and the emitter raises rather
  than guessing a directory
- an **`M_EXPORT_C` on an implementation** — its wrapper is the entry point, and the export
  publishes an unvalidated twin beside it, under a name (`foo_impl`) a binding caller cannot
  tell apart from the validated `foo`. Support routines in the same module (the
  `DM_OUTPUT_FROM` recommend routines) have no wrapper and stay exported, which is why an
  implementation module can itself show up in the bindings — `tox_loess_impl.py` is the
  Python home of `tox_loess_required_workspace`, and holds no implementation at all
- an **`intent(out)` argument the allocating wrapper takes over** — a `<base>_perm` whose
  `<base>` is an argument (it seeds and sorts it itself) or a `DM_OUTPUT_FROM(..., AUTO)` value
  (it calls the producer itself). Either way it would discard what the implementation computes;
  only the expert tier returned it. A `tmp_` work array stays exempt
- an implementation named **`<x>_expert_impl`** or **`<x>_alloc_impl`** — both name a wrapper,
  and the suffix is the generator's to add. `foo_expert_impl` beside `foo_impl` yields two
  procedures called `foo_expert`, and the emitter would strip the suffix and call `foo_impl`
  from the wrong one: wrong code that compiles, because `foo_impl` exists. `foo_alloc_impl`
  yields `foo_alloc`, which every author of this framework reads as the allocating tier while
  being an ordinary second procedure beside the generated `foo` that really is one
- an **`allocatable` local or dummy** in any procedure of the module, implementations and
  helpers alike — the generated wrapper owns the memory. Seen through the declaration, since no
  body is ever read; a `pointer` local aliases and is fine. The dummy is the subtler half: it
  looks like the caller's memory, but whoever fills it allocated it
- a **`use` naming anything but another `_impl` module or a member of
  `Conventions.impl_import_whitelist`**. This is what gives the rule above any reach beyond one
  file: the allocation check reads declarations, so an implementation calling a helper
  *elsewhere* that allocates would pass it, and only the import list can see that. It also
  keeps an implementation from reaching a generated wrapper — a layering inversion, and a
  module cycle within one family. A `use` inside a procedure counts the same as one in the
  module header (`Procedure.uses`), or the rule would sit one indentation level from being
  bypassed.

  The list is checked, not asserted: **every whitelisted module must itself be
  allocation-free**, and the generator verifies it. Without that the bound is only a promise —
  whitelist a module, let it grow an allocation, and every implementation reaching it allocates
  with nothing to say so. It went unchecked for a long time because it would have failed:
  `f42_utils` sat on the list while its `f42_stats` child still had hand-written `_alloc`
  halves, and `tox_conversions` until its string conversions became pointer views. Converting
  the f42 family removed the entry rather than justifying it, which is the move to reach for —
  a module that has to be whitelisted may be a module that should be an implementation

- a **C kind named without `use f42_safeguard`**. `iso_c_binding` gives a C kind the value -1
  where the platform has no counterpart, and a negative kind is a compile error at every
  declaration naming it. `f42_safeguard`'s guards turn that into one error saying *which* kind
  is missing — but only if it compiles first, and a dependency is the only way to order a
  Fortran build. Checked because it has already regressed once: before logicals moved to
  `c_bool` the only modules naming a C kind were the generated C wrappers, every one of which
  had the `use`, and the sweep put C kinds in 39 more compilation units without it. Reads
  argument kinds, which is what the IR carries, so a C kind used only on a module variable is
  still invisible. The generated wrappers are exempt — their `use` list is the emitter's, and
  `emit/fortran_wrapper` writes it.

A `DM_PROLOGUE` is checked too, having had no analysis or validation at all before:

- the named procedure must exist
- it must declare a scalar `logical, intent(out) :: handled` — the wrapper returns early on it
  regardless, so without it that branch reads an undefined value. A default `logical`: the wrapper
  declares its flag so (`emit/fortran_wrapper._prologue_local`) and passes it, and it never
  crosses to C. Any other kind is refused, naming the prologue and the argument: argument
  association does not convert kinds, so a `logical(c_bool)` one passed the generator and then
  failed to compile ("Type mismatch in argument 'handled'")
- a dummy naming an implementation argument is supplied from it; one naming nothing becomes an
  argument of `foo`, which is what the prologue derives *from*. A name **one edit** from an
  implementation argument is refused as a misspelling, since it would otherwise become a new
  argument and leave the prologue and the implementation working from different values
- a dummy that some mode's wrapper does not have — the mode argument, or one scoped to a mode —
  is refused: the prologue runs in all of them
- the implementation must generate an allocating wrapper for it to run in (it does when
  anything is taken over, or when the prologue asks for an argument of its own)
- it may not produce anything the allocations, permutation sorts or recommend calls above it read

At run time, the bindings refuse what the input contract does not accept, before anything is
called and naming the argument (`ToxArgumentError` / a `tox_input_error`, no code; the table is in
`design/language-layers.md`): a bool for a number, a fraction, NaN or out-of-range value for an
integer, an integer or a long double float64 would round for a real (`2**53 + 1`), a bool, a
string or an object for a single complex in Python, any number but 0 and 1 for a logical, a
string, object, datetime or complex for a number, an R factor, Date or `integer64` for a number, a
masked array with masked elements, a non-`str` for a string, a bare `str` for an array of strings,
a NUL, a lone surrogate or invalid UTF-8 (in R, an unmarked string whose bytes are not UTF-8 in
any locale, one marked "bytes", and a latin1 one holding a byte with no character), an NA, a
string longer than its declared length, an unknown or over-long mode, and any other rank. A
refusal that has a fix names it: `.view(np.int32)` for a uint32 bit-mask word, `.astype(str)` for
bytes, `as.integer()` for an `integer64`, `.filled()` for a masked array.

### Documentation that names something

**Every Ford `[[...]]` cross-reference is resolved against the project, and one that names
nothing warns** (`validate.check_doc_links`). It is the failure with no other symptom: Ford
renders a dangling link as the literal text an author typed, and `emit/doc_links` falls back to
plain code on purpose, because an R `\link` to a missing topic fails `R CMD check`. So a link
left pointing at a renamed module keeps generating, keeps compiling, keeps passing, and quietly
stops being a link in four languages at once — which is exactly what converting f42 to `_impl`
modules did to twenty-three of them, unnoticed for a month.

Resolution runs over the **augmented** project, so a link into a generated module is correct:
`[[f42_binary_search_tree(module):build_bst_index]]` names the wrapper a reader can call, and
that is what `doc_links` turns it into. What counts as existing is what Ford would document —
its parse is already filtered by visibility, so a link to a *private* routine warns, as it
dangles in the rendered documentation too.

Generic interfaces and derived types are read for this and nothing else (`ir.entities.Declaration`,
`frontend._declarations`): a binding can express neither, but authors link to both, and an
interface block's own doc comment is prose that belongs to no procedure and no module, so
without it the links inside it would be read nowhere. `_ITEM_TYPES` in `ir/doc.py` was missing
`interface` until 2026-08-12 — a `(interface)` link did not match the pattern at all, so it was
never a link to resolve, to check, or to render.

Errors rather than warnings. The generated documentation *is* the API's documentation in four
languages, and a link is the one part of it that can be checked mechanically -- so it is checked
at the weight of a signature the C layer cannot express.

Two pieces of Ford's syntax the parser has to get exactly right, because getting either wrong
makes a link silently not-a-link:

- **the two type vocabularies are not interchangeable.** A component takes `module`, `type`,
  `submodule`, `program`, `block`, `file`, `namelist` and the procedure/interface kinds, each
  optionally `ext`-prefixed for an entity of another project; an item takes `variable`, `bound`,
  `common`, `constructor`, `final`, `modproc` and the procedure/interface kinds. Ford warns and
  generates nothing for a mismatch. `_COMPONENT_TYPES` and `_ITEM_TYPES` each carried two of the
  other's entries until 2026-08-12.
- **inline code is left verbatim** (Ford 7). ``` `[[my_sub]]` ``` shows the syntax rather than
  following it, so `DocLine.parse` skips any link inside a backtick span -- which also keeps an
  R-style ``` `x[[i]]` ``` in a comment from being read as a link to `i`, and now from failing
  a build over it.

---

## Edge cases handled

Beyond the conventions above, the cases that took care to get right. Several exist *because*
the compiler or a run rejected the first attempt.

- **Null validation order** (`ierr`, then scalars, then arrays) so `c_loc` is never taken of
  a zero-size target; an empty array passes through to the callee's own validation. See
  [`design/c-layer.md`](design/c-layer.md).
- **Nullable optionals** cross as `bind(C)` `OPTIONAL` (interoperable since F2018), a plain
  pointer that is absent when C passes null — no descriptor, no branch.
- **A nullable mode** (optional, no `DM_DEFAULT`) is the one optional that is converted, so it
  cannot be handed straight on: the callee gets a pointer to the mode's integer, nullified and
  associated only when the string was given. The wrapper used to view the string unguarded,
  which segfaulted Python and R when the mode was omitted, and to pass its integer local,
  which is always present (#221).
- **`DM_RESULT_SIZE_IS` counts along the last axis**, at any rank: Python slices
  `[..., :n]`, R takes `utils::head(x, n)` for a vector and `x[, seq_len(n), drop = FALSE]`,
  with one leading comma per further axis, for an array. `head()` on a matrix took rows, the
  padding columns included (#221).
- **`intent(inout)`** is modified in place in Python but copied-and-returned in R, because R
  is copy-on-modify. The R side clones before Fortran touches the data.
- **Characters** carry their length as a leading C extent, and the wrapper takes a zero-copy
  `character(len=n)` pointer view of the caller's buffer rather than converting it. An **array
  of no strings** is viewed too, at the module's `character(1), target` anchor
  (`names.EMPTY_VIEW_ANCHOR`) with shape `[0]`, since `c_loc` may not be given the zero-size
  buffer itself: one `c_f_pointer` of an address chosen by `if`
  (`M_CHECK_CHARACTER_ARRAY_VIEW`, `names.VIEW_ADDRESS`). It used to be refused as
  `ERR_EMPTY_INPUT` before the procedure ran, so no binding could pass a data set with zero
  families (#221). A zero width with strings to view stays refused, as does a single string of
  zero width -- which only C can send, Python and R padding an empty string to one blank. Buffers are
  **blank**-filled, which is Fortran's own convention, so `trim` and comparison work on what the
  view yields and nothing has to scan for a NUL. Numeric outputs are left uninitialised
  (zeroing them is wasteful and dishonest — zeros look like results).
- **`c_bool`** logicals cross as real booleans, and since every stored logical in the library
  is `logical(c_bool)` the wrapper copies nothing — it hands the caller's buffer straight to
  the implementation. Declare yours the same way; a default-kind `logical` dummy would put the
  elementwise copy back.
- **Shape checks** the binding language makes but Fortran cannot: every input axis against what
  it must be, or a wrong answer / an over-read results (see "every input axis is compared" above).
- **Argument-named errors**: `ierr` encodes the offending argument's position, so an error
  names it (`invalid input arguments (argument 'n_dims')`).
- **`M_IMPLICIT_NONE`** (`implicit none (type, external)`) in generated modules, so a typo
  in a called helper is a compile error, not a link-time surprise.
- **An entity's own array-spec wins over `DIMENSION`**: in `dimension(n) :: a, b(m)`, `b` has
  `m` elements (F2018 8.2). Ford keeps listing `dimension(n)` among `b`'s attributes, and
  reading the attribute first sized `b` by `n` in every layer.
- **Assumed shape** (`x(:)`) is validated with the size asked of the array, `size(x, 1)`, where
  the extent `(:)` used to be pasted into the validator call. The C ABI already turned each `:`
  into an `n_<x>_elements` argument (`compute_rdi`'s `dscale`).
- **A module constant in an extent** (`points(N_AXES, n)`, `len=MAX_LEN`) is imported by both
  the Fortran wrapper and the C wrapper, which redeclare every argument. Before, such a wrapper
  compiled only where a range directive happened to name the same constant (#178's
  `MAX_N_BINS`). The C wrapper's import list comes from `CWrapper.size_constants`. Python and
  the R shim cannot name a Fortran constant at all, so they get its evaluated value
  (`CWrapper.with_constant_sizes_evaluated`, applied by each binding emitter's `module`): an
  output `sums(N_AXES)` is allocated with 2, a `len=MAX_LEN` string padded to 8, exactly as if
  the literal had been written. Before, the shim's `Rf_allocVector(REALSXP, N_AXES)` stopped the
  whole library build and Python raised `NameError`.
- **`impure elemental`** is impure. `elemental` implies `pure` unless `impure` says otherwise,
  and Ford lists both; taking `elemental` alone gave the wrapper a `pure` prefix it cannot
  compile with.
- **NA** (R) is rejected only where the check is free — integers, logicals, characters, in the
  shim's conversion pass (`INTEGER_NO_NA` keeps `1:n` O(1)) — never doubles, whose `NA` is a NaN
  Fortran already catches.

The full rationale for each is in [`design/language-layers.md`](design/language-layers.md)
and [`design/c-layer.md`](design/c-layer.md), including the alternatives that were rejected.

---

## How it is tested

Three layers, fastest first:

1. **Unit tests** on the IR, built from hand-made fixtures with no Ford — most of the suite.
2. **Frontend tests** parsing small fixture modules through the real Ford, plus a parse of
   the entire real `src/` that must come out with no diagnostics.
3. **End-to-end tests** that generate from the fixtures, **compile** the output
   (`gfortran -std=f2018`), build a shared library, and **call it** — from generated Python,
   from generated R, and from Fortran itself. This is the only thing that proves the output
   is not just plausible but correct: `fx_sum_matrix` really sums a matrix, and the wrappers
   synthesised from `rank_scores_impl` really behave as two tiers — `rank_scores` allocates
   the work array and seeds and heapsorts the permutation that `rank_scores_expert` takes from
   the caller, both agree with a hand-driven call of the implementation itself, and each
   rejects a value below its documented minimum with the argument position its *own* caller
   sees rather than the implementation's.

   The Fortran one holds *generated* code to `-std=f2018` but compiles the f42 modules it
   links against without it, exactly as fpm does — `f42_math` uses the F2023 `reduce()`
   locality spec and `f42_sort` declares `array(n)` before `n` is typed. Neither is generated,
   so holding them to the generator's conformance bar would only break the test.

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
  touches the frontend or IR. **Read [`emit/README.md`](emit/README.md) first**: it is the
  checklist of every edge case Python and R each had to answer — layout, ownership, absent
  optionals, missing values, partial buffers, error decoding — plus what to wire in and how to
  prove it works. The three existing emitters are the template.
- **A new convention or `DM_` macro**: add the macro to `src/macros.h`, its recognition to
  `ir/directives.py` (or a naming rule to `config.py`), its meaning to `ir/roles.py`, and a
  validation rule if one applies. Nothing about a target language changes.
- **A new type mapping**: extend the tables in `abi/c_abi.py` and the emitters. An unmapped
  kind is already an error, so a gap is loud.

Every convention lives in `config.py`; there are no naming literals scattered through the
code. Every name the generated code introduces beside the arguments -- a local, an import, a
`<argument><suffix>` -- is built from `names.py`: underscored where the language allows it, so it
meets no argument, and otherwise refused where a Fortran wrapper emits it. An emitter that starts
binding a plain name of its own fails `tests/test_names.py`.

---

## Open items

- **ifx**: the F2018 features used (`OPTIONAL` in `bind(C)`, `implicit none (type,
  external)`) are verified with gfortran only. ifx is expected to agree; worth a check.
- **Compile check in CI**: the end-to-end tests need a compiler. They are marked to skip
  without one; wiring them into CI is a configuration step.
