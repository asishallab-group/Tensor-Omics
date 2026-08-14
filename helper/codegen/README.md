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
committed bindings no longer match the sources.

Requirements: Python 3.11+ (`contextlib.chdir`), and `FORD`, `pcpp`, `numpy` (see
`requirements-dev.txt`). The generator reads the Ford settings from `fpm.toml`.

### Build switches the output honours

The generator always emits all of it; these decide what the *build* keeps. Pass them through
`./build.sh --directive=NAME`, which forwards to the preprocessor.

| Directive | Drops | Notes |
|---|---|---|
| `NO_C_BINDING` | the `bind(C)` wrappers, and the R shims with them | they call the `bind(C)` symbols; takes the `f42_safeguard` dependency with it |
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
| `r_binding_dir` | `src/generated/bindings/r` | **output**: the R C `.Call` shims (fpm bundles these into `libtensor-omics.so`) |
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
    __init__.py  library.py  error_handling.py
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
  `synthesize.generated_path_for` is the single place that decides it: the emitter and the
  cleaner both go through it, and the Ford frontend excludes the whole `generated_dir` they
  write into — rather than three literals kept in step by hand, which is how one of them comes
  to write, delete and exclude three different paths.
- **Each target directory is cleaned before writing** (unless `--no-clean`), so a procedure
  that stops being exported leaves no stale wrapper. In the R *package* directory only the
  generated `*.R` are removed — a hand-written `DESCRIPTION` or `NAMESPACE` next to them is
  left alone. The wrapper modules are the one target cleaned file by file rather than by
  directory, because they share `src/generated` with the binding trees; the list comes from
  the sources (`synthesize.generated_wrapper_paths`), so a wrapper whose implementation has
  been *deleted* is not named by it and survives on disk. Harmlessly: Ford excludes
  `generated_dir` as a directory, so nothing stale is ever parsed back in.
- **`--library`** (default `build/libtensor-omics.so`) is not an output path; it is where the
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
| valid ranges, finiteness, masks, distance matrices, optionals | guide §5.1–5.6 |
| work arrays, permutations, recommend sizing | guide §5.7–5.9 |
| modes, per-mode procedures | guide §5.10–5.12 |
| prologues, runtime errors, split families | guide §5.13–5.15 |
| exporting a hand-written procedure | guide §6 |
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
- **A mode crosses as a string.** The binding languages pass the parameter name without its
  prefix, lower-cased (`METHOD_WARD` → `"ward"`), and the C wrapper maps it back to the integer,
  rejecting an unknown one before Fortran is entered. The Fortran wrapper separately checks
  membership against exactly the values the mode table names.
- **`DM_OUTPUT_FROM(..., AUTO)` is called by the language layer**, so the producer must be
  exported — there has to be a wrapper to call. Python imports it *inside* the calling function
  rather than at the top of the file, because two modules may size each other's outputs and a
  module-level import would then be circular. R needs no import: every wrapper is sourced into
  one environment.
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
  list. Each conditional import now has a predicate beside the code that emits it
  (`_converts_an_input`, `_conversion_helpers`), with tests that a module needing none asks
  for none. When auditing this by hand, expand `src/macros.h` first: most apparent unused
  imports are named inside `M_ALLOCATE` and `M_CHECK_NON_NULL`.
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

An implementation module is checked separately, because its procedures are read rather than
exported and so are reached by none of the above. These rules follow the `_impl` suffix
wherever it is written — under `src/tox`, under `src/f42`, anywhere — which is the price of
dropping the fixed directory they used to be scoped to:

- a module **not in a file named after it** — the *module* name is what synthesis triggers on
  and what the generated module is called; the *file* name is what its path is mirrored from,
  and what the cleaner reads without parsing anything. Let the two disagree and `module tox_x`
  is written into some unrelated file — or into none, since a module named `<x>_impl` in a
  file that is not an `_impl` file has no mirrored path at all and the emitter raises rather
  than guessing a directory
- an **`M_EXPORT_C` on an implementation** — its wrapper is the entry point, and the export
  publishes an unvalidated twin beside it, under a name (`foo_impl`) a binding caller cannot
  tell apart from the validated `foo`. Support routines in the same module (the
  `DM_OUTPUT_FROM` recommend routines) have no wrapper and stay exported, which is why an
  implementation module can itself show up in the bindings — `tox_loess_impl.py` is the
  Python home of `tox_loess_required_workspace`, and holds no implementation at all
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

  The list is a curated boundary, not a proof: nothing checks that a member is itself
  allocation-free, and `tox_conversions` is not — two of its procedures return a
  deferred-length string, which Fortran can only express as an `allocatable`. No
  implementation imports it, so nothing is wrong today. `f42_utils` sat on the list while its
  `f42_stats` child still had hand-written `_alloc` halves; converting that family removed the
  entry rather than justifying it, which is the move to reach for — a module that has to be
  whitelisted may be a module that should be an implementation

A `DM_PROLOGUE` is checked too, having had no analysis or validation at all before:

- the named procedure must exist
- it must declare a scalar `logical(c_bool), intent(out) :: handled` — the wrapper returns early on it
  regardless, so without it that branch reads an undefined value
- a dummy naming an implementation argument is supplied from it; one naming nothing becomes an
  argument of `foo`, which is what the prologue derives *from*. A name **one edit** from an
  implementation argument is refused as a misspelling, since it would otherwise become a new
  argument and leave the prologue and the implementation working from different values
- a dummy that some mode's wrapper does not have — the mode argument, or one scoped to a mode —
  is refused: the prologue runs in all of them
- the implementation must generate an allocating wrapper for it to run in (it does when
  anything is taken over, or when the prologue asks for an argument of its own)
- it may not produce anything the allocations, permutation sorts or recommend calls above it read

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
- **`intent(inout)`** is modified in place in Python but copied-and-returned in R, because R
  is copy-on-modify. The R side clones before Fortran touches the data.
- **Characters** carry their length as a leading C extent. Output buffers are zero-filled
  (Fortran fills them only partially; the nulls terminate the strings), numeric outputs are
  not (wasteful and dishonest).
- **`c_bool`** logicals cross as real booleans, and since every stored logical in the library
  is `logical(c_bool)` the wrapper copies nothing — it hands the caller's buffer straight to
  the implementation. Declare yours the same way; a default-kind `logical` dummy would put the
  elementwise copy back.
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
code.

---

## Open items

- **ifx**: the F2018 features used (`OPTIONAL` in `bind(C)`, `implicit none (type,
  external)`) are verified with gfortran only. ifx is expected to agree; worth a check.
- **Compile check in CI**: the end-to-end tests need a compiler. They are marked to skip
  without one; wiring them into CI is a configuration step.
