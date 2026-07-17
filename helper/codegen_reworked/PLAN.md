# Code generator rework — implementation plan

Working branch: `131-codegen-new` (off `122-parallelize-codebase`).
This document is the agreed plan. It is transient: once the generator is in place its
design content moves into `README.md` and this file is deleted.

## 1. Context

The generator on `131-code-gen` works — it runs end-to-end and reproduces its committed
output byte-for-byte. Its *pipeline* (Ford parse → model → per-language serializer) is the
right shape and is kept. Its *serializer mechanism* is not, and is replaced.

### Why not keep the current implementation

`CodeGenerator.__format__` (`helper/codegen/api/utils.py:186`) dispatches to a serializer
held as **global class state** (`CodeGenerator.serializer`), invoking serializer methods with
the *model* object bound as `self`. Emission is therefore keyed on **format-spec strings**
(`"dummy"`, `"arglist"`, `"type_conversion_inputs"`, …). Consequences, in order of severity:

1. **Untestable.** The model cannot be built without a full Ford parse of the whole project,
   and only one serializer can be active process-wide. This is the direct reason there are
   zero tests today.
2. **Silent gaps.** An unhandled spec surfaces as a `ValueError` from deep inside a nested
   f-string, or — worse — as an unbound local. No emitter is checked for completeness.
3. **`self` is a lie.** `Serializer.FordTable(self, spec)` reads `self.header`; `self` is the
   `FordTable`, not the serializer. Every emitter is an unbound visitor in method's clothing.

Live bugs that the design allowed to go unnoticed, all reproduced on `131-code-gen`:

| Location | Bug |
|---|---|
| `codegen/c_wrapper.py:146` | `elif meta["optional_output"] is not None:` — `meta` is undefined in that scope. Any optional-output argument raises `NameError`. `variant` can also be unbound. |
| `codegen/api/utils.py:4` | `eval_expr` calls `np.pi`/`np.arccos` but `utils.py` never imports numpy. Any call raises `NameError`. (The `fortran.py` copy of the same function works — it is duplicated.) |
| `codegen/api/c_wrapper.py:218` | `procedure.fortran_procedure.find_child(f"{procedure.name}_alloc")` searches the *procedure's* children, not the module's. The `_expert_c` naming rule never fires. |

Rewriting is cheaper than retrofitting tests onto a design whose core prevents them.

### Decisions taken (confirmed)

- **Rewrite** under `helper/codegen_reworked/`; `helper/generate_code.py` imports from there.
  The old `helper/codegen/` stays untouched until the rewrite reaches parity, then is removed.
- **No source annotations are touched.** `122-parallelize-codebase` carries no
  `category: C-interface` tags and no `DM_` macros. Those are added by the maintainer later.
  The generator is developed and tested exclusively against **fixture Fortran modules** owned
  by the test suite. The one exception is `src/macros.h`, which must gain the `DM_` *macro
  definitions* — those are generator infrastructure, not annotation (see §4).
- **Tests:** IR unit tests on hand-built fixtures + golden-file tests through the real Ford
  frontend. Compiling generated wrappers with `gfortran` is an opt-in extra (§6.4).

## 2. Target architecture

Strict one-way dependency: `frontend → ir → abi → emit`. Nothing downstream imports Ford;
nothing upstream knows a target language exists.

```
helper/codegen_reworked/
  diagnostics.py        Diagnostic, DiagnosticBag, SourceLocation, CodegenError
  config.py             paths, prefixes/suffixes, macro names — one place, no literals in code
  render/
    writer.py           Writer: line(), blank(), indent() ctx mgr, block()  (no str subclass)
  frontend/
    ford_frontend.py    Ford → ir.Project. The ONLY module importing ford.
    macros.py           pcpp wrapper for src/macros.h (DM_ expansion + regex-escaped variant)
  ir/
    types.py            Intent, FortranType, Dimension, Kind
    doc.py              Doc, DocLine, FordLink, FordTable  (language-neutral doc tree)
    macros.py           DM_* directive parsing → typed Directive objects
    entities.py         Project, Module, Procedure, Argument
    roles.py            semantic pass: temporary / mode / dim / shape / mask-count / …
    errors.py           ErrorCatalogue from tox_errors (ERR_ = error, STAT_ = status)
    validate.py         rule checks → DiagnosticBag
  abi/
    c_abi.py            ir.Procedure → CWrapper IR (the C-ABI decision layer)
    model.py            CWrapper, CWrapperArgument, CWrapperModule
  emit/
    fortran_c.py        CWrapper → src/c_interface/<mod>_c.F90
    python_ctypes.py    CWrapper → python/tensor_omics/<mod>.py   (numpydoc)
    rcpp.py             CWrapper → rcpp/tensor_omics/<mod>.cpp    (roxygen2)
    rcpp_types.py       Tox* SEXP wrapper classes (Rcpp::as interface)
    errors_python.py    tox_errors → python error module
    errors_r.py         tox_errors → R error module
    doc_numpydoc.py     Doc → numpydoc
    doc_roxygen.py      Doc → roxygen2
  cli.py                argparse entry point
  tests/                see §6
  README.md             design + edge cases + dev guide (deliverable)
```

**Emitters are plain classes with real methods and explicit arguments.** No global state, no
format specs, no `__format__`. Each emitter is constructed with its config and a `Writer`, and
exposes intention-revealing methods (`emit_dummy_decl`, `emit_null_checks`, …). Two emitters
can run concurrently; an emitter that forgets a case fails a test, not a user.

**The `abi` layer is where C-specific structure is decided** — synthesising extent/strlen
arguments, appending `ierr`, turning a function's result into an `intent(out)` argument,
rewriting mode integers to characters. Python and Rcpp emitters consume that same `CWrapper`
IR, which is what keeps the three targets consistent by construction. This mirrors what
`C_Wrapper_Modules` does today; it is the one piece of the old design worth preserving.

### What gets ported rather than rewritten

The old serializers encode real, working knowledge. Ported as ordinary functions/tables:

- type maps (Fortran ↔ `c_*` ↔ ctypes ↔ numpy dtype ↔ Rcpp type)
- `tox_conversions` call patterns for logical/character in/out
- mode-table parsing and `MODE_X` → `"x"` string mapping
- Ford link/table parsing (`doc.py` is largely sound; it gains tests and typed directives)
- numpydoc layout

## 3. Source-language contract

Recorded here because it is the generator's input spec and §7 of the README expands it.
Conventions from issue #131 and its comments; anything marked **new** is a proposal to
confirm during implementation.

**Selection.** `category: C-interface` meta tag on the procedure.

**Naming.** `<p>_alloc` → `<p>_c`. `<p>` → `<p>_expert_c` **iff** a `<p>_alloc` exists in the
same module, else `<p>_c`. (Resolved against the module's routines — fixes the bug above.)

**Prefix/suffix conventions.** `tmp_` temporaries; `<arg>_shape` shape arguments;
`n_selected_<arg>` mask counts for `<arg>_mask` / `<arg>_selection_mask`;
`mode`/`method` or `*_mode`/`*_method` mode arguments.

**Doc macros** (`DM_` = doc macro, vs `M_` code macros, `CM_` in-file macros):

| Macro | Meaning | Status |
|---|---|---|
| `DM_DEFAULT(VAL)` | constant-expression default for an optional | exists on `131-code-gen` |
| `DM_REQUIRED_IF_MODE(MODE_ARG, MODULE, MODE_VAR)` | optional required only for a given mode | exists |
| `DM_RESULT_SIZE_IS(ARG)` | `ARG` names the scalar holding the used result count | exists |
| `DM_OPTIONAL_OUTPUT` | output only produced on request | exists |
| `DM_OUTPUT_FROM(ARG, PROC, MODULE, AUTO\|JUST_INFO)` | argument comes from another procedure's output | exists, **unimplemented** |
| `DM_SHAPE_OF(ARG)` | **new** — explicit shape-arg link where naming is insufficient | propose |

`DM_OUTPUT_FROM(..., AUTO)` requires an argument mapping (which of `PROC`'s arguments are fed
from where, and which the caller must still supply). Per the proposal this needs "a following
table". Design a `| Argument | Source |` table analogous to the mode table (§5, step 9).

## 4. `src/macros.h`

`122-parallelize-codebase` has no `DM_` definitions. Port the block from `131-code-gen`'s
`src/macros.h` (7 defines) as a **standalone commit**, unchanged in spirit:
definitions only, no procedure annotations. The fixture modules `#include <src/macros.h>`, so
the tests exercise the real macro definitions rather than a copy.

## 5. Work breakdown

Small commits, in order. Each is independently reviewable; each with code carries its tests.
No pushing.

| # | Commit | Contents |
|---|---|---|
| 1 | plan | this document |
| 2 | scaffolding | package skeleton, `config.py`, `diagnostics.py`, `pyproject`/`requirements-dev` (pytest), `render/writer.py` + tests |
| 3 | macros | `DM_` defines into `src/macros.h`; `frontend/macros.py` (pcpp wrapper, both variants) + tests |
| 4 | ir: types | `Intent`, `Dimension`, `FortranType`, kind handling + tests (incl. `len=*` / `len=:` / `len=n`) |
| 5 | ir: doc | `Doc`/`DocLine`/`FordLink`/`FordTable`, ported and typed + tests (tables, links, malformed input) |
| 6 | ir: directives | `DM_*` → typed `Directive` objects; replaces the regex-match-object dict + tests |
| 7 | ir: entities | `Project`/`Module`/`Procedure`/`Argument`, constructible **without Ford** (the testability keystone) |
| 8 | ir: roles | semantic pass: temporary, mode, dim-of, shape-of, mask-count-of, default, required-if-mode + tests |
| 9 | ir: output_from | `DM_OUTPUT_FROM` arg-mapping table + resolution (currently unimplemented) |
| 10 | ir: errors | `ErrorCatalogue` from `tox_errors`; `ERR_` vs `STAT_` split; arg-position handling |
| 11 | ir: validate | rule checks → diagnostics with source locations; fail-fast CLI reporting |
| 12 | frontend | Ford → IR adapter, `src_dir` injectable so tests parse only fixtures |
| 13 | abi | `c_abi.py`: extent/strlen synthesis, `ierr` injection, function→subroutine, mode→char, optional-as-null, assumed-size for shape-linked args + tests |
| 14 | emit: fortran_c | C wrapper emitter + golden tests |
| 15 | emit: errors | Python + R error modules from the catalogue + tests |
| 16 | emit: python | ctypes wrapper + numpydoc + shape cross-checks + golden tests |
| 17 | emit: rcpp types | `Tox*` SEXP wrappers implementing `Rcpp::as`, validation + preprocessing in construction |
| 18 | emit: rcpp | Rcpp wrapper + roxygen2 + golden tests |
| 19 | cli | `cli.py`, rewire `helper/generate_code.py` |
| 20 | parity | run against fixture corpus covering every edge case in §7; diff-review |
| 21 | README | design, edge-case catalogue, dev guide |
| 22 | cleanup | remove `helper/codegen/`, old `helper_*.py` if superseded |

Order rationale: 2–11 is pure Python with no Ford dependency and no I/O, so it is fully
testable before any generation happens. 12 attaches the parser. 13 is the single place where
"what the C ABI looks like" lives. 14–18 are leaves.

## 6. Test strategy

### 6.1 IR unit tests (fast, no Ford)
Build `Procedure`/`Argument` fixtures directly. Cover roles, directives, defaults, validation
diagnostics, type mapping. This is only possible because of step 7.

### 6.2 Golden tests (real Ford)
`tests/fixtures/fortran/<case>.F90` — small annotated modules, parsed by the real frontend with
`src_dir` pointed at the fixture dir, output diffed against `tests/golden/<case>/…`.
`--update-golden` flag to regenerate deliberately.

### 6.3 Edge-case corpus (fixture modules to write)
Each is a named case with an expected outcome — including the ones that must **fail**:

- naming: `_alloc` present/absent; `_expert_c` collision
- functions with result; result needing conversion
- `intent(in|out|inout)`; `intent`-less dummy
- character `len=*` (synthesised strlen), `len=n`, `len=123`, `len=:` → **error**
- character scalar vs 1-D vs 2-D; 3-D → **error** (unsupported by `tox_conversions`)
- logical scalar/array, in/out/inout
- assumed-shape `(:)` → explicit-size + synthesised extent
- shape argument `<arg>_shape` → assumed size `(*)`; shape arg not `intent(in)` → **error**;
  shape arg without dimension → **error**
- mode argument with valid table; missing table → **error**; malformed table → **error**;
  mode var not `MODE_`/`METHOD_` prefixed → **error**
- optionals: with `DM_DEFAULT`, without default (null), `DM_REQUIRED_IF_MODE`
- `DM_RESULT_SIZE_IS`, `DM_OPTIONAL_OUTPUT`, `DM_OUTPUT_FROM(AUTO|JUST_INFO)`
- `tmp_` temporaries: `intent(out)` and `intent(inout)`
- mask count `n_selected_<arg>`
- procedure without `ierr` → synthesised
- shared dimension names across arguments → cross-check in Python/R
- module with no tagged procedures → no file emitted
- zero-size array + `M_CHECK_NON_NULL` (the standing `c_loc` TODO in `macros.h`, §8)

### 6.4 Optional compile check
`gfortran` 16.1 and `fpm` 0.13 are present. A marked-skip test compiling the generated golden
`_c.F90` against a stub proves standard-compliance. Off by default; enabled locally/CI on demand.

## 7. Open design questions

To resolve during implementation — flagged rather than silently decided.

1. **`logical` ↔ `c_bool`.** The proposal wants `c_bool` "if bool is already part of the C
   standard" to avoid the copy. But `c_bool` is 1 byte while default `logical` is 4, so a copy
   is unavoidable *unless the wrapped routine itself declares `logical(c_bool)`*. Proposal:
   pass through with no copy when the original dummy is `logical(c_bool)`; otherwise convert.
   Emit `logical(c_bool)` (not `integer(c_int)`) in the wrapper signature either way, so
   Python/R pass real booleans. **Needs confirmation** — it changes the existing C ABI.
2. **Character copy avoidance.** For a scalar `character(len=n)` the proposal hopes for a
   pointer rather than a copy. Associating a `character(len=n)` pointer with an incoming
   `c_char` buffer is not straightforwardly standard-conforming. Research item; default to the
   existing copy via `tox_conversions` and document the finding either way.
3. **By reference vs by value.** The proposal hedges for scalars; issue #131 says everything by
   reference, no `value`. Null validation *requires* a pointer. Resolution: **by reference,
   always** — recorded so the hedge isn't reopened.
4. **Per-code R/Python exception types.** "Could have its own Error type if suitable, but not
   necessary." Default: one base `ToxError` + subclasses for the coarse groups already implied
   by `tox_errors`' numeric ranges (I/O 1xx, input 2xx, memory 3xx, internal 9xxx).
5. **`DM_REQUIRED_IF_PRESENT` / `M_DOC_NO_DEFAULT`** appear in issue #131 but not in
   `macros.h`. Confirm whether `DM_REQUIRED_IF_MODE` subsumes them.

## 8. Known upstream issue

`src/macros.h` carries a standing `TODO codegen`: `M_CHECK_NON_NULL(ARG)` calls `c_loc(ARG)`
before the array's declared extent is validated `> 0`, and `c_loc` on a zero-size target is not
standard-conforming. This is a generator concern — the generator emits the call order. Fix
within this work (validate extents before null-checking arrays, or null-check via a scalar
extent argument first) and cover with an edge case in §6.3.
