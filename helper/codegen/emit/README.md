# Adding a target language

One emitter per target lives here. Each turns the same `abi` model into a binding, and each
had to answer the same set of questions — about layout, missing values, ownership, absent
optionals, partial buffers — before it worked.

**This document is that set of questions**, with the answer Python and R each gave. It exists so
that a third language (Lua, MATLAB, Julia, whatever) is a matter of *deciding* each row rather
than rediscovering that the row exists. Most of them were found the hard way; several were found
by a segfault.

Rationale for the existing answers is in [`../design/language-layers.md`](../design/language-layers.md)
(Python and R) and [`../design/c-layer.md`](../design/c-layer.md) (the C ABI beneath them). This
is the working checklist; those are the arguments.

---

## Contents

- [1. What is already decided when your emitter runs](#1-what-is-already-decided-when-your-emitter-runs)
- [2. Two shapes of a language layer](#2-two-shapes-of-a-language-layer)
- [3. The pipeline every wrapper follows](#3-the-pipeline-every-wrapper-follows)
- [4. The edge cases](#4-the-edge-cases)
  - [4.1 Memory and layout](#41-memory-and-layout)
  - [4.2 Ownership and mutation](#42-ownership-and-mutation)
  - [4.3 Arguments the caller never sees](#43-arguments-the-caller-never-sees)
  - [4.4 Optionals and defaults](#44-optionals-and-defaults)
  - [4.5 Types that do not cross cleanly](#45-types-that-do-not-cross-cleanly)
  - [4.6 Outputs](#46-outputs)
  - [4.7 Errors](#47-errors)
  - [4.8 Documentation](#48-documentation)
- [5. What a language layer must not do](#5-what-a-language-layer-must-not-do)
- [6. Wiring a new target in](#6-wiring-a-new-target-in)
- [7. Proving it works](#7-proving-it-works)

---

## 1. What is already decided when your emitter runs

You consume [`abi/model.py`](../abi/model.py). Everything in it is derived once, so the targets
cannot disagree — **never re-derive any of it in an emitter**:

| You are handed | Meaning |
|---|---|
| `CWrapper.name` / `.stripped_name` | the exported C symbol, and the name your binding should use (`_c` dropped, and a hand-written `_alloc` pair already translated) |
| `CArgument.type`, `.dimension`, `.intent` | the type *as C sees it*, already mapped from the Fortran kind |
| `.origin` | `ARGUMENT`, `RESULT`, `EXTENT`, `STRLEN` or `ERROR` — the last three are invented for C and are yours to supply, not the caller's |
| `.conversion` | `NONE`, `LOGICAL`, `CHARACTER` or `MODE` — the value needs work between your language and C |
| `.optional`, `.default`, `.has_default` | nullability, and the default **already evaluated** from `DM_DEFAULT` |
| `.mode` | the accepted values of a mode argument, as a table of name → integer |
| `.sizes` / `.axis` | for an `EXTENT`: which argument and which 0-based axis it measures |
| `.shape_arg` | for an array whose extents travel separately, the argument carrying them |
| `.size_extents` | the extents to multiply for the element count |
| `.source.roles` | the IR roles: the `is_*` predicates (`is_extent`, `is_shape_arg`, `is_mask_count`, `is_temporary`, `is_computed`, `is_derived`, `is_mode`) and the links between arguments (`extent_of`, `shape_of`, `shape_arg`, `mask_count_of`, `count_arg`, `result_size_arg`, `computed_from`) |
| `CWrapper.validation_order` | the order in which pointers may safely be null-checked |
| `.doc` | the parsed documentation, links and all |

If something you need is not there, add it to the ABI layer — not to your emitter. The moment
two emitters answer the same question separately, they start to drift, and
[`../tests/test_emitter_parity.py`](../tests/test_emitter_parity.py) exists because that has
happened repeatedly.

---

## 2. Two shapes of a language layer

| | **Direct FFI** (Python) | **Compiled shim** (R) |
|---|---|---|
| How the host reaches Fortran | `ctypes` loads `libtensor-omics.so` and calls the `_c` symbol | a C `.Call` shim compiled *into* the same `.so`, called from an R wrapper |
| What is generated | one `.py` per module | a `.c` shim **and** an `.R` wrapper per module, plus a registration `init.c` |
| Who converts host values | the host language | the shim (`REAL(x)` → `double*`) |
| Who validates | the host language | the host language — **never the shim** |
| Cost | none beyond the FFI | one host-language call per invocation |

Which one you need depends on whether the host can hand a raw pointer to a C function. LuaJIT's
FFI and MATLAB's `loadlibrary` can, so they look like Python; a MEX file or a CPython C extension
cannot avoid a shim, so they look like R.

If you do need a shim, the split that made R work is: **the host language decides and raises; C
marshals and calls.** The shim receives an opaque host object, so it cannot name the argument in
an error message; the host can, for free, and can raise a typed error the caller can catch.
Validating in C means calling back into the host runtime to construct that error — fragile and
slow. See `../design/language-layers.md`, "validation lives in R, not in the C marshalling layer".

---

## 3. The pipeline every wrapper follows

Both emitters generate the same nine steps in the same order. The order is not stylistic — each
step needs what the one before it settles.

```
1  parameters      which C arguments the caller actually supplies (§4.3)
2  coerce inputs   type check + convert; raise naming the argument (§4.5)
3  producers       call DM_OUTPUT_FROM(AUTO) producers (§4.3)
4  derive extents  from the arrays, shape arguments and masks you were given (§4.3)
5  cross-check     arguments sharing an extent must agree (§4.3)
6  allocate        outputs, with the right initialisation (§4.6)
7  call            the C symbol
8  check ierr      decode and raise; statuses do not raise (§4.7)
9  return          the outputs, in a shape idiomatic for the host (§4.6)
```

Two ordering traps worth stating outright:

- **Producers bracket extent derivation.** An extent a producer *reads* must be settled before
  the producer is called; an extent read off what a producer *makes* can only be settled after.
  So step 4 runs in two halves, around step 3. (`python_ctypes.py:363`)
- **Null checks have a forced order** — error code, then every scalar, then arrays, and an array
  whose shape travels separately comes last. `c_loc` may not be taken of a zero-size target, so
  an array cannot be checked before the extents sizing it are readable. `validation_order` hands
  you this; do not re-sort it.

---

## 4. The edge cases

The checklist proper. Each row is a decision your emitter has to make; the last column is what to
think about for a new host.

### 4.1 Memory and layout

| Case | Python | R | For a new host |
|---|---|---|---|
| **Array order** | C-order by default, so a rank-≥2 input may cost a **transposing copy** | already column-major — free | Fortran is column-major. Julia and MATLAB are too (free); Lua tables and NumPy defaults are not. This decides whether rank-≥2 inputs are cheap or copied |
| **Rank 1** | the two orders coincide — no copy | same | never pay for a transpose at rank 1 |
| **Contiguity** | `np.ascontiguousarray` / `asfortranarray` as the rank requires | R vectors are contiguous | a host with strided or lazy arrays must materialise before passing a pointer |
| **Element count** | product of `size_extents` | same | C sees one flat block; the rank lives only in your wrapper |

### 4.2 Ownership and mutation

| Case | Python | R | For a new host |
|---|---|---|---|
| **`intent(inout)`** | modified **in place**, returns `None` | **copied and returned** (`v <- f(v)`) | follow the host's own convention. R is copy-on-modify, so mutating through a pointer is a bug, not an optimisation |
| **An `intent(inout)` that would need conversion** | **rejected** with a message saying why — converting copies, so the caller would never see the change | coerced, since it is copied anyway | if your host mutates in place, a conversion silently breaks the contract; refuse it |
| **Scalar outputs** | boxed (`byref(c_int(...))`) — Fortran writes through the pointer, and Python cannot rebind the caller's name | the shim allocates and returns | every host needs somewhere for C to write |

### 4.3 Arguments the caller never sees

The single most valuable thing this layer does, and the easiest to get subtly wrong: an argument
both *asked for* and *derived* is silently overwritten. There is a parity test for exactly that.

| Kind | Where it comes from | Note |
|---|---|---|
| **extents** (`n_dims`) | the array's shape | `.sizes` + `.axis` say which array and which axis |
| **character length** (`STRLEN`) | the item size of an input buffer — **only** for an input; for an output the caller must say how long the strings are, since the buffer does not exist yet |
| **shape arguments** (`<arg>_shape`) | the input's shape; for an *output* sized by a shape, the caller supplies the shape because it is what they want produced |
| **mask counts** (`n_selected_x`) | `sum()` / `count()` of the mask |
| **work arrays** (`tmp_`) | allocated by you, never returned, never asked for |
| **`DM_OUTPUT_FROM(AUTO)`** | call the producer's **own generated wrapper**, not its C symbol — validation and error handling come free, and the value you pass is already prepared, so the producer's conversions are no-ops |
| **the error code** | you allocate it, you decode it, the caller never sees it |
| **cross-checks** | two arguments sharing an extent must agree; Fortran cannot check it, and a mismatch is a wrong answer or a segfault. Name **both** arguments in the message |

Two producer subtleties: match producer inputs to consumer arguments **by name** unless the
argument's table renames them (the keyword is the *producer's* name, the value the consumer's);
and if the producer refines a value you handed it — capping it, say — adopt what it returned.
In Python the producer import goes *inside* the function, because two modules may size each
other's outputs and a module-level import would be circular.

### 4.4 Optionals and defaults

| Case | Handling |
|---|---|
| **optional with `DM_DEFAULT`** | **always passed.** The default is evaluated once in the ABI layer, so no two targets can disagree about it. Put it in the signature — a real default (`span = 0.7`), not a null |
| **optional without a default** | nullable: a null pointer is how C says "absent" |
| **null-checking an optional** | never — a null *is* the value. This is why extents and shape arguments may not be optional |
| **host quirk** | `ctypes` rejects `None` for a checked `argtype`, so the argtype has to allow it explicitly. Expect your host's FFI to have an equivalent wrinkle |

### 4.5 Types that do not cross cleanly

| Type | What C sees | What your emitter must do |
|---|---|---|
| **logical** | `c_bool`, one byte | the host's boolean is probably an `int` (R) or an object (Python) — convert, and note the scan is free, so check for missing values while you are there |
| **character in** | a `c_char` buffer with the length as the **leading extent** | encode; the item size *is* the string length |
| **character out** | the same buffer | **zero-fill it** (§4.6) |
| **mode argument** | an integer, but the caller passes a **string** | lower-case the string and let the C wrapper map it; an unknown mode is rejected before Fortran is entered. The accepted set is in `.mode` |
| **missing values** | nothing — Fortran has no `NA` | check where the check is free. R checks integers (`anyNA`, ALTREP-aware so `1:n` is O(1)), logicals and characters (converted anyway), never doubles (`NA_real_` is a NaN payload Fortran already catches). A host with a distinct missing marker needs the same table; one without needs none of it |
| **kinds with no mapping** | — | already an error upstream. Never guess a mapping: a wrong guess compiles and lies |

### 4.6 Outputs

| Case | Handling |
|---|---|
| **numeric output** | allocate **uninitialised** (`np.empty`). Fortran fills it; zeroing is an O(n) write of data about to be overwritten, and where Fortran fills only part, zeros *look like results* while garbage looks like garbage |
| **character output** | allocate **zero-filled**. Fortran fills a character buffer only partially and the nulls terminate the strings; an uninitialised byte past the written data is read back *as part of a string* |
| **result trimming** (`DM_RESULT_SIZE_IS`) | return only the first *n* elements, where *n* is the argument named |
| **serialized array out** | reshape to the shape argument's contents before returning, column-major — the caller wants the n-d array, not a flat buffer plus homework |
| **shape of the return** | Python returns a scalar for one output and a `dict` for several; R returns the value or a list. Follow the host |
| **result mutability** | freeze it if the host can. Python sets `flags.writeable = False` on the allocated buffer after the call — O(1), and it propagates into the reshape and slice views the return builds, which then cannot be unfrozen at all. R has no counterpart and needs none: it is copy-on-modify, so a returned vector cannot be aliased. Freeze only what the wrapper allocated — never an `intent(inout)` array, which is the caller's — and say so in the docs with the host's escape hatch (`.copy()`) |
| **`intent(inout)` in the return** | see §4.2 — it depends on whether you mutate |

### 4.7 Errors

`ierr` never reaches a caller. Decode it and raise; a successful call returns results and nothing
else.

- The code is packed as `M_ERR_ARG_POS_FACTOR * arg_pos + code` (the factor is read from
  `src/macros.h`, never hardcoded).
- `arg_pos` indexes the **Fortran** dummy list at every layer, so pass the argument-name tuple in
  that order and index it.
- Pass the *source* of each derived argument alongside, so the message blames something the
  caller actually wrote: `(argument 'vec1', via 'n_elements')`.
- `ERR_*` codes raise; `STAT_*` codes are outcomes and never do.
- Generate a per-language error module from the catalogue (`errors_python.py`, `errors_r.py` are
  the templates): a decode function, the messages taken from the Fortran documentation, and one
  error type **per group** (IO, input, memory, runtime, internal) rather than per code — thirty
  types nobody would name individually is not a usable API. Prefix them, so a host-language
  `except`/`tryCatch` on the base type cannot catch something unrelated.
- If you use a shim, raise from the **host**, not from C.

### 4.8 Documentation

- Render the host's own docstring dialect from the parsed `Doc` — numpydoc and roxygen2 are the
  two instances (`doc_numpydoc.py`, `doc_roxygen.py`).
- Render **literals** per language: `0.7_real64` → `0.7`, `.true.` → `True` / `TRUE`
  (`doc_literals.py`). A Fortran kind suffix in a doc is noise everywhere else — and in an
  *extent expression* it is a syntax error, so the same stripping applies there.
- Resolve Ford links to **what the reader can call in your language** (`doc_links.py`): a link to
  an implementation becomes the binding that wraps it, a link to a mode parameter becomes
  the mode string.
  Resolve parsed link spans, never by text substitution — a text pass over `[[...]]` mangles R's
  own `x[[i]]`.
- Never rewrite author prose. If a sentence is wrong, it is wrong in the Fortran.
- State provenance: which Fortran procedure this came from, since its argument names are what
  errors report.

---

## 5. What a language layer must not do

- **Re-derive the ABI.** If you compute a C signature, you have already diverged.
- **Validate in a shim.** The host knows the argument names; C does not.
- **Ask for a derived argument.** It will be overwritten, and the caller will never know why.
- **Skip the cross-checks.** They are the thing Fortran cannot do for itself.
- **Invent a naming convention.** `stripped_name` is the name. A generated pair is already
  called `foo` / `foo_expert` in Fortran; a hand-written `foo` / `foo_alloc` pair has been
  translated into those same two names before it reaches you.
- **Hardcode a constant that lives in `src/macros.h`** — the arg-pos factor above all.
- **Mirror `NO_INPUT_VALIDATION`.** That directive drops the *Fortran* wrappers' checks at
  compile time, for a caller who has established their inputs are good. Your layer's own
  validation is a run-time thing and cannot be preprocessed out, so do not try to make it
  conditional: a host caller who wants no checks calls through C or Fortran instead. The
  cross-extent checks in particular (§4) are the one thing Fortran cannot do for itself, so
  they are exactly what should *not* become optional.

---

## 6. Wiring a new target in

| Touch | What for |
|---|---|
| `emit/<lang>.py` | the emitter itself |
| `emit/errors_<lang>.py` | the error module, from the `tox_errors` catalogue |
| `emit/doc_<dialect>.py` | the docstring renderer, if the host has its own dialect |
| `emit/doc_literals.py` | a literal-rendering entry for the host, if `True`/`TRUE` are not enough |
| `emit/doc_links.py` | how a link renders in the host |
| `config.py` `Paths` | `<lang>_out_dir`, plus a `<lang>_binding_dir` if you generate a compiled shim |
| `generate.py` | a `_<lang>_files()` builder, and its branch in `generate()` |
| `generate.py` `_clean` | so a procedure that stops being exported leaves no stale wrapper |
| `cli.py` | add the target to `C_BINDING_TARGETS` (and thus `ALL_TARGETS`) |
| `emit/vscode_snippets.py` | a snippet flavour for the host, if you want editor support |
| `.gitattributes` | mark the output tree `linguist-generated` |
| `README.md`, `design/language-layers.md` | what the layer costs in this host, and every decision you took above |

A shim-based target additionally needs its shim sources compiled into `libtensor-omics.so` (fpm
picks up anything under `src/`, which is why the R shims live in `src/generated/bindings/r`) and
a registration entry point. Note the trap the R shims hit: their host-API symbols are undefined
until the host loads the library, so Python's eager `ctypes` load of the *same* `.so` would fail
on them — they are marked **weak** so the load resolves them to null, while a genuinely missing
symbol still fails loudly.

---

## 7. Proving it works

Three levels, and the third is the one that matters:

1. **Unit tests** over hand-built IR — no Ford, fast, most of the suite.
2. **Parity** — [`../tests/test_emitter_parity.py`](../tests/test_emitter_parity.py). Extend it
   with your target. It asserts what no single-target test can see: that every target asks its
   caller for the *same* arguments, that a shim is called with the arguments it declares, and
   that nothing is both asked for and derived. Run it over the fixtures **and over the real
   `src/`** — every divergence found so far was in `src/`, not in the fixtures.
3. **End to end** — generate from the fixtures, **compile**, and **call** it, the way
   `test_end_to_end.py` (Python) and `test_end_to_end_r.py` (R) do. Assert real behaviour:
   the sum is right, an out-of-range value is rejected with the argument position its own caller
   sees, and the allocating entry point really allocates and sorts. Anything less proves the
   output is plausible, not correct.

The fixtures in [`../tests/fixtures/`](../tests/fixtures/) are the specification: a deliberate,
complete set of the constructs the generator supports, held to a no-diagnostics bar. If your host
needs a construct they do not cover, the fixture set is what grows.
