# Designing the language layers

Why the Python and R bindings are shaped the way they are, and what was rejected on the
way. If you find yourself asking "why didn't they just...", the answer should be here. If
it isn't, that is a gap worth filling.

Companion documents: [`c-layer.md`](c-layer.md) for the Fortran/C wrappers,
[`../README.md`](../README.md) for the generator itself, and
[`../emit/README.md`](../emit/README.md) — the same ground as a **checklist for adding a third
language**, where this document is the argument for what Python and R chose.

---

## The rule the layers follow

**A generated function asks only for what it cannot work out, checks what Fortran cannot,
and returns results rather than error codes.**

Everything below follows from that, plus one observation: **the same rule produces
different code in Python and R, because the host languages differ.** Where the two layers
disagree, it is deliberate, and the reason is recorded.

---

## What each layer costs, and why it differs

| | Python | R |
|---|---|---|
| Array layout | C order by default, so a rank-≥2 `intent(in)` input costs a **transposing copy** | **already column-major** — zero work |
| Passing to Fortran | `ndpointer` over the numpy buffer, no copy when the dtype matches | `REAL(x)` is a plain `double*`, no copy when the type matches |
| `intent(inout)` | **modified in place**, returns `None` | **copied and returned** |
| Missing values | none — `np.nan` is a float, Fortran's `validate_in_range_real` catches it | `NA` for every type, and `NA_integer_` is an ordinary number to Fortran |

The layout row is worth internalising: **R and Fortran agree natively, Python does not.**
A 2-D input in R is free; the same input in Python may cost a full copy. Nothing can be
done about it, and it is not a defect — but it is why the Python layer works harder.

---

## Decision: `intent(inout)` is in-place in Python, by value in R

**Python** modifies the caller's array and returns `None`, the way `list.sort()` does.
An `intent(inout)` array is therefore never converted — a conversion would copy, and the
caller would never see the modification. It is rejected instead, with a message saying so:

```
TypeError: 'vector' is modified in place, so it must already be a numpy array of np.float64
```

**R** copies and returns. R is copy-on-modify: mutating the caller's vector through a
pointer is not an optimisation, it is a bug that shows up as spooky action at a distance
in unrelated code.

```r
v <- fx_normalize(v)     # R
```
```python
fx_normalize(v)          # Python: v is modified
```

**This looks like an inconsistency and is not.** Both are the idiomatic answer in their
language. A generator that forced one convention on both would be wrong in one of them.

*Rejected:* returning the array from Python as well, for symmetry. It would imply a copy
was made and invite `v = fx_normalize(v)`, which reads as though `v` were rebound when it
was not.

---

## Decision: validation lives in R, not in the C marshalling layer

The natural instinct is to validate in C, where the data is. It is the wrong instinct.

The `.Call` shim receives only the `SEXP`; `REAL(x)` / `INTEGER(x)` hand back a bare
pointer that cannot know the argument was called `weights`, so C cannot say so:

```
Not compatible with requested type: [type=character; target=double].
```

R can, for free, and better:

```r
if (!is.double(vector)) {
  if (!is.numeric(vector)) tox_type_error("vector", "numeric", vector)
  storage.mode(vector) <- "double"     # copies only if it must
}
```

`is.double()` is O(1). The message names the argument. `stop()` raises a **classed
condition** (`tox_type_error` / `tox_error` / `error`) that `tryCatch` can select on --
and constructing one from C means calling back into the package environment, which is
fragile and slow.

And once R guarantees the type, `REAL(x)` in the C shim **cannot fail**, so the whole
problem disappears rather than being worked around.

So: **R validates and marshals nothing; C marshals and validates nothing.**

*Rejected:* the `Tox*` C++ wrapper classes (Rcpp `as<>`), as sketched in issue
#131. They exist to produce good type errors, and R produces better ones for less. The
`Tox` prefix survives on the **conditions** (`tox_error`), where it earns its keep.

*Known consequence:* validation is skippable by calling the internal `.Call("<name>_call", ...)`
directly. That is acceptable for a dot-prefixed internal entry point. If the C ever needs to
be defensive on its own — because something else starts calling it — this decision flips,
and the marshalling layer grows validation of its own.

---

## Decision: NA is checked where the check is free

`NA` has no Fortran counterpart, and `NA_integer_` and `NA_LOGICAL` are both `INT_MIN` --
an ordinary number. Fortran would happily compute on `-2147483648`.

The cost, per type:

| Type | Checked | Cost | Why |
|---|---|---|---|
| `double` | no | none | `NA_real_` is a NaN payload; `validate_in_range_real` already catches it |
| `logical` | yes | **none** | R logical is `int`, Fortran wants 1-byte `c_bool`; the conversion scans anyway |
| `character` | yes | **none** | must be copied into a `c_char` buffer anyway |
| `integer` | yes | one `anyNA()` | the only real cost |

`anyNA()` rather than a hand-rolled loop: it is ALTREP-aware, so for a compact sequence
like `1:n` it answers **in O(1) without touching memory**, and it is a tuned C loop
otherwise.

In practice the large arrays are doubles, which are not scanned at all; integers are
indices, counts and shapes.

*Rejected:* checking everything uniformly. It would scan the big arrays for a sentinel
that cannot occur in them.

*Rejected:* checking nothing and letting Fortran validate. `NA_integer_` would be reported
as an out-of-range value rather than as `NA`, if it were reported at all.

*Rejected:* checking lazily, only after Fortran reports an error at argument *k*. Elegant
— zero cost on the happy path — but it only works for routines that validate, and one
that does not would compute silently on `INT_MIN`.

---

## Decision: numeric outputs are `np.empty`, character outputs are `np.zeros`

For a **numeric** output, `np.empty`: Fortran fills it, so zeroing is an O(n) write of
data that is immediately overwritten. It is also more honest — where Fortran fills only
part of a numeric array, zeros *look like results*, while uninitialised memory looks like
the garbage it is. On the error path nothing is returned, so a caller never sees either.

For a **character** output, `np.zeros` — and here the reasoning inverts, which is why the
rule is not uniform. Fortran fills a character buffer only *partially*:
`string_as_c_char_1d` writes the string and **one** null; `string_as_c_char_2d` fills only
as many columns as it has strings. Whatever is read back is terminated by the first null,
so any uninitialised byte past the written data would be read *as part of a string*. The
zeros are the null padding that stops that. `np.empty` here would return trailing garbage
inside the returned strings.

So: **`empty` where Fortran fills the whole buffer, `zeros` where it fills only part.**
Numeric is the former, character the latter.

---

## Decision: errors are raised by the binding language, not by C

`ierr` never reaches a caller. `check_err_code` decodes it and raises, in Python and in R
alike, so a successful call returns its results and nothing else.

The message is the Fortran documentation of the code, and because the argument position is
packed into `ierr` (`M_ERR_ARG_POS_FACTOR*arg_pos + error`) and the wrapper knows its own
argument list, the error can **name** the argument:

```
ToxInputError: invalid input arguments (argument 'n_dims')
```

Status codes (`STAT_*`) are outcomes, not failures, and never raise.

*Rejected:* raising from C in the R layer. C has no way to construct a classed R condition
without calling back into the package environment.

*Consequence:* in R, type errors are raised by R before the call and Fortran errors by R
after it — one place, one set of condition classes.

---

## Decision: shape cross-checks happen in the binding language — and in R that is R, not C

`matrix(n_rows, n_cols)` and `weights(n_cols)` share an extent. They agree **by
declaration only** — Fortran has no way to check that the actual arguments agree, and a
mismatch surfaces as a wrong answer or a segfault.

So the binding language checks it, naming both arguments:

```
ValueError: 'weights' has 3 along axis 0, but 'matrix' implies n_cols == 2
```

This is the single most valuable thing the generated layer does that a hand-written
wrapper usually forgets.

In R the check is in the **R** wrapper, not the C. Everything a user can trip over --
wrong type, wrong shape, an `NA` -- is caught in one place and raised as a **classed
condition** (`tox_type_error`, `tox_shape_error`, `tox_error`), so `tryCatch` can select
on it. C raising these would mean calling back into the package environment, which is
fragile. The split is clean: **R decides and raises; C marshals and calls.** C
re-derives the extents it needs for the call from the R objects it is handed (`Rf_length(x)`
is free), rather than being passed them.

---

## Decision: derived arguments are never parameters

Extents, shape arguments (`<arg>_shape`), mask counts (`n_selected_<arg>`) and work arrays
(`tmp_`) all come from the real inputs. Asking for them would be asking twice, and would
let the caller disagree with themselves.

```python
fx_normalize(vector)          # not fx_normalize(vector, n_dims)
```

---

## Decision: an optional with a default is required in C

From issue #131: the binding languages know the default and pass it, so C always
receives a value. Only an optional with **no** default is nullable.

This is what keeps the C wrapper flat. Every nullable optional is a branch, and the
default has to be applied somewhere regardless — applying it in Python and R keeps it in
one place and out of the wrapper.

The default must therefore be a **constant expression**, evaluable at generation time. It
is evaluated once, in the ABI layer, so Python and R cannot disagree about it.

---

## Decision: R keeps a hand-generated wrapper over the C `.Call`

A raw `.Call` entry is not a usable R function. The R wrapper over it buys three things:

1. a real signature — `span = 0.1`, not `span = NULL`
2. roxygen2 on an actual R function, plus validation in the caller's own language
3. the clean user-facing name (`fx_normalize`, not `.Call("fx_normalize_call", ...)`)

The cost is one R function call per invocation, which is nothing against the Fortran work.

---

## Decision: DM_OUTPUT_FROM(AUTO) calls the producer's generated wrapper

An argument documented `DM_OUTPUT_FROM(count, producer, module, AUTO)` is obtained by
calling `producer` and taking its `count` output, so the caller never supplies it. The
binding languages call the producer's **own generated wrapper**, not the C function --
so its validation, error checking and result handling all come for free, and the value
passed in is the consumer's already-prepared input, making the producer's conversions
no-ops rather than a second copy.

The producer's inputs are matched to the consumer's arguments **by name**. That is all the
real cases need: `mask_chunk_count(n_genes, count)` is called from a consumer that also has
`n_genes`. A producer input with no same-named consumer argument is an **error** for now; a
markdown table for renamed inputs (the same shape as the mode table) is the eventual
override, deferred until a case needs it.

Python computes it inline, since the whole wrapper is one function. R computes it in the R
wrapper and passes it to the C `.Call` shim, because the shim cannot call an R wrapper. So
an AUTO argument is a *parameter* of the `.Call` shim (R fills it) but never a parameter of
the user-facing R or Python function. Two constraints, enforced with
diagnostics: the producer must be exported (else there is no wrapper), and in the same
module (a cross-module call would need an import).

## Open

- **AUTO across modules**, and **AUTO with renamed producer inputs** (the markdown table),
  are the two deferred pieces of `DM_OUTPUT_FROM`. Both error clearly until implemented.
- **ifx** — `optional` in `bind(C)` and `implicit none (type, external)` are verified with
  gfortran only. Both are F2018; ifx is expected to agree.
