# Designing the language layers

Why the Python and R interfaces are shaped the way they are, and what was rejected on the
way. If you find yourself asking "why didn't they just...", the answer should be here. If
it isn't, that is a gap worth filling.

Companion documents: [`c-layer.md`](c-layer.md) for the Fortran/C wrappers,
[`../README.md`](../README.md) for the generator itself.

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

## Decision: validation lives in R, not in the C++ layer

The natural instinct is to validate in C++, where the data is. It is the wrong instinct.

`Rcpp::as<T>(SEXP)` receives only the SEXP. It cannot know the argument was called
`weights`, so it cannot say so:

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
and constructing one from C++ means calling back into the package environment, which is
fragile and slow.

And once R guarantees the type, `as<NumericVector>` in C++ **cannot fail**, so the whole
problem disappears rather than being worked around.

So: **R validates and marshals nothing; C++ marshals and validates nothing.**

*Rejected:* the `Tox*` C++ wrapper classes with an `as<>` interface, as sketched in issue
#131. They exist to produce good type errors, and R produces better ones for less. The
`Tox` prefix survives on the **conditions** (`tox_error`), where it earns its keep.

*Known consequence:* validation is skippable by calling the internal `.<name>_rcpp`
directly. That is acceptable for a dot-prefixed internal. If the C++ ever needs to be
defensive on its own — because something else starts calling it — this decision flips, and
the wrapper classes come back.

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

## Decision: errors are raised by the interfacing language, not by C

`ierr` never reaches a caller. `check_err_code` decodes it and raises, in Python and in R
alike, so a successful call returns its results and nothing else.

The message is the Fortran documentation of the code, and because the argument position is
packed into `ierr` (`M_ERR_ARG_POS_FACTOR*arg_pos + error`) and the wrapper knows its own
argument list, the error can **name** the argument:

```
ToxInputError: invalid input arguments (argument 'n_dims')
```

Status codes (`STAT_*`) are outcomes, not failures, and never raise.

*Rejected:* raising from C++ in the R layer. `Rcpp::stop` cannot construct a classed
condition without calling back into the package environment.

*Consequence:* in R, type errors are raised by R before the call and Fortran errors by R
after it — one place, one set of condition classes.

---

## Decision: shape cross-checks happen in the interfacing language — and in R that is R, not C++

`matrix(n_rows, n_cols)` and `weights(n_cols)` share an extent. They agree **by
declaration only** — Fortran has no way to check that the actual arguments agree, and a
mismatch surfaces as a wrong answer or a segfault.

So the interfacing language checks it, naming both arguments:

```
ValueError: 'weights' has 3 along axis 0, but 'matrix' implies n_cols == 2
```

This is the single most valuable thing the generated layer does that a hand-written
wrapper usually forgets.

In R the check is in the **R** wrapper, not the C++. Everything a user can trip over --
wrong type, wrong shape, an `NA` -- is caught in one place and raised as a **classed
condition** (`tox_type_error`, `tox_shape_error`, `tox_error`), so `tryCatch` can select
on it. C++ raising these would mean calling back into the package environment, which is
fragile. The split is clean: **R decides and raises; C++ marshals and calls.** C++
re-derives the extents it needs for the call from the R objects it is handed (`x.size()`
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

From issue #131: the interfacing languages know the default and pass it, so C always
receives a value. Only an optional with **no** default is nullable.

This is what keeps the C wrapper flat. Every nullable optional is a branch, and the
default has to be applied somewhere regardless — applying it in Python and R keeps it in
one place and out of the wrapper.

The default must therefore be a **constant expression**, evaluable at generation time. It
is evaluated once, in the ABI layer, so Python and R cannot disagree about it.

---

## Decision: R keeps a hand-generated wrapper over the Rcpp function

`// [[Rcpp::export]]` already generates an R function. The extra layer buys three things
it cannot:

1. a real signature — `span = 0.1`, not `span = NULL`
2. roxygen2 on an actual R function rather than passed through `RcppExports.R`
3. the clean user-facing name (`fx_normalize`, not `.fx_normalize_rcpp`)

The cost is one R function call per invocation, which is nothing against the Fortran work.

---

## Open

- **`DM_OUTPUT_FROM(..., AUTO)`** — calling another procedure to obtain an argument (work
  array sizes) is not implemented. It needs an argument-mapping syntax that renders as
  real prose in Ford; see the README's contract section.
- **ifx** — `optional` in `bind(C)` and `implicit none (type, external)` are verified with
  gfortran only. Both are F2018; ifx is expected to agree.
