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

## Decision: what an argument accepts -- one table, both languages

Python and R accept the same things for the same Fortran argument, and refuse the rest by name
(#221). Before, each binding accepted whatever its host's conversion did: `np.ascontiguousarray`
truncated `1.7` to 1 and wrapped `2**32 + 5` to 5 (past `DM_MAX`, which then saw 5), turned
`'no'` and `0.5` into a True mask element and `None` into False, parsed `'1.5'` as a number and
counted a `datetime64` in days since 1970; `str()` turned `None`, a list or bytes into a file
name; R's `as.integer` truncated as well, and a matrix passed for a vector was flattened.

**The rule.** For every target that is not a real: look at the input's own **kind** first, then
accept it only if converting it to the Fortran type and back gives **exactly** the same values.
The kind is an allowlist, checked before anything converts: numpy's `dtype.kind` in Python
(`i`, `u`, `f` for a number, `b` too for a logical; never `O`, `U`, `S`, `M`, `m`, `c`), and
`is.numeric` / `is.logical` / `is.character` in R -- never `typeof`, which calls a factor an
integer (its codes are not its values), and never `storage.mode`. The round trip is one
comparison that catches a fraction, a NaN, an infinity and an out-of-range value at once, and it
compares after the range is known, so the cast is never the thing that decides. A single value
is the 0-d case of the same rule.

| Fortran | Accepted | Refused, naming the argument |
|---|---|---|
| `integer(int32)` | Python and numpy integers, 0-d arrays; floats whose every value is whole and in range; R integer, and double with whole values (`2` is a double in R: 96 calls in the suites pass one); `-0.0` as 0 | a bool or R logical; a fraction, NaN/NA, an infinity, out of range; strings, objects, datetimes, complex; an R factor, Date, `integer64`; a masked array with masked elements |
| `real(real64)` | Python and numpy integers and floats; R integer and double; NaN/Inf as `DM_ALLOW_NAN` says | a bool or R logical; strings (never parsed), complex, datetimes, objects; an R Date; a masked array with masked elements |
| `logical` | bools; any number whose every value is exactly 0 or 1 (12 Python call sites pass 0/1 integers, and `c(1, 0, 1)` is R's own literal for a mask) | any other number, NaN/NA/None, strings, objects |
| a string | a Python `str`, or an `os.PathLike` whose path is a `str`; an R character vector of length 1, not NA, marked UTF-8 or latin1, or unmarked with UTF-8 bytes | bytes, None, a sequence, a number; length other than 1; NA; a NUL, a lone surrogate, invalid UTF-8; an R string marked "bytes", or latin1 with a byte that is no character |
| an array of strings | a list, tuple, ndarray or Series of `str`; an R character vector, or a factor (its labels) | a bare `str`; a set or dict, which has no order; any element that is not a `str`; NA, a NUL, a surrogate, invalid UTF-8; a value longer than a declared length |
| a mode | a valid string matching one of the modes whole (case aside, as before), after its trailing blanks are dropped like any string's (`'ward '` is `ward`) | anything else, a leading blank included; the message lists the modes |
| rank | numpy `ndim` == the declared rank; R `length(dim(x))`, where a plain vector, a 1-D array and a table are rank 1, and an all-numeric data.frame is a matrix (as a DataFrame is in Python) | any other rank: nothing is flattened |

**Trailing blanks are Fortran semantics, accepted and documented** (FES, 2026-10-26). A string
crosses blank-padded, so `'a '` and `'a'` are the same value to Fortran and the blank never
comes back; like 1-based indexing, this is what the library means, and every generated string
parameter says so (`emit.doc_contract`). A string longer than a *declared* length is refused:
the buffer would cut it short, which turned `delimiter='::'` into `':'` and a mode
`'weightedXYZ'` into `'weighted'`. A mode is therefore matched whole in the binding, which also
lets the refusal name the argument and list the modes, where the C wrapper's named neither.
Strings reach Fortran as UTF-8 in R too, and come back marked UTF-8.

**An R string crosses by what R marks it as, never by the session's locale** (FES 2026-10-29).
Marked UTF-8, it is passed as it stands; marked latin1, it is translated (so it equals its UTF-8
twin); unmarked -- what `readLines` and `rawToChar` give -- it is passed as it stands if its
bytes are well-formed UTF-8, in any locale. Everything else is refused by name: an unmarked
string that is not UTF-8, one marked "bytes", which has no encoding, and a latin1 one holding a
byte Windows-1252 (R's reading of latin1) leaves undefined. The first design translated every
string the way `enc2utf8` does, and that is locale-dependent: in the C locale (`Rscript` under
cron, a container with `LANG` unset) R takes every non-ASCII native byte for one it cannot
convert and writes it out as escape text, so `"\u00e9t\u00e9"` read from a UTF-8 file reached
Fortran as the ASCII `<c3><a9>t<c3><a9>` -- valid UTF-8, so nothing refused it -- where main had
passed the bytes unchanged. *Rejected:* translating an unmarked string from the locale's charset
and refusing it where that fails: in the C locale it would refuse every non-ASCII string read
from a file, the commonest case, which main handled. *Consequence:* in a latin1 locale, an
unmarked latin1 string is refused, with a hint to mark it (`Encoding(x) <- "latin1"`); R has run
in UTF-8 by default on Windows since 4.2, as on macOS and the Linux distributions.

Some refusals say what to do instead: a uint32 bit-mask word of 2^31 or more passes as
`.view(np.int32)`; a bytes array decodes with `.astype(str)`; an `integer64` converts with
`as.integer()` or `as.double()`; a masked array fills with `.filled()`.

Edges that are documented rather than fixed: a Python list `[True, 0, 3]` is an int64 array
before any check sees the bool in it; `-2**31` is an int32 to Python but `NA_integer_` to R,
which refuses it; a `Decimal` or `Fraction` is an object, and so refused for a real it used to
reach through `__float__`, as is a Python int too large for int64; a double above 2^53 is whole
by construction, so neither language can tell that a typed value was rounded before it arrived.

**Where each check runs.** Python's whole contract is one generated module,
`tensor_omics/_contract.py` (`emit.python_contract`), which every wrapper calls; it is free for an
input that already has the Fortran type (no pass at all) and costs one comparison where a
conversion runs anyway. R checks what an argument *is* -- kind, rank, and a single value in full,
which R itself uses as an extent -- in `tox_validate.R`; the elements of an integer, logical or
character array are checked by the C shim while it converts them (`tox_marshal.h`), in the
pass it makes anyway. That is the change to the decision below.

*Rejected:* converting first and checking the result. A post-conversion check cannot see what
was converted from, and numpy's unsafe casts and parsing accept far more than numbers.
*Rejected:* a float range test (`x.max() <= 2**63 - 1`): `float(2**63 - 1)` is `2**63`, so it
accepts a value the cast then wraps; the round trip does not. *Rejected:* refusing R doubles for
an integer, and numbers for a logical: the suites rely on both, and the languages would then
disagree on their most natural input. *Rejected:* refusing trailing blanks, which would make
fixed-width data need an explicit strip -- FES chose to document them.

---

## Decision: R checks what an argument is, the C shim checks every element

The `.Call` shim receives only the `SEXP`, and `REAL(x)` / `INTEGER(x)` hand back a bare
pointer that cannot know the argument was called `weights`. So the first design had R do all of
it, and C none: R can say which argument, for free --

```r
if (!is.double(vector)) {
  if (!is.numeric(vector)) tox_type_error("vector", "numeric", vector)
  storage.mode(vector) <- "double"     # copies only if it must
}
```

-- and constructing a classed condition in C seemed to mean calling back into the package
environment. That held while R's checks were O(1). The contract (above) needs every element of
an integer and a logical array read, and in R that costs more than the call it guards: at 1e7
elements the round trip `as.integer(x); anyNA(y) || !all(y == x)` took 137-156 ms and `x !=
trunc(x)` 380 ms, against 70 ms for the whole of `validate_gene_to_family_mapping` on native
integers. In C, in the loop that converts anyway, the same check costs what `as.integer` alone
did (46 ms; `attack-contract`).

So the split is by *what* is checked, not by language:

- **R** checks the kind and the rank of every argument, coerces a double (`storage.mode`), and
  checks a single value in full, because R uses those itself -- as extents in the shape checks,
  as inputs of a producer -- before the call.
- **The C shim** checks every element of an integer, logical or character array while it
  converts it for Fortran (`tox_int_in`, `tox_bool_in`, `tox_char_in`). The generated `.R` reads
  no element of an array.
- **Both raise the same conditions.** C is handed the argument's name as a literal, builds the
  condition R's `.tox_raise` builds -- the same classes, fields and NA message, emitted from
  `emit.errors_r` -- and signals it through base R's `stop()`, which needs nothing of the
  environment the wrappers were sourced into. It is the refusal path only, and an error longjmp
  leaks nothing: every buffer is `R_alloc`ed.

And once R has checked the kind, `REAL(x)` in the shim **cannot fail**, as before.

*Rejected:* the `Tox*` C++ wrapper classes (Rcpp `as<>`), as sketched in issue
#131. They exist to produce good type errors, and R produces better ones for less. The
`Tox` prefix survives on the **conditions** (`tox_error`), where it earns its keep.

*Rejected:* keeping the element checks in R with cheaper idioms (`min`/`max` for a 0/1
integer, 37 ms). There is none for a double that must be whole, and two places to check
elements would be two rules to keep in step.

*Known consequence:* the R-level checks are skippable by calling the internal
`.Call("<name>_call", ...)` directly. The shim refuses a value of the wrong R type rather than
reading it as another (`tox_refuse_type`), but a wrong rank or shape would reach Fortran. That
is acceptable for a dot-prefixed internal entry point.

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
| `integer` | yes | one pass in the shim, or none | a double is converted anyway, and checked in that pass; an integer vector is scanned for `NA_INTEGER` unless R already knows it holds none |

`INTEGER_NO_NA` first: it is ALTREP-aware, so for a compact sequence like `1:n` it answers
**in O(1) without touching memory**, as `anyNA()` did in R.

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
rule is not uniform. Fortran writes through a pointer view of this buffer, so it blank-fills
the whole width of every element it *assigns* — but it may assign none at all, or only the
leading elements: an early `ierr` return, or an array it fills partly. The bindings strip
trailing blanks on the way out, which turns an unwritten element into whatever bytes were
already there. The zeros make that `""`, because NumPy drops the trailing nulls of an `S`
item and R's `tox_char_alloc` blank-fills for the same reason. `np.empty` here would return
trailing garbage inside the returned strings.

Not `np.full(..., b" ")`, which looks like the obvious spelling and is not: NumPy fills each
item with one blank and nulls the rest.

So: **`empty` where Fortran fills the whole buffer, `zeros` where it fills only part.**
Numeric is the former, character the latter.

---

## Decision: a returned array is read-only in Python, and unfrozen in R

A result is a value, not a scratch buffer. Python marks every array it allocated
`flags.writeable = False` after the call, so modifying one in place raises where the mistake
is rather than wherever the stale array is read next. It costs O(1) and no copy, and the flag
propagates into the views the return builds -- the reshape of a serialized output, the slice
of a partially filled one -- which then cannot be unfrozen at all.

Freezing happens **after** the call: `ctypes` is handed the raw data pointer, so a flag set
earlier would not stop Fortran writing, only misdescribe it.

Three things stay writeable, because there is nothing to protect: an `intent(inout)` array is
the *caller's*, modified in place, and freezing it would break exactly what it is for; scalars
come back as immutable Python numbers; a character buffer is never returned, only the strings
decoded out of it, which are fresh per call.

**R gets no counterpart, and that is not an omission.** R has no read-only vector, and it does
not need one: copy-on-modify means a returned vector cannot be aliased by anything else, so the
class of bug this guards against does not exist there. The three near-misses were all rejected.
`lockBinding` locks a *name in an environment*, not a value. `MARK_NOT_MUTABLE` in the shim
makes later modifications copy **silently** -- a hidden allocation rather than an error, which
teaches the caller nothing. An S3 class with `[<-` methods that `stop()` would raise, but costs
`Ops`, printing and every idiom besides. Same shape as `intent(inout)` (in-place in Python,
copy-return in R): the idiomatic answer in each host, not one convention forced on both.

*Consequence:* `result += 1`, `np.add(..., out=result)` and anything demanding a writeable
buffer (`torch.from_numpy`) now fail on a result. `.copy()` is the escape hatch, and the
generated `Returns` section names it.

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

*Rejected:* decoding `ierr` in the R shim's C. It would mean raising the code's own condition
class from C, from a table the R error module already holds; the shim raises only the binding's
own refusals of an element (above), which need no code at all.

*Consequence:* in R, a refusal is raised before the call -- by R, or by the shim for an array
element -- and a Fortran error by R after it, all of one set of condition classes.

**A refusal of the binding's own carries no tox code, and is an input error all the same**
(#221). A wrong shape, type or value is found before anything is called, so there is no code to
report, and main's convention is to report none: Python's `assert_error` without a code and R's
mean exactly that. But the caller made the same mistake as one the library reports itself --
the C wrapper's `ERR_SIZE_MISMATCH` for an axis only Fortran can evaluate (`c-layer.md`) -- and
had to catch the two apart: a `ValueError` in one case, a `ToxInputError` in the other. So Python raises `ToxArgumentError(ToxInputError, TypeError,
ValueError)` with `code` None: an `except ToxInputError` catches both, and an `except ValueError`
written against the old refusal (the numpy idiom) still works. R's `tox_type_error`,
`tox_shape_error` and `tox_na_error` gained `tox_input_error` between themselves and `tox_error`,
with `code` NA, which R's `assert_error` already reads as "no tox code".

*Rejected:* giving the refusal a code (`ERR_SIZE_MISMATCH` for a shape). It would claim the
library reported something it was never asked, and turn six Python tests that assert "no tox
code" into failures for no gain.

**The package exports the error module whole** (#221 C, FES). `from tensor_omics import
ERR_NAN_INF, ToxInputError` works, as do every `ERR_*` and `STAT_*` code, `ARG_POS_FACTOR`, every
exception class and `check_err_code`, all listed in `__all__`: Python is where R already was, whose
`error_handling.R` is sourced into the global environment. Before, a caller who wanted to compare
a code or catch by class had to know the submodule (`tensor_omics.error_handling`), while the
package itself offered `ToxError` and `check_err_code` alone. The names come from one helper,
`errors_python.public_names`, that the error module is emitted from too, so #220's `STAT_*` codes
follow with no edit. A name the top level would then bind twice -- a procedure or module called
like one of these or like a submodule, ignoring case -- is a generator error
(`python_ctypes.check_package_names`).

*Rejected:* exporting them without listing them in `__all__`. They would be importable and still
invisible to `import *`, documentation tools and linters: a half state. `import *` gains only
`ERR_`, `STAT_` and `Tox` names, and `ARG_POS_FACTOR` and `check_err_code`, which is acceptable.

---

## Decision: shape checks happen in the binding language — and in R that is R, not C

`matrix(n_rows, n_cols)` and `weights(n_cols)` share an extent. They agree **by
declaration only** — Fortran has no way to check that the actual arguments agree, and a
mismatch surfaces as a wrong answer or a read past the buffer.

So the binding language checks it, by one rule (#221, `emit.shape_checks`): **every input axis
that names an extent is compared with that extent's value, except the axis the value was read
from, and an axis fixed at a constant with the constant.** The value is read off the first input
that must be passed, counted from a mask, computed by a producer or passed by the caller; an
optional array is compared when given. It used to be "two inputs sharing an extent agree", which
left out an optional array (a 1x1 `neighbor_mask` read past), an input sized by a mask's count
(compared with nothing) and a constant axis (`fields(n, 2, m)` read a field without its target).
One module decides and both languages write it out, and Python reads an extent off the very axis
the rule leaves out (`extent_provider`), so they cannot drift. An axis the declaration computes
(`sum(reps_per_tissue)`) is not theirs: they pass its extent to C, and the C wrapper compares it
(`c-layer.md`). The message names the arrays and shows the declared shape:

```
ToxArgumentError: 'weights' has 3 along axis 0, but 'matrix' implies n_cols == 2; 'weights' is declared with shape (n_cols)
```

This is the single most valuable thing the generated layer does that a hand-written
wrapper usually forgets.

In R the shape check is in the **R** wrapper, not the C: it needs no element, only extents,
which R has for free, and it is raised as a **classed condition** (`tox_type_error`,
`tox_shape_error`, `tox_input_error`, `tox_error`), so `tryCatch` can select on it. The shim
raises the same classes for an element it refuses (above). C re-derives the extents it needs for
the call from the R objects it is handed (`Rf_length(x)` is free), rather than being passed them.
A mask that sizes an input is counted in R as the shim counts it, its exact TRUEs and 1s, since
its elements are checked only later, in C.

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

It is evaluated with **Fortran's arithmetic**, because the value a binding passes has to be the
one a Fortran caller gets when it omits the argument. Python's arithmetic differs in exactly the
places a default reaches (#221, B9): `round` rounds a half to even (`nint(2.5)` gave 2, Fortran
3), `**` with a negative exponent gives a float (`2**-1` gave 0.5, Fortran's integer result is
0), `int(a / b)` goes through a float and loses integers above 2^53, a negative base raised to a
real power gives a complex number where Fortran refuses, and an unbounded `int` never
overflows. So every integer carries its kind (`constants.Integer`), overflowing it is an error,
division and `mod` truncate exactly, and `nint` compares the exact distance to the floor with a
half rather than computing `floor(x + 0.5)` -- whose addition rounds `0.49999999999999994` to 1.
The evaluated default must also fit the argument's own kind: `ctypes.c_int` wraps an
out-of-range value silently, so a default of 3000000000 on an `integer(int32)` is refused.
A sign binds as Fortran's grammar says, to a whole add-operand below `*` and `**`: `-2**2` is
-4, where the first grammar computed `(-2)**2` and handed the bindings +4.

Real kinds are **not** tracked: every real literal is read as real64. An unkinded `0.95` is a
real32 0.949999988 to the compiler and 0.95 here. Reading it the compiler's way would change
main's `DM_DEFAULT(0.95)` (`tox_data_integration_jsd_impl`) away from the `0.95_real64` the
implementation itself falls back to, and refusing unkinded reals would break the build; so the
rule is documented instead -- write `_real64`, and the two readings agree.

---

## Decision: an extent the binding computes is parsed once, in a grammar every layer agrees on

Python and the R shim allocate outputs and work arrays, and pad strings, before the call, so they
evaluate the Fortran extent themselves (`n_points - 1`, `max(0_int32, n_timepoints - 1)`,
`size(gene_to_fam)`). Two regex translators did it, one per language, and each passed through
what it did not recognise (#221 B11): `/` reached Python as float division, `mod` and `**` reached
the shim's C, which did not compile, and the shims are bundled into the library.

`abi/extents.py` parses an extent into a tree, and each binding writes it out from that tree. The
grammar is what means the same number in Fortran, Python and C -- names, module constants, integer
literals with or without a kind, `+ - * /`, parentheses, `max`, `min`, `size(a[, k])` -- and
anything else is refused before emission (`check_binding_extents`), pointing at `DM_OUTPUT_FROM`.
Integer `/` is the one operator that needs writing out: Fortran and C truncate toward zero,
Python's `//` floors, so `-7/2` is -3 against -4. Python gets the exact
`(a // b if (a >= 0) == (b > 0) else -(-a // b))`, never `int(a / b)`, which goes through a float.

*Rejected:* refusing `/` along with `mod`, which the first design proposed. Writing it out is one
line, and a size that halves is common enough to want. *Rejected:* a generated "extent oracle"
procedure per wrapper, which would let Fortran compute every size: a second call and a new symbol
per procedure, for a grammar no output on main leaves.

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

The producer's inputs are matched to the consumer's arguments **by name**, which is what most
cases need: `mask_chunk_count(n_genes, count)` is called from a consumer that also has
`n_genes`. Where the two spell a quantity differently, or where a producer input is a
*constant* the consumer has no argument for, the consumer's argument carries a markdown table
mapping them (the same shape as the mode table) -- `tox_loess_required_workspace` is supplied
`n_dim = 1_int32` and `save_factorization = .false.` that way. A producer input that is neither
name-matched nor in the table is an **error** naming it.

Python computes it inline, since the whole wrapper is one function. R computes it in the R
wrapper and passes it to the C `.Call` shim, because the shim cannot call an R wrapper. So
an AUTO argument is a *parameter* of the `.Call` shim (R fills it) but never a parameter of
the user-facing R or Python function. The producer must be exported, else there is no wrapper
to call; it may live in **any** module, and Python imports it inside the calling function
rather than at the top of the file, because two modules may size each other's outputs and a
module-level import would then be circular.

## Decision: a default is rendered in the vocabulary of the layer that reads it

`DM_DEFAULT(MODE_ROBUST)` quotes the author's Fortran back into the documentation, and in
Fortran a mode is an integer. Everywhere a binding reads that documentation the mode is a
*string*, so the sentence named the one value such a caller must not pass — a Python
docstring said `mode : str, one of 'plain' | 'robust', optional, default 'robust'` and then,
two lines down, ``The default value is `1`.``

The fix is not "render the mode string" but the more general rule that already governs
`doc_literals`: **a literal is written the way the layer that reads it writes it.** So the
default now reads `'robust'` in the Python docstring, `"robust"` in the roxygen block — each
matching the quoting its own type line already uses for the accepted modes — and `'robust'`
in the C wrapper, whose dummy is `character(len=1, kind=c_char)`.

**The generated Fortran wrapper is deliberately left alone.** Its dummy really is
`integer(int32)`, so `1_int32` is the correct thing to print there, and rewriting it would be
the same bug pointing the other way. That is what makes the negative control in
`test_generate` worth as much as the three positive ones. The split falls out of the layering
rather than being asserted per emitter: the three string layers all work from `CArgument`,
which knows the mode table and carries the translated default, while the Fortran wrapper
works from the IR `Argument`, which does not.

`render_mode_default` anchors on the whole sentence `DM_DEFAULT` expands to, which is what
keeps it off the author's own prose — the wording is the macro's, so rewriting the literal
inside it is the generator editing its own output. **Rejected — rewriting the `Doc` block
upstream in `c_abi`**, which would have fixed all three layers in one place: it also fixes
them to one spelling, and the whole reason `doc_literals` exists is that `.true.` is `True`
in Python and `TRUE` in R.

**The gap that let this ship** was in the fixtures, not the emitters: `fx_modes` had two mode
arguments and neither was optional, so "a mode with a default" was generated nowhere in the
test suite. `link_method` now carries `DM_DEFAULT(METHOD_WARD)`.

---

## Decision: an argument's name must work in every layer, and is refused where it would not

A dummy's name is published verbatim: a Python parameter, an R parameter, a parameter of the C
`.Call` shim and a dummy of every generated wrapper. Nothing checked it, and a name one layer
reserves broke more than its own procedure -- a parameter called `lambda` is a SyntaxError, and the
package `__init__` imports every module, so `import tensor_omics` failed; `next` failed
`r/load_tensor_omics.R` for every wrapper; `default` failed the C build of the whole library; `np`
hid numpy inside the function. So `validate` refuses the words each layer reserves
(`names.language_word`), for an argument and for a procedure's own name, which is a Python and an
R function.

The names the generated code binds itself are a different problem, and are designed out where the
language lets them be. Everything the R shim's C and a generated Python module bind beside the
arguments starts with `_` (`_v_<name>`, `_nprot`, `_np`, `_check_err_code`, `_error`); no Fortran
name can, so no argument can meet one. The functions the shim's header defines (`tox_int_in`,
`tox_char_in`, ...) cannot be underscored -- a leading `_` at file scope is reserved to C's
implementation -- so they are listed (`names.R_SHIM_HELPERS`, held to the emitted header by a
test) and an argument of one of their names is refused: as a shim parameter it would hide the
helper, and the library would not build. A Fortran wrapper's own locals cannot start with `_`, so
they stay `<name>_f`, `<name>_mode_f`, `<name>_mode_p`, `<name>_address`, `<name>_value` and
`handled` (and the C wrapper module's `empty_view_anchor`), and an argument is refused where the
wrapper emits the same name -- read off the wrapper as emitted, together with what it calls,
imports, uses as a kind and finds declared in its module (`emit/name_clashes.py`).

A procedure's or a module's name meets one more set of names at the Python package's top level:
the error module's exports and the package's submodules, beside which `tensor_omics` binds it. A
clash there is refused too, ignoring case (`python_ctypes.check_package_names`; the errors
decision above).

**Rejected -- renaming instead.** The published name is the API, and a binding that quietly calls
an argument `lambda_` is one more thing a reader of the Fortran cannot predict.

**Rejected -- refusing Fortran keywords and R function names too.** Fortran reserves no word, and
an R argument cannot shadow a function call (R skips non-function bindings when it looks a call
up), while main publishes `value`, `mode`, `vector`, `array`, `expression`, `factor`, `points` and
`weights`. Refusing them would refuse working API for a clash that does not happen.

**Rejected -- refusing every `<other>_<suffix>` and every name a wrapper might use** (WP2's first
form). It refused `p_value` beside `p` and `window_size` beside `window`, ordinary names whose local
is never emitted, and a blanket list of intrinsics would refuse `count`, which main publishes in
`mask_chunk_count`. Refusing where the wrapper really emits the name is both narrower and complete.

**Rejected -- a list of what each emitter might reference**, kept beside the emitter. It would be a
second description of the emitted code that drifts from the first; reading the rendered wrapper
cannot.

---

## Open

- **ifx** — `optional` in `bind(C)` and `implicit none (type, external)` are verified with
  gfortran only. Both are F2018; ifx is expected to agree.
