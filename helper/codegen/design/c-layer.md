# Designing the C wrapper layer

Why the generated Fortran `bind(C)` wrappers look the way they do, and what was rejected.

Companion documents: [`language-layers.md`](language-layers.md) for Python and R,
[`../README.md`](../README.md) for the generator itself.

---

## The rule the wrappers follow

**A wrapper presents a plain C ABI — pointers and nothing clever — validates what it can
without dereferencing anything it should not, converts only what C cannot express, and
reports through `ierr`.**

Everything below follows from that.

---

## The shape of a wrapper

```fortran
subroutine fx_serialized_c(data, data_shape, n_data_shape_elements, ierr) &
        bind(C, name="fx_serialized_c")
    use fx_edges, only: fx_serialized

    integer(c_int), intent(in), target :: n_data_shape_elements
    real(c_double), dimension(*), intent(in), target :: data
    integer(c_int), dimension(n_data_shape_elements), intent(in), target :: data_shape
    integer(c_int), intent(out), target :: ierr

    M_CHECK_IERR_NON_NULL
    call set_ok(ierr)
    M_CHECK_NON_NULL(n_data_shape_elements)
    M_CHECK_ARRAY_NON_NULL(data_shape, n_data_shape_elements)
    M_CHECK_ARRAY_NON_NULL(data, product(data_shape))

    call fx_serialized(data = data(1:product(data_shape)), ...)
end subroutine
```

In order: check `ierr` itself, set it ok, null-check in an order `c_loc` tolerates,
convert, call, convert back.

---

## Decision: the wrappers live in their own `<module>_c` modules

Not in the source files they wrap. From issue #131:

1. `f42_safeguard` is needed only by C-facing code, so only these modules use it
2. the whole C binding can be compiled out with one directive — the wrappers are wrapped
   in `#ifndef NO_C_BINDING`, and with them goes the `f42_safeguard` dependency
3. the source files are not bloated with wrapper definitions
4. the generator never touches hand-written source

---

## Decision: the R C shims are bundled into the library, omittable like the C binding

The R binding is pure C (`.Call`) shims — generated into `src/generated/bindings/r/*.c` — that
marshal R objects and call the `bind(C)` wrappers. They are **compiled by fpm into the one
`libtensor-omics.so`** (fpm already scans `src/` recursively and rebuilds only changed
files), so there is a single artifact and no separate R build step; the R loader just
`dyn.load`s it and `.Call`s the entry points by name (resolved by dynamic symbol lookup —
the `.so`'s hyphenated name means the registration `R_init_*` never auto-fires, which is
fine).

The shims are guarded by `#if !defined(NO_R_BINDING) && !defined(NO_C_BINDING)`, mirroring
the Fortran wrappers' `#ifndef NO_C_BINDING`. So:

- `./build.sh --directive=NO_R_BINDING` drops the R layer (the shims compile to empty
  objects that need no R headers) while keeping the C ABI for Python/direct-C use;
- `NO_C_BINDING` implies no R layer too — the shims call the `bind(C)` symbols, which are
  gone — so the guard tests both;
- if the R layer is wanted but R is not installed, the build **auto-disables it with a
  warning** (it needs R's headers, via `R CMD config --cppflags`, only when included).

A third switch, `NO_INPUT_VALIDATION`, drops the *generated wrappers'* input checks rather than
any binding layer; it is the same kind of whole-build decision and is designed in
[`impl-layer.md`](impl-layer.md). It deliberately leaves the null checks here alone: those
prevent a segfault rather than reject a bad value.

The R headers reach the C compiler through fpm's `--c-flag` (its `--flag` is Fortran only);
the `.so` links with R's symbols left undefined, resolved when R loads it, so it carries no
`libR` dependency.

**The same `.so` is loaded by Python too** (`ctypes`), where there is no R to resolve those
symbols. So *every* R API symbol the shims reference is marked `#pragma weak` (in
`tox_marshal.h` and `init.c`, listed in `c_call._WEAK_R_SYMBOLS`): the ordinary eager load
resolves them to null — Python never calls the R code that would use them — while a
*genuinely* missing symbol (a build regression) still fails loudly at load. Weakening only
the *data* symbols and loading `RTLD_LAZY` was rejected: it works, but defers real
missing-symbol errors to first call. Adding an R symbol to the shims means adding it to the
list; `test_ctypes_loads_the_r_bundled_library` catches an omission (its eager load fails).

*Build caveat:* fpm content-hashes C sources and does not track their dependency on
`tox_marshal.h`, so a change to that header alone will not trigger a recompile of the
`.c` shims — force a clean build (a `--directive=...` does, or clear `build/<compiler>_*`).

---

## Decision: the null-check order is `ierr`, scalars, arrays

This is the fix for the standing `TODO codegen` in `src/macros.h`.

`c_loc` **may not be given a zero-size target**. A wrapper that null-checks
`M_CHECK_NON_NULL(array)` before the extent sizing that array is known to be readable is
not standard-conforming — and a caller passing a legitimately empty array hits exactly
that.

So:

1. **`ierr` first.** It is the only channel for reporting anything. If its own pointer is
   null there is no way to signal failure at all, so `M_CHECK_IERR_NON_NULL` returns
   silently rather than dereferencing it to say so.
2. **`call set_ok(ierr)`.** A procedure that declares no `ierr` of its own will never set
   it, and the wrapper would return whatever was on the stack.
3. **Every scalar**, which includes every extent and every string length.
4. **The arrays**, each guarded by the element count that step 3 has just made safe to
   read:

   ```fortran
   #define M_CHECK_ARRAY_NON_NULL(ARG, N) if ((N) > 0) then; M_CHECK_NON_NULL(ARG); end if
   ```

Within step 4, an array whose shape travels in a separate argument waits for that
argument — its count is `product(data_shape)`, which means reading it.

**Consequence, and it is deliberate:** an empty array passes through untouched. The
callee's own `validate_dimension_size` decides whether empty is an error for that routine,
which is where that policy already lives.

**This is why extents and shape arguments may never be optional** — the wrapper reads them
before it may take `c_loc` of the arrays they size. `ir/validate.py` enforces it.

---

## Decision: an argument position always numbers the Fortran dummy list

`ierr` packs the position of the argument it blames. The C wrapper's argument list is *not*
the Fortran one — it inserts an `n_<name>_elements[_dim_k]` per assumed-shape extent and a
`<name>_strlen` per `len=*` character, so `compute_rdi_expert_c` carries `n_dscale_elements`
at C position 5 that the Fortran procedure does not have. The obvious conclusion is that the C
layer should renumber. It should not.

**The position numbers the wrapped Fortran procedure's dummy arguments, at every layer. The
size and string-length arguments the C wrapper inserts are not counted. A direct C caller reads
positions against the Fortran signature.**

Renumbering in C would break the two languages that consume the position. Python and R both
build their argument-name tuple from the *Fortran* procedure (`emit/python_ctypes.py`,
`emit/r_wrapper.py`) and index it with `arg_pos` when they turn an `ierr` into an exception or
a condition — `_COMPUTE_RDI_EXPERT_ARGUMENTS` is the ten Fortran names, with no
`n_dscale_elements` among them. Renumber, and every position from the insertion point onward
names the wrong argument in both. Renumber *and* fix both emitters, and their error messages
start naming arguments no Python or R caller ever passes.

The C layer packs no position of its own, so there is nothing for it to reconcile: neither
`M_CHECK_NON_NULL` nor the unmatched-mode `case default` supplies `arg_pos`. And it never drops
a Fortran argument it can map — an argument it cannot map is a hard generator diagnostic, not a
silently renumbered signature.

The *label* a binding prints is translated where it can be: an `arg_pos` naming an argument the
caller never wrote — an extent read off an array, a `<x>_shape`, an `n_selected_<x>` — is
reported as `(argument 'vec1', via 'n_elements')`, the caller's word first and the Fortran one
after. The position itself is untouched; only the name shown changes, and an argument with no
caller-visible source keeps its Fortran name alone.

`map_err_arg_pos` in `tox_errors` stays for hand-written code where two dummy lists really are
known to correspond. What a *generated* wrapper does with an inherited position is the opposite
and is decided elsewhere: see "The argument position an `ierr` carries" in
[`impl-layer.md`](impl-layer.md).

## Decision: nullable optionals are `OPTIONAL`, not pointers

An `OPTIONAL` dummy of a `bind(C)` procedure is absent exactly when C passes a null
pointer (TS 29113, folded into F2018). The C prototype stays a plain pointer:

```c
void fx_nullable_c(..., const int *ortholog_set, ...);
```

And the wrapper hands the argument **straight to the callee**, because an optional
associated with an absent optional is itself absent. Presence propagates with no branch:

```fortran
call fx_nullable(ortholog_set = ortholog_set, ...)
```

*Rejected:* a nested if-else chain over every combination of present optionals — O(2ⁿ),
and what the previous generator attempted.

*Rejected:* a `POINTER` dummy plus `c_f_pointer`. F2018 does permit a pointer dummy in
`bind(C)`, but it is passed as a **CFI descriptor** — the prototype would become
`CFI_cdesc_t*` and the plain-pointer ABI would be gone. It would have compiled and been
wrong.

*Note:* `OPTIONAL` was not interoperable before TS 29113 (2012). On an F2008 draft this
design would not have been available.

---

## Decision: logicals cross as `c_bool`

C passes a real `bool`, so Python and R hand over a boolean rather than a 0/1 integer.

**And since 2026-08-14 the sources are `logical(c_bool)` too**, so there is no copy at all:
the wrapper hands the caller's buffer straight to the implementation. It used to declare an
automatic array of default `logical` and convert elementwise, because `c_bool` is one byte and
a default `logical` is four. 31 such array temporaries across 16 wrappers are now one — a
scalar, and only because `mask_check_state` is a `logical function` whose *result* stays
default-kind, where there is no memory to win and a condition accepts any kind anyway.

The decision was memory, not speed: no boundary buffer, and a quarter of the bytes for every
mask in the library. `misc/bench/logical-kinds/RESULTS.md` has the measurements, and they point the same way
— at the sizes where the memory argument bites, `c_bool` is *faster* for the same reason it is
smaller.

*Rejected:* passing `integer(c_int)` 0/1 and converting with `c_int_as_logical`, as the
previous generator did. It makes every caller in every language do bool→int→logical.

**The lesson the switch taught, which cost four defects to learn: "does Fortran need a copy"
and "does R need a copy" are different questions, and the code had one flag for both.**
`CArgument.needs_conversion` answers the first, and it goes false the moment a dummy is
`logical(c_bool)`. The R emitter had keyed four sites on it — because while every logical was
default-kind the two questions had the same answer. R's logical vector is `LGLSXP`, an array
of 4-byte `int`, whatever kind the Fortran declares, so R *always* rebuilds the bytes. Three
sites then emitted `REAL()` on a logical SEXP, which R rejects at run time, and one crashed the
generator outright. `c_call._r_marshals` is now the single definition all of them ask.

Adding `LGLSXP` to the `R_SEXPTYPE` table would have silenced the crash and been wrong: the
pointer path would then have handed `LOGICAL(sexp)` — an `int*` — to an `unsigned char*`
parameter, and corrupted quietly instead of failing loudly.

Two smaller ones from the same change, both in the Fortran wrapper emitter: kind constants have
to be imported from the module that owns them (`c_bool` is `iso_c_binding`, not
`iso_fortran_env`), and a literal the *emitter itself* writes needs its kind too —
`allow_nan=.true.` against a `logical(c_bool)` dummy is a compile error, because argument
association does not convert kinds. The second means the import list cannot be built from the
arguments alone; `_emitted_literal_kinds` exists for the kinds no argument carries.

---

## Decision: a converted local is allocated, never automatic

A wrapper that converts an argument needs somewhere to convert into, and the obvious spelling
is an automatic array:

```fortran
logical, dimension(n_vectors) :: vectors_selection_mask_f   ! what this used to emit
```

That is stack storage sized by whatever C passed. At 10 000 000 elements a logical mask is
40 MB, and **ifx segfaults on exactly this construct** against a default 8 MB stack; gfortran
survives only because it quietly rehouses large automatics on the heap. It is reachable from
the published API, not a theoretical size: `serialize_logical` takes `arr(n_elements)` of the
caller's choosing, and `is_outlier(n_genes)` is an output copy of the same shape.

Every such local is now `allocatable` and reaches the heap through `M_ALLOCATE`, which is the
same macro the generated Fortran wrappers already use for their work arrays. That fixes the
storage on every compiler, and it turns an allocation failure into `ERR_ALLOC_FAIL` returned
to the caller instead of a crash inside a wrapper whose whole job is to report errors. The
bare `allocate(...)` calls this replaced had no `stat=` at all, so they aborted too.

A scalar `logical` local stays automatic: four bytes, and not a size the caller chooses. A
character local is not allocated at all — it is a pointer view of the caller's buffer, see
*a fixed-length character view* below — so as of that change the real project's C layer
contains no `M_ALLOCATE` at all, and `test_generate.py` asserts exactly that.

*Rejected:* `-heap-arrays` in the ifx profile. It fixes the crash, but per compiler, by a flag
anyone can drop, and it leaves the failure unchecked — an allocation that fails still aborts
rather than returning `ERR_ALLOC_FAIL`. The allocation belongs in the emitted code, where it
can be reasoned about, not in a build flag.

The measured cost of the copies these locals exist for is in `misc/bench/logical-kinds/RESULTS.md`; it is
small, and the stack is the reason this changed, not the speed.

---

## Decision: a string array crosses as a fixed-width char matrix, not `char**`

A Fortran `character(len=n), dimension(m)` is **one contiguous block** of `n*m` bytes. C's
idiomatic `char*[]` is an array of *pointers* to scattered blocks. There is no zero-copy
bridge between those two memory models, so the wire format is the shape both languages can
express contiguously: a column-major `strlen x n_strings` matrix of
`character(len=1, kind=c_char)`, one column per string, plus `strlen` as its own argument
because a `char*` carries no length.

Each binding pays a different price for that, and neither pays much. NumPy's fixed-width `S`
dtype *is* the layout, byte for byte, so Python passes `arr.itemsize` as the `strlen` and
copies nothing. R's `STRSXP` is an array of pointers to scattered `CHARSXP`s, so `tox_char_in`
genuinely packs it into a `len x n` buffer.

The cost is padding to the longest string, which is why a ragged array of mostly-short IDs
carries the longest one's width throughout. That is inherent to the format; the alternative is
an offsets array, which is a different and much larger design.

**Fixed-width padded blocks are the normal C representation for bulk string data**, not an
oddity of this project: HDF5 fixed-length string datatypes, Arrow's `FixedSizeBinary`, tar and
FITS headers, `sockaddr_un.sun_path`, and NumPy's own `S`/`U` dtypes are all exactly this. For
*arrays* of strings it is essentially the only zero-copy option in C too.

---

## Decision: a string crosses as a fixed-length pointer view, blank-padded, and is never converted

### Multibyte characters: the width is right, the encoding flag is not

Checked 2026-08-14, because a fixed-width format truncates the moment the width is measured in
different units from the payload.

**It is not.** Both bindings measure bytes, which is the same unit the buffer is filled in, so a
UTF-8 string is never truncated by the sizing. R's `tox_max_strlen` uses `LENGTH()` on a
`CHARSXP`, which is the byte length and not the character count (`nchar()` would have been the
bug). Python builds the buffer with `.encode()` — UTF-8 bytes — into a NumPy `S` dtype whose
`itemsize` is likewise bytes. Fortran agrees: `character(len=n)` of kind `c_char` is n bytes.

**What is lost is the encoding flag, not the content.** `tox_char_out` rebuilds R strings with
`Rf_mkCharLen`, which marks the result native/unknown. A UTF-8 string handed in comes back
byte-identical but with its `Encoding()` reset, which is harmless in a UTF-8 locale and wrong
outside one. `Rf_mkCharLenCE(p, m, CE_UTF8)` is the fix, and it needs a decision rather than a
patch: the shim would be asserting an encoding for bytes Fortran may have built itself, so the
honest version marks UTF-8 only where the input was UTF-8, or the API declares that it always
is.

**For the pointer view specifically, blank padding is safe.** No UTF-8 continuation byte is
`0x20` — they are all `0x80`-`0xBF` — so `trim` and `len_trim` cannot split a character, and
neither could the NUL scan they replaced. What can split one is any *byte-index* slice of a
`character(len=n)` (`s(1:k)`), but that is inherent to Fortran character handling and no worse
under a pointer view than it was before.


The wrapper used to convert the matrix into `character(len=:), allocatable, dimension(:)` with
`c_char_2d_as_string`, which scanned each column for `c_null_char` and copied. That conversion
is now a pointer remap costing nothing: `character(len=n), dimension(m)` and
`character(len=1), dimension(n, m)` have **identical layout** -- contiguous, column-major,
string *i* at bytes `[(i-1)n+1 .. i*n]` -- so the view is pure reinterpretation.

Probed on 2026-08-14 against gfortran 16.1 and ifx 2026.1, at default, `-std=f2018`/`-stand
f18` and `-std=f2023`/`-stand f23`:

| form | gfortran | ifx |
|---|---|---|
| `c_f_pointer` to `character(len=n), pointer` (scalar) | OK, no diagnostics | OK, no diagnostics |
| `c_f_pointer` to `character(len=n), pointer :: v(:)` from an assumed-size `dimension(strlen,*)` dummy, `intent(in)` and `intent(out)` | OK, no diagnostics | OK, no diagnostics |
| `c_f_pointer` to `character(len=:), pointer` (deferred) | compiles, **empty string** | compiles, **SIGSEGV** |
| `c_f_strpointer` (F2023) | **absent from `iso_c_binding`** | works, but warns #8893 even at `-stand f23` |

**Never remap onto a deferred-length pointer.** It compiles clean on both compilers at every
standard level, with no warning, and then produces two different wrong answers. If this shape
is ever written it must be caught by review, because no compiler will catch it.

`c_f_strpointer` is the route F2023 added for exactly this, and it is not usable: the current
gfortran release does not have it at all, and ifx does not consider it standard even under its
own F2023 flag.

**On conformance.** By the letter of F2018 a `character(len=n)` scalar with `n > 1` is not an
interoperable *type*, so `c_f_pointer` onto one is arguably outside the standard. The reason is
not NUL-termination -- C has no string type at all, and `char[n]` and `character(len=n)` have
exactly the same layout. It is that C's type system has no length parameter, so the
correspondence is one Fortran object to *n* C objects, and interoperability is defined between
types. That is why the len>1 allowance exists only as a **dummy-argument** rule for `bind(C)`
procedures, which `c_f_pointer` cannot reach. There is no ABI question here, only a formal one,
which is why every compiler does the obvious thing.

The fully conforming fallback, if that ever matters, is for the implementation to take the 2-D
`character(len=1)` buffer directly -- no pointer, no remap, and no conversion either, at the
cost of clumsier string handling inside the kernel.

### What adopting it required

**The wire format is blank-padded, both ways.** The conversion trimmed at the first NUL; a
remap yields the full padded width, so whatever pads a short string *is* part of the value.
Fortran's own convention is *blank* padding -- `trim`, `len_trim` and string comparison all
work on blanks, not NULs -- so the bindings blank-pad instead: `memset(buf, ' ', total)` in
`tox_char_in` and `tox_char_alloc`, `bytes.ljust(width)` on the NumPy side, and
`tox_char_out` / `.rstrip(' ')` trimming blanks rather than scanning for a NUL. The trade is
that a value genuinely ending in blanks stops round-tripping, which for gene and family IDs
is no loss.

**For arrays that is behaviour-preserving.** `c_char_2d_as_string` already allocated
`character(len=strlen) :: str_out(n)` and assigned each string into it, which blank-pads to
the full width. So implementations were already receiving blank-padded fixed-width arrays of
exactly this width, and every consumer already trims. For *scalars* the length changes from
trimmed to blank-padded, which is safe because 17 of the 21 boundary character arguments are
filenames and `FILE=` removes trailing blanks (F2018 12.5.6.10) -- safer padded than
NUL-terminated -- and because nothing in `src/` ever reads `len()` of a boundary scalar.

**The published C ABI is byte-identical and its semantics are not.** The symbol set and the
prototypes do not move, but a direct C caller must now blank-pad rather than NUL-terminate:
`char buf[16] = "genes.tsv"` used to work and now yields a filename with embedded NULs, and a
`method` buffer holding `"ward\0\0\0\0"` now matches no `case` and returns
`ERR_INVALID_INPUT`. Nothing in-tree does that -- Python and R are the only callers -- but it
is the one thing `NO_R_BINDING`'s direct-C audience has to be told.

**A nullable optional needs a target and an explicit `nullify`.** `c_loc` requires TARGET,
which an optional was deliberately not given; the eight optional characters of
`save_tox_data_c` now carry `optional, target`, which is legal on a `bind(C)` dummy and costs
the ABI nothing. The absent case is a *disassociated* pointer, which F2018 15.5.2.12 makes an
absent optional dummy -- the same rule the unallocated allocatable used. It must be a
`nullify` statement and never `=> null()`: an initialiser gives the local an implicit SAVE,
which is both thread-unsafe and formally impossible on an automatic-length object, and both
compilers accept the spelling silently.

**Zero-size buffers are fine, and were probed rather than guarded.** `M_CHECK_ARRAY_NON_NULL`
lets an empty array through unchecked, so the remap is the first place the wrapper takes an
unguarded `c_loc` -- and `read_tox_data_into` on a partial archive really does reach it with
`family_id_len == 0, n_family_ids == 0`, where R's `R_alloc(0, 1)` hands back NULL. Probed
2026-08-14: `strlen=0, n=0` and `strlen=0, n=2`, both over a NULL buffer, give a well-behaved
zero-size/zero-length view on gfortran 16.1 and ifx 2026.1 under `-fcheck=all` / `-check all`.
A guard would not have helped anyway: a non-optional dummy cannot receive a disassociated
pointer. This is the same class of formal-only objection as the `len>1` point above.

**Only `Conversion.MODE` still converts.** It is a fixed-width `select case` lookup, and
`c_char_1d_as_string` on a blank-padded buffer simply returns the full width, which
`select case` matches because a Fortran character comparison blank-pads the shorter operand.
So the mode path was left exactly as it was. That, and one hand-written caller, are why
`c_char_1d_as_string` survives: `get_zip_entry_name` maps a fixed 4096-byte window over
libzip's `zip_get_name` and has no bound to remap against, so the NUL scan is the only thing
that can find the length. The other three -- `c_char_2d_as_string` and
`string_as_c_char_{1d,2d}` -- are gone, with their 42 call sites.

The probe, reduced to the case that decides it:

```fortran
subroutine take_in(strlen, n_strings, arr)
    integer(c_int), intent(in) :: strlen, n_strings
    character(kind=c_char, len=1), dimension(strlen, *), intent(in), target :: arr
    character(len=strlen), pointer :: view(:)
    integer :: i

    call c_f_pointer(c_loc(arr), view, [n_strings])
    do i = 1, n_strings
        print '(a,i0,3a)', 'in', i, '=[', trim(view(i)), ']'
    end do
end subroutine take_in
```

Re-run it against any new compiler before relying on this.

---

## Decision: kinds map through an explicit table

A kind with no entry in the table is an **error**, not a guess. Guessing produces a
wrapper that compiles and lies — which is exactly what the previous generator did for
`complex`, emitting `real(c_double_complex)`: a *real* of a complex kind. It compiles,
because `c_double_complex` and `c_double` are both 8 bytes.

---

## Decision: generated modules use `M_IMPLICIT_NONE`

A bare `implicit none` constrains only **variables**. A call to a procedure that does not
exist still compiles, as an implicit external, and fails at link time with a message
naming a symbol rather than a line.

Generated code is exactly where that must not be possible — a typo there is the
generator's fault and should surface at its output. `M_IMPLICIT_NONE` expands to
`implicit none (type, external)` (F2018), which makes it a compile error at the call site.

Worth adopting in hand-written modules for the same reason.

---

## Decision: an array with a separate shape is assumed size, and sliced on call

Serialised arrays arrive flat with their shape in a `<arg>_shape` argument, so the wrapper
declares `dimension(*)`. But an assumed-size actual **cannot be passed to an
assumed-shape dummy**, so the call slices it to its real extent:

```fortran
call fx_serialized(data = data(1:product(data_shape)), ...)
```

`product`, not `size`: `size(data_shape)` is the **rank**. (The previous generator emitted
`data(1:size(data_shape))`, which would have passed a 1000-element array as 2 elements.)

---

## Decision: declarations are ordered so extents come first

Referring to a symbol typed further down the specification part is a **GNU extension**,
not standard Fortran — `gfortran -std=f2018` rejects `dimension(n)` above `integer :: n`.
So anything named in someone's extents is hoisted; everything else keeps the author's
order.

---

## Decision: a function is exposed as a subroutine

C receives the result as an `intent(out)` argument, placed **before** the error code, so
that adding a result to a subroutine does not move `ierr` for existing callers.

---

## How this is checked

Not by eye. The test suite generates the wrappers from the fixture modules, compiles them
with `gfortran -std=f2018`, builds a shared library, and calls it from the generated
Python. See `tests/test_end_to_end.py`.

Three of the decisions above exist *because* the compiler rejected the first attempt.
