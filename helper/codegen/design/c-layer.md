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

1. `safeguard` is needed only by C-facing code, so only these modules use it
2. the whole C interface can be compiled out with one directive — the wrappers are wrapped
   in `#ifndef NO_C_INTERFACE`, and with them goes the `safeguard` dependency
3. the source files are not bloated with wrapper definitions
4. the generator never touches hand-written source

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

The copy into a default `logical` stays — `c_bool` is one byte, a default `logical` is
four — but it is one copy rather than a conversion at both ends, and it is a plain
intrinsic assignment (`flag_f = flag`) rather than a call into `tox_conversions`. A dummy
already declared `logical(c_bool)` needs no copy at all.

*Rejected:* passing `integer(c_int)` 0/1 and converting with `c_int_as_logical`, as the
previous generator did. It makes every caller in every language do bool→int→logical.

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
