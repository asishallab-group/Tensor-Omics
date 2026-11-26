# Fixture modules

These are the generator's **specification**. They are not part of the library and are
never compiled into it.

The real sources carry no `category: C-binding` tags yet, so there is nothing to
generate from. Rather than wait, the generator is developed and tested against these:
small, complete specimens of every construct it claims to support. When the real sources
are annotated, these stay — a fixture is a case someone chose deliberately, which is
worth more as a regression test than whatever the library happens to contain today.

## The `fx_` prefix

`fx` is short for **fixture**. Everything here carries it — modules, procedures,
parameters — for two reasons:

- `fx_basics` next to `tox_normalization` is unmistakable. Nobody should have to wonder
  whether a module is library code or scaffolding.
- `grep -r fx_` finds every piece of test scaffolding in one pass.

It has no meaning to the generator: nothing keys off it, and a real module could be called
`fx_anything` without consequence.

## The rule these files must satisfy

**They parse to a clean IR with no diagnostics at all** — no errors, *and no warnings*.
`test_frontend.py::TestParsingProducesCleanIR` enforces it.

That bar exists because they are what the emitters are tested against. A fixture that
warns is a fixture nobody reads carefully.

The exception is deliberate: `fx_internal` breaks the export contract on purpose (a
deferred-length string, a mode argument with no table) to prove that an **internal**
procedure is held to none of it. It carries no `category: C-binding`, so nothing is
generated from it, and it must produce no diagnostics either.

## What is where

| File | Covers |
|---|---|
| `src/fx_basics.F90` | the ordinary cases: extents, a shared extent, optionals with defaults, temporaries, mode and method arguments, a function, one with its own `ierr` (its C wrapper takes the result before it), an internal procedure |
| `src/fx_edges.F90` | the awkward ones: every character length form, a separately-travelling shape, a mask and its count, a matrix result trimmed by a count, `c_bool`, a hand-written two-tier pair, `fx_cluster` / `fx_cluster_expert`, published as written, `DM_OUTPUT_FROM`, a nullable optional, a procedure with no `ierr`, two arrays in one declaration (`dimension(n) :: a, b(m)`), an extent spelled in another case than its declaration, a kinded logical default |
| `impl/fx_ranks_impl.F90` | an implementation with a documented minimum, a work array and a permutation, so both generated wrappers are meaningful -- compiled and run by `test_end_to_end_fortran.py` |
| `impl/fx_prologue_impl.F90` | an implementation with a prologue that asks for an argument of its own, beside the implementation's own `ierr` -- compiled and run by `test_end_to_end_fortran.py`, through Python and R by `test_end_to_end_impl.py` (the prologue's argument before `ierr`), and read back by `test_names.py` for the locals a prologue wrapper declares |
| `impl/fx_shapes_impl.F90` | an assumed-shape input the wrapper validates, an extent and a character length named by a module constant, an `impure elemental` implementation -- compiled and run by `test_end_to_end_fortran.py`, C wrappers included, and through Python and R by `test_end_to_end_impl.py`, R shims compiled into the library |
| `impl/fx_functions_impl.F90` | function implementations: a scalar, an array and a `logical(c_bool)` result, an `elemental` function, one with its own `ierr` and one with a work array, every tier publishing the result right before `ierr`; a `complex(real64)` result, and subroutines with a `complex(real64)` output, input, optional input and `intent(inout)` argument, each of which Python boxes in two doubles -- compiled and run by `test_end_to_end_fortran.py` (by position), through C, Python and R by `test_end_to_end_impl.py` |
| `impl/fx_modes_impl.F90` | an optional mode without a default, which C may pass as null, and an output saying whether `present()` saw it -- compiled and run by `test_end_to_end_fortran.py`, through Python and R by `test_end_to_end_impl.py` (given and omitted), and deleted there to show that regenerating leaves nothing of it behind |
| `impl/fx_extents_impl.F90` | the extents a binding computes and the ones it checks (#221): an output sized by an integer division Python must truncate as Fortran does, a constant input axis, an optional array sharing an extent, an input sized by a mask's count, and input axes the declaration computes (`n + 1`, `sum(group_sizes)`, one of them optional), which the C wrapper compares; an extent with a `DM_EQUALS` (#203) -- compiled and run by `test_end_to_end_fortran.py`, through Python and R by `test_end_to_end_impl.py` |
| `impl/fx_positions_impl.F90` | implementations and a prologue that fail at a position the test sets, in their own dummy lists: a permutation, a work array, a recommend-sized value (from a producer in the same module), a fixed mode, another mode's optional, a prologue-produced local and a function result ahead of kept arguments -- compiled and run by `test_end_to_end_fortran.py`, every position in turn, to show each wrapper returns the position in its own list |
| `impl/fx_producers_impl.F90`, `impl/fx_sized_impl.F90` | `DM_OUTPUT_FROM` producers that are generated wrappers, and their consumers: a work array the Fortran wrapper sizes by calling the producer and a result Python and R size by calling its binding, each with a renamed input; two producers, the first of which refuses its input; a producer in its consumer's own module; a producer extent the consumer calls by another name; an output nothing takes; a producer taking a shape the bindings read off an array; a `JUST_INFO` size, whose link names the wrapper; an array a producer fills, with a permutation sorted by it and a work array sized by another producer; one call giving an array and a size, with any number of columns, and declared before the array the caller passes; a value only a bound names; an integer array a producer fills; optional producer-input suppliers, one with a `DM_DEFAULT` under another name and one passed through absent; one call giving data and a size that sizes a result, or an input the caller passes, so the call is kept whole (rule K); producer inputs left out with `OMITTED`, an optional without a default and one with, on both routes, and omitted in one table while supplied in another (rule O); a mode argument supplying a producer's mode input that accepts one more mode, and a constant naming one of its mode parameters, on both routes (rule T); an optional of the consumer's whose default differs from the producer's, and optionals without a default, an integer and a mode, left out of the binding's call (P1) -- through Python and R by `test_end_to_end_impl.py` (each producer error names the caller's argument; the expert takes the data and computes the sizes), compiled and run by `test_end_to_end_fortran.py` |
| `carrier/fx_carrier_impl.F90` | a value a producer computes that sizes only an input, which the Python and R expert read off it, while the plain tier computes it and checks the input against it; and a lone wrapper whose computed capacity sizes an `intent(inout)` array -- built and run through Python and R by `test_end_to_end_impl.py`, its parameters traced by `test_parameters_reach_the_call.py` |

The implementation fixture sits in `impl/`, not beside the other two in `src/`, because `src/` is
parsed as one project whose module list is asserted exactly (`test_frontend.py`) and globbed
whole (`test_emit_fortran_c.py`). A third module there would break both for no gain.
`carrier/` is a directory of its own, built into a library of its own: it is the case a plain
tier's R shim once declared twice (the value R computed, and a second one read off the input
that carries it), and kept apart it shows that shim compiling without the rest of `impl/`.

## Adding one

Some of these have real bodies rather than empty ones. That is not decoration:
`test_end_to_end.py` builds them into a shared library and checks the answers through the
generated Python, so `fx_sum_matrix` really has to sum a matrix.

When adding a case:

1. give it a body if an end-to-end test would be meaningful, otherwise leave it empty
2. keep it to the one construct it is there to demonstrate
3. run the suite — it must stay diagnostic-free
4. if it is an edge case, say so in its `summary:`; the fixture is documentation too
