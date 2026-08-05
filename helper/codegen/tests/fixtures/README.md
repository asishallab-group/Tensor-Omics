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
| `src/fx_basics.F90` | the ordinary cases: extents, a shared extent, optionals with defaults, temporaries, mode and method arguments, a function, an internal procedure |
| `src/fx_edges.F90` | the awkward ones: every character length form, a separately-travelling shape, a mask and its count, `c_bool`, an alloc/expert pair, `DM_OUTPUT_FROM`, a nullable optional, a procedure with no `ierr` |
| `kernel/fx_ranks_kernel.F90` | a kernel with a documented minimum, a work array and a permutation, so both generated wrappers are meaningful -- compiled and run by `test_end_to_end_fortran.py` |

The kernel fixture sits in `kernel/`, not beside the other two in `src/`, because `src/` is
parsed as one project whose module list is asserted exactly (`test_frontend.py`) and globbed
whole (`test_emit_fortran_c.py`). A third module there would break both for no gain.

## Adding one

Some of these have real bodies rather than empty ones. That is not decoration:
`test_end_to_end.py` builds them into a shared library and checks the answers through the
generated Python, so `fx_sum_matrix` really has to sum a matrix.

When adding a case:

1. give it a body if an end-to-end test would be meaningful, otherwise leave it empty
2. keep it to the one construct it is there to demonstrate
3. run the suite — it must stay diagnostic-free
4. if it is an edge case, say so in its `summary:`; the fixture is documentation too
