# Tensor Omics

## Fortran coding guidelines

All Fortran code in this repo (core Scientific Kernel, C wrappers, helpers) must strictly
follow `./misc/Fortran_Coding_Guides.tex` (the "F42" standard). Read that file before
writing or reviewing any Fortran. Highlights:

- `real(real64)` (`iso_fortran_env`) everywhere, literals with `_real64` suffix.
- Every argument gets an explicit `intent`; use `pure` wherever possible.
- Optional arguments only in core Fortran routines, never in `bind(C)`/R-facing wrappers.
- Three-layer architecture: pure SK core -> `_C` `bind(C)` wrapper -> Python (`ctypes`)/R (`Rcpp`).
- `_C` wrappers: minimal (wrap only, no computation), `iso_c_binding` types only, explicit-dimension
  arrays (no assumed-shape), `target` attribute on all args, null-checked via `src/macros.h`
  (`M_USE_NULL_VALIDATION`, `M_CHECK_IERR_NON_NULL`, `M_CHECK_NON_NULL`), always an `ierr` output.
- Error handling via `tox_errors` (`set_ok`/`set_err`/`is_err`, standard error codes), errors
  propagate up unmodified.
- No `GOTO` — use `exit`/`cycle`/`select case`; manual cleanup subroutines instead of `final ::`.
- SK routines: no I/O, no `allocate`/`deallocate`, deterministic, allocatable+contiguous args.
- DataTables are column-major: prefer column access (`arr(:,j)`, no copy) over row access
  (`arr(i,:)`, triggers a copy).
- FORD doc comments (`!>` description, `!|` per-argument, `!!` examples/warnings).
- OpenMP: `!$omp simd` only in the leaf (innermost) loop doing arithmetic; `!$omp parallel do`
  only at the top-level outer loop; don't combine as `parallel do simd` without explicit review.
  Alignment via `precompiler_constants.F90` / `DEFAULT_ALIGNMENT` macro, never hardcoded.
- Mandatory workflow for new functionality: core Fortran -> C wrapper -> Python wrapper -> R
  wrapper -> docs -> tests (Fortran tests the algorithm; Python/R tests only call-ability).
- Structure-of-Arrays (SoA), never Array-of-Structures, for every data representation.
- Every function does one thing well, with a signature designed to be rigorously unit
  testable in isolation (single responsibility, no hidden coupling to caller state).

## Shape Truthful Clustering (STC) integration plan

STC (Shape Truthful Clustering) is a renormalization-group-inspired ensemble-growth
clustering method. It is planned to replace LoManLe's Steps 1-5b (adaptive neighborhood
growth + greedy anchor selection). Full design docs: `misc/STC_for_LoManLe.md` (converged
decisions and rationale for the LoManLe integration) and
`misc/STC_current_algorithm_draft.md` (the general algorithm, abstract API, and the
density-based and tangent-space instantiations).

Modules:

- `src/shape_truthful_clustering.F90`: new module. Implement STC here first, standalone,
  strictly to the F42 guidelines above (SoA, pure SK kernels, `_alloc` wrappers doing all
  allocation outside the pure kernels, `src/macros.h` validation macros, single-responsibility
  testable functions).
- `src/lomanle.F90`: existing module. Only integrate STC into LoManLe's pipeline after STC
  itself is implemented and tested in isolation -- do not interleave the two.
