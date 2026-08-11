# Tensor Omics

## Fortran coding guidelines

All Fortran code in this repo must strictly follow `./misc/Fortran_Coding_Guides.tex` (the
"F42" standard) for general Scientific Kernel (SK) philosophy and coding style, **and**
`codegen_guide.md` for naming/annotation conventions and the kernel/wrapper split -- read
both before writing or reviewing any Fortran. `codegen_guide.md` supersedes F42 section
10.3.2's older `_alloc -> unsuffixed -> _helper` naming triad; F42 still governs everything
else (purity, types, control flow, parallelization, data layout).

F42 highlights that still apply everywhere:

- `real(real64)` (`iso_fortran_env`) everywhere, literals with `_real64` suffix.
- Every argument gets an explicit `intent`; use `pure` wherever possible.
- Error handling via `tox_errors` (`set_ok`/`set_err`/`is_err`, standard error codes), errors
  propagate up unmodified.
- No `GOTO` -- use `exit`/`cycle`/`select case`; manual cleanup subroutines instead of `final ::`.
- DataTables are column-major: prefer column access (`arr(:,j)`, no copy) over row access
  (`arr(i,:)`, triggers a copy).
- FORD doc comments (`!>` description, `!|` per-argument, `!!` examples/warnings).
- OpenMP: `!$omp simd` only in the leaf (innermost) loop doing arithmetic; `!$omp parallel do`
  only at the top-level outer loop; don't combine as `parallel do simd` without explicit review.
  Alignment via `precompiler_constants.F90` / `DEFAULT_ALIGNMENT` macro, never hardcoded.
- Structure-of-Arrays (SoA), never Array-of-Structures, for every data representation.
- Every function does one thing well, with a signature designed to be rigorously unit
  testable in isolation (single responsibility, no hidden coupling to caller state).

## Code generation

TensorOmics generates its public API from hand-written kernels. Read `codegen_guide.md`
before writing or reviewing any kernel -- it is the authoritative, exhaustive reference;
this section is only a pointer, not a summary to keep in sync by hand.

The load-bearing rule: **only write kernels (`src/kernel/`) and other hand-written,
non-generated trees (`src/f42/`, `src/data/`)**. Never hand-edit `src/generated/`,
`python/tensor_omics/`, `r/tensor_omics/`, or `snippets/` -- the generator deletes and
rewrites all four on every run, so anything written there by hand survives exactly until the
next build.

- A kernel is a `pure` procedure named `<name>_kernel`, in a module named
  `tox_<family>_kernel` under `src/kernel/`. No validation, no `ierr` (except a genuine
  runtime failure no input check could foresee -- an external library reporting failure,
  e.g.), no allocation anywhere in the module.
- Every argument's `!!` doc, plus its `DM_MIN`/`DM_MAX`/`DM_SENTINEL`/`DM_ALLOW_NAN`/
  `DM_ALLOW_INFINITE`/`DM_DEFAULT`/`DM_OUTPUT_FROM`/`DM_PROLOGUE` annotations where
  applicable, is the *entire* contract the generator has to work from -- there is no other
  place to state a constraint. Undocumented or under-annotated arguments generate an
  unvalidated or wrongly-validated public API.
- The generator writes the validating wrapper, the allocating wrapper (`_alloc`, published
  as the plain name), the expert tier (published as `_expert`), and the C/Python/R bindings.
- IO, infrastructure, and anything that is not a numeric pipeline kernel is hand-written and
  exported via `M_EXPORT_C` instead (`src/data/`, `src/f42/`) -- see `codegen_guide.md` Part II.

## Unit Test Writing

Create very stressing and atomic test cases for each procedure we created in
`test/`. Use the assertions from @test/asserts.F90. Run the tests with `bash
test_runner.sh --skip-kinds-test`. A new `test/mod_test_foo.F90` file will
contain the tests, uses asserts, test_case from test_suite, tox_errors and
`int32`/`real64` from `iso_fortran_env` and the to-be-tested module. As entry
point for `test/run_tests.F90`, create a `get_all_tests_foo` function, matching
the `get_all_interface` from `test/run_tests.F90`. Also register that one
there. The tests should be very stressing. The code generator generates the
Python/R functions for our framework into @r/tensor_omics,
@python/tensor_omics. The tests we add to @test/ also go to @python/test and
@r/tests. In Python use `assert` statements and the
`{python,r}/test_helpers.*::run_all_tests`. In R use for assertions the ones in
`r/test_helpers.R`. For both languages don't use third-party test frameworks,
only `assert` statements and test helpers. No `stop()` in any R test file, just
assertions.

**Every crash or error is a regression test.** Whenever something -- a build, a test run, an
experiment script, manual exploration -- crashes or raises an error you did not expect, that
is a unit test you are missing, not just a bug to patch and move past. Before considering the
issue resolved: (1) write a test in `test/` that reproduces the failure condition (never one
that reproduces an actual crash/undefined-behavior itself -- see the note below), (2) confirm
it fails against the current code, (3) fix the underlying code, (4) confirm the same test now
passes, (5) mirror it to `python/test/` and `r/tests` per the convention above. If the
failure is genuinely not reproducible from Python/R (e.g. it depended on a raw Fortran call
omitting an optional argument -- something no binding does, since both always resolve
defaults explicitly before crossing the ABI boundary), say so in a comment at the mirror site
instead of forcing an artificial test; see
`test/mod_test_shape_truthful_clustering_seeding.F90`'s
`test_seeds_omitted_k_density_is_clamped` and its Python/R siblings for a worked example of
both halves of this (the Fortran-only crash regression, and the adjacent behavior that *is*
reachable from Python/R, covered there instead).

If the failure was an actual crash (a segfault, corrupted memory, undefined behavior) rather
than a clean, validated error: do not write a test that re-triggers the crash to prove it is
fixed -- a test that sometimes segfaults is worse than no test, since a crash can take the
whole test binary down with it, mask every later test in the same run, and behave differently
across platforms/builds/compiler flags. Fix the code first, so the failure mode becomes a
clean, deterministic result (a correct value, or a proper `ierr`) -- only then write the
regression test, asserting *that* deterministic result.

## Shape Truthful Clustering (STC)

STC (Shape Truthful Clustering) is a renormalization-group-inspired ensemble-growth
clustering method, planned to replace LoManLe's Steps 1-5b (adaptive neighborhood growth +
greedy anchor selection). `misc/mod_STC.md` is the sole, authoritative implementation spec
(algorithm, stop conditions, output shapes, naming). `misc/STC_for_LoManLe.md` and
`misc/STC_current_algorithm_draft.md` remain as earlier-draft/LoManLe-integration-rationale
references only -- implement against `mod_STC.md`.

Kernel modules, one per major step, under `src/kernel/shape_truthful_clustering/`:

- `tox_shape_truthful_clustering_kernel.F90` -- parent; holds `ensemble_identification`'s
  own kernel(s) directly (the family's natural top-level entry point) *and* `use`s the six
  children below -- a deliberate deviation from `codegen_guide.md` section 5.15's own
  `tox_data_integration_kernel` example, where the parent holds no procedures of its own.
  `ensemble_identification`/`_merged` also report `low_confidence_mask`/
  `ensemble_low_confidence_masks`, one column per seed regardless of `stop_reason` -- see
  `misc/mod_STC.md`, "Ensemble identification", "Output".
- `tox_shape_truthful_clustering_seeding_kernel.F90` -- `density_labels`, `seeds`. Also
  `use`s its sibling `tox_shape_truthful_clustering_ensemble_growing_kernel` directly, to
  reuse `calc_ensemble_growth_radius_kernel` for `seeds`' own coverage radius rather than a
  second, separately-implemented computation -- the first sibling-to-sibling dependency in
  this family.
- `tox_shape_truthful_clustering_ensemble_growing_kernel.F90` -- `calc_ensemble_growth_radius`,
  `grow_ensemble`.
- `tox_shape_truthful_clustering_observable_kernel.F90` -- `observable`, `normal_error`,
  `tangent_scales`.
- `tox_shape_truthful_clustering_accept_kernel.F90` -- `accept_ensemble`.
- `tox_shape_truthful_clustering_reconciliation_kernel.F90` -- `ensemble_reconciliation`.
- `tox_shape_truthful_clustering_parameter_estimation_kernel.F90` -- `sample_estimator_anchors`,
  `grow_estimator_anchor_clouds`, `estimate_stc_parameters`. A separate, optional pipeline
  step (not run automatically by `seeds`/`ensemble_identification`) that proposes starting
  values for `k_min`/`k_density`/`density_quantile`/`alpha_max`/`G_max`/`d_max` directly from
  the data -- see `misc/mod_STC.md`, "Estimate parameters from data", for the full algorithm
  and its deliberate simplicity trade-offs. Not yet wired into
  `misc/STC-experiments/run_stc.py`.

STC's own output layer, hand-written (not generated) but not kernels either:

- `src/tox_stc_html_assets.F90` -- generated (via `helper/embed_stc_html_assets.py`, not the
  main codegen pipeline; do not hand-edit) compile-time Fortran string constants for the
  vendored D3 bundle and the interactive report template.
- `src/tox_stc_json.F90` -- `M_EXPORT_C` module: `serialize_stc_results_as_json`,
  `write_stc_interactive_html_report`. Takes STC's raw pipeline arrays directly.
- `src/tox_stc_csv.F90` -- `M_EXPORT_C` module: `serialize_stc_points_as_csv`,
  `serialize_stc_ensemble_overlap_as_csv`, `serialize_stc_super_ensembles_as_tsv` --
  plain-text companions to the JSON/HTML report, for post-processing in Python/R/etc.
- `C-layer/` -- the C CLI (`stc_cli`, GNU argp + libcsv), see `C-layer/README.md` and
  `misc/mod_STC.md`'s own "Command line interface (CLI) in C" section. Intended, eventually,
  to replace `misc/STC-experiments/run_stc.py`/`plot_stc.R`, contingent on feature parity;
  every individual STC function stays callable from Python/R directly regardless.

`src/lomanle.F90`: existing module. Only integrate STC into LoManLe's pipeline after STC
itself is implemented and tested in isolation -- do not interleave the two.
