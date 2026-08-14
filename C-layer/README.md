# stc_cli

A command line interface to Shape Truthful Clustering (STC), see `misc/mod_STC.md`'s own
"Command line interface (CLI) in C" section for the design this implements. Reads an
all-real-number CSV, runs STC's pipeline through the generated `_c` bindings (this directory
never re-implements any of STC itself), and writes an interactive HTML/D3 report plus the
same results in a few plain-text shapes for post-processing in Python/R/etc.

Has replaced `misc/STC-experiments/run_stc.py`/`plot_stc.R`/`export_json.py`/
`render_interactive.py`/`run_stc_pair.py` (all removed) -- a compiled CLI over hand-run
scripts. Every individual STC function stays callable from Python/R directly regardless
(`python/tensor_omics`/`r/tensor_omics`, generated the same way this CLI's own calls are).

## Files

- `csv_table.h`/`csv_table.c` -- a small, reusable, libcsv-backed reader for an
  all-real-number CSV, producing a flat buffer already laid out as Fortran's own
  column-major `vectors(n_dimensions, n_records)` expects (handed to the `_c` bindings below
  with zero copy). Not STC-specific; nothing here mentions STC at all.
- `stc_cli.c` -- the argp-based CLI itself: parameter parsing/validation, the k-d
  tree/seeding/ensemble-identification/reconciliation/estimation pipeline, and the five
  output writers. Pure wiring -- every actual computation happens in the Fortran kernels this
  file calls into.
- `test_stc_cli.sh` -- a smoke test (not part of `test_runner.sh`, which only exercises the
  Fortran/Python/R suites): builds the project, runs the CLI end to end against a small
  synthetic fixture (manual parameters, `--estimate-parameters`, and two validation-error
  paths), and checks the expected output files and exit codes.

## Prerequisites

- `libcsv` (LGPL-2.1-or-later), as a system package -- linked dynamically and scoped to just
  this executable (see `fpm.toml`'s own comment on why: `libtensor-omics.so`, what Python/R
  load, never needs it). On Fedora: `dnf install libcsv-devel`. On Debian/Ubuntu:
  `apt install libcsv-dev`.
- `argp` -- part of glibc; nothing extra to install on Linux.

## Building and running

```bash
./build.sh                      # builds stc_cli alongside the rest of the project
find build -name stc_cli        # locate the built executable
LD_LIBRARY_PATH=build build/<compiler-hash>/app/stc_cli --help
```

## Usage

```bash
LD_LIBRARY_PATH=build <path-to-stc_cli> \
    --input data.csv --header --n-records 500 --output-dir results/ \
    --k-min 30 --k-density 30 --chordal-dist-max-as-prcnt-of-range 0.3 --d-max 2 \
    --g-max 3.0 --rmse-change-max 3.0 --o 20
```

writes into `results/`:

| File | Contents |
|---|---|
| `report.html` | Self-contained interactive D3 report (vendored D3 + template + this run's JSON, all baked into one file -- see `src/tox_stc_html_assets.F90`) |
| `results.json` | The same data as `report.html` embeds, as a standalone JSON file |
| `points.csv` | One row per input vector: which ensemble(s)/super-ensemble(s) it belongs to, by the input CSV's own row number -- the file to join back to your own data in pandas/R |
| `ensemble_overlap_coefficients.csv` | The full pairwise ensemble Overlap Coefficient matrix (only pairs with a nonempty intersection) |
| `super_ensembles.tsv` | One line per super-ensemble, gene-family-file style: `<id>` TAB `<comma-separated member ensemble ids>` |

Every STC/estimation parameter is a flag; run `stc_cli --help` for the full, current list with
its defaults (this file does not duplicate that list, to avoid it drifting out of sync).
Broadly:

- `--chordal-dist-max-as-prcnt-of-range`, `--d-max`, `--g-max`, `--rmse-change-max`, `--o` have
  no kernel-side default and are always required (unless `--estimate-parameters` is given, for
  the first three of those five).
- Every other STC parameter (`--k-min`, `--k-density`, `--bandwidth-percentile`,
  `--exclusion-radius-percentile`, `--f-max`, `--a`, `--reconciliation-mode`,
  `--min-overlap-coefficient`, `--max-group-size`) is optional, defaulting to the same value
  the underlying kernel documents.
- `--reconciliation-exclude-stop-reasons` (comma-separated Stop Condition names, e.g.
  `rejected_immediately`, or `rejected_immediately,rejected_after_stable`), `--filter-d-min`/
  `--filter-d-max` (an ensemble's *final* intrinsic dimension, both inclusive -- distinct from
  `--d-max`, which bounds dimension *drift* during growth, not the final value), and
  `--filter-var-explained-min` (an ensemble's final classical variance explained,
  `sum(tangent eigenvalues) / (sum(tangent eigenvalues) + normal_error)`) each independently
  exclude ensembles from reconciliation entirely -- never merged into a super-ensemble, never a
  pair in `ensemble_overlap_coefficients.csv`/`results.json`'s `overlap_coefficient_matrix` --
  while leaving them fully visible everywhere else (points, the report's main plot, the
  ensembles list, each now also carrying `reconciliation_eligible`/`excluded_by` so a report can
  show which criterion, if any, excluded it). All four default to no filtering. See
  `misc/mod_STC.md`'s "Filtering ensembles before merging".

### `--estimate-parameters`

Runs `estimate_stc_parameters` first and applies its own
`k_min`/`k_density`/`chordal_dist_max_as_prcnt_of_range`/`G_max`/`d_max` as this run's actual
parameters (its own `density_quantile` output has no run parameter to apply to; it is only
ever reported, as `estimated_density_quantile` in the JSON). Supplying any of
`--k-min`/`--k-density`/`--chordal-dist-max-as-prcnt-of-range`/`--g-max`/`--d-max` together
with `--estimate-parameters` is a validation error -- pick estimation or manual values, never
both. `--n-anchors`/`--seed-max-set-size`/`--first-quartile-percentile` tune the estimation
step itself and are only meaningful in this mode.

```bash
LD_LIBRARY_PATH=build <path-to-stc_cli> \
    --input data.csv --header --n-records 500 --output-dir results/ \
    --estimate-parameters --rmse-change-max 3.0 --o 20
```

Whichever mode is used, the six estimated values are always reported in `results.json`'s
`params` object under an `estimated_` prefix, whether or not they were actually applied to
the run.
