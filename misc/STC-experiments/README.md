# STC Experiments

Runnable, visual experiments for Shape Truthful Clustering (STC) -- the same synthetic
point-cloud datasets, generated the same way, that branch `smoothing`'s LoManLe experiments
(`run_lomanle_tests.sh`, `python/generate_smoothing_datasets.py`,
`r/plot_lomanle_spheres.R`) used to visualize a learned manifold skeleton. Here we run STC's
own pipeline over the same data and visualize what STC produces instead: **clusters
(ensembles), their tangent spaces, and their intersections** -- not a stitched manifold.

This is not a Jupyter notebook, on purpose: no notebook tech stack, just plain scripts you
run from a terminal, plus this file walking through what each one does and why. Read it
top to bottom once, then use it as a reference.

## Why this exists

Every kernel in this family (`seeds`, `ensemble_identification`, `ensemble_reconciliation`,
...) already has thorough unit tests -- see `test/mod_test_shape_truthful_clustering*.F90`
and their mirrored `python/test/`/`r/test/` suites. Those confirm each piece is *correct* on
small, hand-constructed, exactly-verifiable fixtures. What they cannot show is what the whole
pipeline actually *looks like* running on realistic, messy, 500-point data -- how many seeds
`seeds` finds on a noisy S-curve, whether `ensemble_identification`'s stop conditions trigger
sensibly, whether `ensemble_reconciliation` finds the intersections a human would expect by
eye. That is what this directory is for: a fast, visual sanity check and a place to explore
how STC's parameters trade off against each other, on the same test data your data scientist
was already using to explore LoManLe.

## Prerequisites

- The project builds (`./build.sh` from the repository root -- `run_stc_experiments.sh` does
  this for you).
- Python 3 with `numpy` and `pandas` (already used by `python/test/`).
- R with the `ggplot2` package (already used by nothing else in this repo, but commonly
  preinstalled; `install.packages("ggplot2")` if it is missing). Deliberately nothing else --
  see "Why only ggplot2" below.

## Quick start

```bash
# from the repository root, or from misc/STC-experiments/ -- both work
misc/STC-experiments/run_stc_experiments.sh misc/STC-experiments/data/2d/s_curve_2d_noise_low.csv
```

The first run builds the project and generates every dataset (a few seconds), then runs STC
on the one file you named with default parameters, and writes a PDF report to
`misc/STC-experiments/results/plots/`. Open it.

To sweep parameters, pass comma-separated lists (see `run_stc_experiments.sh --help`-style
usage message, printed when you run it with no arguments, for the exact positions):

```bash
misc/STC-experiments/run_stc_experiments.sh misc/STC-experiments/data/2d/kinked_curve_2d_noise_medium.csv \
    15,30 30,60 1,2 3.0,8.0 merge_jsi,merge_any 0.1
#   ^k_min ^alpha_max_deg ^d_max ^g_max ^reconciliation_mode ^min_jsi
```

produces one PDF per combination (2×2×2×2×2×1 = 16 here). To run every generated dataset at
once with one parameter combination:

```bash
misc/STC-experiments/run_stc_experiments.sh all
```

This can take a while (there are ~120 generated CSVs) -- start with one file while you are
still exploring parameters, and only sweep `all` once you have settled on values worth
surveying broadly.

## The pipeline, file by file

```
generate_datasets.py  -->  data/2d/*.csv, data/3d/*.csv
                                    |
                                    v
run_stc.py  (per CSV, per parameter combination)
    seeds -> ensemble_identification_merged -> ensemble_reconciliation
                                    |
                                    v
              results/data/<prefix>_{points,membership,ensembles,
                                      super_ensembles,super_ensembles_jsi,params}.csv
                                    |
                                    v
plot_stc.R  -->  results/plots/<prefix>.pdf

run_stc_experiments.sh orchestrates all three, sweeping parameter lists like
branch `smoothing`'s run_lomanle_tests.sh did for LoManLe.
```

### `generate_datasets.py`

A direct adaptation of branch `smoothing`'s `python/generate_smoothing_datasets.py`: same
seed (42), same 500 points, same three noise levels (low/medium/high = sigma 0.02/0.08/0.2),
same ten generator families (linear, cubic, exponential, s_curve, circular_arc,
kinked_curve, bifurcation_2way, bifurcation_3way, heteroscedastic, mixed_generators). The
only functional change: the original had a `DIM` global you edited and reran for; this
writes both the 2D and 3D variant of every dataset in one run, into `data/2d/` and `data/3d/`
respectively, instead of `results/data/2d,3d/`. Every CSV has columns `x1, ..., xk,
y_original` -- STC does not distinguish a "response" column from a "feature" column (unlike
LoManLe, which fits a manifold *through* the data), so `run_stc.py` just treats every numeric
column as one ambient dimension.

Deterministic: delete `data/` and rerun any time to regenerate byte-identical CSVs.

### `run_stc.py`

Runs `seeds` -> `ensemble_identification_merged` -> `ensemble_reconciliation` (the same three
calls `misc/mod_STC.md`'s "Tangent Space Variant" describes) on one CSV via the
`python/tensor_omics` bindings directly -- there is no need for a compiled Fortran test
driver the way LoManLe's experiments needed `test_aux/test_lomanle.f90`, since STC's whole
pipeline is already exposed to Python.

```
python3 run_stc.py <input.csv> [options]
```

Run `python3 run_stc.py --help` for the full option list. The most important ones, and what
they mean algorithmically (see `misc/mod_STC.md` for the full definitions):

| Option | STC concept | Default |
|---|---|---|
| `--density-quantile` | `seeds`' density-radius percentile | 0.15 |
| `--k-min` | growth radius: k-NN pool size for the median-distance radius | 30 |
| `--alpha-max-deg` | `accept_ensemble`'s max principal angle, in degrees (radians internally) | 30 |
| `--d-max` | `accept_ensemble`'s max tolerated change in intrinsic dimension | 1 |
| `--g-max` | `accept_ensemble`'s max tolerated `\|log(G_tp1/G_t)\|` | 3.0 |
| `--f-max` | Stop Condition 1's ensemble-size-fraction ceiling | 0.95 |
| `--a` | Stop Condition 2's "stably accepted" threshold | 2 |
| `--o` | trailing observable-history window depth | 10 |
| `--reconciliation-mode` | `report` / `merge_jsi` / `merge_any` | `merge_jsi` |
| `--min-jsi` | minimum JSI for `merge_jsi` | 0.1 |
| `--max-group-size` | cap on ensembles per super-ensemble | `min(1024, n_ensembles)` |

It writes six files per run (`<prefix>_points.csv`, `..._membership.csv`,
`..._ensembles.csv`, `..._super_ensembles.csv`, `..._super_ensembles_jsi.csv`,
`..._params.{json,txt}`) -- see "Output CSV schemas" below.

### `plot_stc.R`

```
Rscript plot_stc.R <prefix>
```

reads the six files above and writes `<prefix>.pdf`, an 7-8 page report (some pages are
skipped when there is nothing to show, e.g. no ensemble reaches intrinsic dimension 2, or no
super-ensemble was found at the current reconciliation threshold):

1. **Input point cloud and seeds** -- the raw data, with `seeds`' own output marked.
2. **Clusters (ensembles)** -- every point colored by which ensemble(s) it belongs to; grey
   = never joined any ensemble.
3. **Raw intersection density** -- every point colored by *how many* ensembles it belongs to,
   before reconciliation groups anything. This is the "ground truth" overlap Ensemble
   Reconciliation works from.
4. **Growth radii** -- one circle per ensemble, centered on its seed, radius =
   `calc_ensemble_growth_radius`'s output for that seed.
5. **Tangent space, first principal direction** -- one segment per ensemble, centered on its
   final `mu`, direction = `U`'s first column, length = the first `tangent_scales` entry
   (recovered from the retained singular value, see "Deriving `tangent_scales` without an
   extra binding call" below).
6. **Tangent space, second principal direction** -- same, for ensembles whose estimated
   intrinsic dimension is >= 2 (skipped if none are).
7. **Ensemble Reconciliation: super-ensembles** -- `ensemble_reconciliation`'s own output:
   filled points are ensembles that got grouped, colored by super-ensemble; open circles are
   ensembles that did not qualify at the current mode/threshold.
8. **Ensemble Reconciliation: JSI along each chain** -- one segment per consecutive pair
   within a super-ensemble's column (see `misc/mod_STC.md`'s `super_ensembles_JSI`
   definition), colored by JSI -- this is a spanning path through the group, not every
   pairwise JSI, exactly matching what the kernel itself reports.

#### Why only `ggplot2`

`r/plot_lomanle_spheres.R` also used `ggforce` (for `geom_circle`), the standalone `viridis`
package, and optionally `patchwork`. None of those three are installed in the environment
this was written in. Rather than depend on installs that may not be there:
`scale_color_viridis_c`/`_d` are built into `ggplot2` itself since 3.0, so the color scale
needs nothing extra; `geom_circle` is replaced by a ~10-line polygon helper
(`circle_df`/`circles_df` at the top of the script); `patchwork`'s only use in the original
was optional side-by-side panels, dropped here in favor of one page per plot (which is what
most of the original's own pages already were). If `ggforce`/`viridis`/`patchwork` are
available in your environment and you would rather use them, swapping them back in is a
small, self-contained change.

#### Deriving `tangent_scales` without an extra binding call

`ensemble_identification_merged` retains **singular values** (`ensemble_S_history`), not
eigenvalues or `tangent_scales` -- see `misc/mod_STC.md`'s "Output" section for why
(recomputable on demand from `S_history` + `k_history`, so nothing is stored redundantly).
`run_stc.py` recovers the segment length directly, the same formula
`tox_shape_truthful_clustering_observable_kernel`'s own `tangent_scales_kernel` uses:
`tangent_scale_j = sqrt(eigenvalue_j) = S_j / sqrt(k - 1)`. This is simple enough to inline
rather than round-trip through the `tangent_scales` Python binding.

### `run_stc_experiments.sh`

Orchestrates the three scripts above: builds the project, generates datasets if `data/` does
not exist yet, then loops `run_stc.py` + `plot_stc.R` over every combination of the parameter
lists you give it (or every dataset, with `all`) -- the same "one CSV or `all`,
comma-separated parameter lists, nested loop" shape as branch `smoothing`'s
`run_lomanle_tests.sh`. Run it with no arguments for the usage message.

## Output CSV schemas

For reuse outside `plot_stc.R` (a different plotting tool, a notebook, ad-hoc analysis).
`<prefix>` is whatever `--out-prefix` (or `run_stc_experiments.sh`'s own naming) chose.

- **`<prefix>_points.csv`** (wide, one row per input point): `point_id`, one column per
  ambient dimension (named after the input CSV's own numeric columns, e.g. `x1`,
  `y_original`), `n_ensembles` (how many ensembles this point belongs to).
- **`<prefix>_membership.csv`** (long, one row per point that belongs to at least one
  ensemble, repeated once per ensemble it belongs to): `point_id`, `ensemble_id`, `is_seed`
  (1 if `point_id` was this ensemble's own seed).
- **`<prefix>_ensembles.csv`** (wide, one row per ensemble -- including ones that ended up
  empty, e.g. Stop Condition 1, so row count always equals the seed count `seeds` found):
  `ensemble_id`, `seed_point_id`, `stop_reason` (`max_size` / `rejected_after_stable` /
  `rejected_immediately` / `fixed_point` / `error`), `growth_radius`, `size`, one `mu_<dim>`
  column per ambient dimension (the ensemble's final center), `d` (final estimated intrinsic
  dimension), `G` (final spectral gap), and `u1_<dim>`/`s1`, `u2_<dim>`/`s2` (the first two
  tangent directions and their scales -- `s2` is 0 when `d < 2`).
- **`<prefix>_super_ensembles.csv`** (long, one row per (super-ensemble, member ensemble)
  pair): `group_id`, `ensemble_id`.
- **`<prefix>_super_ensembles_jsi.csv`** (long, one row per consecutive pair within a
  super-ensemble's column): `group_id`, `ensemble_id_from`, `ensemble_id_to`, `jsi`.
- **`<prefix>_params.json`** / **`<prefix>_params.txt`**: every parameter `run_stc.py` ran
  with, plus `n_vectors`/`n_dimensions`/`n_ensembles`. Two formats for the same content --
  `.txt` (plain `key=value` lines) is what `plot_stc.R` reads, so it needs no JSON-parsing
  package; `.json` is there for anything that wants to parse it programmatically.

## Known limitations

- **3D datasets are plotted using only their first two coordinates.** The clustering itself
  runs in the full ambient dimension (3D data is genuinely clustered in 3D), but `plot_stc.R`
  only ever draws `dim_names[1]` vs. `dim_names[2]` -- there is no 3D rendering or
  projection here yet. A worthwhile follow-up, not attempted in this pass.
- **`data/` and `results/` are gitignored**, not committed -- unlike branch `smoothing`,
  which committed its generated `results/data/*.csv`. Both are fully reproducible
  (`generate_datasets.py` is seeded; `run_stc.py`/`plot_stc.R` are deterministic given the
  same inputs and parameters), so nothing is lost by regenerating on demand, and the
  repository stays smaller for it -- there are ~120 generated dataset CSVs alone.
- **No exhaustive default sweep has been run and committed.** This directory is a tool for
  exploring STC's behavior, not a benchmark suite with checked-in reference output; running
  `all` yourself, with whatever parameters you are currently investigating, is the intended
  workflow.

## Where the actual algorithm lives

This directory only *drives and visualizes* STC -- for the algorithm itself:

- `misc/mod_STC.md` -- the authoritative specification.
- `src/kernel/shape_truthful_clustering/` -- the implementation (`tox_shape_truthful_clustering_kernel.F90`
  for `ensemble_identification`(`_merged`), the sibling `_seeding_kernel.F90`,
  `_ensemble_growing_kernel.F90`, `_observable_kernel.F90`, `_accept_kernel.F90`,
  `_reconciliation_kernel.F90` modules for everything it orchestrates).
- `test/mod_test_shape_truthful_clustering*.F90` (and their `python/test/`/`r/test/`
  mirrors) -- the unit tests these experiments complement, not replace.
