# STC Experiments

Runnable, visual experiments for Shape Truthful Clustering (STC) -- the same synthetic
point-cloud datasets, generated the same way, that branch `smoothing`'s LoManLe experiments
(`run_lomanle_tests.sh`, `python/generate_smoothing_datasets.py`,
`r/plot_lomanle_spheres.R`) used to visualize a learned manifold skeleton. Here we run STC's
own pipeline over the same data and visualize what STC produces instead: **clusters
(ensembles), their tangent spaces, and their intersections** -- not a stitched manifold.

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
  this for you), which now also builds `stc_cli` (see `C-layer/README.md`), including its
  `libcsv` system dependency (`dnf install libcsv-devel` on Fedora).
- Python 3, only to generate the synthetic datasets (`generate_datasets.py`; `numpy`).
- Nothing else. `stc_cli` writes the interactive HTML report itself (D3 v7 is vendored
  in-repo, `vendor/d3.v7.min.js`, and baked into the executable at compile time -- see
  `src/tox_stc_html_assets.F90`); no R, no separate plotting step, no server, no npm.

## Quick start

```bash
# from the repository root, or from misc/STC-experiments/ -- both work
misc/STC-experiments/run_stc_experiments.sh misc/STC-experiments/data/2d/s_curve_2d_noise_low.csv
```

The first run builds the project and generates every dataset (a few seconds), then runs STC
on the one file you named with default parameters, and writes five output files (an
interactive HTML report among them) to their own directory under
`misc/STC-experiments/results/`. Open `report.html`.

To sweep parameters, pass comma-separated lists (see `run_stc_experiments.sh`'s own usage
message, printed when you run it with no arguments, for the exact positions):

```bash
misc/STC-experiments/run_stc_experiments.sh misc/STC-experiments/data/2d/kinked_curve_2d_noise_medium.csv \
    15,30 0.5,0.7 1,2 3.0,8.0 merge_overlap_coefficient,merge_any 0.9
#   ^k_min ^chordal_dist_max_as_prcnt_of_range ^d_max ^g_max ^reconciliation_mode ^min_overlap_coefficient
```

produces one output directory per combination (2×2×2×2×2×1 = 16 here). Two more optional,
trailing lists sweep `seeds`' own parameters -- `k_density` (its neighborhood size,
independent of `k_min`) and `exclusion_radius_percentile` (its coverage/exclusion radius,
independent of the growth-phase radius, see "Known limitations" below):

```bash
misc/STC-experiments/run_stc_experiments.sh misc/STC-experiments/data/2d/kinked_curve_2d_noise_medium.csv \
    30 0.5 1 3.0 merge_overlap_coefficient 0.9 15,30 25,50
#   ^k_min ^chordal ^d_max ^g_max ^mode ^min_overlap_coefficient ^k_density ^exclusion_radius_percentile
```

Each output directory's name includes every one of these values (e.g.
`..._k30_kd30_..._erp50.0`), so results from different combinations never collide or get
mixed up.

To run every generated dataset at once with one parameter combination:

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
                    stc_cli  (per CSV, per parameter combination)
    seeds -> ensemble_identification_merged -> ensemble_reconciliation
                                    |
                                    v
       results/<prefix>/{report.html, results.json, points.csv,
                          ensemble_overlap_coefficients.csv, super_ensembles.tsv}
```

`run_stc.py`, `plot_stc.R`, `export_json.py`, `render_interactive.py`, and `run_stc_pair.py`
-- the Python/R scripts that used to drive this pipeline and produce PDF/interactive-HTML
reports from separate CSV intermediates -- have been removed. `stc_cli` (`C-layer/`, GNU
argp + libcsv) now does the whole thing directly from Fortran, in one call, with no CSV/JSON
intermediate: see `C-layer/README.md` for its own full option list and output-file
descriptions, and `misc/mod_STC.md`'s "Command line interface (CLI) in C" section for the
design. Every STC parameter this directory's own scripts used to expose is still exposed,
via `stc_cli`'s own flags; run `stc_cli --help` (via `run_stc_experiments.sh` locating the
built binary, or directly -- see `C-layer/README.md`'s own "Building and running") for the
authoritative, current list.

`generate_datasets.py` is unchanged and still needed: it produces the CSV inputs
`stc_cli` (via `run_stc_experiments.sh`) reads.

### `generate_datasets.py`

A direct adaptation of branch `smoothing`'s `python/generate_smoothing_datasets.py`: same
seed (42), same 500 points, same three noise levels (low/medium/high = sigma 0.02/0.08/0.2),
same ten generator families (linear, cubic, exponential, s_curve, circular_arc,
kinked_curve, bifurcation_2way, bifurcation_3way, heteroscedastic, mixed_generators). Writes
both the 2D and 3D variant of every dataset, into `data/2d/` and `data/3d/` respectively.
Every CSV has columns `x1, ..., xk, y_original` -- STC does not distinguish a "response"
column from a "feature" column (unlike LoManLe, which fits a manifold *through* the data),
so `stc_cli` just treats every numeric column as one ambient dimension.

Deterministic: delete `data/` and rerun any time to regenerate byte-identical CSVs.

### `run_stc_experiments.sh`

Orchestrates the pipeline above: builds the project (which builds `stc_cli` too), generates
datasets if `data/` does not exist yet, then loops `stc_cli` over every combination of the
parameter lists you give it (or every dataset, with `all`) -- the same "one CSV or `all`,
comma-separated parameter lists, nested loop" shape as branch `smoothing`'s
`run_lomanle_tests.sh`, and as this script's own prior version. Two parameters `stc_cli` has
no default for at all (`--rmse-change-max`, `--o`) are fixed at this script's own top
(matching the removed `run_stc.py`'s former defaults, `|log(1.5)|` and `10`) rather than
exposed as more sweep dimensions -- call `stc_cli` directly if you need to vary them. Run
this script with no arguments for the usage message.

### No stitched manifolds

The interactive HTML report never draws a super-ensemble's members as a connected, ordered
path/skeleton -- `ensemble_reconciliation` only ever reports *which* ensembles intersect
enough to group, not in what order or with what geometry they connect. `super_ensembles`'
member order within a group is ensemble-*discovery* order (highest density first), not
spatial order, so drawing it as a chain would zigzag across the data with no relationship to
the actual manifold shape. Turning a set of intersecting ensembles into one stitched manifold
is LoManLe's job, downstream of STC, not something STC's own report should imply it has
already done.

## Output file schemas

For reuse outside the interactive HTML report (a notebook, ad-hoc analysis, another plotting
tool). See `C-layer/README.md` for the full table; briefly:

- **`points.csv`** (one row per input point, by the input CSV's own row number): which
  ensemble(s)/super-ensemble(s)/low-confidence-ensemble(s) it belongs to, and which ensemble
  (if any) it is the seed of.
- **`ensemble_overlap_coefficients.csv`**: the full pairwise ensemble Overlap Coefficient
  matrix (`|intersection| / min(|A|,|B|)`), only pairs with a nonempty intersection.
- **`super_ensembles.tsv`**: one line per super-ensemble, gene-family-file style (`<id>` TAB
  `<comma-separated member ensemble ids>`).
- **`results.json`**: everything above, plus per-ensemble geometry (`stop_reason`,
  `growth_radius`, `size`, `d`, `G`, `mu`, tangent directions `u1`/`s1`/`u2`/`s2`, tangent-line
  endpoints `line_start`/`line_end`, `t_final` -- the last accepted growth iteration --
  `observable_history` -- one `{iteration, g, rmse, drift}` entry per retained accepted growth
  iteration, `drift` being the consecutive-iteration chordal distance, see `misc/mod_STC.md`'s
  "Ensemble Observable Plots" -- and `final_chordal_distance`, the chordal distance actually
  tested when the ensemble's real last growth step was accepted, present only when
  reconstructable), `reconciliation_eligible`/`excluded_by` (whether this ensemble was eligible
  to be merged into a super-ensemble, and if not, which of `stop_condition`/`dimension`/
  `variance_explained` excluded it -- see `misc/mod_STC.md`'s "Filtering ensembles before
  merging"; an ineligible ensemble still has every other field above, it just never appears in
  `super_ensembles` or `ensemble_overlap_coefficients.csv`), each point's per-ensemble
  `residual_length` (distance off that ensemble's tangent subspace), and every parameter the run
  used -- the same data `report.html` embeds, as a standalone file.

## Known limitations

- **The interactive HTML report only plots the first 2 ambient dimensions.** `results.json`'s
  `points[].coords` still has every ambient coordinate, so a 3D-aware extension to
  `interactive_template.html` is possible; not attempted so far.
- **No PDF report and no 3D report at all**, unlike this directory's earlier, now-removed
  `plot_stc.R` (which drew a `scatterplot3d`-based static 3D PDF for 3+ ambient dimensions).
  `stc_cli` runs correctly on 3D+ data -- clustering itself is dimension-agnostic -- only its
  *visualization* is currently 2D-only, per the point above.
- **`data/` and `results/` are gitignored**, not committed -- unlike branch `smoothing`,
  which committed its generated `results/data/*.csv`. Both are fully reproducible
  (`generate_datasets.py` is seeded; `stc_cli` is deterministic given the same inputs and
  parameters), so nothing is lost by regenerating on demand, and the repository stays smaller
  for it -- there are ~120 generated dataset CSVs alone.
- **No exhaustive default sweep has been run and committed.** This directory is a tool for
  exploring STC's behavior, not a benchmark suite with checked-in reference output; running
  `all` yourself, with whatever parameters you are currently investigating, is the intended
  workflow.
- **Seed-exclusion suppression is purely geometric and does not know about curvature.**
  `seeds`' coverage/exclusion radius (tunable via `--exclusion-radius-percentile`, see "Quick
  start" above) can suppress seed placement across a region -- typically a curvature extremum
  (a peak, trough, or kink on a wavy manifold) -- that the suppressing seed's own later growth
  never actually reaches, because `accept_ensemble`'s curvature-based stop conditions
  (`chordal_dist_max_as_prcnt_of_range`/`g_max`/`d_max`/`rmse_change_max`) can halt growth well
  short of that same geometric radius.
  The result is points with no seed of their own and no membership in any grown ensemble
  either. Shrinking `exclusion_radius_percentile` reduces how much territory each seed
  suppresses, but does not close the gap on its own -- it only shrinks it. Fixing this
  properly needs seed selection to know about actual grown-ensemble membership, not just
  geometric coverage; not attempted in this pass. `points.csv`'s
  `low_confidence_ensembles` column reports, rather than closes, this gap: it shows which of
  these orphaned points still have a low-confidence fallback available (a nearby seed's
  iteration-1 mask) versus which have none at all (empty list).

## Where the actual algorithm lives

This directory only *drives and visualizes* STC -- for the algorithm itself:

- `misc/mod_STC.md` -- the authoritative specification.
- `src/tox/shape_truthful_clustering/` -- the implementation (`tox_shape_truthful_clustering_impl.F90`
  for `ensemble_identification`(`_merged`), the sibling `_seeding_impl.F90`,
  `_ensemble_growing_impl.F90`, `_observable_impl.F90`, `_accept_impl.F90`,
  `_parameter_estimation_impl.F90` (`estimate_stc_parameters` and its helpers),
  `_reconciliation_impl.F90` modules for everything it orchestrates), and
  `src/tox_stc_json.F90`/`src/tox_stc_csv.F90` for the output writers `stc_cli` calls.
- `C-layer/` -- `stc_cli` itself; see `C-layer/README.md`.
- `test/mod_test_shape_truthful_clustering*.F90` (and their `python/test/`/`r/test/`
  mirrors) -- the unit tests these experiments complement, not replace.
