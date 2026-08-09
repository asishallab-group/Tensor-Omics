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
- Optional, for 3D reports: the `scatterplot3d` package (pure R, no compiled system
  dependencies -- `install.packages("scatterplot3d")`). Without it, 3D+ datasets still get
  their normal 2D report, just no `<prefix>_3d.pdf` -- see "3D reports" below.
- The `gridExtra` package (pure R, no compiled system dependencies --
  `install.packages("gridExtra")`), for the parameters-table page.
- Nothing extra for the interactive HTML report: D3 v7 is vendored in-repo
  (`vendor/d3.v7.min.js`), and the output is a single self-contained `.html` file you open
  directly in a browser -- no server, no npm, no build step at view time.

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

produces one PDF per combination (2×2×2×2×2×1 = 16 here). Two more optional, trailing lists
sweep `seeds`' own parameters -- `k_density` (its neighborhood size, independent of `k_min`)
and `exclusion_radius_percentile` (its coverage/exclusion radius, independent of the
growth-phase radius, see "Known limitations" below):

```bash
misc/STC-experiments/run_stc_experiments.sh misc/STC-experiments/data/2d/kinked_curve_2d_noise_medium.csv \
    30 30 1 3.0 merge_jsi 0.1 15,30 25,50
#   ^k_min ^alpha ^d_max ^g_max ^mode ^min_jsi ^k_density ^exclusion_radius_percentile
```

Before this was wired up, `k_density` and `exclusion_radius_percentile` were not
forwardable at all -- every run silently used `seeds`' hardcoded defaults (`k_density=30`,
`exclusion_radius_percentile=50.0`, the median) regardless of what a results filename's
leading `k...` label said (that label only ever moved `k_min`, the growth-phase parameter).
Filenames produced from here on include both explicitly (`..._k30_kd30_..._erp50.0`) so this
cannot silently happen again.

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
run_stc.py  (per CSV, per parameter combination)
    seeds -> ensemble_identification_merged -> ensemble_reconciliation
             (-> estimate_stc_parameters, informational, with --estimate-parameters)
                                    |
                                    v
              results/data/<prefix>_{points,membership,low_confidence_membership,
                                      ensembles,super_ensembles,super_ensembles_jsi,
                                      ensemble_jsi_matrix,params}.csv
                                    |
                          +---------+---------+
                          v                   v
        plot_stc.R  -->  <prefix>.pdf   export_json.py --> <prefix>.json
        (+ <prefix>_3d.pdf                                       |
         for 3+ ambient dims)                                    v
                                              render_interactive.py --> <prefix>_interactive.html

run_stc_experiments.sh orchestrates the PDF and interactive-HTML paths for a single set of
parameters, sweeping parameter lists like branch `smoothing`'s run_lomanle_tests.sh did for
LoManLe. run_stc_pair.py instead runs run_stc.py twice on the same input -- once with the
parameters you give it, once with estimate_stc_parameters' own proposal actually applied --
and produces both PDF+interactive-HTML reports for each, side by side (see "run_stc_pair.py"
below).
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
| `--k` | shortcut: sets both `--k-min` and `--k-density` at once (either still overrides individually) | -- |
| `--k-min` | growth radius: k-NN pool size for the median-distance radius | 30 |
| `--k-density` | seeding: k-NN pool size for both `density_labels`' adaptive bandwidth and `seeds`' coverage radius | 30 |
| `--alpha-max-deg` | `accept_ensemble`'s max principal angle, in degrees (radians internally) | 30 |
| `--d-max` | `accept_ensemble`'s max tolerated change in intrinsic dimension | 1 |
| `--g-max` | `accept_ensemble`'s max tolerated `\|log(G_tp1/G_t)\|` | 3.0 |
| `--f-max` | Stop Condition 1's ensemble-size-fraction ceiling | 0.95 |
| `--a` | Stop Condition 2's "stably accepted" threshold | 2 |
| `--o` | trailing observable-history window depth | 10 |
| `--reconciliation-mode` | `report` / `merge_jsi` / `merge_any` | `merge_jsi` |
| `--min-jsi` | minimum JSI for `merge_jsi` | 0.1 |
| `--max-group-size` | cap on ensembles per super-ensemble | `min(1024, n_ensembles)` |
| `--exclusion-radius-percentile` | `seeds`' own coverage/exclusion radius quantile | 50.0 |
| `--bandwidth-percentile` | `density_labels`' local KDE bandwidth quantile | 68.27 |
| `--estimate-parameters` | also run `estimate_stc_parameters` and report its estimates (informational only, never auto-applied -- see `misc/mod_STC.md`, "Estimate parameters from data") | off |
| `--n-anchors` / `--seed-max-set-size` / `--first-quartile-percentile` | `estimate_stc_parameters`'s own parameters, only used with `--estimate-parameters` | 5 / 5.0 / 25.0 |

It writes eight files per run (`<prefix>_points.csv`, `..._membership.csv`,
`..._low_confidence_membership.csv`, `..._ensembles.csv`, `..._super_ensembles.csv`,
`..._super_ensembles_jsi.csv`, `..._ensemble_jsi_matrix.csv`, `..._params.{json,txt}`) --
see "Output CSV schemas" below.

### `plot_stc.R`

```
Rscript plot_stc.R <prefix>
```

reads the files above and writes `<prefix>.pdf`, a 7-10 page report (some pages are skipped
when there is nothing to show, e.g. no ensemble reaches intrinsic dimension 2, no
super-ensemble was found at the current reconciliation threshold, or `--estimate-parameters`
was not passed):

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
   ensembles that did not qualify at the current mode/threshold. Grouping only -- deliberately
   *not* drawn as an ordered chain/skeleton, see "No stitched manifolds" below.
7b. **JSI heatmap** -- every non-empty ensemble pair's Jaccard similarity, from
    `<prefix>_ensemble_jsi_matrix.csv`, as a `geom_tile()` heatmap (viridis fill, diagonal
    forced to 1.0). Unlike Page 7, this includes pairs that never reached the reconciliation
    threshold, so it is the "before thresholding" view Page 7 is a slice of.
8. **Low-confidence fallback coverage** -- every point colored by one of three categories: in
   a retained ensemble; orphaned by the retained pipeline but reachable by some seed's
   iteration-1 `low_confidence_mask` (a fallback LoManLe may choose to use, see
   `misc/mod_STC.md`, "Ensemble identification", "Output"); or genuinely uncovered by
   anything.
9. **Parameters table** -- only when the run used `--estimate-parameters`: two side-by-side
   `gridExtra::tableGrob` tables (not a plot), the input parameters this run actually used and
   `estimate_stc_parameters`' own proposal (`k_min`/`k_density`/`density_quantile`/
   `alpha_max`/`G_max`/`d_max`), for a direct visual diff.

### No stitched manifolds

`plot_stc.R` never draws a super-ensemble's members as a connected, ordered path/skeleton --
`ensemble_reconciliation` only ever reports *which* ensembles intersect enough to group, not
in what order or with what geometry they connect. An earlier version of this script did draw
such a chain (via `super_ensembles_JSI`'s per-column consecutive-pair JSI values), and it was
actively misleading: `super_ensembles`' member order within a group is ensemble-*discovery*
order (highest density first), not spatial order, so the "chain" zigzagged across the data
with no relationship to the actual manifold shape. Turning a set of intersecting ensembles
into one stitched manifold is LoManLe's job, downstream of STC, not something STC's own
report should imply it has already done.

### 3D reports

For CSVs with 3 or more ambient dimensions, `plot_stc.R` also writes `<prefix>_3d.pdf`: the
same input/clusters/low-confidence-coverage/tangent-direction pages as the 2D report, but
actually using all 3 (of the first 3) dimensions, via the `scatterplot3d` package. Before
this, 3D data was silently only ever plotted in its first two dimensions -- the third
ambient axis was computed and clustered on correctly, just never drawn.

The original intent (matching branch `smoothing`'s `r/plot_lomanle_spheres.R`,
`plot_3d_report`) was an **interactive** `plotly` + `htmlwidgets` HTML export, rotatable in a
browser. That was not achievable in the environment this was written in: `plotly`,
`htmlwidgets`, and `rgl` all transitively need compiled system libraries (`libcurl`,
`libssl`, `libuv`) whose `-devel` headers were not installed, installing them needs `sudo`,
and that environment's own `conda` was broken (permission denied on the `conda` executable
itself) -- all confirmed directly, not assumed. `scatterplot3d` is pure R with no compiled
system dependency, so it installs and runs anywhere R does; the trade-off is a fixed viewing
angle per page instead of a rotatable one. If `plotly`/`htmlwidgets` become installable in
your environment (`sudo dnf install libcurl-devel openssl-devel libuv-devel` on
Fedora/RHEL, then `install.packages(c("plotly","htmlwidgets"))`), the 3D report section at
the bottom of `plot_stc.R` is the one to replace with the original interactive approach.

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

Orchestrates the pipeline above: builds the project, generates datasets if `data/` does
not exist yet, then loops `run_stc.py` + `plot_stc.R` + `render_interactive.py` over every
combination of the parameter lists you give it (or every dataset, with `all`) -- the same
"one CSV or `all`, comma-separated parameter lists, nested loop" shape as branch
`smoothing`'s `run_lomanle_tests.sh`. `--estimate-parameters` is always passed to
`run_stc.py`, so Page 9 (parameters table) is always present. Run it with no arguments for
the usage message.

### The interactive HTML report

The 2D pages of `plot_stc.R`'s PDF report (Pages 1-7b: point cloud/seeds, clusters, overlap
density, growth radii, both tangent directions, super-ensembles, JSI heatmap) also exist as
one interactive HTML page you can open in a browser, `<prefix>_interactive.html` -- built
from three pieces:

- **`export_json.py <prefix>`** -- reads the same `<prefix>_*.csv`/`_params.txt` files
  `plot_stc.R` does and consolidates them into one `<prefix>.json`: per-point
  `ensembles`/`low_confidence_ensembles`/`seed_of` lists precomputed (not left for the
  browser to join at view time), per-ensemble `super_ensemble_id`, the full JSI matrix, and
  `params`.
- **`interactive_template.html`** -- the reusable part; a self-contained D3 v7 page with
  `__STC_D3_JS__`/`__STC_DATA__` placeholders. Controls panel (point color mode, which single
  ensemble to highlight, layer checkboxes for growth-radius circles/tangent segments/JSI
  heatmap) plus a hover tooltip showing every relevant fact about a point (seed status,
  which ensemble(s) and super-ensemble(s) it belongs to, low-confidence coverage). Replaces
  the old approach of drawing every information layer at once and letting them overplot each
  other -- exactly one ensemble-point relation is highlighted at a time, chosen from the
  dropdown. Edit this file to change the visualization; nothing else needs to change per
  report.
- **`render_interactive.py <prefix>`** -- combines the template, the vendored D3 bundle
  (`vendor/d3.v7.min.js`), and `<prefix>.json` into one self-contained `<prefix>_interactive.html`
  (calling `export_json.py` first if the JSON does not exist yet). The data is inlined at
  generation time rather than `fetch()`-ed at view time, because a plain `file://` page's
  `fetch()` of a sibling JSON file is blocked by browsers' CORS rules for local files -- this
  way the output is one file, directly double-clickable, no local server needed.

`run_stc_experiments.sh` and `run_stc_pair.py` both call `render_interactive.py`
automatically; call it by hand (`python3 render_interactive.py <prefix>`) after a bare
`run_stc.py` invocation if you want the interactive report for a run that didn't go through
either wrapper.

### `run_stc_pair.py`

```
python3 run_stc_pair.py <input.csv> [run_stc.py options for the "original" run]
```

Runs `run_stc.py` twice on the same input: once with the parameters you give it (the
"original" run, with `--estimate-parameters` always added so its estimate is available),
once with that estimate's `k_min`/`k_density`/`alpha_max_deg`/`g_max`/`d_max` actually
applied as real run parameters (the "estimated" run) -- not just reported, unlike
`--estimate-parameters` on its own, which never feeds back into the same run.
`density_quantile` has no direct CLI flag (see the options table's "no equivalent" note) and
is skipped; `reconciliation_mode`/`min_jsi`/`max_group_size` are held fixed across both runs,
since the estimator does not propose values for them. Writes
`<out-prefix>_original_*`/`<out-prefix>_estimated_*` (PDF, 3D PDF where applicable, JSON,
interactive HTML) side by side, so the two can be compared directly -- this is how "run two
experiments per dataset case: original input parameters, and the estimator's own proposal
actually applied" is done.

## Output CSV schemas

For reuse outside `plot_stc.R` (a different plotting tool, a notebook, ad-hoc analysis).
`<prefix>` is whatever `--out-prefix` (or `run_stc_experiments.sh`'s own naming) chose.

- **`<prefix>_points.csv`** (wide, one row per input point): `point_id`, one column per
  ambient dimension (named after the input CSV's own numeric columns, e.g. `x1`,
  `y_original`), `n_ensembles` (how many ensembles this point belongs to),
  `n_low_confidence_ensembles` (how many seeds' iteration-1 low-confidence mask reaches this
  point, see `low_confidence_membership.csv` below).
- **`<prefix>_membership.csv`** (long, one row per point that belongs to at least one
  ensemble, repeated once per ensemble it belongs to): `point_id`, `ensemble_id`, `is_seed`
  (1 if `point_id` was this ensemble's own seed).
- **`<prefix>_low_confidence_membership.csv`** (long, one row per point covered by some
  seed's iteration-1 bootstrap mask, repeated once per such seed): `point_id`, `ensemble_id`.
  Reported for every seed regardless of `stop_reason`, including ones Stop Condition 1 later
  discarded entirely -- see `misc/mod_STC.md`, "Ensemble identification", "Output",
  `low_confidence_mask`.
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
- **`<prefix>_ensemble_jsi_matrix.csv`** (long, one row per unordered pair of non-empty
  ensembles): `ensemble_id_1`, `ensemble_id_2`, `jsi`. The *full* pairwise matrix -- including
  pairs that never qualified for reconciliation at all (`jsi=0`), unlike
  `super_ensembles_jsi.csv`, which only ever reports consecutive pairs within an already-formed
  group. Feeds the JSI heatmap (Page 7b, see below).
- **`<prefix>_params.json`** / **`<prefix>_params.txt`**: every parameter `run_stc.py` ran
  with, plus `n_vectors`/`n_dimensions`/`n_ensembles`, plus (only with
  `--estimate-parameters`) `estimated_k_min`/`estimated_k_density`/
  `estimated_density_quantile`/`estimated_alpha_max_deg`/`estimated_g_max`/`estimated_d_max`.
  Two formats for the same content -- `.txt` (plain `key=value` lines) is what `plot_stc.R`
  reads, so it needs no JSON-parsing package; `.json` is there for anything that wants to
  parse it programmatically.

## Known limitations

- **The interactive HTML report only plots the first 2 ambient dimensions.** Unlike
  `plot_stc.R`'s own separate `<prefix>_3d.pdf` (see "3D reports" below),
  `interactive_template.html` draws one 2D scatter (`DATA.dim_names[0]`/`[1]`) regardless of
  how many ambient dimensions the dataset has -- there is no interactive 3D report. All
  ambient coordinates are still present in `<prefix>.json`'s `points[].coords`, so a 3D-aware
  extension to the template is possible; not attempted in this pass.
- **The 3D report is static, not interactive.** See "3D reports" above -- `scatterplot3d`
  gives a genuine, correct 3D plot (all 3 axes actually used, not a 2D projection), but each
  page is one fixed viewing angle rather than something you can rotate in a browser the way
  the original `plotly`-based plan intended. Datasets with 4+ ambient dimensions only ever
  plot their first 3.
- **`data/` and `results/` are gitignored**, not committed -- unlike branch `smoothing`,
  which committed its generated `results/data/*.csv`. Both are fully reproducible
  (`generate_datasets.py` is seeded; `run_stc.py`/`plot_stc.R` are deterministic given the
  same inputs and parameters), so nothing is lost by regenerating on demand, and the
  repository stays smaller for it -- there are ~120 generated dataset CSVs alone.
- **No exhaustive default sweep has been run and committed.** This directory is a tool for
  exploring STC's behavior, not a benchmark suite with checked-in reference output; running
  `all` yourself, with whatever parameters you are currently investigating, is the intended
  workflow.
- **Seed-exclusion suppression is purely geometric and does not know about curvature.**
  `seeds`' coverage/exclusion radius (tunable via `exclusion_radius_percentile`, see "Quick
  start" above) can suppress seed placement across a region -- typically a curvature extremum
  (a peak, trough, or kink on a wavy manifold) -- that the suppressing seed's own later growth
  never actually reaches, because `accept_ensemble`'s curvature-based stop conditions
  (`alpha_max_deg`/`g_max`/`d_max`) can halt growth well short of that same geometric radius.
  The result is points with no seed of their own and no membership in any grown ensemble
  either. Shrinking `exclusion_radius_percentile` reduces how much territory each seed
  suppresses, but does not close the gap on its own -- it only shrinks it. Fixing this
  properly needs seed selection to know about actual grown-ensemble membership, not just
  geometric coverage; not attempted in this pass. The "Low-confidence fallback coverage"
  page (see `plot_stc.R` above) reports, rather than closes, this gap: it shows which of
  these orphaned points still have a low-confidence fallback available (a nearby seed's
  iteration-1 mask) versus which have none at all.

## Where the actual algorithm lives

This directory only *drives and visualizes* STC -- for the algorithm itself:

- `misc/mod_STC.md` -- the authoritative specification.
- `src/kernel/shape_truthful_clustering/` -- the implementation (`tox_shape_truthful_clustering_kernel.F90`
  for `ensemble_identification`(`_merged`), the sibling `_seeding_kernel.F90`,
  `_ensemble_growing_kernel.F90`, `_observable_kernel.F90`, `_accept_kernel.F90`,
  `_parameter_estimation_kernel.F90` (`estimate_stc_parameters` and its helpers),
  `_reconciliation_kernel.F90` modules for everything it orchestrates).
- `test/mod_test_shape_truthful_clustering*.F90` (and their `python/test/`/`r/test/`
  mirrors) -- the unit tests these experiments complement, not replace.
