# Shatter Specifications

Engineering notes for `tox_stc`, generated from `src/tox/tox_stc_impl.F90`.
The methodology is described in `misc/Tensor_Omics_Methods/Shatter_Clustering.tex`.

## Generated API

Prepare the data with `calculate_density_radius`, `calculate_labels_as_density`, and
`identify_ensemble_seeds`. Seed radii are indexed by input point; gather them in seed-index
order before calling the per-seed growth entry point.

- `obtain_ensembles_global_radius` grows every seed with the scalar radius `r`.
- `obtain_ensembles_seed_radius` uses one radius per seed.
- Both return one raw mask per seed and a stop reason for each seed. The raw ensemble count
  is the number of seed indices supplied, or the number of columns of `ensemble_matrix`;
  it is not a separate output.
- `grow_ensemble` performs one unconditional density-compatible surface step, taking separate
  input and output masks.
- `merge_ensembles` combines raw masks. Python and R return only surviving columns;
  Fortran and C also return the surviving count.

The standalone `compute_ensemble_observable` and `accept_ensemble` routines share one window:
the column count of the supplied history matrix. Both require at least two columns.
The `obtain_ensembles_*` window parameter `t_observables` must also be at least two. Its
workspace width is `max(2, min(t_observables, n_vectors))`. Use 10 for the benchmark setting.
A zero ambient median absolute deviation is an error when there are seeds to grow.

The generated wrappers validate input arguments. The implementation reports zero-MAD runtime
errors against `density_labels` using its own argument position. Per codegen_guide section 5.14,
generated wrappers clear implementation argument positions because their signatures differ.

Regenerate all bindings with `python3 helper/generate_code.py`; verify with `--check`.

## Growth and workspace layout

Let N be the number of vectors, D their dimension, S the number of seeds,
T = min(n_tiles, max(1, S)), and W = max(2, min(t_observables, N)).
The allocating API owns workspaces; expert callers supply the same shapes.
Inputs `vectors(D, N)`, density labels, dimension order, tree indices and seed indices are
shared read-only. Ambient density statistics are calculated once before processing the seeds.

`obtain_ensembles_impl` distributes seeds across T independent concurrent tiles. Each tile
processes its seeds serially with `grow_single_seed_helper`. The tiles reuse these workspaces:

| Workspace | Shape | Role |
|---|---|---|
| `tmp_stack` | 3 x 64 x T | One K-D traversal stack per tile |
| `tmp_vicinity_mask` | N x T | Result of one neighborhood query |
| `tmp_surface_mask` | N x T | Cached union of accepted members' neighborhoods |
| `tmp_perm` | N x T | Sorting scratch, reused to hold the next frontier |
| `tmp_abs_diff` | N x T | Density sorting/deviation scratch |
| `tmp_observables` | 5 x W x T | Sliding observable history |
| `tmp_current_mask` | N x T | Accepted membership for each tile |
| `ensemble_matrix` | N x S | Raw membership output; also holds trial membership during growth |

Each seed starts with its own point. Growth queries only the frontier of newly accepted members
and adds their neighborhoods to the cached surface. `build_ensemble_trial_helper` preserves all
accepted members and adds density-compatible candidates. The observable update records the trial,
then acceptance compares consecutive heterogeneity values. Rejected trial membership is discarded;
accepted membership becomes the next current mask. No backup-mask array or N-by-N vicinity matrix
is allocated. At exit, the raw output contains the last accepted mask.

The five observables are arithmetic-mean density, harmonic-mean density, their ratio,
ensemble size, and candidate count. A full history shifts left by one column on each update.
Growth workspace scales as O(NT + WT), in addition to O(NS) raw output and O(DN) input storage.
Logical masks use `logical(c_bool)`; their storage size is platform dependent.

## Reconciliation

`merge_ensembles_impl` starts with raw masks, uses union-find parents and active-column flags
to merge overlapping ensembles, then compacts surviving columns and clears the unused tail.
The returned merged count may be smaller than S and remains part of the merge API.

# Per-seed growth radius

Growth can use one radius per seed instead of one radius shared by all seeds. `identify_ensemble_seeds`
already computes, for each selected seed, the chosen percentile (`coverage_quant`, default `0.5`)
of that seed's `k_seeding` nearest-neighbor distances, and uses it to mark the region the seed
covers. It also returns that value in `seed_radii(N)` -- one real array indexed by *vector* index,
holding the radius at each seed's own index and `0` everywhere else. Indexing by vector rather
than by selection rank keeps it order-independent: a caller that builds `seed_indices` in any
order gathers the growth radii as `seed_radii(seed_indices)`.

The growth mode is chosen by entry point, not by a runtime flag. `growth_mode` on
`obtain_ensembles_impl` has a procedure column, so the generator emits one entry point per mode
and fixes the mode inside each:

- `obtain_ensembles_global_radius` takes the scalar `r` and grows every seed at that radius.
- `obtain_ensembles_seed_radius` takes `seed_radii(S)`, aligned with `seed_indices`, and grows
  seed $i$ -- and every later expansion of its ensemble -- at `seed_radii(i)`.

Each entry point requires its own radius argument and does not have the other's, so neither can
be called with the wrong one. Inside the implementation the mode only chooses each seed's radius --
`seed_radii(i)` or `r` -- and growth is identical from there, so the global mode is exactly the
constant-array case of the per-seed one, with no radius array allocated for it.

Cost: $8S$ bytes, and nothing else -- no extra queries, since the radii fall out of work seeding
already does.

## What this buys, and what it costs

- Dense regions grow tightly and sparse regions generously, instead of one radius being
  simultaneously too coarse for the former and too tight for the latter.
- Growth is no longer symmetric: seed $i$ reaches point $p$ at $r_i$ while seed $j$ does not reach
  it at $r_j < |p - x_j|$. Membership therefore depends on which seed a region is grown from, and
  the reconciliation stage (`merge_ensembles`) is what puts the asymmetric raw ensembles
  back together.
- A seed sitting in a sparse pocket inside an otherwise dense region gets a large radius and can
  overshoot into neighboring structure; the density-compatibility test (`alpha_mad`) is the only
  thing holding it back, and it is scaled by the *ambient* MAD, not a local one.
- The radii are estimated from `k_seeding` neighbors only, so they are noisy for small
  `k_seeding` -- the percentile choice trades that noise against how aggressively the radius
  tracks the local scale.
- With `obtain_ensembles_seed_radius`, a single global `r` is no longer a meaningful knob to tune
  or report: growth behavior is a function of `k_seeding` and `coverage_quant` instead.
