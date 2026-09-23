# Shatter Specifications

Working notes on the current Shatter (`tox_shatter_cluster_data`, `src/tox_shatter_cluster_data.F90`)
implementation, kept separate from the methodology writeup
(`misc/Tensor_Omics_Methods/Shatter_Clustering.tex`) so it can track engineering detail --
memory layout, complexity, optimization -- without polluting the math spec.

Notation used throughout: $N$ = `n_vectors` (ambient point count), $D$ = `n_dimensions`,
$S$ = `n_seeds`, $W$ = `t_observables` (observable history window, default 10). Byte sizes assume
gfortran's actual default kinds: `real64` = 8 bytes, `int32` = 4 bytes, default `logical` = 4 bytes
(the last one is compiler-defined, not guaranteed by the Fortran standard, but it is what gfortran
actually uses and every size below assumes it).

# Memory Layout: Single-Seed Ensemble Growth

`grow_single_seed_helper` is the unit of work: it grows exactly one seed into one ensemble via a
`do while` loop that alternates `grow_ensemble_helper` (surface expansion) and
`accept_ensemble_helper` (stability check on the observable history). Everything below is the
data this one call touches, per call.

## Shared, read-only inputs (not owned by the growth call)

These are computed once for the whole dataset and passed down by reference into every seed's
growth, seed-selection call, and the eventual merge step:

- `vectors(D, N)`, `real64` -- the ambient point cloud, $8DN$ bytes.
- `density_labels(N)`, `real64` -- one density scalar per point (neighbor count within the
  density radius, from `calculate_labels_as_density_helper`), $8N$ bytes.
- `dimension_order(D)`, `int32` -- the K-D tree's per-depth split-axis sequence, $4D$ bytes.
- `kd_indices(N)`, `int32` -- the K-D tree's index permutation (implicit balanced tree: a
  contiguous `[left, right]` range plus a depth is a subtree; no explicit node/pointer structure
  exists at all), $4N$ bytes.

## Per-seed state carried across growth iterations

- `current_mask(N)`, `logical` -- which of the $N$ points are currently ensemble members. This is
  the entire "ensemble" data structure: there is no compact member list, only a dense $O(N)$ mask
  regardless of how small the ensemble actually is.
- `backup_mask(N)`, `logical` -- a full copy of `current_mask` taken at the top of every iteration,
  restored verbatim if the iteration's growth step is later rejected by `accept_ensemble_helper`.
  This is the only rollback mechanism; there is no incremental undo.
- `tmp_observables(5, W)`, `real64` -- the sliding-window observable history
  (`CM_OBSERVABLE_COUNT = 5`: arithmetic-mean density, harmonic-mean density, their ratio as a
  heterogeneity score, active ensemble size, candidate count added this step). Once `iter > W`
  (default $W=10$), every write shifts all $W-1$ prior columns left by one
  (`observables(:, 1:W-1) = observables(:, 2:W)`) before writing the new column -- an $O(5W)$
  copy on every iteration past the tenth, forever, for a window whose only actual read is the
  last two columns (current vs. previous, inside `accept_ensemble_helper`).

## Per-iteration scratch, allocated once per seed and reused destructively every iteration

- `tmp_stack(3, 64, N)`, `int32` -- the K-D tree traversal stack
  (`KD_STACK_ENTRY_SIZE=3` x `KD_TRAVERSAL_STACK_DEPTH=64`), one full copy *per point*, i.e.
  $768N$ bytes, even though only one traversal runs at a time inside a given seed's growth (the
  extra $N$ multiplicity exists solely so `grow_ensemble_helper`'s `do concurrent` over points can
  give each iteration of that inner loop its own private stack slice
  `tmp_stack(:, :, i_vec)` -- see below).
- `tmp_vicinity_mask(N, N)`, `logical` -- a full $N \times N$ dense matrix, $4N^2$ bytes. Column
  `i_vec` holds the boolean vicinity-query result *from* member `i_vec` *to* every other point.
  Two things about this buffer matter for its actual memory behavior:
  - It is **entirely overwritten every single growth iteration**, active columns and inactive
    columns alike (`vicinity_vectors_helper` itself starts every call with
    `vicinity_mask = .false.`, an $O(N)$ memset, and `grow_ensemble_helper` explicitly zeroes the
    columns of non-members too). Nothing about its contents ever persists or is reused across
    iterations -- it is logically a per-iteration scratchpad shaped like a persistent matrix.
  - Its $N^2$ shape is sized for the *global* point count, not the ensemble's current member
    count or its (typically much smaller) neighbor count; a 1000-point dataset allocates
    4 MB for this buffer alone, per seed, on every call.
- `tmp_perm(N)`, `int32` and `tmp_abs_diff(N)`, `real64` -- generic sort/percentile scratch,
  $4N + 8N$ bytes. Notably **reused for two logically unrelated computations within one call** to
  `grow_ensemble_helper`: first to sort all $N$ points' densities for the ambient median/MAD
  (steps 1-3), then repurposed (only the first `k = count(ensemble_mask)` slots, `tmp_abs_diff(1:k)`
  / `tmp_perm(1:k)`) to sort the *active members'* densities for the ensemble-center median
  (step 4). Same physical storage, two different logical arrays, back to back -- correct today
  only because the second use starts by overwriting every slot it reads.

## Output

- `out_mask(N)`, `logical` -- the accepted final membership mask, a copy of `current_mask` at
  loop exit.

## Batched layout: growing all seeds at once (`obtain_ensembles_helper`)

`obtain_ensembles_helper` parallelizes the above over all $S$ seeds with a `do concurrent` over
`i_seed`. Since `do concurrent` iterations must not alias, **every per-seed buffer above gets an
extra trailing dimension of size $S$**, allocated once for the whole batch:

- `tmp_stack(3, 64, N, S)`, `int32` -- $768NS$ bytes.
- `tmp_vicinity_mask(N, N, S)`, `logical` -- $4N^2S$ bytes. This is the single largest allocation
  in the entire pipeline by a wide margin, and it scales as $O(N^2 S)$: for $N=1000$, $S=50$ this
  is 200 MB, entirely of scratch that is discarded (recomputed from scratch) every growth
  iteration of every seed.
  - Also worth noting: since every seed's column-`i_vec` slice only holds a KD-tree query result
    for one point (values entirely determined by `vectors`, `r`, and the tree), independent seeds'
    slices along the third dimension can and do overlap in *value* even though they are stored in
    physically disjoint memory -- there is no sharing of identical vicinity-query results across
    seeds whose ensembles happen to pass through the same point at the same radius.
- `tmp_perm(N, S)`, `int32`, `tmp_abs_diff(N, S)`, `real64` -- $4NS + 8NS$ bytes.
- `tmp_observables(5, W, S)`, `real64` -- $40WS$ bytes.
- `tmp_current_mask(N, S)` / `tmp_backup_mask(N, S)`, `logical` -- $4NS$ bytes each.
- `ensemble_matrix(N, S)`, `logical` -- the actual output, $4NS$ bytes; everything else in this
  list is scratch that exists only to make the outer loop over seeds race-free.

## Reconciliation-stage layout (`merge_ensembles_helper`)

Runs after all $S$ raw ensembles exist. State: `merged_masks(N, S)` (starts as a copy of the raw
masks, $4NS$ bytes) and `tmp_active_flag(S)` (which columns are still "live" after being absorbed
into another). The algorithm itself is a repeat-to-fixpoint double loop over active column pairs
$(i, j)$, each pairwise test costing `count(merged_masks(:,i) .and. merged_masks(:,j))`, i.e. an
$O(N)$ mask AND-and-count -- so one full pass over all pairs is $O(S^2 N)$, repeated until no
merge happens in a full pass (worst case $O(S)$ passes, i.e. $O(S^3 N)$ total).

# Dominant terms, at a glance

| Structure | Shape | Scaling | Persists across iterations? |
|---|---|---|---|
| `vicinity_mask` (batched) | $N \times N \times S$ | $O(N^2S)$ | No -- fully rewritten every iteration |
| `tmp_stack` (batched) | $3 \times 64 \times N \times S$ | $O(NS)$ | No -- pure per-call scratch |
| `current_mask` / `ensemble_matrix` | $N \times S$ | $O(NS)$ | Yes -- the actual state |
| `merge_ensembles` working set | $N \times S$ | $O(NS)$ space, $O(S^3N)$ time | Yes, mutated in place |
| `tmp_observables` | $5 \times W \times S$ | $O(WS)$ | Yes, but only last 2 columns are ever read |

$N^2S$ dominates memory; the same recompute-from-scratch-every-iteration pattern behind it also
dominates time, since every growth iteration re-queries the KD-tree for every current member
instead of only the members added in the previous step.

# Per-seed growth radius

Growth no longer shares one scalar radius across all seeds. `identify_ensemble_seeds_helper`
already computes, for each selected seed, the chosen percentile (`coverage_quant`, default
`CM_SEEDING_COVERAGE_PERCENTILE = 0.5`) of that seed's `k_seeding` nearest-neighbor distances, and
uses it to mark the region the seed covers. That same value is now also returned, in
`seed_radii(N)` -- one shared real array indexed by *vector* index, holding the radius at each
seed's own index and `0` everywhere else. Indexing by vector rather than by selection rank keeps
it order-independent: a caller that builds `seed_indices` in any order maps through it as
`seed_radii(seed_indices)`.

`obtain_ensembles_alloc` takes that array as the optional `seed_radii(S)`, aligned with
`seed_indices`, and it takes precedence over the uniform `r`. Below the allocating layer the array
is the only contract: `obtain_ensembles` and `obtain_ensembles_helper` take `seed_radii(S)` instead
of a scalar, and `grow_single_seed_helper` -- which already grows exactly one seed -- receives
`seed_radii(i_seed)` as its `r`. So seed $i$ and every later expansion of the ensemble around it
query the K-D tree at the radius fitted to seed $i$'s own neighborhood. When `seed_radii` is absent
the allocating layer broadcasts the scalar `r` (or the computed density radius) into the array, so
the previous uniform-radius behavior is exactly what a constant array produces.

Cost: $8S$ bytes, and nothing else -- no extra queries, since the radii fall out of work seeding
already does.

## What this buys, and what it costs

- Dense regions grow tightly and sparse regions generously, instead of one radius being
  simultaneously too coarse for the former and too tight for the latter.
- Growth is no longer symmetric: seed $i$ reaches point $p$ at $r_i$ while seed $j$ does not reach
  it at $r_j < |p - x_j|$. Membership therefore depends on which seed a region is grown from, and
  the reconciliation stage (`merge_ensembles_helper`) is what puts the asymmetric raw ensembles
  back together.
- A seed sitting in a sparse pocket inside an otherwise dense region gets a large radius and can
  overshoot into neighboring structure; the density-compatibility test (`alpha_mad`) is the only
  thing holding it back, and it is scaled by the *ambient* MAD, not a local one.
- The radii are estimated from `k_seeding` neighbors only, so they are noisy for small
  `k_seeding` -- the percentile choice trades that noise against how aggressively the radius
  tracks the local scale.
- A single global `r` is no longer a meaningful knob to tune or report: growth behavior is now a
  function of `k_seeding` and `coverage_quant` instead.
