# Idea Box

Notes captured from an end-of-day design discussion (2026-08-05) between Asis and Claude, on STC and LoManLe more broadly: hoped-for novelty, performance goals, and a recursive surface-learning extension. None of this is implemented or validated yet -- these are ideas to revisit, not decisions.

## Novelty -- STC and LoManLe

Hope, on both sides: STC's combination (adaptive per-seed local-scale radius, RG-style joint stopping rule on a tangent/dimension/gap tuple, planned noise-tube-gated reconciliation) and LoManLe's combination (overlapping tangent charts, inverse-variance stitching, MST backbone with structural roles) are not aware of an exact match in the published literature.

Caveat, stated plainly: "not aware of prior art" is not the same as "confirmed novel." The honest next step is a direct benchmark against the closest published relatives once both are implemented -- ElPiGraph, SimplePPT, Slingshot, and Atlas/IAN/LTSA (already tracked in `misc/lomanle_literature_resume.md`) -- rather than assuming novelty from absence of a known match.

## Performance goals

Fortran core with a robust C ABI, exposed to R and Python (and potentially other languages), following the project's existing three-layer architecture (pure SK core -> `_C` `bind(C)` wrapper -> Python/R). Coarrays are a good structural fit for STC specifically, since seeds grow independently -- the same embarrassingly-parallel shape that currently maps to `!$omp parallel do` across seeds could instead map to coarrays for distributed-memory scaling across nodes, if that is ever needed beyond single-node core-count scaling.

Caveat: coarray support quality varies significantly by compiler (gfortran's is functional but less mature/performant than Intel's ifx). `build_utils.sh` already supports gfortran/ifx/nvfortran; decide which coarray runtime is actually being targeted before writing code that depends on it.

## Candidate selection in atomic growth steps

Concern raised: STC's density-instantiation atomic growth (add all qualifying candidates simultaneously, no ordering) can take one large, uncontrolled jump per iteration.

Idea: jackknife the candidates -- leave-one-candidate-out perturbation of the observable, to catch a candidate that would destabilize the fit before adding it. Principled, but costly: naively `O(|candidates|)` extra observable evaluations per growth step.

Agreed cheaper alternative: sort candidates by some criterion (e.g. ascending distance to the current tangent plane) and add them in batches of `n_batch`, rather than either fully atomic (`n_batch = |candidates|`) or fully sequential (`n_batch = 1`). Sorting by geometric proximity is a cheap proxy for what jackknifing actually measures (each candidate's effect on observable stability), not equivalent to it, but should catch the worst case (a single step bridging two branches) reasonably well.

Suggested calibration: run both extremes (`n_batch=1` and `n_batch=|candidates|`) on real test data and observe where results stop changing as the batch size shrinks -- picks `n_batch` empirically rather than by guessing.

## Recursive surface learning: feeding LoManLe's residuals back into LoManLe

Idea: within each accepted ensemble, take the points in the top 90th percentile of (orthogonal) distance to the fitted local manifold -- the worst-fit points, i.e. the ones farthest from the tangent-plane approximation. Pool these across all ensembles into a "surface cloud," and run LoManLe again on that cloud, this time targeting a higher intrinsic dimension (`manifold_dim > 1`) to learn the secondary surface structure that the primary pass treated as noise/residual.

This is conceptually a multi-resolution decomposition in the same renormalization-group spirit as STC itself: each pass extracts the dominant scale-invariant structure at its level, and the residual becomes the next level's input, rather than being discarded.

### This directly motivates LoManLe's unimplemented d>1 topology

Step 10's stitching already projects correctly onto a `manifold_dim`-dimensional tangent plane for any `manifold_dim`, but Phase XII's backbone construction is 1-D by construction regardless of `manifold_dim` -- it threads one ordering coordinate along a single tangent direction, with no second axis or face structure. This has been an open, unforced item. The recursive-surface-learning idea gives it a concrete forcing use case (the second pass explicitly needs a real surface topology, not another 1-D chain) and is a reason to reprioritize it.

### Open structural question: does the pooled residual cloud actually cohere?

Whether the top-90th-percentile points from *different* ensembles form one connected, analyzable surface sample, or are scattered/unrelated per-ensemble outlier tails, depends entirely on whether the residual structure is itself consistent along the primary manifold. Intuition: a tube around a curve -- cross-sectional "outliers" at every point along the curve, pooled, genuinely trace the tube's surface, because the secondary structure has constant character everywhere. If the true residual is just isotropic noise with no such consistency, the second pass should, if the stability criteria are well-calibrated, simply fail to grow any stable ensemble on the residual cloud -- correctly reporting "no coherent secondary structure" rather than hallucinating one.

### Confounds to rule out before trusting this on real data

- **Boundary artifacts.** Local covariance is already documented as becoming asymmetric near manifold boundaries (a stated LoManLe limitation). An ensemble sitting at an edge will show elevated residuals for a purely geometric reason, unrelated to any genuine secondary structure.
- **Concentration of measure in high ambient dimension.** A normal-space residual norm is a sum over `D - d` roughly-independent components; by concentration of measure, that norm's distribution can be fairly narrow. "Top 90th percentile" may just sample the tail of ordinary noise rather than isolating something real, and this gets worse as ambient dimension `D` grows.

### Suggested validation plan

- **Positive control:** synthetic data with a *known*, deliberately injected secondary structure (e.g. a curve with a controlled orthogonal tube/wiggle). Confirm the second LoManLe pass recovers exactly that structure.
- **Negative control:** synthetic primary trajectory plus genuinely unstructured, isotropic residual noise (no injected secondary structure). Confirm the second pass correctly refuses to converge on anything -- i.e. does not hallucinate structure that isn't there.

Both controls mirror how bifurcation detection was already validated (`bifurcation_2way`/`bifurcation_3way` synthetic datasets) -- the same discipline should apply here before trusting this on real biological data, where there is no ground truth to check against.
