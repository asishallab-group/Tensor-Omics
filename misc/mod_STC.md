---
title: |
    | Shape Truthful Clustering (STC)
    | Algorithm Definitions
author: Prof. Dr. Asis Hallab
date: \today
geometry: margin=2.0cm
numbersections: true
---

In this document we find the pseudocode definitions of the functions for the
Shape Truthful Clustering (STC) variant tangent space.

# Data Definitions

* Input multivariate vectors as an 2D Fortran real array.
* Logical mask defining the seeds.

# Numerical Linear Algebra: SVD vs. Eigendecomposition

This module never forms a covariance or Gram matrix explicitly and
eigendecomposes it. Forming $C=YY^\top$ (or $M^\top M$ for a comparison
between two bases) squares the condition number of the underlying data
relative to computing singular values of $Y$ (or $M$) directly, which
disproportionately degrades precision in the *small* eigenvalues -- exactly
the normal-space ones that `normal_error` and `tangent_scales` depend on.
This is not a question of which eigensolver to use: even LAPACK's most
robust symmetric eigensolver cannot recover precision already lost by
forming the squared matrix in the first place. The fix is to never form it:
compute singular values and vectors of the actual centered data (or
comparison) matrix directly.

Two distinct problem shapes occur in this module, with two different
recommended routines:

* **Full spectrum** (all $D$ singular values/vectors of a $D\times k$
  centered ensemble matrix $Y_\mathcal{E}$, needed by `observable`,
  `normal_error`, and `tangent_scales`): use `dgesdd` (divide-and-conquer
  SVD), economy mode (`JOBZ='S'`). It is the fastest LAPACK SVD routine for
  computing full singular vectors, and is standard/portable across any
  LAPACK/OpenBLAS build. It requires a genuine workspace *query* (call once
  with `LWORK=-1`, read the optimal size back from `WORK(1)`) rather than a
  closed-form minimum, plus an integer workspace array (`IWORK`, size
  $8\cdot\min(D,k)$) -- size both in the `_alloc` layer.
* **Small $d\times d$ matrix** (comparing two tangent bases, e.g.
  $M=U_{\mathcal{E}_r}^\top U_{\mathcal{E}_{t+1}}$, for principal angles in
  `accept_ensemble`): use `dgesvd`, not `dgesdd`. At $d$'s typical size
  (the intrinsic dimension, usually 1-5), divide-and-conquer's speed
  advantage does not materialize and its extra `IWORK` bookkeeping is pure
  overhead for no benefit; `dgesvd`'s simpler single-workspace-array
  bookkeeping wins. Its singular values are $\cos\theta_j$ directly -- no
  squaring, no `sqrt` needed back out. `accept_ensemble` now performs up to
  $o+1$ such comparisons per growth step, one per retained reference basis
  (see "SKG `accept_ensemble`" below), not just one -- still all small
  $d\times d$ problems, so the routine choice is unaffected, only the count.

Do not use `dsyev`, `dsyevd`, or `dsyevr` anywhere in this module. They
solve "eigendecompose a given symmetric matrix as accurately as possible,"
which is not the problem this module has -- it can choose not to form that
symmetric matrix at all.

# Functions / Subroutines

## Coding / Naming convention

Follow `misc/Fortran_Coding_Guides.pdf` (the F42 standard) for general
Scientific Kernel (SK) philosophy (section 10.2: purity, no I/O, no
`allocate`/`deallocate` in the numerical core, determinism), and
`codegen_guide.md` for the naming scheme and how each procedure is written --
it supersedes F42 section 10.3.2's older `_alloc -> unsuffixed -> _helper`
suffix triad. Concretely:

* a kernel is a `pure` procedure named `<name>_kernel`, in a module named
  `tox_<family>_kernel` under `src/kernel/`: no validation, no `ierr` (except
  a genuine runtime failure no input check could foresee, e.g. an SVD that
  fails to converge -- `codegen_guide.md` section 5.14), no allocation
  anywhere in the module;
* every argument's constraints are stated as `!!`/`DM_*` documentation
  (`DM_MIN`, `DM_MAX`, `DM_SENTINEL`, `DM_DEFAULT`, `DM_OUTPUT_FROM`,
  `DM_PROLOGUE`, ...) -- that documentation is the generator's only input,
  and the *sole* place validation for that argument is stated;
* the generator derives the validating wrapper, the allocating wrapper
  (`_alloc`, published as the plain name), the expert tier (published as
  `_expert`), and the C/Python/R bindings from the kernel and its
  annotations -- none of those are hand-written.

Related routines following this pattern are called SK-groups (SKGs). See
"Implementation Modules" below for where each SKG's kernel lives.

## Implementation Modules

Each major step gets its own kernel module under
`src/kernel/shape_truthful_clustering/`, following `codegen_guide.md` section
5.15 ("a family too big for one file"):

| Module file | SKGs |
|---|---|
| `tox_shape_truthful_clustering_kernel.F90` (parent) | `ensemble_identification` (its own kernel(s), directly in the parent -- see note below) |
| `tox_shape_truthful_clustering_seeding_kernel.F90` | `density_labels`, `seeds` |
| `tox_shape_truthful_clustering_ensemble_growing_kernel.F90` | `calc_ensemble_growth_radius`, `grow_ensemble` |
| `tox_shape_truthful_clustering_observable_kernel.F90` | `observable`, `normal_error`, `tangent_scales` |
| `tox_shape_truthful_clustering_accept_kernel.F90` | `accept_ensemble` |
| `tox_shape_truthful_clustering_reconciliation_kernel.F90` | `ensemble_reconciliation` |
| `tox_shape_truthful_clustering_parameter_estimation_kernel.F90` | `sample_estimator_anchors`, `grow_estimator_anchor_clouds`, `estimate_stc_parameters` |

The parent module holds `ensemble_identification`'s own kernel(s) directly,
in addition to `use`-ing the five children -- a deliberate deviation from
section 5.15's own `tox_data_integration_kernel` example, whose parent holds
no procedures of its own. `ensemble_identification` is this family's natural
top-level entry point (it orchestrates every other SKG), unlike
`data_integration`'s parent, which really is just a bag of siblings with no
natural "top" member.

## Seeding

Seeding requires two dedicated SKGs: one to assign every vector a local
density label, the other to greedily select seeds from those labels.

### SKG `density_labels`

Receives an optional $k_{\text{density}}$ (default 30 -- the same default
`calc_ensemble_growth_radius`'s own $k_{\min}$ uses, see "Local Radius
Identification" below, but a genuinely independent argument: `density_labels`
and `seeds` are called before `ensemble_identification` in the pipeline, so
there is no single call through which one value could reach both) and an
optional $p_{\text{bw}}$, `bandwidth_percentile` (default 68.27). For each
vector $x_i$:

1. Find its $k_{\text{density}}$ nearest neighbors (excluding itself), giving
   distances $d_{i,1}, \ldots, d_{i,k}$.
2. $b_i = \text{percentile}_{p_{\text{bw}}}\big(d_{i,1}, \ldots,
   d_{i,k}\big)$ -- the $p_{\text{bw}}$-th percentile of those distances
   (linear interpolation, via `calc_percentile_helper`), floored at a small
   epsilon (coincident neighbors at distance 0 would otherwise give $b_i = 0$,
   a zero-bandwidth Gaussian below).
3. $\text{density}(x_i) = \dfrac{1}{b_i^D}\sum_{j=1}^{k}
   \exp\!\left(-\dfrac{d_{i,j}^2}{2\, b_i^2}\right)$, where $D$ is the
   ambient dimensionality.

This is an adaptive-bandwidth kernel density estimate: unlike a single
dataset-wide radius, $b_i$ shrinks in dense regions and grows in sparse ones,
so the resulting density labels reflect *local*, not global, density. The
$1/b_i^D$ normalization is not optional polish -- without it, $d_{i,j}/b_i$ is
scale-invariant, so a cluster and the same cluster uniformly stretched out
score identically, which defeats the purpose of a *density* label.

$p_{\text{bw}}$ is deliberately exposed as a tunable, undisguised heuristic,
not presented as a calibrated standard deviation. The 68.27 default is
anchored to the familiar 1-D-Gaussian fact "68.27% of the mass lies within
one SD" -- but that fact is 1-D-specific: for distances (norms) in $D$
dimensions under an isotropic-Gaussian-scatter assumption, distance$/\sigma
\sim \chi_D$, and the correct percentile for "one SD" shrinks with $D$
(39.35% at $D=2$, 19.88% at $D=3$, ...). Rigorously correcting for this would
require each point's *local intrinsic dimension*, not the ambient $D$ -- and
that is only available later, from `observable`'s per-ensemble SVD, not at
the seeding stage. Using ambient $D$ instead would actively over-correct for
STC's manifold-concentrated data (real local dimension is typically far
below $D$). So $p_{\text{bw}}$ stays a plain tunable parameter, not a
disguised statistical estimate.

### SKG `seeds`

Input:

* data_vectors - 2D real array

output:

* logical mask defining the selected seed points.

### Algorithm

First, assign each $x_i \in \text{data\_vectors}$ a density label via
`density_labels`.

Sort the density labels descending. Start with the highest-density unvisited
vector and mark it a seed. Its coverage radius is the *same computation*
`calc_ensemble_growth_radius` uses for its own growth radius -- a percentile
(`exclusion_radius_percentile`, default 50.0, the median) of the distances to
its own $k_{\text{density}}$ nearest neighbors -- called directly on the
newly-selected seed with $k_{\text{density}}$ in place of $k_{\min}$, not a
separate radius computed some other way. Mask every vector within that radius
as visited. Continue with the next highest-density unvisited vector until
none remain.

Tying the coverage radius to the same local-scale notion growth itself uses,
rather than a single dataset-wide radius as an earlier draft of this
algorithm did, is deliberate: a fixed global radius can suppress seed
placement across a region much larger than what that seed's own ensemble
will ever actually grow into, leaving points "covered" by seed-exclusion but
never reached by any grown ensemble -- see `misc/STC-experiments/README.md`
for a concrete demonstration of the resulting coverage gaps.

## Ensemble identification

In parallel grow ensembles around each seed vector. Use OMP paradigms if and
only if `do concurrent` parallelism is not supported due to external library
calls (gfortran).

The SKG is called `ensemble_identification` and is an iteration wrapper for the
below steps.

Consider each seed an ensemble of size one. Growing is done iteratively until
stop conditions are reached. Each iteration has an ensemble $\mathcal{E}_t$, an
2D array of $o$ observables $\mathcal{O}_{t-o+1}, \mathcal{O}_{t-o+2}, \dots,
\mathcal{O}_{t}$.

### First growth step

A seed is a single point: its `observable` would be the singular value
decomposition of a one-column centered matrix, which is degenerate -- no
meaningful $U$, $d$, or $G$ exists to compare against. `ensemble_identification`
therefore performs its first `grow_ensemble` + `observable` call
unconditionally, without invoking `accept_ensemble`; there is nothing yet to
compare the result against. `accept_ensemble` is invoked starting from the
second growth iteration onward, comparing iteration 2's observable against
iteration 1's.

This bootstrap rule lives here, in the iteration wrapper, and not in `seeds`,
`grow_ensemble`, `observable`, or `accept_ensemble` themselves. `seeds`'
purpose is identifying starting points, nothing else -- and it runs *before*
`calc_ensemble_growth_radius`, so it has no growth radius available to grow
by even if it wanted to. Keeping `grow_ensemble`, `observable`, and
`accept_ensemble` fully general and iteration-unaware keeps each
independently testable in isolation. `ensemble_identification` is the one
place that already knows which iteration it is, without introducing any new
state elsewhere.

### Stop Conditions

Growth for a single seed's ensemble stops under exactly one of the following
four conditions. `ensemble_identification` reports which one arose via an
integer result-code output, `stop_reason` (see "Output" below) -- a scheme
distinct from, and orthogonal to, this codebase's `ierr` error-code
convention: all four conditions below are valid, non-error algorithmic
outcomes, so `ierr` stays 0/success in every case.

1. **Maximum ensemble size reached.** The ensemble has grown to at least
   $\boldsymbol{f_{\text{max}}}\cdot N$ members, where $N$ is the total
   number of input vectors and $f_{\text{max}}$ is an optional input
   argument defaulting to 0.95. Growth reaching this size is treated as
   evidence against the hypothesis of intrinsic local structure (the
   ensemble has essentially absorbed the whole dataset), so **no ensemble is
   returned** for this seed -- only the result-code flagging this outcome.
2. **Rejected after being stably accepted.** The ensemble has been accepted
   at least $\boldsymbol{a}$ times (an optional input argument, default 2)
   and the current growth step is rejected by `accept_ensemble`. The last
   *accepted* ensemble and its observable history are returned, together
   with the result-code for this condition.
3. **Rejected immediately.** The first time `accept_ensemble` is actually
   invoked (comparing iteration 2 against iteration 1 -- see "First growth
   step" above) returns `.false.`, before the ensemble ever reached $a$
   accepted iterations. The bootstrap ensemble (iteration 1, always accepted
   unconditionally) is returned as the last accepted result, together with a
   result-code distinct from condition 2's, so callers can tell a
   well-supported stop from one resting on a single, unconfirmed growth
   step.
4. **No further growth possible.** `grow_ensemble` returns the same
   membership it was given (no candidate vectors lie within the growth
   radius of any current member) -- a natural fixed point. The last accepted
   ensemble and its observable history are returned, together with the
   result-code for this condition. This is the expected, "clean" way for
   growth to end.

Conditions 2-4 all return the last *accepted* ensemble; only condition 1
returns no ensemble at all. Evaluating condition 2 requires
`ensemble_identification` to track a running count of how many growth
iterations have been accepted in total for the current seed -- a separate,
cumulative counter, independent of the trailing $o$-window history below,
since that window can be shorter than the true total once growth runs longer
than $o$ iterations.

### Output

`ensemble_identification` returns, per seed:

* `final_ensemble_mask(n_vectors)`, logical -- the accepted ensemble's
  membership, using the same logical-mask convention as `seeds`. All
  `.false.` when `stop_reason` signals condition 1 (maximum ensemble size
  reached).
* `stop_reason`, integer -- the result-code for whichever of the four "Stop
  Conditions" above applied.
* The trailing observable history, the last (up to) $o$ iterations. Each
  $\mathcal{O}_t$ is not a flat vector -- $U$ is a $D\times d_t$ matrix whose
  width $d_t$ can itself change between iterations, which is exactly what
  $\delta_d$ checks. Storing the **full** $D\times D$ decomposition per
  retained iteration, as already decided for `observable` above, sidesteps
  the variable-width problem: slice $U(:,1{:}d_t)$ out on demand instead of
  needing variable-shape storage. Concretely:
  * `U_history(D, D, o)`, real;
  * `S_history(D, o)`, real -- the singular values; eigenvalues,
    `normal_error`, and `tangent_scales` at any retained iteration are
    derived on demand from these plus `k_history` below, via the
    `normal_error`/`tangent_scales` SKGs, rather than stored redundantly;
  * `d_history(o)`, integer;
  * `G_history(o)`, real;
  * `mu_history(D, o)`, real;
  * `k_history(o)`, integer -- ensemble size per retained iteration, needed
    to recover eigenvalues from singular values,
    $\lambda_j = s_j^2/(k-1)$;
  * `accepted_history(o)`, logical -- whether the growth iteration retained
    in the corresponding column was accepted. Iteration 1 (the bootstrap
    step) is recorded as `.true.` by convention, even though
    `accept_ensemble` is never actually invoked for it (see "First growth
    step" above).
* `U_first(D, D)`, real, and `d_first`, integer -- the tangent+normal basis
  and intrinsic dimension of iteration 1 (the bootstrap step), retained for
  the ensemble's entire growth and never evicted, independently of the
  trailing $o$-window above. Distinct storage from `U_history`/`d_history`
  rather than widening those to $o+1$ or pinning one of their existing
  columns, so the trailing window keeps meaning exactly what it always has
  ($o$ genuinely trailing iterations) and no existing caller's
  interpretation of `U_history`/`d_history` changes. Set once, at iteration
  1, alongside `low_confidence_mask` below (same call, no extra SVD). Used
  by `accept_ensemble`'s tangent-space-drift criterion (see "SKG
  `accept_ensemble`" below) to bound cumulative, not just step-to-step,
  rotation -- see `misc/STC-experiments/README.md`'s "P5" discussion for the
  failure mode this closes: an ensemble that "walks" arbitrarily far in
  tangent-space orientation via many individually-small-enough steps, none
  of which alone would trip a step-to-step-only check.
* Optionally, `member_added_at_step(n_vectors)`, integer -- 0 (or a
  documented sentinel) for non-members, the seed's own designation for the
  seed itself, and the growth-iteration index at which each subsequent
  member joined otherwise. Answers "which members were added when" and "did
  this ensemble oscillate during growth" directly.
* Always collected, the same reasoning as `member_added_at_step` above (a
  byproduct of computation already happening unconditionally, not a
  separate cost to gate behind its own flag): `low_confidence_mask(n_vectors)`,
  logical -- the membership from this seed's iteration 1 (the unconditional
  bootstrap `grow_ensemble` +
  `observable` call, see "First growth step" above), reported *regardless*
  of what `stop_reason` the seed's full growth eventually reaches, including
  seeds discarded entirely under condition 1 (maximum ensemble size
  reached), for which `final_ensemble_mask` above is all-`.false.`.
  Iteration 1 already has a genuine, non-degenerate `observable`/SVD behind
  it (unlike the single-point seed itself), so it is real ensemble data, not
  a guess -- but it has never survived even one `accept_ensemble` comparison,
  hence "low confidence": callers should treat it as a fallback, not a
  result on par with `final_ensemble_mask`. This is deliberately reported
  unconditionally per seed rather than the kernel itself deciding which
  points are "orphaned": that decision belongs to whoever consumes
  `ensemble_masks` (e.g. LoManLe), by diffing a point against the *retained*
  ensembles -- keeping `ensemble_identification` itself free of any opinion
  about what its caller does with an uncovered point.

#### Merged output

Suppose, Ensemble Identification found $N_{\mathcal{E}}$ ensembles (clusters).

The above per ensemble outputs are then collected into a structure of arrays
representation of the output, where each of the per ensemble outputs is
collected in an array with an extra trailing dimension of size
$N_{\mathcal{E}}$, each index $i_{\mathcal{E}}$ along that dimension
representing the data for ensemble number $i_{\mathcal{E}}$.

Thus the Ensemble Identification module's main function returns the following
two, three, or four dimensional Fortran arrays: each per-seed output above
simply gains one extra trailing dimension of size $N_{\mathcal{E}}$ once
merged. `U_history(D, D, o)` is already three dimensional per seed, so its
merged form, `ensemble_U_history`, is Tensor Omics' first four dimensional
data structure: `ensemble_U_history(D, D, o, N_{\mathcal{E}})`. This is
exactly where $o$ earns its keep: $D$ appears twice and $N_{\mathcal{E}}$
once more in that product, so keeping $o$ small is what keeps this array's
memory footprint tractable.

* `logical :: ensemble_masks`
* `real :: ensemble_U_history`
* `real :: ensemble_S_history`
* `integer :: ensemble_d_history`
* `real :: ensemble_G_history`
* `real :: ensemble_mu_history`
* `integer :: ensemble_k_history`
* `logical :: ensemble_accepted_history`
* `real :: ensemble_U_first` -- the merged form of `U_first` above,
  `ensemble_U_first(D, D, N_{\mathcal{E}})`.
* `integer :: ensemble_d_first` -- the merged form of `d_first` above,
  `ensemble_d_first(N_{\mathcal{E}})`.
* `integer :: ensemble_stop_reason` -- the result-code for whichever of the
  four "Stop Conditions" above applied to that ensemble, see "Output" above.
* Always collected (a byproduct of the per-seed growth loop, not a separate
  cost worth gating -- see `member_added_at_step` above): `integer ::
  ensemble_member_added_at_step` - The integer values indicate at what
  iteration members were added. The row index in an ensemble's
  `member_added_at_step` column vector refer to the `.true.` row indices in the
  respective ensemble's `ensemble_masks` column vector.
* Always collected, same reasoning: `logical :: ensemble_low_confidence_masks`
  -- the merged form of `low_confidence_mask` above, one column per seed
  regardless of `ensemble_stop_reason`.

### Local Radius Identification

For growing ensembles we need a locally adapted radius with which at each step
candidates for the addition are identified. Those candidates are vectors within
the local ensemble specific radius.

The SKG should be called `calc_ensemble_growth_radius` or similar.

Use a fixed-count kNN pool and compute a percentile of the distances among a
seed's own $k_{\min}$ nearest neighbors, `radius_percentile` (0 to 100,
default 50.0 -- the median, this SKG's original, non-parameterized behavior).
Make $k_{\min}$ an optional argument with default value 30.

For each seed we store the growth radius in a 1D real array called
`ensemble_growth_radii`.

`seeds` reuses this exact SKG for its own coverage radius (see "Seeding"
above), called with $k_{\text{density}}$ in place of $k_{\min}$ -- not a
second, separately-implemented percentile-of-k-NN-distances computation.
`seeds` exposes its own call's `radius_percentile` as `exclusion_radius_percentile`,
independently of whatever percentile any other caller of this SKG uses for the
actual growth-phase radius: a fixed, dataset-wide-in-spirit median exclusion
radius was observed empirically to suppress seed placement across curvature
extrema (peaks, troughs, kinks on a wavy manifold) that a seed's own growth
never actually reaches -- see `misc/STC-experiments/README.md`. Shrinking
`exclusion_radius_percentile` below 50 trades more, smaller ensembles for less
over-eager suppression; it does not change the growth-phase radius itself.

### Tangent Space Variant

This section describes the implementation used for the tangent space variant of
STC. 

Each iteration of `ensemble_identification` runs:

* `grow_ensemble`
* `observable`
* `accept`

#### SKG `grow_ensemble`

Identify the candidate points within the input ensemble's growth-radius $r_{\mathcal{E}}$:
$$
\mathcal{C}_{\mathcal{E}_{t+1}} = \{ x_k \in \text{data\_vectors} ~|~ ||x_k - x_i|| \leq r_{\mathcal{E}} ~\exists x_i \in \mathcal{E} \}
$$

#### SKG `observable`

Calculate the tuple $\mathcal{O}_{\mathcal{E}_{t+1}} = ( U_{\mathcal{E}_{t+1}},
d_{\mathcal{E}_{t+1}}, G_{\mathcal{E}_{t+1}}, \boldsymbol\mu_{\mathcal{E}_{t+1}},
\text{normal\_error}_{\mathcal{E}_{t+1}},
\text{tangent\_scales}_{\mathcal{E}_{t+1}})$ for the input ensemble
$\mathcal{E}_{t+1}$, where

**(1)** $U_{\mathcal{E}_{t+1}}$, the tangent basis, is obtained from the
**singular value decomposition** of the Ensemble's centered member vectors
-- not from an eigendecomposition of their covariance (see "Numerical
Linear Algebra" above for why). To obtain it, compute the neighborhood
center

$$
\boldsymbol{\mu}_{\mathcal{E}_{t+1}}
=
\frac{1}{|\mathcal{E}_{t+1}|}
\sum_{x_j \in \mathcal{E}_{t+1}}\mathbf{x}_j,
$$

construct the centered matrix

$$
Y_{\mathcal{E}_{t+1}}=
[
\mathbf{x}_{1}-\boldsymbol{\mu}_{\mathcal{E}_{t+1}},
\ldots,
\mathbf{x}_{\mathcal{N}}-\boldsymbol{\mu}_{\mathcal{E}_{t+1}}
],
\text{for}~x_1, \ldots, x_{\mathcal{N}} \in \mathcal{E}_{t+1},
$$

and its economy-mode singular value decomposition, using LAPACK `dgesdd`
(`JOBZ='S'`), over the **full** set of $D$ singular values/vectors, not
only the top-$d$ tangent ones:

$$
Y_{\mathcal{E}_{t+1}} = U_{\mathcal{E}_{t+1}} S_{\mathcal{E}_{t+1}} V_{\mathcal{E}_{t+1}}^\top.
$$

The retained "leftover" pairs are what make center, normal error, and
tangent scales free below, and are also exactly what a later noise-tube
covariance $\Sigma_\perp$ (needed once ensembles are reconciled) would be
built from.

The corresponding eigenvalues are recovered from the singular values as

$$
\lambda^{\mathcal{E}_{t+1}}_{j} = \frac{\left(s^{\mathcal{E}_{t+1}}_j\right)^2}{k_{\mathcal{E}_{t+1}}-1},
\qquad
\lambda^{\mathcal{E}_{t+1}}_{1}\geq
\lambda^{\mathcal{E}_{t+1}}_{2}\geq\cdots\geq
\lambda^{\mathcal{E}_{t+1}}_{D}\geq0.
$$

LAPACK SVD routines already return singular values in *descending* order
(unlike `dsyev`'s ascending eigenvalues), so this matches the document's own
convention directly: the tangent basis and its eigenvalues are the *first*
$d$ columns/entries of `dgesdd`'s output, and the normal-space basis and
eigenvalues are the remaining $D-d$.

**(1b)** $\boldsymbol\mu_{\mathcal{E}_{t+1}}$, the ensemble center, is simply
the mean vector already computed above to form
$Y_{\mathcal{E}_{t+1}}$ -- no further work needed.

**(1c)** $\text{normal\_error}_{\mathcal{E}_{t+1}}$ and
$\text{tangent\_scales}_{\mathcal{E}_{t+1}}$ are obtained from the retained
eigenvalues via the `normal_error` and `tangent_scales` SKGs defined below.

**(2)** The Spectral Gap is calculated as follows.

For every candidate intrinsic dimension
$$
r = 1, \ldots, D-1
$$
compute its Spectral Gap
$$
G(r) = \frac{\lambda^{\mathcal{E}_{t+1}}_{r}}{\lambda^{\mathcal{E}_{t+1}}_{r+1} + \epsilon},
$$
where $\epsilon$ is a minimal positive real number to avoid division by zero.

**(3)** and the estimated intrinsic dimension
$$
d_{\mathcal{E}_{t+1}} = \text{arg}~ \operatorname*{max}_{r} G(r)
$$

and set the current ensemble's spectral gap to
$$
G_{\mathcal{E}_{t+1}} = G(d_{\mathcal{E}_{t+1}})
$$

Note that an optional argument $o$ with default value of e.g. 10 indicates how
many observables of the last $o$ iterations are stored. Thus, each ensemble
$\mathcal{E}_{t+1}$ has a two dimensional real array of observables, one per
column: $$ \Omega_{\mathcal{E}_{t+1}} = [ \mathcal{O}_{\mathcal{E}_{t-o+2}},
\mathcal{O}_{\mathcal{E}_{t-o+3}}, \ldots, \mathcal{O}_{\mathcal{E}_{t+1}} ] $$

In addition to this trailing window, $\mathcal{O}_{\mathcal{E}_1}$ (iteration
1, the bootstrap observable) is *always* retained as `U_first`/`d_first`
(see "Output" above), independently of $o$'s FIFO eviction. Together,
$\Omega_{\mathcal{E}_{t+1}}$'s trailing window and $\mathcal{O}_{\mathcal{E}_1}$
form the reference set `accept_ensemble`'s tangent-space-drift criterion
compares the new candidate against -- see "SKG `accept_ensemble`" below. For
$t+1 \leq o$, $\mathcal{O}_{\mathcal{E}_1}$ already sits in the trailing
window too ($U_{\text{history}}(:,:,1)$); comparing the candidate against it
twice (once via `U_first`, once via the window) is harmless -- it does not
change a $\max$ -- and simpler than special-casing the overlap away.

#### SKG `normal_error`

Given intrinsic dimension $d_{\mathcal{E}}$ and the eigenvalues
$\lambda^{\mathcal{E}}_{1}\geq\cdots\geq\lambda^{\mathcal{E}}_{D}$ from
`observable`, the mean squared residual of the ensemble's members off its
$d_{\mathcal{E}}$-dimensional tangent subspace is

$$
\text{normal\_error}_{\mathcal{E}} = \sum_{j=d_{\mathcal{E}}+1}^{D} \lambda^{\mathcal{E}}_{j}.
$$

No pass over the ensemble's member vectors is required; the sum is already
implied by the singular value decomposition computed in `observable`.

#### SKG `tangent_scales`

Given the same $d_{\mathcal{E}}$ and eigenvalues, the extent along each
tangent direction is

$$
\text{tangent\_scales}_{\mathcal{E}}(j) = \sqrt{\lambda^{\mathcal{E}}_{j}}, \qquad j = 1,\ldots,d_{\mathcal{E}}.
$$

#### SKG `accept_ensemble`

This routine returns a Boolean indicating whether the grown ensemble
$\mathcal{E}_{t+1}$ is accepted or not.

A grown ensemble ($\mathcal{E}_{t+1}$) is assessed against a **reference
set** $\mathcal{R}$ of previously accepted observables, not merely against
$\mathcal{O}_{\mathcal{E}_t}$ alone:
$$
\mathcal{R} = \{\mathcal{O}_{\mathcal{E}_1}\} \cup
\{\mathcal{O}_{\mathcal{E}_{t-o+2}}, \ldots, \mathcal{O}_{\mathcal{E}_t}\}
$$
i.e. the bootstrap observable (`U_first`/`d_first`, always retained, see
"Output" above) plus the trailing $o$-window. Comparing against this whole
set, not just $\mathcal{O}_{\mathcal{E}_t}$, is what closes the "no
cumulative-rotation budget" gap: a step-to-step-only check cannot see an
ensemble that walks arbitrarily far in tangent-space orientation via many
individually-small steps, since each individual step still passes; anchoring
every candidate against the ensemble's own original state (and its recent
history) bounds the *cumulative* drift instead.

$\mathcal{E}_{t+1}$ is accepted if all four of the following hold:

**(1) Tangent-space drift.** For every reference $\mathcal{O}_{\mathcal{E}_r}
\in \mathcal{R}$ with $d_r = d_{\mathcal{E}_{t+1}}$ (references with a
different $d_r$ are skipped for this criterion -- there is no shared tangent
dimension to compare angles over, judged instead by criterion (2)), compute
the principal angles $\theta_1, \ldots, \theta_{d_r}$ between $U_r$ and
$U_{\mathcal{E}_{t+1}}$ via `dgesvd` on $M = U_r^\top U_{\mathcal{E}_{t+1}}$,
whose singular values are $\cos\theta_k$ directly (see "Numerical Linear
Algebra" above), then the **chordal distance**
$$
\mathcal{C}(U_r, U_{\mathcal{E}_{t+1}}) = \sqrt{\sum_{k=1}^{d_r} \sin^2\theta_k}
$$
(the standard Grassmannian chordal/projection-Frobenius distance -- Edelman,
Arias & Smith 1998; related to the Davis-Kahan $\sin\Theta$ theorem). Accept
requires
$$
\max_{r\,:\,d_r = d_{\mathcal{E}_{t+1}}} \mathcal{C}(U_r, U_{\mathcal{E}_{t+1}})
\leq \boldsymbol{\text{chordal\_dist\_max\_as\_prcnt\_of\_range}} \cdot
\sqrt{d_{\mathcal{E}_{t+1}}}
$$
(vacuously satisfied if no reference shares $d_{\mathcal{E}_{t+1}}$).
`chordal_dist_max_as_prcnt_of_range` (optional argument, 0 to 1, default 0.5)
replaces the earlier `alpha_max` (a raw angle in radians): $\mathcal{C}$'s
range is $[0, \sqrt{d}]$, so this parameter states the tolerated fraction of
that range directly, without an intermediate angle-to-sine conversion, and
without needing a separate $\sqrt{d}$-dependent knob per intrinsic dimension
-- the $\sqrt{d}$ scaling is handled internally. It reduces to a familiar
special case at $d=1$: $\mathcal{C} = |\sin\theta|$, so
`chordal_dist_max_as_prcnt_of_range` $= \sin(\alpha_{\max})$ recovers exactly
the old single-angle boundary (default 0.5 $\approx \sin(30^\circ)$, a literal
carry-over of the old default under this substitution -- provisional, like
every other default in this family, pending empirical re-tuning once this
criterion runs on real data).

Only the candidate is compared against each reference -- **not** every pair
within $\mathcal{R}$ itself ($\binom{|\mathcal{R}|}{2}$ pairs). This is
sufficient, not a shortcut: every reference in $\mathcal{R}$ already passed
this same criterion against every *other* reference present in $\mathcal{R}$
at the growth step it was itself accepted (an inductive invariant -- FIFO
eviction only ever removes references from $\mathcal{R}$, which cannot break
a "already mutually within tolerance" property among the ones that remain).
So checking the one new candidate against the (at most) $o+1$ existing
references costs $O(o)$ small SVDs per growth step, not $O(o^2)$.

**(2) Intrinsic-dimension drift.** Both
$$
d_{\text{to\_first}} = |d_{\mathcal{E}_{t+1}} - d_{\mathcal{E}_1}| \leq \boldsymbol{d_{\text{max}}}
\quad\text{and}\quad
d_{\text{to\_last}} = |d_{\mathcal{E}_{t+1}} - d_{\mathcal{E}_t}| \leq \boldsymbol{d_{\text{max}}}
$$
must hold -- i.e. $\max(d_{\text{to\_first}}, d_{\text{to\_last}}) \leq
d_{\text{max}}$. Compared only against the first and the immediately
preceding iteration (not the full $o$-window, unlike criterion (1)) --
cheap, and closes the same first-vs-cumulative gap for $d$ that criterion (1)
closes for orientation, without needing $d$ history beyond what
`d_first`/the previous iteration's $d$ already provide.

**(3) Spectral-gap drift.**
$$
\left|\log \frac{G_{\mathcal{E}_{t+1}}}{G_{\mathcal{E}_t}}\right| \leq \boldsymbol{G_{\text{max}}}
$$
Consecutive only ($t$ vs. $t+1$), unchanged from before. $G$ is a magnitude
(the ratio of two adjacent eigenvalues), not an orientation -- the "walk
arbitrarily far while every single step looks fine" failure mode criterion
(1) fixes is specific to direction, not scale, so this criterion is left as
is for now.

**(4) Residual (RMSE) drift.** Let
$\text{RMSE}_{\mathcal{E}} = \sqrt{\text{normal\_error}_{\mathcal{E}} + \epsilon}$ (the
root of the `normal_error` SKG's output, itself already free -- see
"SKG `normal_error`" above; no new SVD or storage, `normal_error` at any
retained iteration is already derivable on demand from `S_history`/
`k_history`), where $\epsilon$ is the same minimal-positive-real convention
`observable`'s own spectral-gap formula already uses (`epsilon(1.0_real64)`
in Fortran) -- unlike $G$, which is a *ratio* already protected by its own
$+\epsilon$ denominator, `normal_error` is a raw sum of eigenvalues and can
be genuinely, exactly zero (a perfectly flat/collinear ensemble, no noise at
all in the normal directions); without this guard,
$\log(\text{RMSE}_{t+1}/\text{RMSE}_t)$ is $\log(0/0)$, undefined, the first
time growth is ever perfectly noise-free. Accept requires
$$
\left|\log \frac{\text{RMSE}_{\mathcal{E}_{t+1}}}{\text{RMSE}_{\mathcal{E}_t}}\right| \leq \boldsymbol{\text{RMSE\_change\_max}}
$$
(optional argument, default $|\log(1.5)| \approx 0.405$), consecutive only,
same structural form as criterion (3). $G$ (a ratio between adjacent
eigenvalues) and RMSE (the absolute magnitude of *all* residual eigenvalues
combined) are genuinely different signals, not a restatement of the same
information under a different name: an ensemble can hold $G$ essentially
constant across growth (its dimensionality estimate stays well-supported)
while the absolute off-tangent-subspace spread still creeps upward --
criterion (3) alone cannot see that, since it only ever compares $\lambda_d$
against its immediate neighbor $\lambda_{d+1}$, never the total normal-space
spread. This criterion is what catches that.

## Ensemble Reconciliation

This step implements a final polishing, quality-filtering, and conflict
resolution step of the identified ensembles. In the tangent-space-variant
intersections between ensembles are identified and reported. Note that by
default ensembles are not merged based on intersections, because the
intersections are needed by our Manifold Learning algorithm to stith together
local linear manifolds.

This module first identifies ensembles that intersect. It then offers three
different modes of how to process intersections:
* (1) just report them
* (2) merge with a minimal set Overlap Coefficient (SOC; default 0.9)
* (3) merge all, i.e. if the intersection has size >= 1, merge.

The module does not alter the result produced by "Ensemble Identification". It
creates a sparse representation of an adjacency matrix of edges between
ensembles indicating which have intersections (mode 1), or have been merged
(modes 2 and 3).

Detecting an intersection at all already requires computing
$|\mathcal{E}_i \cap \mathcal{E}_j|$, needed by every mode, and each
ensemble's own size $|\mathcal{E}_i|$ is cheap to precompute once per
ensemble ($O(N_{\mathcal{E}})$ total, not once per pair). The **Overlap
Coefficient** (Szymkiewicz-Simpson coefficient) itself,
$$
\text{OC}(\mathcal{E}_i, \mathcal{E}_j) = \frac{|\mathcal{E}_i \cap \mathcal{E}_j|}{\min(|\mathcal{E}_i|, |\mathcal{E}_j|)},
$$
is then a single extra $O(1)$ step per pair -- even cheaper than the Jaccard
Similarity Index this replaces, since the denominator is just the smaller of
the two already-precomputed sizes, not a union requiring both sizes and the
intersection combined. Overlap Coefficient (unlike Jaccard) is exactly the
statistic for "is one ensemble, or most of it, contained within the other,"
which is the actual question a true or partial subset relationship poses:
two ensembles where the smaller is fully contained in the larger score
$\text{OC}=1$ regardless of how much larger the other one is, whereas their
JSI would be driven arbitrarily low by the size difference alone -- exactly
the containment relation JSI fails to capture. Reporting it is made
available regardless of mode, via an optional user flag, rather than tied to
mode 2 specifically. Modes 1 and 3 do not need it for their own decision, so
**the implementation must guard the computation behind that flag** (an `if
(report_overlap_coefficient) ...`, not an unconditional computation): the
added cost is small, but it is not exactly zero, and nothing here should be
computed that nobody asked for.

The output assumes a maximum number of ensembles that can belong to one
merged group, to enable allocation of the sparse adjacency matrix
representation. By default this max number is $\min(1024, N_{\mathcal{E}})$.
Thus the output can be represented as the two dimensional
`integer(int32) :: super_ensembles` array in which each column contains one
"super-ensemble", i.e. a group of ensembles that have intersections (mode 1)
or have been merged (modes 2 and 3); its rows hold that group's member
ensembles, padded with `0` (an invalid ensemble number, since ensembles are
1-indexed) beyond the group's actual size. The content of each cell are the
column indices of module Ensemble Identification's `ensemble_masks` output,
thus identifying the ensembles directly via their indices. If Overlap
Coefficient reporting was requested, an optional two dimensional real
`super_ensembles_overlap_coefficient` array is returned additionally. It has
the same number of columns as `super_ensembles` and one row less. Each
column in `super_ensembles_overlap_coefficient` represents the
super-ensemble in `super_ensembles` directly. Within each column $l$, each
cell index $c_i$ indicates the Overlap Coefficient between the ensemble
number stored in `super_ensembles(c_i, l)` and `super_ensembles(c_i + 1, l)`.

## Estimate parameters from data

The crucial parameters are `k_min`, `k_density`, `density_quantile`,
`chordal_dist_max_as_prcnt_of_range`, `G_max`, `d_max` (`RMSE_change_max`,
`accept_ensemble`'s fourth criterion, is not estimated by this SKG -- see
"SKG `estimate_stc_parameters`" below). Grid-searching them, or the bootstrap-resampling
scheme considered and rejected during this section's design (per-seed
resampled SVDs to measure noise-only wobble in angle/gap), both cost far more
than the rest of the pipeline combined. This SKG estimates all six from one
cheap, coarse pass instead: grow a handful of anchors into rough local
neighborhoods using the same primitives (`density_labels`, `observable`) the
real pipeline already has, then read the parameters off simple summary
statistics of that pass. It is a **separate, optional** step -- it neither
runs automatically as part of `seeds` or `ensemble_identification`, nor
requires that either has already run; a caller may run it, inspect the
proposed values, adjust them, and only then invoke the real pipeline, or skip
it entirely and supply `k_min`/`chordal_dist_max_as_prcnt_of_range`/etc. by
hand as always.

### SKG `sample_estimator_anchors`

Receives `density_labels`' own output (see "Seeding" above -- this SKG does
not recompute density itself) and an optional $n_{\text{anchors}}$ (default
5). Sort the vectors by density label and pick $n_{\text{anchors}}$ of them,
by nearest-rank (not interpolated -- these must be genuine point indices),
at the percentiles $\frac{100}{n_{\text{anchors}}}, \frac{200}{n_{\text{anchors}}},
\ldots, 100$ of that sorted order -- e.g. $n_{\text{anchors}}=5$ gives
20/40/60/80/100%ile, generalizing the original "one point at each 20%ile"
choice to a tunable anchor count. We call the resulting points the estimator
anchors (EAs).

### SKG `grow_estimator_anchor_clouds`

Each EA starts as its own single-point cloud. Growth proceeds in rounds:
every EA proposes the single closest not-yet-claimed point to *any* member of
its own current cloud (not just to the anchor itself -- the cloud grows as a
region, the same way `grow_ensemble` grows an ensemble); among all pending
proposals, the closest one is claimed by its proposing EA, ties broken by
which EA is nearer. This is a multi-source region-growing process --
equivalent to a multi-source Prim's/Dijkstra expansion, one shared frontier
across all EAs, distance-ordered claims -- not $n_{\text{anchors}}$
independent kNN queries: a point equidistant-ish between two EAs must go to
whichever is *currently* closer, which can only be answered by comparing
proposals across EAs at claim time, not by asking each EA in isolation.

Implemented as a brute-force rescan, not via the k-d tree: every round, for
every EA, scan all not-yet-claimed points against all of that EA's current
cloud members directly on `vectors`, no tree traversal. `f42_kd_tree`'s own
primitives answer "nearest point(s) to one query point," not "nearest
unclaimed point to any member of a growing region," so using it here would
mean re-deriving that query shape from scratch -- more implementation
complexity for a saving that does not matter at this SKG's intended scale.
With $n_{\text{anchors}}$ small (default 5) and total growth capped small by
`seed_max_set_size` (default 5% of $N$) precisely so that clouds stay local
(see below), the whole rescan is $O(\text{rounds} \times n_{\text{anchors}}
\times N \times \text{cloud size})$, all of them small factors by
construction -- trivial in absolute terms for the datasets this SKG targets.
This is a deliberate simplicity-over-asymptotic-elegance trade, consistent
with this whole SKG's "minimal complexity" mandate; it is not the right
choice if `seed_max_set_size` is ever pushed towards its unbounded (100%)
extreme on a large $N$, where the original per-EA-partitions-everything
draft's own cost concern (see "Growth stops..." below) returns.

Growth stops when either no unclaimed point remains reachable, or the total
number of claimed points across all EAs reaches an optional
`seed_max_set_size` (0 to 100, default 5.0) percent of $N$. The default is
deliberately small: growing until the *entire* dataset is partitioned among
$n_{\text{anchors}}$ anchors (the original, unbounded draft of this
algorithm) gives each EA a cloud averaging $N/n_{\text{anchors}}$ points --
for $n_{\text{anchors}}=5$, a fifth of the whole dataset, easily spanning a
kink or bend and no longer "local" in any sense `observable`'s SVD could
meaningfully summarize. Bounding total growth to a small percentage of $N$
keeps every EA's cloud a genuinely local neighborhood, at the cost of not
necessarily reaching every point -- acceptable here, since covering the
dataset is `seeds`' job, not this SKG's.

### SKG `estimate_stc_parameters`

The orchestrator. Runs `density_labels`, then `sample_estimator_anchors`,
then `grow_estimator_anchor_clouds`, then `observable` once per EA on its
final cloud (reused directly -- this SKG does not reimplement the SVD). From
the $n_{\text{anchors}}$ per-EA results $(k_i, \bar{d}_i^{\text{median}},
d_i, G_i)_{i=1}^{n_{\text{anchors}}}$, where $k_i$ is EA $i$'s final cloud
size and $\bar{d}_i^{\text{median}}$ the median distance from EA $i$'s anchor
to its own cloud members:

* $k_{\min} \leftarrow \text{median}_i(k_i)$
* $k_{\text{density}} \leftarrow k_{\min}$ -- reused, not separately
  estimated, the same relationship `density_labels`' own default already has
  to `calc_ensemble_growth_radius`'s $k_{\min}$ (see "Seeding" above).
* $\text{density\_quantile} \leftarrow \text{median}_i(\bar{d}_i^{\text{median}})$
  -- a literal radius (data units), used directly wherever a radius is
  needed, not converted into a 0-100 percentile: `seeds`/
  `calc_ensemble_growth_radius`'s own `exclusion_radius_percentile`/
  `radius_percentile` already exist to *choose* a percentile of a kNN pool;
  this SKG instead hands back an already-measured absolute scale, which a
  caller may use directly as a growth or exclusion radius.
* For every pair $i<j$ among the EAs with $d_{ij} = \min(d_i, d_j) > 0$,
  compute the principal angles $\theta_1,\ldots,\theta_{d_{ij}}$ between
  $U_i$ and $U_j$ over their shared rank $d_{ij}$ -- the same
  `dgesvd`-on-$U_i^\top U_j$ machinery `accept_ensemble` uses (see
  "Numerical Linear Algebra" above), but *not* gated on $d_i=d_j$ the way
  `accept_ensemble` gates its own criterion (1) (there, mismatched-$d$
  references are vacuously skipped, judged instead by criterion (2)); here,
  with only $\binom{n_{\text{anchors}}}{2}$ pairs total (10 at the default
  $n_{\text{anchors}}=5$), skipping mismatched pairs risks too few (or zero)
  samples to estimate a quartile from at all. Comparing over the shared rank
  always produces a value, at the cost of diverging slightly from
  `accept_ensemble`'s own convention -- a deliberate trade for this SKG's
  specific purpose, not an oversight. Pairs with $d_{ij}=0$ (at least one EA
  has no meaningful tangent direction at all) are excluded from the sample
  entirely, the same "nothing to compare" convention as `accept_ensemble`'s
  own $d_r > 0$ guard.
* Chordal distance per pair, normalized by its own $\sqrt{d_{ij}}$ so pairs of
  differing rank are directly comparable -- exactly
  `chordal_dist_max_as_prcnt_of_range`'s own definition (see "SKG
  `accept_ensemble`" above), applied here to EA pairs instead of
  candidate-vs-reference pairs:
  $$
  c_{ij} = \mathcal{C}(U_i, U_j) / \sqrt{d_{ij}} = \sqrt{\textstyle\sum_{k=1}^{d_{ij}} \sin^2\theta_k} \big/ \sqrt{d_{ij}}
  $$
* `chordal_dist_max_as_prcnt_of_range` $\leftarrow Q_{p}\big(\{c_{ij}\}_{i<j}\big)$
* $G_{\max} \leftarrow Q_{p}\big(\{|\log(G_i/G_j)|\}_{i<j}\big)$
* $d_{\max} \leftarrow Q_{p}\big(\{|d_i - d_j|\}_{i<j}\big)$

where $Q_p$ is the $p$-th percentile (linear interpolation, via
`calc_percentile_helper`) and $p$, `first_quartile_percentile`, is itself an
optional argument (default 25.0, the first quartile) -- not hardcoded,
following the same "expose the heuristic, do not disguise it" precedent as
`bandwidth_percentile`/`exclusion_radius_percentile` above.

This whole procedure is a heuristic starting point, not a converged answer:
$n_{\text{anchors}}=5$ (10 pairs) is a small sample, anchors are chosen by
density quantile alone (no spatial-spread guarantee -- on a dataset with
strong density heterogeneity, all anchors could land in similar regions), and
`grow_estimator_anchor_clouds` has no curvature-awareness of its own, only
the `seed_max_set_size` size cap. Treat its output the way any of this
family's other tunable defaults are treated: a reasonable value to start
from and refine, not a guarantee.

# Command line interface (CLI) in C

Implemented in `C-layer/` (`csv_table.h`/`csv_table.c`, `stc_cli.c`; see `C-layer/README.md`
for usage) -- not `src/f42/`, since it is C, not Fortran, and not part of the generator's own
kernel/wrapper split.

Use GNU argp to create a C command line interface that leverages `libcsv` to
parse an input table with callbacks to create a 2D real Fortran array of input
data vectors, then call the Fortran `_c` API to run our Shape Truthful
Clustering from the command line interface.

## Real 2D array from a CSV with libcsv

Write a F42 module that uses C's libcsv to parse a CSV and creates a transposed
2D array from the input data. In this first version assume all values are real
numbers. Also note our `src/f42/f42_safeguard.F90` ensuring that C types match
Fortran types.

The module expects the following arguments:
* file, i.e. the complete file path to the input CSV
* separator, i.e. the char used to separate columns
* header, i.e. a flag indicating whether a header line is present

The module implements the following parsing approach:
* Query libcsv for the number of columns, e.g. by parsing the first line, then
  instantiate a real 2D array with nrows equal to the number of columns. If
  possible use some file information to get the number of lines and instantiate
  the the read 2D array with nrows and nlines, i.e. the input CSV will be
  stored transposed in Fortran, because Fortan is column major. 
  * If not possible to obtain the number of lines, makes this a required input
    argument.
* Then use libcsv's callback mechanism to populate the above output 2D array.
  Do all of this still in C.
* Once parsing has passed provide a pointer and dimension information to pass
  the 2D array with zero copy to Fortran. Our codebase has many examples of how
  to do this with zero copy.

## Arguments of CLI

The command line interface will require arguments:
* file, i.e. the complete file path to the input CSV
* separator, i.e. the char used to separate columns
* header, i.e. a flag indicating whether a header line is present
* the arguments required by STC, like k_min, etc.
  * include optional arguments with default values, i.e. the user can provide
    values but is not required to.
* output-dir, required argument of a directory-path to which to write our
  output files, including the HTML/D3 plot.
