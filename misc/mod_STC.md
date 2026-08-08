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
  $M=U_{\mathcal{E}_t}^\top U_{\mathcal{E}_{t+1}}$, for principal angles in
  `accept_ensemble`): use `dgesvd`, not `dgesdd`. At $d$'s typical size
  (the intrinsic dimension, usually 1-5), divide-and-conquer's speed
  advantage does not materialize and its extra `IWORK` bookkeeping is pure
  overhead for no benefit; `dgesvd`'s simpler single-workspace-array
  bookkeeping wins. Its singular values are $\cos\theta_j$ directly -- no
  squaring, no `sqrt` needed back out.

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
* Optionally, `member_added_at_step(n_vectors)`, integer -- 0 (or a
  documented sentinel) for non-members, the seed's own designation for the
  seed itself, and the growth-iteration index at which each subsequent
  member joined otherwise. Answers "which members were added when" and "did
  this ensemble oscillate during growth" directly.

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
* `integer :: ensemble_stop_reason` -- the result-code for whichever of the
  four "Stop Conditions" above applied to that ensemble, see "Output" above.
* _optional_ (user flag decides whether this is collected and returned)
  `integer :: ensemble_member_added_at_step` - The integer values indicate at
  what iteration members were added. The row index in an ensemble's
  `member_added_at_step` column vector refer to the `.true.` row indices in the
  respective ensemble's `ensemble_masks` column vector.

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

A grown ensemble ($\mathcal{E}_{t+1}$) is assessed based on its recent
trajectory of observables $\Omega_{\mathcal{E}_{t+1}}$. 

Acceptance is calculated comparing $\mathcal{O}_{\mathcal{E}_{t+1}}$ with
$\mathcal{O}_{\mathcal{E}_{t}}$, and $\mathcal{E}_{t+1}$ is accepted if:

* principal angles $\alpha_i \leq \boldsymbol{\alpha_{\text{max}}}, i =
  1,\ldots,\min(d_{\mathcal{E}_{t+1}}, d_{\mathcal{E}_t})$ between
  $U_{\mathcal{E}_{t+1}}$ and $U_{\mathcal{E}_t}$ -- obtained via `dgesvd` on
  $M=U_{\mathcal{E}_t}^\top U_{\mathcal{E}_{t+1}}$, whose singular values are
  $\cos\alpha_i$ directly (see "Numerical Linear Algebra" above); check
  $d_{\mathcal{E}_{t+1}}=d_{\mathcal{E}_t}$ first and short-circuit to reject
  on mismatch, since that is cheaper than the SVD and the ensembles do not
  share a common tangent dimension to compare angles over in that case, and
* $| d_{t+1} - d_{t} | \leq \boldsymbol{d_{\text{max}}}$, and
* $|log \frac{G_{\mathcal{E}_{t+1}}}{G_{\mathcal{E}_{t}}} | \leq \boldsymbol{G_{\text{max}}}$

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
* (2) merge with a minimal Jaccard Similarity Index (JSI; default 0.1)
* (3) merge all, i.e. if the intersection has size >= 1, merge.

The module does not alter the result produced by "Ensemble Identification". It
creates a sparse representation of an adjacency matrix of edges between
ensembles indicating which have intersections (mode 1), or have been merged
(modes 2 and 3).

Detecting an intersection at all already requires computing
$|\mathcal{E}_i \cap \mathcal{E}_j|$, needed by every mode, and each
ensemble's own size $|\mathcal{E}_i|$ is cheap to precompute once per
ensemble ($O(N_{\mathcal{E}})$ total, not once per pair). The JSI itself,
$|\mathcal{E}_i \cap \mathcal{E}_j| / |\mathcal{E}_i \cup \mathcal{E}_j|$, is
then a single extra $O(1)$ step per pair, via
$|\mathcal{E}_i \cup \mathcal{E}_j| = |\mathcal{E}_i| + |\mathcal{E}_j| -
|\mathcal{E}_i \cap \mathcal{E}_j|$ -- so reporting it is made available
regardless of mode, via an optional user flag, rather than tied to mode 2
specifically. Modes 1 and 3 do not need the JSI for their own decision, so
**the implementation must guard the computation behind that flag** (an `if
(report_jsi) ...`, not an unconditional computation): the added cost is
small, but it is not exactly zero, and nothing here should be computed that
nobody asked for.

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
thus identifying the ensembles directly via their indices. If JSI reporting
was requested, an optional two dimensional real `super_ensembles_JSI` array
is returned additionally. It has the same number of columns as
`super_ensembles` and one row less. Each column in `super_ensembles_JSI`
represents the super-ensemble in `super_ensembles` directly. Within each
column $l$, each cell index $c_i$ indicates the JSI between the ensemble
number stored in `super_ensembles(c_i, l)` and `super_ensembles(c_i + 1, l)`.

## Estimate parameters from data

The crucial parameters are `k_min`, `k_density`, `density_quantile`,
`alpha_max`, `G_max`, `d_max`.

Once we have the density labels, we can sample one point at each 20%ile of the
sorted densities. We call those the estimator anchors (EAs). Now we iterate,
adding one point at a time to each of the EAs in parallel. Points are only
added if they are not yet covered. If more than one EA would add the same point
the one with less distance to that point takes precedence. Stop growing if
there are no points to add. Now we estimate the above parameters as the median
of the following measures which we obtain for each EA:
* k_min <- the number of neighbors in the EAs grown cloud
* k_density <- use k_min
* density_quantile <- median of distances to an EA's neighbors
* execute SVD on each EA's neighbor set, then:
  * alpha_max as the first quartile of all pairwise principal angles between
    all EA's neigborhoods
  * G_max the first quartile of all pairwise EA's G
  * d_max the first quartile of all pairwise EA's d_max

Importantly, this `first quartile` must be an optional parameter, too.
