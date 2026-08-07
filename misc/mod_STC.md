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

Follow `misc/Fortran_Coding_Guides.pdf` (the F42 standard) strictly, in
particular section 10.2 ("Scientific Kernel (SK)") for SK purity/allocation
rules and section 10.3.2 ("Naming Convention", p. 22) for the required suffix
scheme:

* no suffix: pure SK routine (the numerical core; no I/O, no
  `allocate`/`deallocate`, deterministic);
* `_alloc`: allocation helper that computes required workspace sizes,
  allocates them, and calls the validated entry point below;
* `_C` / `_R`: C-/R-facing API wrappers;
* `_helper`: general-purpose orchestration helper not bound to a specific API.

The guide's suffix list has no separate "validate" suffix; input validation
(via `src/macros.h`'s precompiler macros) lives in the *unsuffixed* routine,
which validates its arguments and then calls the pure `_helper` core. This
three-tier `_alloc` -> unsuffixed validated entry -> `_helper` pattern is
already used on this branch in `src/tox_shatter_cluster_data.F90`
(`calculate_density_radius_alloc` -> `calculate_density_radius` ->
`calculate_density_radius_helper`); follow that same shape here rather than
introducing a new naming scheme.

Related routines following this pattern are called SK-groups (SKGs).

## Seeding

Seeding requires two dedicated SKGs. One to find the radius used to measure
local density centered on a vector $\vec{x}_i$, the other to identify seeds.

The SKG `calculate_density_radius` to identify the radius receives a percentile
as optional input (default 15%). The algorithm finds the mean vector
$x = \frac{1}{N} \cdot \sum_{i=1}^{N} x_i$ with 
$x_i \in \text{data\_vectors}$
Next, all distances are computed to this mean vector. Finally, the argument
percentile distance is returned as the radius to be used in density label
calculation.

A `seeds` SKG is used with input:

* data_vectors - 2D real array
output:
* logical mask defining the selected seed points.

### Algorithm

First, assign each $x_i \in \text{data\_vectors}$ a density label. Use a or
reuse the function `density_labels` for this that populates a 2D real array of
length `n-columns(data_vectors)`. 

Sort the density labels descending. Start with the first and identify all input
data vectors within its radius. Mask those as visited. Continue with the next
highest non visited density label until none remain.

## Cluster identification

In parallel grow ensembles around each seed vector. Use OMP paradigms if and
only if `do concurrent` parallelism is not supported due to external library
calls (gfortran).

The SKG is called `cluster_identification` and is an iteration wrapper for the
below steps.

Consider each seed an ensemble of size one. Growing is done iteratively until
stop conditions are reached. Each iteration has an ensemble $\mathcal{E}_t$, an
2D array of $o$ observables $\mathcal{O}_{t-o+1}, \mathcal{O}_{t-o+2}, \dots,
\mathcal{O}_{t}$.

### First growth step

A seed is a single point: its `observable` would be the singular value
decomposition of a one-column centered matrix, which is degenerate -- no
meaningful $U$, $d$, or $G$ exists to compare against. `cluster_identification`
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
independently testable in isolation. `cluster_identification` is the one
place that already knows which iteration it is, without introducing any new
state elsewhere.

### Stop Conditions

Growth for a single seed's ensemble stops under exactly one of the following
four conditions. `cluster_identification` reports which one arose via an
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
`cluster_identification` to track a running count of how many growth
iterations have been accepted in total for the current seed -- a separate,
cumulative counter, independent of the trailing $o$-window history below,
since that window can be shorter than the true total once growth runs longer
than $o$ iterations.

### Output

`cluster_identification` returns, per seed:

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

### Local Radius Identification

For growing ensembles we need a locally adapted radius with which at each step
candidates for the addition are identified. Those candidates are vectors within
the local ensemble specific radius.

The SKG should be called `calc_ensemble_growth_radius` or similar.

Use a fixed-count kNN pool and compute the median distance among a seed's own
$k_{\min}$ nearest neighbors. Make $k_{\min}$ an optional argument with default
value 30.

For each seed we store the growth radius in a 1D real array called
`ensemble_growth_radii`.

### Tangent Space Variant

This section describes the implementation used for the tangent space variant of
STC. 

Each iteration of `cluster_identification` runs:

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

