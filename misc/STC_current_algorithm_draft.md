---
title: |
    | Shape Truthful Clustering (STC) — Algorithm and Abstract API
author: Asis Hallab
date: \today
geometry: margin=2.0cm
numbersections: true
---

## Motivation

Shape Truthful Clustering (STC) is a clustering procedure inspired by the renormalization-group idea in statistical physics: repeatedly coarse-grain a system and track how a macroscopic observable behaves under that coarse-graining. A scale-invariant region is one whose observable stops changing (beyond some tolerance) as the region grows further — a fixed point of the coarse-graining operator. STC applies this idea to a point cloud: an ensemble (candidate cluster) is grown from a seed, an observable is computed after every growth step, and growth stops once the observable is judged stable. The accepted ensemble is a "shape truthful" cluster: its description is faithful to its own extent, not an artifact of an arbitrarily chosen scale.

This document specifies STC as an abstract algorithm with a fixed five-step shape and a pluggable observable, then gives two concrete instantiations: a density-based one (the first working prototype) and a local-tangent-space one, which infers piecewise-linear manifold structure rather than density clusters.

## Abstract API

STC is defined by five conceptual steps. Each step is pluggable: concrete instantiations differ in what an "observable" is and how growth and acceptance are defined, and therefore differ in argument types and counts. For this reason STC is *not* implemented as a single Fortran generic interface; each instantiation provides its own concrete subroutines following the same five-step shape, rather than dispatching through one shared interface block.

1. **Seed identification.** Identify seed vectors from which ensembles are grown. Each seed starts an ensemble of size one.
2. **Growth.** Expand an existing ensemble by adding candidate vectors. The default mechanism is a fixed surface-radius expansion via a k-d-tree, but other growth rules are anticipated — for example, adding vectors that are ambient-distant but similar under some other criterion (metadata, direction, text).
3. **Observable computation.** Compute an observable for the grown ensemble. The observable is a vector-valued summary of the ensemble's current macroscopic state. Retain either the last few or all iterations' observables, so that change across iterations can be evaluated.
4. **Accept or reject.** Compare the current observable against one or more previous iterations and decide whether the growth step is accepted. If rejected, the ensemble from the previous iteration is the final result for that seed.
5. **Reconciliation.** Once every seed's ensemble has converged, identify relationships between ensembles (typically their pairwise intersections). Whether and how to merge on top of that is instantiation-specific: some applications want overlapping ensembles merged, others only want the intersections reported and merging left to a downstream consumer.

```text
seed_identification()
    -> grow()
    -> compute_observable()
    -> accept()
    -> (repeat 2-4 until rejected, per seed, seeds independent)
    -> reconcile()
```

Instantiations anticipated for this API include: a **density** observable (candidate count per volume; described in full below, the first implemented case), an **angular / directional similarity** observable (coherence of orientation among ensemble members), a **textual / metadata similarity** observable (agreement among symbolic annotations attached to each vector, rather than the vectors' ambient coordinates), and a **local tangent space** observable (local linear manifold geometry; described in full below).

## Definitions common to all instantiations

Let $\mathcal V=\{\mathbf v_1,\ldots,\mathbf v_N\}\subset\mathbb R^D$ be the ambient vectors. An ensemble after growth iteration $t$ is a subset $\mathcal E_t\subseteq\mathcal V$. An observable is a function $O_t=f(\mathcal E_t)$ returning a fixed-length real vector (or tuple of such). A difference operator $\mathcal D(O_t,O_{t+1})$ compares consecutive observables; the accept criterion is a threshold test on $\mathcal D$'s output.

## Density-based instantiation (first prototype)

### Ambient density label

$$
\rho_i=\sum_{j=1}^{N}\mathbf 1\!\left(d(\mathbf v_i,\mathbf v_j)\le r\right),
$$

where $r$ is a fixed neighborhood radius (default: the 15th percentile of distances from every vector to the global mean vector).

### Seed identification

Choose high-density seeds (default: top 15 percent by $\rho_i$). Each seed initializes one ensemble.

### Growth

Find all vectors adjacent to the ensemble's current surface (via the k-d-tree). A candidate is accepted into the candidate set $\mathcal C_t$ iff

$$
\left|\rho_i-\operatorname{median}(\rho_{\mathcal E_t})\right|
\le
\alpha_{\mathrm{MAD}}\cdot \mathrm{MAD}_{\mathrm{ambient}},
\qquad
\mathrm{MAD}_{\mathrm{ambient}}=\operatorname{median}_i\left(\left|\rho_i-\operatorname{median}(\rho)\right|\right),
$$

i.e. candidates must have a density label within $\alpha_{\mathrm{MAD}}$ median absolute deviations (over the whole ambient set) of the ensemble's own current median density. $\mathrm{MAD}$ is used instead of the standard deviation for robustness to outliers and skewed density distributions. All qualifying candidates are added at once (atomic growth):

$$
\mathcal E_{t+1}=\mathcal E_t\cup\mathcal C_t.
$$

No sequential ordering of candidates is performed within a growth step.

### Observable

$$
O_t=\left(\rho_{\mathrm{arith}},\ \rho_{\mathrm{harm}},\ H_\rho,\ |\mathcal E_t|,\ |\mathcal C_t|\right),
\qquad
H_\rho=\frac{\rho_{\mathrm{arith}}}{\rho_{\mathrm{harm}}},
$$

the arithmetic and harmonic mean density over the ensemble, their ratio (density heterogeneity), the ensemble size, and the candidate-set size. $H_\rho$ is redundant given the first two entries but is stored for convenience and debugging. Observable history is stored as a real matrix with one column per iteration.

### Accept / reject

Reject growth (return the previous ensemble $\mathcal E_t$ as final) if

$$
\left|\log_2\!\left(\frac{H_{\rho,t+1}}{H_{\rho,t}}\right)\right|\ge \alpha_{\mathrm{accept}}
$$

(default $\alpha_{\mathrm{accept}}=1$, a two-fold change in density heterogeneity), or if the maximum iteration count is reached. Otherwise continue growing.

### Output

A logical membership matrix (one column per ensemble, one row per ambient vector), the observable trajectories, and, optionally, a post-processing merge of ensembles that share members.

## Local tangent-space instantiation

### Motivation

Rather than a density-homogeneous region, this instantiation grows a region over which a single affine tangent plane is a faithful first-order approximation of an underlying manifold, and monitors that approximation's own stability instead of a density statistic. Growth is accepted once the region's tangent subspace, inferred intrinsic dimension, and spectral separation all stop changing appreciably from one iteration to the next.

### Seed identification

Sequential, density-ranked selection, **excluding vectors already covered** by an already-accepted ensemble. Unlike the density instantiation's independent top-percentile seeding, seeds here must not be chosen freely in parallel: without an exclusion rule, a dense region would spawn many near-duplicate, heavily overlapping ensembles, since reconciliation (below) does not deduplicate.

### Growth and local scale

Fixed-radius surface expansion via the k-d-tree, as in the density instantiation, but the radius is derived **per seed** rather than once globally. Each seed's own local scale is the median distance to its $k_{\min}$ nearest neighbors:

$$
\text{local\_scale}=
\begin{cases}
\text{dist. to the } \frac{k_{\min}+1}{2}\text{-th nearest neighbor}, & k_{\min}\text{ odd}\\[4pt]
\tfrac12\left(\text{dist. to the }\tfrac{k_{\min}}{2}\text{-th}+\text{dist. to the }\left(\tfrac{k_{\min}}{2}+1\right)\text{-th nearest neighbor}\right), & k_{\min}\text{ even}.
\end{cases}
$$

Using a single dataset-wide radius here would reintroduce the classical problem this instantiation exists to avoid: sampling density, curvature, and noise vary across a manifold, so a fixed scale is systematically wrong somewhere. Deriving the radius per seed keeps growth locally adaptive.

### Observable

For an ensemble $\mathcal E_t$ with member coordinates centered on their mean $\boldsymbol\mu_t$, form the covariance $C_t$ and its eigendecomposition $C_t=U\Lambda U^\top$, eigenvalues $\lambda_1\ge\cdots\ge\lambda_D\ge0$. The observable is

$$
O_t=(U_t,\ d_t,\ G_t),
$$

with $U_t\in\mathbb R^{D\times d_t}$ the top-$d_t$ eigenvectors (the tangent basis), $d_t$ the inferred intrinsic dimension, and

$$
G_t=\frac{\lambda_{d_t}}{\lambda_{d_t+1}+\varepsilon}
$$

the spectral-gap statistic separating tangent variance from normal variance. The same eigendecomposition additionally yields, at no extra cost, the ensemble's center $\boldsymbol\mu_t$, the mean squared residual off the tangent subspace ("normal error"), the per-tangent-direction extent, and the normal-space covariance

$$
\Sigma_\perp = U^\perp \Lambda^\perp (U^\perp)^\top,
$$

built from the $D-d_t$ eigenvector/eigenvalue pairs *not* selected as tangent directions — these are already computed while forming $U_t$ and are otherwise discarded. $\Sigma_\perp$ characterizes the ensemble's local noise thickness, anisotropically, and is retained for use by reconciliation (below) and by any downstream consumer needing to test whether two ensembles' manifolds plausibly meet within their observed noise.

### Difference operator

Given tangent bases $U_t,U_{t+1}\in\mathbb R^{D\times d}$ (equal dimension across the compared iterations), form $M=U_t^\top U_{t+1}$ and its SVD $M=P\Sigma_M Q^\top$. The singular values $\sigma_j$ are the cosines of the principal angles between the two subspaces, $\theta_j=\arccos(\sigma_j)$. Define

$$
\delta_U=\max_j\theta_j,
\qquad
\delta_d=|d_{t+1}-d_t|,
\qquad
\delta_G=\left|\log\frac{G_{t+1}}{G_t}\right|.
$$

$\delta_U$ operates on the subspace as a whole, not on individual eigenvectors, since a rotation within a $d$-dimensional tangent space represents the same tangent space.

### Accept / reject

Accept the grown ensemble, and continue growing, while

$$
\delta_U<\tau_U,
\qquad
\delta_d=0,
\qquad
\delta_G<\tau_G.
$$

Reject as soon as any one condition fails, or the maximum iteration count is reached; the previous iteration's ensemble is then the final result for that seed. Requiring $\delta_d=0$ means the ensemble's inferred intrinsic dimension, not only its tangent direction, must itself be stable before growth is trusted.

### Reconciliation and pairwise compatibility

Reconciliation reports pairwise intersections between accepted ensembles; it does not merge them, leaving that decision to the downstream consumer. Because every ensemble already carries $\Sigma_\perp$ from its observable, a downstream consumer can test whether two intersecting ensembles $a,b$ plausibly represent the same underlying manifold, rather than two manifolds that merely happen to pass close to one another, via a Mahalanobis-type statistic:

$$
D_{ab}^2=\delta_{ab}^\top\, \Sigma_{\perp,ab}^{+}\, \delta_{ab},
\qquad
\Sigma_{\perp,ab}=\Sigma_{\perp,a}+\Sigma_{\perp,b},
$$

where $\delta_{ab}$ is the separation between the two ensembles' nearest tangent-plane points and $+$ denotes the Moore-Penrose pseudoinverse ($\Sigma_\perp$ is generally rank-deficient in ambient coordinates). Two ensembles are compatible if $D_{ab}^2\le\tau_{\text{noise}}$.

Tangent-angle divergence between two otherwise-compatible ensembles should be used to *classify* the relationship (small angle: continuation of the same manifold; large or divergent angle: a branch point) and not to reject it. A true bifurcation can have an arbitrarily large, even orthogonal, angle between its branches; requiring small angle as a condition for compatibility would incorrectly discard genuine junctions.

### Pseudocode

```text
for each seed (independently, in parallel if desired):
    E = {seed}
    O_prev = observable(E)
    repeat:
        C = candidates_within_radius(E, local_scale(seed))
        E_next = E union C
        O_next = observable(E_next)
        d = difference(O_prev, O_next)
        if d.delta_U < tau_U and d.delta_d == 0 and d.delta_G < tau_G:
            E = E_next
            O_prev = O_next
        else:
            break                      # E is the accepted ensemble for this seed
    until max_iterations reached

reconcile(all accepted ensembles)      # pairwise intersections only, no merge
```

## Extensions

Additional observables are anticipated for further instantiations of the same abstract API: entropy, free energy, angular coherence, and semantic (metadata/text) descriptors. In general, an iteration's observable vector $\Omega_t$ (density, tangent, or otherwise) traces a trajectory $\Omega_1,\ldots,\Omega_T$ that itself forms a point cloud in observable space. Future work may analyze that trajectory directly — via its covariance, bootstrap confidence regions, Mahalanobis distances, Hotelling's $T^2$, or manifold learning applied to the trajectory itself — to discover which observables are scale-invariant and discriminative for a given application, rather than fixing that choice in advance.
