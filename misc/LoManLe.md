# LoManLe — Local Manifold Learning

## 1. Objective

Given a noisy finite point cloud

$$
X=\{\mathbf{x}_1,\ldots,\mathbf{x}_n\},
\qquad
\mathbf{x}_i\in\mathbb{R}^{D},
$$

LoManLe estimates a low-dimensional manifold

$$
\mathcal{M}\subset\mathbb{R}^{D},
\qquad d\ll D,
$$

where $D$ is the ambient dimension and $d$ is the intrinsic manifold dimension.

The central approximation is that a sufficiently small neighborhood of a smooth manifold can be represented by its local tangent space. LoManLe estimates these tangent spaces from noisy observations using local covariance eigendecomposition / singular value decomposition (SVD), constructs overlapping local charts, projects observations onto them, stitches compatible charts, and optionally iteratively refines the resulting manifold.

Conceptually,

$$
\boxed{
\text{adaptive neighborhoods}
\rightarrow
\text{local tangent models}
\rightarrow
\text{overlapping atlas}
\rightarrow
\text{projection}
\rightarrow
\text{stitching}
\rightarrow
\text{refinement}
}
$$

The local tangent model is the first-order approximation of an underlying smooth manifold.

---

## 2. Input and output

### Input

A matrix

$$
X\in\mathbb{R}^{D\times n},
$$

with one column $\mathbf{x}_i$ per observation.

Required parameters:

- $k_{\min}$: minimum neighborhood size.
- $k_{\max}$: maximum allowed neighborhood size, potentially $n/4$.
- $g>1$: geometric neighborhood growth factor, currently approximately $1.25$.
- $d$: desired intrinsic manifold dimension, unless inferred locally.
- overlap requirement for atlas construction.
- numerical criteria for stable local tangent estimation.
- maximum number of refinement iterations.
- convergence/stopping criterion for iterative refinement.

Optional:

- infer $d$ locally from the singular/eigenvalue spectrum;
- smoothing of adaptive kernel bandwidths;
- noise-aware stitching;
- bifurcation/multifurcation detection.

### Output

At minimum:

1. projected manifold points

$$
\hat X=\{\hat{\mathbf{x}}_1,\ldots,\hat{\mathbf{x}}_n\};
$$

2. atlas anchors and their memberships;

3. local tangent bases;

4. local intrinsic dimensions if inferred;

5. atlas intersection/connectivity graph;

6. explicit connectivity of the learned manifold where available:
   - edges for $d=1$;
   - higher-dimensional mesh/simplicial representation for $d>1$, once implemented;

7. optionally local residual/noise models and reconstruction diagnostics.

---

## 3. Mathematical foundation

Assume locally that an unknown $d$-dimensional manifold $\mathcal M$ can be parameterized by

$$
f:\mathbb{R}^{d}\rightarrow\mathbb{R}^{D}.
$$

Around intrinsic coordinate $\mathbf{u}_0$,

$$
f(\mathbf{u}_0+\Delta\mathbf{u})
=
f(\mathbf{u}_0)
+
J_f(\mathbf{u}_0)\Delta\mathbf{u}
+
O(\|\Delta\mathbf{u}\|^2).
$$

Therefore, locally,

$$
\mathcal M
\approx
\mathbf{c}
+
\operatorname{span}(U_d),
$$

where $\mathbf c$ is a local center and

$$
U_d=
[\mathbf{u}_1,\ldots,\mathbf{u}_d]
$$

contains an orthonormal basis for the estimated tangent space.

LoManLe estimates $U_d$ from the dominant local variance directions.

Thus each chart is a **data-derived first-order local approximation** of the manifold.

---

## 4. Phase I — Spatial index

Construct a k-d tree or another nearest-neighbor index over $X$.

It provides

$$
\operatorname{kNN}(\mathbf{x}_i,k)
$$

for arbitrary $k$.

The index is used for adaptive neighborhood construction, density estimation and local geometric operations.

In sufficiently high ambient dimension, the nearest-neighbor implementation itself may need replacement because k-d trees deteriorate with increasing $D$.

---

## 5. Phase II — Adaptive neighborhood estimation

A fixed $k$ is undesirable because sampling density, curvature and noise vary across the manifold.

For every observation $\mathbf{x}_i$, begin with

$$
k_i=k_{\min}.
$$

Obtain its neighborhood

$$
N_i(k_i)=\operatorname{kNN}(\mathbf{x}_i,k_i).
$$

Compute the neighborhood center

$$
\boldsymbol{\mu}_i
=
\frac{1}{|N_i|}
\sum_{j\in N_i}\mathbf{x}_j.
$$

Construct the centered matrix

$$
Y_i=
[
\mathbf{x}_{j_1}-\boldsymbol{\mu}_i,
\ldots,
\mathbf{x}_{j_k}-\boldsymbol{\mu}_i
].
$$

Equivalently compute the local covariance matrix

$$
C_i=
\frac{1}{k_i-1}Y_iY_i^\top.
$$

Its eigendecomposition is

$$
C_i
=
U_i\Lambda_iU_i^\top,
$$

with

$$
\lambda_{i,1}\geq
\lambda_{i,2}\geq\cdots\geq
\lambda_{i,D}\geq0.
$$

This is equivalent, for the required purpose, to obtaining the local principal directions through SVD of $Y_i$.

---

## 6. Phase III — Grow neighborhoods until tangent estimation is stable

For each $\mathbf{x}_i$, repeatedly increase the neighborhood geometrically:

$$
k_i^{(t+1)}
=
\min
\left(
\left\lceil gk_i^{(t)}\right\rceil,
k_{\max}
\right).
$$

At every scale, estimate the local eigensystem.

The neighborhood is accepted at the **smallest scale** at which the local manifold estimate is sufficiently stable.

This early stopping is essential: increasing $k$ reduces variance but eventually introduces curvature bias and may mix nearby branches.

Hence:

$$
\boxed{\text{choose the smallest neighborhood giving a stable tangent estimate}}
$$

rather than the largest statistically convenient neighborhood.

### Tangent stability

For fixed intrinsic dimension $d$, compare tangent spaces between consecutive neighborhood sizes.

Let

$$
U_d^{(t)},U_d^{(t+1)}
$$

be the two tangent bases.

Use their principal angles

$$
0\leq\theta_1\leq\cdots\leq\theta_d\leq\frac{\pi}{2}.
$$

For example, require

$$
\theta_{\max}
=
\max_j\theta_j
\leq \tau_\theta
$$

for one or more consecutive growth steps.

The criterion should operate on the **subspace**, not individual eigenvectors, because rotations within a $d$-dimensional tangent space represent the same tangent space.

### Spectral separation

Additionally require evidence that tangent variance is distinguishable from normal variance, for example through

$$
G_i(d)
=
\frac{\lambda_{i,d}}
{\lambda_{i,d+1}+\epsilon}.
$$

Require

$$
G_i(d)\geq\tau_G.
$$

The precise stability criterion and thresholds are algorithm parameters and must be reported.

If no acceptable neighborhood is found before $k_{\max}$, flag that location as geometrically unresolved rather than silently treating the largest neighborhood as reliable.

---

## 7. Optional local intrinsic dimensionality estimation

Instead of prescribing $d$, inspect

$$
\lambda_1\geq\lambda_2\geq\cdots\geq\lambda_D.
$$

Candidate dimensionalities correspond to stable spectral separations such as

$$
G(r)
=
\frac{\lambda_r}
{\lambda_{r+1}+\epsilon}.
$$

Select a dimension $d_i$ only if the inferred dimension and associated tangent subspace remain stable while the neighborhood grows.

Thus different charts may have

$$
d_i\neq d_j.
$$

This permits locally changing intrinsic dimensionality, although stitching charts of different dimensions requires explicit topology rules.

---

## 8. Phase IV — Local density / reliability labels

Each point receives a local density or reliability measure derived from its adaptive neighborhood.

For a Gaussian kernel,

$$
w_{ij}
=
\exp
\left(
-\frac{\|\mathbf{x}_j-\mathbf{x}_i\|^2}
{2\sigma_i^2}
\right).
$$

A simple local density estimate is

$$
\rho_i
=
\sum_{j\in N_i}w_{ij}.
$$

Alternatively a weighted-distance measure can be retained if that is the implementation used.

The adaptive bandwidth $\sigma_i$ is derived from the local k-nearest-neighbor scale.

Because raw k-nearest-neighbor bandwidths can vary abruptly from point to point, the bandwidth field may itself be smoothed over neighboring observations before subsequent use.

Density/reliability is used primarily for **atlas construction and weighting**, not as a definition of the manifold itself.

---

## 9. Phase V — Greedy atlas construction

A chart is defined by an anchor $a$, its adaptive neighborhood $N_a$, center $\boldsymbol\mu_a$, tangent basis $U_a$, and optionally local intrinsic dimension $d_a$.

Sort candidate anchors by decreasing reliability/density:

$$
\rho_{(1)}\geq\rho_{(2)}\geq\cdots\geq\rho_{(n)}.
$$

Select the first high-density candidate as an anchor.

Subsequent anchors are chosen such that:

1. they add previously insufficiently covered observations;
2. neighboring charts retain a prescribed overlap;
3. already adequately represented regions do not generate redundant anchors.

The intended result is

$$
\mathcal A=\{A_1,\ldots,A_m\}
$$

covering the sampled manifold with overlapping local charts.

This anchor-selection procedure is presently an **algorithmic heuristic**, rather than part of the local tangent-space estimator itself.

---

## 10. Phase VI — Recompute definitive chart geometry

For each selected anchor $A_a$, gather its final adaptive neighborhood

$$
N_a.
$$

Compute its final center

$$
\boldsymbol{\mu}_a
=
\frac{1}{|N_a|}
\sum_{i\in N_a}\mathbf{x}_i
$$

and tangent basis

$$
U_a\in\mathbb R^{D\times d_a}.
$$

The local manifold chart is

$$
M_a
=
\left\{
\boldsymbol{\mu}_a+U_a\mathbf{z}
:
\mathbf{z}\in\mathbb R^{d_a}
\right\},
$$

restricted to the local region represented by the chart.

---

## 11. Phase VII — Project observations onto local manifolds

For observation $\mathbf{x}_i\in N_a$, define

$$
\mathbf{v}_{ia}
=
\mathbf{x}_i-\boldsymbol{\mu}_a.
$$

Its tangent coordinates are

$$
\mathbf{z}_{ia}
=
U_a^\top\mathbf{v}_{ia}.
$$

Its projection onto chart $a$ is

$$
\boxed{
\hat{\mathbf{x}}_{ia}
=
\boldsymbol{\mu}_a+
U_aU_a^\top
(\mathbf{x}_i-\boldsymbol{\mu}_a)
}
$$

and its residual is

$$
\mathbf{r}_{ia}
=
\mathbf{x}_i-\hat{\mathbf{x}}_{ia}.
$$

By construction,

$$
U_a^\top\mathbf{r}_{ia}=0.
$$

The scalar orthogonal reconstruction error is

$$
r_{ia}
=
\|\mathbf{r}_{ia}\|_2.
$$

Importantly, **the same original observation identifier $i$ is retained for every chart projection**. Thus projections of the same observation into overlapping charts provide natural correspondences for stitching.

---

## 12. Phase VIII — Build the chart-overlap graph

Construct the membership relation

$$
P_{ia}
=
\begin{cases}
1,&\mathbf{x}_i\in N_a,\\
0,&\text{otherwise}.
\end{cases}
$$

Two charts are potential neighbors if

$$
N_a\cap N_b\neq\varnothing.
$$

Define the preliminary chart graph

$$
G=(V,E),
$$

where

$$
V=\{A_1,\ldots,A_m\}
$$

and

$$
(A_a,A_b)\in E
\iff
N_a\cap N_b\neq\varnothing.
$$

Sparse relations such as

- point $\rightarrow$ overlap edges,
- edge $\rightarrow$ shared points

are stored using **Compressed Sparse Row (CSR)** structures.

Connected overlap regions can then be identified using **Breadth-First Search (BFS)**.

However:

$$
\boxed{\text{neighborhood overlap does not imply manifold connectivity}}
$$

because overlapping spheres can represent:

- consecutive pieces of one manifold;
- a true bifurcation or multifurcation;
- intersecting manifolds;
- separate manifold branches passing close to each other.

Therefore sphere overlap creates only a **candidate stitching relation**.

---

## 13. Phase IX — Estimate the local normal/noise structure

For each chart $a$, use residual vectors

$$
\mathbf r_{ia}
$$

to estimate the distribution of observations normal to the local manifold.

Let

$$
U_a^\perp
$$

span the orthogonal complement of $U_a$.

Normal coordinates are

$$
\mathbf q_{ia}
=
(U_a^\perp)^\top\mathbf r_{ia}.
$$

Estimate their covariance

$$
\Sigma_{\perp,a}
=
\operatorname{Cov}
\{\mathbf q_{ia}:i\in N_a\}.
$$

This represents the local anisotropic noise/thickness around the manifold.

A scalar tube is the simpler special case

$$
R_a
=
Q_q
\left(
\{\|\mathbf r_{ia}\|\}
\right),
$$

where $Q_q$ is a robust local quantile such as the 95th percentile.

The covariance representation is preferable when sufficient data are available because normal variation need not be isotropic.

---

## 14. Phase X — Decide whether candidate charts should actually be stitched

For every candidate chart pair $(a,b)$ from the overlap graph, determine their closest relevant manifold locations

$$
\mathbf m_a\in M_a,\qquad
\mathbf m_b\in M_b.
$$

Define

$$
\boldsymbol\delta_{ab}
=
\mathbf m_b-\mathbf m_a.
$$

The charts are eligible for stitching only if their local manifold models approach each other within the noise/thickness supported by the original observations.

In the scalar approximation,

$$
\|\boldsymbol\delta_{ab}\|
\leq
R_a+R_b
$$

provides a simple tube-overlap condition.

With anisotropic residual models, evaluate the separation relative to the combined normal uncertainty. Conceptually,

$$
\Sigma_{\perp,ab}
=
\Sigma_{\perp,a}
+
\Sigma_{\perp,b},
$$

after expressing both covariance structures in a compatible ambient-space representation.

Then evaluate a generalized Mahalanobis-type separation

$$
D_{ab}^2
=
\boldsymbol\delta_{ab}^{\top}
\Sigma_{\perp,ab}^{+}
\boldsymbol\delta_{ab},
$$

where $+$ denotes the Moore-Penrose pseudoinverse because normal covariance is generally rank deficient in ambient coordinates.

Require

$$
D_{ab}^2\leq\tau_{\text{noise}}.
$$

Thus:

> Two charts that merely pass nearby are kept separate if the distance between their manifolds exceeds the locally observed residual/noise structure.

---

## 15. Tangent compatibility

Noise-tube overlap alone is insufficient.

For charts $a,b$, compute principal angles between their tangent spaces.

For equal dimension $d$,

$$
\cos\theta_j
=
s_j(U_a^\top U_b),
$$

where $s_j$ are the singular values.

These angles distinguish possible relationships.

### Continuation

Small principal angles indicate that two charts likely represent consecutive pieces of the same smooth manifold.

### Furcation

Charts whose noise-supported manifold regions meet but whose tangent structures diverge may represent a bifurcation or multifurcation.

Such a furcation should be introduced only when:

1. the local noise regions support physical contact;
2. original observations populate the junction region — the **non-gap requirement**;
3. tangent geometry supports multiple outgoing manifold directions — the **angle requirement**.

### Near miss

If anchor spheres overlap but manifold noise regions do not, do **not** stitch.

Therefore the conceptual rule is

$$
\boxed{
\text{stitch}
=
\text{local support}
\land
\text{noise-region compatibility}
\land
\text{topological/tangent compatibility}.
}
$$

The exact thresholds for distinguishing continuation, furcation and crossing remain parameters requiring empirical validation.

---

## 16. Phase XI — Stitch projections from overlapping charts

Suppose observation $i$ belongs to charts

$$
A_{a_1},\ldots,A_{a_s}.
$$

It then has chart-specific projections

$$
\hat{\mathbf{x}}_{ia_1},
\ldots,
\hat{\mathbf{x}}_{ia_s}.
$$

Because all projections retain original point identifier $i$, no correspondence search is necessary.

The simplest stitched position is

$$
\hat{\mathbf{x}}_i
=
\frac{1}{s}
\sum_{j=1}^{s}
\hat{\mathbf{x}}_{ia_j}.
$$

Preferably use chart reliability weights,

$$
\hat{\mathbf{x}}_i
=
\frac{
\sum_j w_{ia_j}\hat{\mathbf{x}}_{ia_j}
}{
\sum_jw_{ia_j}
}.
$$

Weights may incorporate:

- local density;
- projection residual;
- tangent stability;
- local noise variance;
- distance from the chart center.

An inverse-variance form is

$$
w_{ia}
\propto
\frac{1}
{\sigma_{ia}^2+\epsilon}.
$$

Only chart projections belonging to a geometrically accepted common manifold component are combined.

---

## 17. Phase XII — Construct explicit manifold topology

The stitched projected points

$$
\hat X
$$

are samples from the estimated manifold, but are not by themselves an explicit manifold representation.

For intrinsic dimension $d=1$, construct a graph connecting neighboring manifold points while respecting the accepted chart topology.

This yields a piecewise-linear principal curve/skeleton.

For $d=2$, the corresponding representation requires faces, e.g. local triangulation, constrained by chart membership and accepted chart adjacency.

For general $d$, the analogous representation is a local simplicial complex or atlas representation.

Crucially, nearest Euclidean projected points must **not** automatically be connected, because folds can make geodesically distant manifold regions close in ambient space.

Connectivity should therefore derive primarily from the accepted chart graph and local manifold coordinates.

The current $d=1$ connectivity algorithm exists; generalization to $d>1$ remains an explicit development task.

---

## 18. Phase XIII — Iterative manifold refinement

The stitched manifold can retain small discontinuities, roughness and chart-boundary artifacts.

Apply iterative manifold learning/smoothing to the stitched projected points.

The important principle is that refinement should reduce noise without destroying the geometry recovered by the atlas.

For each manifold point $\hat{\mathbf x}_i$, obtain an adaptive local neighborhood and Gaussian weights

$$
w_{ij}
=
\exp
\left(
-\frac{\|\hat{\mathbf{x}}_j-\hat{\mathbf{x}}_i\|^2}
{2\sigma_i^2}
\right).
$$

The bandwidth field $\sigma_i$ may itself be smoothed to avoid abrupt changes between neighboring kernels.

The update produces a new manifold estimate

$$
\hat X^{(t+1)}
=
F(\hat X^{(t)}).
$$

The previously developed iterative ManLe/ANWIL-style procedure can provide $F$.

A critical implementation choice is whether the current point contributes to its own update. Excluding it entirely can produce excessive displacement in sparsely supported regions; therefore this choice should be treated explicitly and validated.

Likewise, hard truncation at exactly the $k$-th neighbor can introduce discontinuities as observations enter or leave adjacent kernels. Full Gaussian support, efficient radius truncation at negligible weights, or smoothly varying neighborhoods should be considered.

---

## 19. Stop refinement before oversmoothing

At each iteration measure at least:

- roughness $r_t$;
- coverage $c_t$;
- reconstruction error, e.g. root mean squared error $e_t$.

A geometric mean of these quantities was experimentally found unsuitable because it continued refinement until the maximum iteration count.

A weighted arithmetic objective performs better:

$$
L_t
=
w_r\,\tilde r_t
+
w_c\,\tilde c_t
+
w_e\,\tilde e_t,
$$

where quantities must be normalized to comparable scales and their orientation chosen so that smaller $L$ consistently means better.

Experiments indicate useful relative emphasis around

$$
w_r:w_c\approx2:3
$$

or

$$
3:2,
$$

depending on the dataset.

Refinement stops when further iteration no longer improves the objective according to the defined stopping rule.

The exact objective and stop condition are empirical hyperparameters and should not yet be treated as universal constants.

---

## 20. Complete algorithm

```text
LOMANLE(X, parameters)

    # ----------------------------------------------------------
    # A. Local geometry
    # ----------------------------------------------------------

    build nearest-neighbor index for X

    for every observation i:

        k ← k_min

        repeat:

            N ← kNN(x_i, k)

            estimate local center mu
            estimate local covariance / SVD
            estimate tangent subspace U_d

            if intrinsic dimension is adaptive:
                infer candidate d from spectral structure

            compare tangent structure with previous scale

            if tangent structure is stable
               AND tangent/normal spectral separation is sufficient:
                   accept neighborhood
                   break

            k ← ceil(g * k)

        until k >= k_max

        store:
            neighborhood
            local scale
            density/reliability
            tangent stability
            intrinsic dimension
            spectral information


    optionally smooth the local bandwidth/scale field


    # ----------------------------------------------------------
    # B. Atlas construction
    # ----------------------------------------------------------

    rank candidate anchors by density/reliability

    greedily select anchors such that:
        manifold coverage increases
        required overlap is maintained
        unnecessary redundant charts are avoided

    for each selected anchor:
        recompute definitive neighborhood
        compute center mu_a
        compute tangent basis U_a
        project member observations onto local manifold
        retain original observation identifiers
        compute residual vectors


    # ----------------------------------------------------------
    # C. Candidate topology
    # ----------------------------------------------------------

    construct point × anchor membership

    create candidate edge (a,b)
        whenever charts share observations

    store sparse incidence structures using CSR

    identify connected overlap regions using BFS


    # ----------------------------------------------------------
    # D. Geometric topology
    # ----------------------------------------------------------

    for every candidate chart edge (a,b):

        estimate residual/noise structure of chart a
        estimate residual/noise structure of chart b

        determine closest relevant manifold regions

        if manifold regions do not overlap within
           their supported noise regions:
               reject edge
               continue

        analyze tangent-space compatibility

        if geometry supports continuation:
               classify as continuation edge

        else if:
               observations support the junction
               AND tangent geometry supports branching:
               classify as furcation edge

        else:
               reject or mark unresolved


    # ----------------------------------------------------------
    # E. Stitching
    # ----------------------------------------------------------

    for every original observation represented
    in multiple compatible charts:

        collect its chart-specific projections

        combine projections using
        reliability / inverse-variance weights

        obtain one stitched manifold point


    # ----------------------------------------------------------
    # F. Explicit topology
    # ----------------------------------------------------------

    if d == 1:
        connect manifold points into edges using
        local coordinates + accepted chart topology

    if d > 1:
        construct local faces/simplices from chart topology
        and local manifold coordinates


    # ----------------------------------------------------------
    # G. Refinement
    # ----------------------------------------------------------

    repeat:

        perform adaptive manifold smoothing

        evaluate:
            roughness
            coverage
            reconstruction error

        compute weighted objective

        if stop criterion reached:
            break

    until max_iterations


    return manifold points,
           manifold topology,
           atlas,
           local tangent models,
           residual/noise models,
           diagnostics
```

---

## 21. What is fundamental LoManLe and what is currently heuristic

It is useful to make this distinction explicit.

### Geometric core

The following has a direct mathematical interpretation:

$$
\boxed{
\text{local neighborhood}
\rightarrow
\text{local SVD}
\rightarrow
\text{tangent subspace}
\rightarrow
\text{orthogonal projection}
}
$$

For a sufficiently smooth underlying manifold and sufficiently local neighborhoods, this estimates its first-order local geometry.

### Adaptive estimation

Growing neighborhoods until tangent structure becomes stable addresses the variance-versus-curvature trade-off:

$$
\text{too small}
\rightarrow
\text{unstable/noisy tangent}
$$

$$
\text{too large}
\rightarrow
\text{curvature/branch mixing}.
$$

The smallest stable neighborhood is therefore the desired operating point.

### Atlas construction

Anchor selection is presently heuristic. Its purpose is computational and representational: obtain sufficient overlapping local models without fitting a chart around every observation.

### Stitching/topology

This is currently the least mathematically settled component.

The developing principle is:

$$
\boxed{
\text{chart overlap proposes connectivity;
data-supported local geometry decides connectivity.}
}
$$

Residual/noise structure, tangent compatibility and observed support around a junction provide the evidence.

### Iterative refinement

Refinement is also algorithmic rather than part of the first-order tangent approximation. Its purpose is removal of finite-sample and stitching artifacts without destroying the recovered geometry.

---

## 22. Principal limitations that must remain explicit

LoManLe does **not** escape the curse of dimensionality.

Its most important expected failure modes are:

**High ambient dimension $D$.** Euclidean nearest-neighbor relationships may deteriorate, and accumulated noise across many normal dimensions can obscure tangent structure.

**High intrinsic dimension $d$.** Reliable estimation of a $d$-dimensional tangent space requires increasingly many local observations.

**Sparse sampling.** No local algorithm can recover geometry that is unsupported by observations.

**High curvature.** A neighborhood large enough for stable variance estimation may already be too large for a first-order tangent approximation.

**Nearby folds/branches.** k-nearest-neighbor neighborhoods can combine observations that are close in ambient Euclidean distance but far apart geodesically.

**Junctions.** A single tangent space may be undefined or misleading at a true bifurcation/multifurcation.

**Boundaries.** Local covariance becomes asymmetric near manifold boundaries and requires explicit recognition.

**Variable density.** Neighborhood size and reliability must adapt locally.

These are not merely implementation issues; they define the regime in which the method must be experimentally characterized.

---

## 23. Compact definition

LoManLe can therefore be defined as:

> **An adaptive local manifold estimator that reconstructs a noisy low-dimensional manifold embedded in multivariate space by finding the smallest neighborhoods supporting stable local tangent-space estimates, representing those estimates as overlapping affine charts, orthogonally projecting observations onto the charts, inferring global topology from data-supported chart relationships, stitching corresponding projections, and optionally refining the resulting manifold while controlling reconstruction quality and roughness.**

Or mathematically, its core operation is simply

$$
\boxed{
\mathbf{x}
\longmapsto
\boldsymbol{\mu}_a
+
U_aU_a^\top
(\mathbf{x}-\boldsymbol{\mu}_a)
}
$$

where both the appropriate local scale and the local tangent basis $U_a$ are learned from the noisy point cloud.

Everything around that operation answers three questions:

$$
\boxed{
\begin{aligned}
&\text{At what scale is this local approximation trustworthy?}\\
&\text{Which local approximations belong to the same manifold?}\\
&\text{How should those approximations be assembled globally?}
\end{aligned}}
$$

Those are, at present, the defining problems LoManLe is trying to solve.
