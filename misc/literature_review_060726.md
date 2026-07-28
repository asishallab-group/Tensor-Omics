Each paper is briefly described and compared against what we do in lomanle. For all of them, we write down what we could leverage in the future.

These articles were reviewed during the week of July 6, 2026, so the comparison with the code corresponds to that date and not necessarily to the final code.

# 1. Before we start: what does “learning a manifold” mean?

Let's imagine we have many points scattered in a space.

Even though the points may have many coordinates, they could be organized around a much simpler shape:

* a curved line;
* a spiral;
* a folded surface;
* several branches;
* a shape with bifurcations;
* several surfaces that cross each other.

The problem is that we only have the points. We don't know beforehand:

* which are the true neighbors;
* what the direction of the shape is;
* how many local dimensions it has;
* where there are curves, branches, or intersections;
* how to join the small local regions.

The papers reviewed solve different parts of that problem.

Some want to **unroll** the manifold and give low-dimensional coordinates. Others want to **reconstruct or smooth it**. Others simply want to discover good neighborhoods or separate regions with different dimensions.

Our current LoManLe is closer to a combination of:

1. adaptive neighborhoods;
2. local PCA;
3. construction of a sphere atlas;
4. projection onto tangent axes;
5. merging of the projections in the overlaps.

---

# 2. What our LoManLe currently does

Before comparing the papers, it's worth describing our code exactly.

## 2.1 General flow

Our implementation roughly follows this process:

```text
Original points
       ↓
k-d tree construction
       ↓
k-NN neighbors for each point
       ↓
Adaptive increase of k
       ↓
Local PCA and spectral gap computation
       ↓
Local radius and density
       ↓
Anchor selection by density
       ↓
Construction of partially overlapping spheres
       ↓
New PCA on each anchor-sphere
       ↓
Detection of spheres that share points
       ↓
Grouping of intersection regions via BFS
       ↓
Projection onto each anchor's principal direction
       ↓
Weighted average in zones covered by several anchors
       ↓
Final skeleton
```

## 2.2 Adaptive neighborhoods

For each point we start with `k_min`.

We compute its neighbors and do PCA. Then we look at the gap between:

* the last eigenvalue we consider tangential;
* the first eigenvalue we consider normal.

Essentially we ask:

> “Can I already clearly see an important direction separate from the noise or the normal directions?”

If the gap isn't large enough, we increase `k` by about 25% and recompute the PCA.

This is a good idea: we don't force every point to use exactly the same number of neighbors.

However, the adaptation is based exclusively on the **spectral gap**. We don't explicitly check whether the neighborhood:

* crosses into another branch;
* jumps across a gap;
* joins two parts that are close in Euclidean distance but far apart on the manifold;
* mixes points from a bifurcation.

## 2.3 Density

Our density is roughly:

[
\rho_i =
\frac{\sum_j \exp\left(-d_{ij}^2/(2\sigma_i^2)\right)}
{\sigma_i^{d}}
]

where (d) is `manifold_dim`.

This combines:

* how many points are nearby;
* how close they are;
* the approximate volume of the neighborhood.

We then sort the points from highest to lowest density. Dense points get priority for becoming anchors.

## 2.4 Atlas construction

The first anchor is the highest-density point.

Subsequent anchors must have an overlap within the interval:

```text
o_min ≤ overlap proportion ≤ o_max
```

The idea is:

* if there's no overlap at all, the charts can't connect;
* if they overlap too much, they're redundant.

When no point satisfies this condition, our code opens a new component by taking the discovered point with the highest density.

That lets us handle disconnected data and, potentially, branches. However, it doesn't mean the code has proven that a bifurcation exists. It only means:

> “I couldn't find another anchor that satisfied the current overlap rules.”

It could be:

* another real component;
* a branch;
* low density;
* a gap in the sampling;
* overlap parameters that are too strict.

## 2.5 Intersections between charts

We build a matrix:

```text
point_in_anchor(point, anchor)
```

Then we define an edge between two anchors when they share at least one point.

We then create Compressed Sparse Row (CSR) structures:

* point → edges;
* edge → points.

Finally we run Breadth First Search (BFS) to group the connected overlap regions.

Computationally, that part is well thought out: CSR avoids scanning every edge for every point.

But we need to distinguish two concepts:

### What we currently detect

We detect:

> regions where the spheres of two or more anchors overlap.

### What we still don't necessarily detect

We don't geometrically identify:

* a real intersection between manifolds;
* a bifurcation;
* two consecutive charts of the same curve;
* two different branches that happen to be close by accident.

All of these cases can produce overlapping spheres.

## 2.6 Stitching

For a point covered by a single anchor:

1. we take the anchor's center;
2. we take the first principal direction;
3. we project the point onto that line.

For a point covered by several anchors:

1. we project it separately onto each anchor's first direction;
2. we reconstruct one position per anchor;
3. we average the positions using:

[
w_a = \rho_a^4
]

The power of four makes the densest anchor dominate strongly.

This produces a less abrupt transition, but it's a heuristic rule. The exponent of four isn't derived from any proven geometric property.

---

# 3. Paper 1: Local Tangent Space Alignment — LTSA

## What problem does it try to solve?

LTSA seeks to transform high-dimensional points into low-dimensional coordinates without destroying their local structure.

Example:

* we have a sheet rolled up in 3D;
* each small region looks like a plane;
* LTSA tries to unroll it and represent it in 2D.

The paper sets out two related goals:

1. obtaining low-dimensional global coordinates;
2. reconstructing a principal manifold that passes through the center of the data.

## Main idea

The intuition is:

> Each neighborhood knows a small piece of the answer, but they all use different coordinate systems. They need to be rotated and aligned to build a global map.

It's like having many small photographs of a territory:

* each photograph correctly shows its own area;
* but they're rotated and shifted;
* LTSA computes how to arrange them to form a complete map.

## Algorithm explained simply

### Step 1: choose neighbors

For each point (x_i), (k) neighbors are selected.

### Step 2: center the neighborhood

The local mean is subtracted to study its shape regardless of absolute position.

### Step 3: local PCA/SVD

PCA is done on each neighborhood.

The first (d) directions form an approximation of the local tangent plane.

### Step 4: create local coordinates

Each neighbor is expressed within that tangent plane.

This way, the paper obtains small local coordinates for each region.

### Step 5: align all local systems

This is LTSA's distinctive part.

Global coordinates (T) are sought such that, when each neighborhood is viewed within (T), it resembles its local coordinates as closely as possible.

This becomes a global eigenvalue problem.

### Step 6: obtain the embedding

The relevant eigenvectors give the low-dimensional coordinates.

### Optional step 7: reconstruct

Once it has the global coordinates, it can learn a function that maps back from the small space to the original space.

## What it produces

It mainly produces:

* a global coordinate of dimension (d) for each point;
* optionally, a reconstruction of the manifold.

## Similarities with LoManLe

Both:

* select local neighborhoods;
* do local PCA/SVD;
* estimate tangent spaces;
* depend heavily on the neighbors correctly representing the geometry;
* try to combine local information to describe something global.

The paper explicitly acknowledges that the neighborhood size must balance density, noise, and curvature: too few neighbors can produce a poor tangent, while too many introduce bias.

## Fundamental differences

### LTSA aligns coordinates; we average positions

LTSA solves a global optimization so that all local systems are compatible.

Our stitching does:

```text
anchor 1 projection
anchor 2 projection
...
weighted average of positions
```

That doesn't explicitly align the tangent bases with each other.

### LTSA produces a global embedding

LoManLe keeps the coordinates in the original space and produces a smoothed skeleton.

### LTSA uses one neighborhood per point

We first compute information per point, but afterward we select only certain points as anchors.

### LTSA assumes a regular manifold with no self-intersections

The main formulation assumes a regular, non-self-intersecting manifold.

Our current goal includes bifurcations and intersections, so that assumption is too restrictive for us.

## What we could leverage

The most interesting thing would be to incorporate a **chart-compatibility** criterion.

Instead of joining two anchors just because they share points, we could measure:

> Maybe for stitching?

1. angle between their tangent spaces;
2. transformation needed to align their coordinates;
3. residual error after alignment;
4. consistency of the projections of the shared points.

Two charts should only be connected as a continuation of the same branch if their local coordinates can be aligned with little error.

---

# 4. Paper 2: Data Segmentation Based on Local Intrinsic Dimension — Hidalgo

## What problem does it try to solve?

This paper doesn't try to reconstruct a curve or a surface.

Its question is:

> “Do all regions of the data have the same dimension?”

For example:

* one part might look like a line;
* another like a surface;
* another like a volume.

Most methods assume a constant intrinsic dimension. This paper shows that assumption can be wrong and proposes segmenting the data by their local dimension.

## Main idea

For each point, only the distances to its first two neighbors are looked at:

* its first neighbor;
* its second neighbor.

It's defined:

[
\mu_i = \frac{r_{i,2}}{r_{i,1}}
]

The distribution of this ratio contains information about the local dimension.

Intuitively:

* on a line there are few directions where neighbors can appear;
* on a surface there are more;
* in a volume there are even more.

Because of this, the pattern of distances changes with dimension.

## Algorithm explained simply

### Step 1: find the first neighbors

For each point it's computed:

* (r_1): distance to the nearest neighbor;
* (r_2): distance to the second neighbor.

### Step 2: compute the ratio

[
\mu = r_2/r_1
]

### Step 3: assume several possible populations

The method assumes there are (K) regions or manifolds, and that each has a dimension (d_k).

Each group generates a distinct Pareto distribution of the (\mu) values.

### Step 4: favor local groupings

It's not enough to group points by having similar ratios.

A preference is also added for neighboring points to belong to the same group.

In other words:

> “If two points are close, they're more likely to be part of the same dimensional region.”

### Step 5: Gibbs sampling

The algorithm keeps alternating between:

* assigning points to groups;
* estimating the dimension of each group;
* the proportion of points in each group.

The process is repeated many times until a posterior distribution is obtained.

### Step 6: uncertain points

If a point doesn't have a sufficiently high probability of belonging to a group, it's flagged as uncertain.

The method is called **Hidalgo**. It uses (K), a number of neighbors for local homogeneity, and a homogeneity probability as parameters.

## What it produces

It produces:

* an estimated dimension for each region;
* a region label for each point;
* a membership probability;
* points with uncertain assignment.

It does not produce:

* a smoothed curve;
* a reconstructed surface;
* stitching;
* low-dimensional coordinates.

## Similarities with LoManLe

Both use:

* local information;
* neighbors;
* the possibility that distinct regions exist;
* point labels;
* a notion of intrinsic dimension.

## Differences

Our code receives `manifold_dim` as a global input.

Hidalgo allows:

```text
local dimension = 1 in one region
local dimension = 2 in another
local dimension = 3 in another
```

Our spectral gap currently answers:

> “Can I properly identify the dimension I was given?”

It doesn't answer:

> “What is the correct dimension at this point?”

## What it could give us

It could be useful before building the anchors.

For example:

```text
Estimate local dimension
        ↓
Separate 1D, 2D, etc. zones
        ↓
Build charts compatible with that dimension
```

It could also help distinguish:

* a 1D branch;
* an intersection where two directions locally appear;
* a 2D surface;
* noise that artificially inflates the dimension.

But it has an important limitation: a bifurcation of several curves can locally look like a higher dimension, even though it's really a union of several 1D manifolds. By itself, it doesn't solve the topology.

---

# 5. Paper 3: Manifold Approximation by Moving Least-Squares Projection — MMLS

## What problem does it try to solve?

MMLS wants to build a smooth surface from noisy points.

It doesn't necessarily try to unroll the data. It wants to define an operation:

[
P(r) = \text{projection of } r \text{ onto the approximated manifold}
]

The result is a smooth manifold onto which points can be projected and operations performed.

## Main idea

Let's imagine we lay a small flat sheet over a cloud of points. That sheet roughly describes the local orientation.

Then, on top of the sheet, we fit a polynomial curve or surface that better reproduces the local shape.

That is:

```text
first: local plane
then: curved correction on top of that plane
```

## Algorithm explained simply

To project a point (r):

### Step 1: find a local origin (q)

A point close to the center of the relevant region is found.

### Step 2: find a local affine space (H)

A plane of dimension (d) representing the nearby points is computed.

It's not necessarily a plain PCA centered at (r). The position (q) and the local space are jointly optimized using distance-dependent weights.

### Step 3: project neighbors onto (H)

Each neighbor gets coordinates within the local plane.

### Step 4: fit a polynomial

A polynomial function is fitted:

[
p:H\rightarrow \mathbb{R}^n
]

that reconstructs the original manifold from the local coordinates.

### Step 5: evaluate the polynomial

The final projection of (r) is the value of the polynomial at its local position.

The paper summarizes precisely these two big steps: finding a local affine space and then a polynomial approximation of the manifold over that coordinate system.

## What it produces

It produces a smooth projection operator.

It can be used for:

* denoising;
* surface approximation;
* projection of new points;
* evaluation of functions over the manifold.

## Similarities with LoManLe

This is one of the papers closest to our projection stage.

Both:

* work in the original space;
* build local approximations;
* use local PCA/SVD;
* project points;
* seek to reduce noise;
* don't necessarily need a global embedding.

## Fundamental differences

### Our local model is linear

Our projection uses:

[
c + u_1u_1^T(x-c)
]

MMLS can use a polynomial of degree 2, 3, etc.

Therefore, MMLS can represent curvature within the same chart.

### Our atlas depends on discrete anchors

MMLS can build a local approximation specific to each query point.

### MMLS optimizes the center and the plane

In our code, the tangential center is the centroid of the points in the sphere, but the sphere is defined around the anchor's original position.

In MMLS, the local origin and the affine space are an explicit part of the fitting problem.

### MMLS is designed for smooth manifolds

It's not designed primarily for bifurcations or unions of manifolds.

A bifurcation isn't locally equivalent to a single Euclidean space at the junction point, so it breaks the classical smooth-manifold assumption.

## A very relevant detail

In its experiments with a hemisphere, global PCA shows a large error, while first- and second-degree MMLS greatly improve the approximation; degree two represents the local geometry of the sphere almost exactly.

This shows a concrete limitation of our projection onto PC1:

> When the chart contains visible curvature, collapsing it onto a straight line can deform the shape.

## What we could leverage

A natural improvement for LoManLe would be to allow two types of chart:

### Linear chart

For nearly straight regions:

```text
center + PCA
```

### Polynomial chart

For curved regions:

```text
tangential coordinate t
reconstructed position = c + a₁t + a₂t²
```

We don't need to implement the full MMLS theory right away. A local quadratic fit could already reduce cuts and breaks between segments.

---

# 6. Paper 4: Learning a Manifold as an Atlas

## What problem does it try to solve?

This paper says:

> “A manifold doesn't have to unroll completely onto a single plane. Better to represent it as what it mathematically is: a collection of overlapping charts.”

This allows representing:

* closed loops;
* closed surfaces;
* periodic motions;
* manifolds that don't fit correctly into a single global coordinate system.

The method learns the charts and the assignment of points to them simultaneously.

## Main idea

Let's imagine a globe.

We can't make a perfect flat map of the whole Earth. Instead we use several maps:

* map of America;
* map of Europe;
* map of Asia;
* etc.

The maps overlap and each one describes its region well.

That's an atlas.

## Algorithm explained simply

### Step 1: propose local charts

Each chart is an affine subspace of dimension (d), similar to a local PCA.

### Step 2: assign points to charts

Each point must be assigned to the interior of some chart.

It can also belong to additional charts to generate overlaps.

### Step 3: measure the cost

A chart is good when:

* it explains the assigned points well;
* it doesn't have too much complexity;
* it preserves neighborhoods;
* it covers a coherent region.

### Step 4: alternate two optimizations

The method alternates between:

1. improving the local planes via PCA;
2. improving the discrete point assignment via graph cuts or α-expansion.

### Step 5: remove or create charts as needed

The optimization decides which charts are worth keeping.

## What it produces

It produces:

* a set of charts;
* the parameters of each chart;
* point assignments;
* overlaps;
* a dual graph between charts.

It doesn't need to unroll the whole manifold.

Points can even be classified by searching for neighbors within each chart. The paper's figure shows precisely how a query can have a correct neighbor within the chart even when the original Euclidean neighbor is misleading.

## Similarities with LoManLe

This is probably the paper conceptually closest to our LoManLe's **overall architecture**.

Both:

* represent the structure through local charts;
* allow overlap;
* use affine/PCA approximations;
* avoid relying exclusively on a global embedding;
* can represent closed structures;
* create relationships between charts through shared points.

Our set of anchors and spheres is, in effect, a simple form of atlas.

## Fundamental differences

### The paper optimizes the assignment

Our atlas is built through a greedy process:

1. sort by density;
2. take the best allowed anchor;
3. mark its points as covered;
4. repeat.

The paper alternates between chart parameters and assignments to reduce an objective function.

Because of this, a bad initial decision can be corrected in the paper. In our implementation it usually stays fixed.

### Charts can adapt more flexibly

Our charts are spheres defined by a radius.

In Atlas, the set of points assigned to a chart doesn't have to match a Euclidean sphere exactly.

### It doesn't identify all overlaps as equivalent

Our edge between anchors is created if they share at least one point.

The paper uses the assignments and topological relationships as part of a more controlled optimization.

### Our stitching collapses onto PC1

Atlas keeps a chart of dimension (d). It doesn't force everything to end up as a line.

## What we could leverage

The clearest improvement would be to include a refinement phase after the greedy step:

```text
1. Build initial anchors
2. Assign each point to the chart that best explains it
3. Recompute each chart
4. Remove redundant charts
5. Reassign points
6. Repeat until stable
```

The “chart that best explains it” could minimize:

[
\text{normal error}
+
\lambda,\text{distance to center}
+
\gamma,\text{penalty for breaking connectivity}
]

That would make our atlas much less dependent on the initial density ordering.

---

# 7. Paper 5: IAN — Iterated Adaptive Neighborhoods

## What problem does it try to solve?

IAN focuses on something that happens before almost every manifold learning algorithm:

> “Who are each point's real neighbors?”

Choosing a fixed `k` can cause two kinds of errors:

### `k` too small

* breaks the manifold;
* creates gaps;
* disconnects regions.

### `k` too large

* joins nearby branches;
* crosses concavities;
* creates short-circuits;
* fills gaps that should stay empty.

The paper proposes different neighborhoods for each point, adapted to geometry, dimension, and sampling.

## Main idea

IAN keeps two representations:

1. a discrete graph of possible neighbors;
2. a continuous Gaussian kernel with a different scale per point.

It then forces both representations to be consistent.

When a neighbor makes the local volume look strange, that connection is removed.

## Algorithm explained simply

### Step 1: build a Gabriel graph

Two points are connected if the sphere having the segment between them as diameter contains no other point.

Intuitively:

> If there's a third point clearly “between” both, the endpoints are probably not immediate neighbors.

This gives a conservative initial graph.

### Step 2: assign a scale (\sigma_i)

Each point gets its own Gaussian kernel width.

### Step 3: optimize the scales globally

A linear program looks for scales that are small but large enough to cover the necessary connections.

### Step 4: compare discrete and continuous volume

The algorithm compares:

* how many neighbors the graph suggests;
* how much volume the continuous kernel implies.

### Step 5: detect suspicious connections

If a neighbor forces the neighborhood to inflate too much, it's considered a possible geometric outlier.

### Step 6: remove the edge

The suspicious connection is pruned.

### Step 7: repeat

The scales are re-estimated and this continues until the graph no longer changes.

The algorithm is, therefore, an **iterative sparsification of the graph**.

## What it produces

It produces:

* an unweighted adaptive graph;
* a weighted kernel;
* a local scale per point.

Those results can later be fed into:

* Isomap;
* diffusion maps;
* geodesic estimation;
* dimension estimation;
* other manifold learning algorithms.

## Similarities with LoManLe

Both:

* look for adaptive neighborhoods;
* use a local scale;
* consider that a single `k` isn't enough;
* have Gaussian kernels or densities;
* want to respect changes in density.

## Differences

### LoManLe increases k until it finds a spectral gap

IAN also asks:

> “Does this neighbor respect the topology and the local volume?”

Our spectral gap can be large even when the neighborhood contains points from two branches. For example, two nearby parallel lines could still produce a clear principal direction.

### IAN removes edges

Our neighborhood is always formed by the Euclidean-nearest points.

IAN can say:

> “Even though this point is close in Euclidean distance, geometrically it shouldn't be a neighbor.”

### IAN seeks global scale consistency

Our adaptation is carried out independently for each point.

## What we could leverage

IAN is probably the most useful paper for improving our problem of incorrect connections.

We don't necessarily need to implement all of its linear programming right away. We could introduce simpler filters over our neighbors:

1. **Gabriel condition**: exclude connections whose diameter contains other points.
2. **Tangent consistency**: exclude neighbors with an excessively normal displacement relative to the chart.
3. **Mutual neighbors**: require that (i) consider (j) a neighbor and vice versa.
4. **Local scale test**: reject neighbors much farther away than the typical local scale.
5. **Reach proxy**: avoid connections that cross empty regions or concavities.

This should happen before the PCA, because an incorrect neighborhood produces:

```text
incorrect PCA
→ incorrect radius
→ incorrect anchor
→ incorrect stitching
```

---

# 8. Paper 6: Locally Defined Principal Curves and Surfaces — SCMS

## What problem does it try to solve?

This paper wants to find the ridges of the data density.

A principal curve is a line that passes through “the center” of a distribution.

The paper's idea is to define that center using:

* the gradient of the density;
* the Hessian of the density.

This unifies:

* clustering;
* principal curves;
* principal surfaces;
* manifold learning.

## Main idea

Let's imagine a mountain range.

The density of the data is like the height of the terrain.

* mean shift takes a point up to the top of a mountain;
* SCMS doesn't always take it to a single peak;
* it takes it to the ridge of the mountain range.

On a principal curve, we want to move toward the center only in the directions normal to the curve, but not along the curve.

## Algorithm explained simply

### Step 1: estimate a density

KDE or a mixture of Gaussians is used.

### Step 2: compute gradient and Hessian

The gradient indicates where the density increases.

The Hessian indicates how the density curves and what its principal directions are.

### Step 3: separate tangential and normal directions

The Hessian's directions let us determine:

* where the ridge extends;
* what the directions perpendicular to it are.

### Step 4: constrained movement

A mean-shift-type step is computed.

But only the part normal to the principal surface is kept.

### Step 5: iterate

The point moves repeatedly until the normal component of the gradient is practically zero.

At that point it has reached the ridge.

## What it produces

It produces:

* principal curves;
* principal surfaces;
* projections of points onto density ridges;
* potentially, principal graphs with junctions.

## Similarities with LoManLe

This is the paper closest to the **output we're looking for** when we want a central skeleton.

Both:

* make points move toward the center of the cloud;
* use local directional information;
* depend on density;
* produce a curve or central structure;
* can project many points onto a lower-dimensional structure.

## Fundamental differences

### SCMS moves points by gradient

LoManLe projects directly onto PCA lines defined by anchors.

### SCMS uses density as a mathematical object

Our density is used mainly to:

* order anchors;
* weight the stitching.

We don't compute the gradient or Hessian of that density.

### SCMS defines the curve as a ridge

Our skeleton is defined through the union of local projections.

### SCMS can better tolerate bifurcations and intersections

The paper notes that density ridges can intersect, so its definition can produce principal graphs without extra rules for each case.

This doesn't mean every bifurcation automatically comes out perfect: the KDE and its bandwidth matter a lot. But conceptually it's better prepared than a single smooth curve.

## Important limitation

KDE becomes hard in high dimensions because of the curse of dimensionality. The paper itself makes clear it doesn't intend to replace graph-based methods for very high-dimensional spaces.

## What we could leverage

We don't necessarily need to implement full KDE and Hessians.

A version compatible with LoManLe could be:

1. use our local neighbors;
2. compute a local weighted mean;
3. get the displacement toward that mean;
4. remove the tangential part;
5. move the point only in the normal direction;
6. repeat.

It would be a kind of **tangent-restricted mean shift**:

[
x_{\text{new}}
================

x + (I-UU^T)(\mu_w-x)
]

where:

* (U) contains the tangential directions;
* (\mu_w) is the local weighted mean;
* (I-UU^T) keeps only the normal movement.

This closely resembles the spirit of SCMS and could be integrated with our anchors.

---

# 9. Global comparison

| Method         | Main goal                                   | Neighborhoods             | Local model                      | Global joining                        | Output                       |
| -------------- | -------------------------------------------- | ------------------------- | --------------------------------- | -------------------------------------- | ----------------------------- |
| LTSA           | Unroll the manifold                          | k-NN                      | PCA tangent                       | Global alignment via eigenproblem      | Low-dimensional coordinates   |
| Hidalgo        | Separate regions by dimension                | First neighbors           | Statistical distance model        | Bayesian mixture + homogeneity         | Labels and dimensions         |
| MMLS           | Approximate and project onto a smooth manifold | Weighted neighborhood     | Plane + polynomial                | Doesn't need discrete stitching        | Projection operator           |
| Atlas          | Represent via overlapping charts             | Local graph                | Affine subspaces                  | Optimization of charts and assignments | Atlas                         |
| IAN            | Find correct neighbors                       | Adaptive Gabriel graph     | Multiscale kernel                 | Discrete/continuous consistency        | Graph and scales              |
| SCMS           | Find density ridges                          | Local kernel               | Gradient + Hessian                | Convergence toward ridges              | Principal curve/surface       |
| Current LoManLe | Build local skeleton                        | Gap-adaptive k-NN          | Local PCA                         | Overlap + averaging                    | Skeleton in original space    |

---

# 10. Which are closest to each part of LoManLe

## Neighborhood selection

The most relevant paper is:

**IAN**

Because it directly attacks the risk of joining points that are close in Euclidean distance but don't belong to the same geometric neighborhood.

## Tangent estimation

The closest one is:

**LTSA**

Both use local PCA/SVD to approximate tangents.

## Anchor and chart construction

The closest one is:

**Learning a Manifold as an Atlas**

Our anchor-sphere procedure is a greedy, geometric version of an atlas.

## Projection and denoising

The closest one is:

**MMLS**

Especially if we want to replace local lines with curved models.

## Central skeleton

The closest one is:

**SCMS**

Because it formally defines a central density structure.

## Variable dimension

The relevant one is:

**Hidalgo**

Because our code currently assumes a predefined global dimension.

---

# 11. Which parts of LoManLe look like our own

The full combination we have doesn't correspond exactly to any of the papers.

In particular, our code combines:

* adaptive `k` selection via spectral gap;
* density normalized by radius;
* anchors chosen by density;
* coverage controlled by overlap proportion;
* graph of intersections between spheres;
* BFS via CSR;
* stitching weighted by density raised to the fourth power.

Each ingredient has relatives in the literature, but the specific sequence and the stitching criterion are our own implementation's.

---

# 12. Concrete problems we see in the current code

## 12.1 Sphere overlap isn't the same as geometric intersection

An edge is created by:

```fortran
any(point_in_anchor(:,i) .and. point_in_anchor(:,j))
```

This includes:

* consecutive charts;
* redundant charts;
* nearby branches;
* real crossings;
* accidental overlaps.

Before joining them, they need to be classified.

## 12.2 Missing tangential compatibility

Two anchors should only be considered a continuation of the same branch when:

* their tangents are compatible;
* the line between their centers is aligned with the tangents;
* the projections of the shared points are coherent.

## 12.3 Averaging can cut through curves

Averaging two local lines can produce a position halfway between both that doesn't belong to either one.

This explains why, in some configurations, joints appear to get cut or incorrectly cross a bifurcation.

## 12.4 `density**4` can create abrupt changes

Even though it's described as a smooth transition, such a high power makes small density differences produce strong dominance.

When the dominant anchor switches from one to another, the output can change abruptly.

## 12.5 Coverage center and projection center aren't the same

A sphere's membership is computed around:

```fortran
work_coords(:, anchor)
```

but the projection uses:

```fortran
anchor_centers(:, anchor)
```

which is the recomputed centroid.

If both centers drift apart, the chart used to decide membership doesn't exactly match the chart used for projection.

## 12.6 The orphans pass can grow a sphere too much

Growing the nearest anchor until it covers an orphan guarantees coverage, but it can:

* cross a gap;
* touch another branch;
* alter the topology;
* create new artificial overlaps.

The PCA is recomputed afterward, which helps, but doesn't eliminate the geometric problem.

---

# 13. What architecture we'd recommend for the next version

We wouldn't replace all of LoManLe. Our base structure is good. The best evolution would be:

*Note: this section was originally written as a roadmap of future ideas. After that version, `lomanle.F90` went through a backbone/edge-construction rewrite, so each phase now includes how much of the original idea actually held up.*

## Phase 1: geometrically safe neighborhoods

Take ideas from IAN:

```text
initial k-NN
→ mutual-neighbor filter
→ Gabriel filter
→ tangent residual filter
→ final neighborhood
```

**What was actually implemented:** not the Gabriel graph or IAN's mutual-neighbor filter. Instead, the adaptive neighborhood growth (Steps 1-3) was rewritten to score each candidate size with a *quality score* (capped spectral gap + tangential stability − normalized normal error − normalized radius) and keep the best size seen so far, not the last one evaluated. Growth stops as soon as the tangent basis becomes unstable for several consecutive steps (`patience`, which self-calibrates based on how noisy each point's smallest neighborhood already looks). This targets the same problem Phase 1 was worried about — growing `k` silently crossing into another branch — but through a different mechanism: tangential stability and normal error built directly into the growth loop, instead of a separate graph filter.

## Phase 2: charts

Keep our anchors, but refine them like Atlas:

```text
initial greedy selection
→ point reassignment
→ recompute center and tangent
→ remove redundant charts
→ repeat
```

**What was actually implemented:** partially, and in a different way than proposed. Anchor selection (Step 5) went from a single greedy pass to a loop that keeps adding anchors until every point is covered, and the `o_min`/`o_max` constraint is now checked *before* accepting an anchor instead of after — so a redundant anchor (whose sphere would already be almost entirely covered) simply never gets accepted, and the old separate "Redundancy Cleanup" pass disappeared entirely, absorbed into the acceptance rule itself. After the atlas is built, each anchor's center and tangent basis are indeed recomputed over its final membership (Step 6), as this phase proposed. What's missing: there's no pass that reassigns points to a different anchor after the fact, and no removal of redundant charts beyond the upfront constraint — it's still a lighter version than Atlas's full alternating optimization.

## Phase 3: classify relationships between charts

For each overlapping pair, compute:

### A. Tangential angle

[
\theta_{ij}
===========

\arccos(|u_i^Tu_j|)
]

### B. Center-center alignment

[
a_i=
\frac{|u_i^T(c_j-c_i)|}{|c_j-c_i|}
]

### C. Normal error of shared points

### D. Density and overlap size

Then classify:

```text
continuation
bifurcation
crossing
accidental overlap
redundancy
```

**What was actually implemented:** this explicit classifier doesn't exist. An edge between two anchors is still created simply because they share at least one point (Step 8) — the same limitation already flagged in section 12.2. Something related did change though: the `density**4` weight in stitching (section 12.5) was replaced by inverse variance (`1/normal_error` of each anchor, see Phase 4), so an anchor that fits its own neighborhood poorly carries less weight regardless of how dense it is. It's not the relationship classification this phase proposed, but it addresses the same underlying concern (that density alone shouldn't dominate stitching).

## Phase 4: topological stitching

Don't automatically average all the anchors covering a point.

First determine which branch it belongs to.

In an ordinary region:

```text
blend of compatible charts
```

In a bifurcation:

```text
keep several branches
+ create a junction node
```

In a crossing with no connection:

```text
don't join the branches
```

**What was actually implemented:** in a fairly different way than proposed, but aiming at the same goal. A new subroutine, `build_skeleton_edges`, was added that builds a backbone graph over the ANCHORS (not over individual points): two anchors are edge candidates when some point already uses them as its primary+secondary anchor (the same pair that already drives Step 10's stitching). Over that candidate graph a minimum spanning tree is computed (Kruskal + union-find), which automatically discards redundant or spurious overlaps and guarantees an acyclic backbone with no separate loop detection needed. Each anchor's degree in that MST classifies it directly: endpoint (degree ≤1), pass-through (degree 2), or bifurcation/junction (degree ≥3) — the same idea of "creating a junction node" at a bifurcation, only now it emerges from the graph instead of being inferred heuristically. Every branch of the tree is walked and each point gets one continuous coordinate along the whole branch (instead of resetting to zero at every chart boundary), so points near a boundary between anchors end up correctly interleaved instead of producing duplicated or zig-zagging fragments.

What's still not implemented: the explicit continuation/bifurcation/crossing/accidental/redundancy classification from Phase 3 doesn't exist, so the MST doesn't distinguish a "crossing with no connection" from a real bifurcation — it simply connects whatever the candidate graph allows.

**Something that wasn't in any phase of this plan:** the whole pipeline (KD-tree, densities, atlas, stitching) is now wrapped in an outer convergence loop (`lomanle_compute` repeatedly calling `lomanle_pass`) that iterates over the resulting coordinates until the maximum displacement of any point falls below a tolerance relative to the dataset (a fraction of the median nearest-neighbor distance), instead of running just once.
