# Introduction

This file contains the instructions to reproduce the results of the smoothing experiments.

So far, the following methods have been tested:

- Loess (Location: `src/tox_loess.F90`)
- Nadaraya–Watson (Location: `src/knn_smoothing_nadaraya_watson.F90`)
- ANWIL (Location: `src/anwil.F90`)
- ManLe (Location: `src/manle_module.F90`)
- AManLe (Location: `src/manle_module.F90`)

At the moment, both ManLe and AManLe are implemented in the same module.

*Important: Since experiments are still ongoing, there may be allocations inside the subroutines. This will be fixed in the future.* 

Please refer to the method details in [Fortran_Coding_Guides.pdf](https://gitlab.rlp.net/a.hallab/tensor-omics/-/blob/Smoothing-Algorithm-Descriptions/misc/Tensor_Omics_Methods.pdf?ref_type=heads).

Now we're testing another approach called LoManle. This tests are independent from the previous ones so please read carefully how to run these tests.

# Project compilation

To reproduce the tests or modify the code and try it out, make sure you always compile the project.  
Apply your changes inside the `src/` folder and compile as usual using:

```bash
./build.sh
````

To run the test suite for all modules in the project, execute:

```bash
./test_runner.sh
```

# Data generation

We generate datasets to experiment with and compare the methods using a Python script. This script generates CSV files in `results/data/2d` for different datasets with low, medium, and high noise levels. You can find the script at:

```
python/generate_smoothing_datasets.py
```

You can also find a script to plot those datasets in:

```
python/plot_datasets.py
```

If you want to generate datasets with more dimensions, you need to modify the `DIM` variable in both scripts.

# Smoothing-specific tests

To generate the smoothing-specific tests, we use a Fortran program that reads the previously generated datasets, applies smoothing using all methods, and writes CSV files for each dataset containing all relevant output from each method, with a header of the form:

```
x_original,y_original,x_loess,y_loess,x_anwil,y_anwil,x_anwil_iterative,y_anwil_iterative,x_nw,y_nw,y_nw_knn,x_manle,y_manle,x_manle_svd,y_manle_svd,x_amanle,y_amanle
```

The source file is located at:

```
test_aux/smooth_all_methods.f90
```

## Automatic compiling

### Loess, Anwil, Nadaraya Watson, Manle, Amanle

To generate smoothing-specific tests, we use a Fortran program that applies various smoothing methods (Loess, ANWIL, NW, ManLe, AManLe) to the generated datasets.

Use the script `./run_smoothing_tests.sh` to compile the project, run the smoothing across a grid of hyperparameters, and generate comparison plots. Now we are generating plots for ANWIL, Manle, AMANLE and local sigma smoothing only. Please modify the script to generate all plots or the local sigma plots separatly (uncomment the desired script).

Please consider that for the latest experiments we used complete nadaraya watson for local sigma smoothing. If you want to use knn, please modify the call `smooth_vectors_gaussian_adaptive_nw` in `anwil` module to use `smooth_type` 2.

**Usage:**
```bash
./run_smoothing_tests.sh <all|dataset_name> <k_list> <span_list> <iters_list> <kernel_list> <sigma_k_list> <method_id> <score_list> <w_rough> <w_rmse> <w_cov>
```

**Example**
```bash
./run_smoothing_tests.sh results/data/2d/circular_arc_noise_medium.csv 30 0.30 10 1 30 0 1 3 3 4
```

**Parameters:**
* `all|dataset_name`: If `all` is used, all datasets in `results/data/2d/`  directory are going to generate the smoothed results. You can also only generate the results for one single dataset passing the name of the file.
* `k_list`: Comma-separated list of neighbors (e.g., `10,20`).
* `span_list`: Comma-separated list of spans for Loess (e.g., `0.30,0.50`).
* `iters_list`: Comma-separated list of max iterations.
* `kernel_list`: `1` for Gaussian, `2` for Tricubic.
* `sigma_k_list`: Neighbors used for local sigma smoothing in ANWIL.
* `method_id`: ID of the method to execute (0=All, 1=Loess, ..., 7=AManLe). See below.
* `score_list`: Calculate smoothing score using 1 = arithmetic mean, 2 = geometric mean
* `w_rough`: Comma-separated weights for Roughness
* `w_rmse`: Comma-separated weights for RMSE
* `w_cov`: Comma-separated weights for Coverage.


Where `method_id`:
- 0 to execute all methods
- 1 to execute loess
- 2 to execute isotropic anwil 
- 5 to execute nadaraya watson
- 6 to execute manle
- 7 to execute amanle 



The smoothed results are located under `results/data/2d/`.

The plots are located in `results/plots/` with the corresponding parameters that you used in the call.

You can find some results experiments in `results/plots/smoothing/`

### Automated Manifold Learning Pipeline (`LoManLe`)

# LoManLe Stitching: Current Implementation and Open Questions

## 1. General Idea

LoManLe is designed as a local manifold learning approach. Instead of fitting one global manifold, we build a local atlas composed of several local tangent approximations.

Each anchor defines:

* a local sphere of influence,
* a local density estimate,
* a local tangent direction,
* and a set of points covered by that sphere.

The stitching step is needed because neighboring anchor spheres overlap. Points inside these overlaps can be projected onto more than one local manifold, so we need a rule to combine those local projections into a single updated position.

### Provenance

The overall design -- a local atlas of overlapping tangent-plane anchors, stitched together where charts overlap, rather than fitting one global embedding -- draws on ideas from the local/piecewise manifold-learning literature, IAN in particular, alongside more general local-PCA / local-tangent-space approaches. LoManLe is not a direct reimplementation of any single published method; the adaptive-radius growth (§4), the inverse-variance stitching rule (§15), and the anchor-graph backbone construction (§16) are our own design choices, arrived at through the iteration documented throughout this file.

### Testing so far

The algorithm has been run on both 2D and 3D ambient data. Datasets are organized by ambient dimension under `results/data/2d/` and `results/data/3d/` (bifurcation and curve/surface shapes, at low/medium/high noise levels). Rendered results for both live in `results/plots/`. See "How the plots are generated" under §18 for exactly what each report contains and why 3D datasets get both.

### Status: experimental, and the API reflects that

`lomanle_compute_alloc` currently returns many intermediate/diagnostic work arrays directly to the caller -- `densities`, `gap_values`, `k_selected`, `stability_values`, `growth_stopped_complex`, both the iteration-1 and the final snapshot of nearly every quantity, and more -- well beyond what a "just give me the skeleton" API would expose. This is intentional for now: we are still actively exploring what these intermediate quantities look like across datasets, and having them all available made every plot in §16 possible without re-running anything. It does mean the public interface is wider than a finished module's would be, and it will need to be trimmed down once the algorithm design has settled.

### Last updated

The Fortran source (`src/lomanle.F90`) was last revised the week of 2026-07-20, to bring it in line with the project's coding guidelines (`misc/Fortran_Coding_Guides.tex`): consistent naming conventions, the Helper/Main/Alloc allocation pattern, `do concurrent` parallelization for the independent per-point/per-anchor loops, and `tox_errors`-based input validation on the public entry point. 

---

## 2. Iterative Structure

The public routine `lomanle_compute` runs `lomanle_pass` iteratively.

At each iteration:

1. The algorithm receives the current coordinates.
2. It builds the local atlas again.
3. It performs stitching.
4. The stitched coordinates become the input for the next iteration.
5. The maximum displacement of any point is computed.
6. If the maximum displacement is below `conv_tol`, the algorithm stops.

```text
initial coordinates
        ↓
LoManLe pass 1
        ↓
stitched coordinates
        ↓
LoManLe pass 2
        ↓
stitched coordinates
        ↓
...
        ↓
converged skeleton
```

The current defaults (set in `test_aux/test_lomanle.f90`, not hardcoded inside `lomanle_compute` itself) are:

```fortran
max_iterations    = 50
relative_conv_tol = 0.01
```

The stopping criterion is based on the largest point movement between two consecutive iterations:

```text
max_disp = max distance between old point position and stitched position
```

This means the algorithm stops when the local manifold estimates and the stitched coordinates become stable.

### `conv_tol` is relative to the dataset, not an absolute number

`relative_conv_tol` (e.g. `0.01`) is not itself the tolerance used in the `max_disp < conv_tol` check. `lomanle_compute` first computes the median nearest-neighbor distance in the *original* input coordinates, and only then derives:

```text
conv_tol = relative_conv_tol * median_nearest_neighbor_distance
```

The reason: an absolute tolerance like `1.0e-3` only means something relative to whatever units and scale a particular dataset happens to use — the same `1.0e-3` could be far too loose for data living in `[-1, 1]` and far too strict for data living in the thousands. Tying `conv_tol` to a fraction of the data's own typical point spacing instead makes it self-calibrating: the same `relative_conv_tol = 0.01` behaves sensibly whether the input coordinates are tiny or huge, without needing to be re-tuned per dataset.

---

## 3. Step 0: KD-tree Construction

The first step builds a KD-tree over the current coordinates.

The KD-tree is used to query local nearest neighbors efficiently. This is important because LoManLe repeatedly needs local neighborhoods around each point.

Without the KD-tree, nearest-neighbor searches would require scanning all points for every point. The KD-tree reduces this cost and makes the local neighborhood construction more efficient.

---

## 4. Steps 1–3: Adaptive Local Radius and Tangent Estimation

For every point, the algorithm starts with `k_min` neighbors and grows the neighborhood adaptively, in steps of `k_curr = nint(k_curr * 1.25)`, up to a hard cap `k_limit = n_points / 4`. This part of the description above is still accurate. What changed is *how the algorithm decides a neighborhood is good enough to stop growing, and which one it keeps* — the original rule (grow until the spectral gap crosses a threshold, or the size cap is hit) turned out to have real failure modes. The current implementation is a direct response to those.

### What was wrong with "grow until the gap is big enough"

* **A bad gap always meant "add more points"**, but a bad spectral gap can also mean the neighborhood just crossed into a different branch, an intersection, or a region of high curvature — cases where growing further makes the estimate *worse*, not better. A single criterion (the gap) could not tell these two situations apart.
* **The k-d tree query included the point itself** (at distance 0) and returned neighbors unsorted, both of which make radius/jump comparisons between consecutive growth steps meaningless unless corrected first.
* **If the gap never got good enough, the code accepted whatever neighborhood it had grown to last** (the largest one) — which, precisely because it never satisfied the gap threshold, was plausibly the worst neighborhood evaluated, not the best.

### The current adaptive-growth rule

For every growth step (`k_curr`), the algorithm now still queries neighbors and computes the spectral gap as before, but also:

1. Excludes the query point itself and sorts the returned neighbors by distance, so the radius and any distance comparisons are meaningful.
2. Computes a **local scale** (median distance among the first `k_min` neighbors) and a **normal reconstruction error** (mean squared residual of the neighbors off the tangent subspace through the local center).
3. Computes **tangent stability**: how similar the current step's tangent basis is to the previous step's (first compared against the previous *outer iteration*'s tangent, then step-to-step). For `manifold_dim == 1` this is a plain dot product; for `manifold_dim > 1`, the smallest singular value of the two bases' cross product (the largest principal angle between the subspaces).
4. Combines all of this into a single **quality score**:

```text
quality = min(gap / g_threshold, 2.0)
        + stability
        - normal_error / local_scale^2
        - radius / local_scale
```

5. **Keeps the best-scoring neighborhood seen so far**, not the last one evaluated — so if the gap never crosses the threshold, the point still gets the neighborhood that best balanced gap, stability, fit, and radius, instead of whatever the largest attempted neighborhood happened to be.
6. Stops growing once stability drops below `stability_threshold` for `patience` consecutive steps in a row (a single bad reading can just be small-sample noise settling down; a real branch or curvature transition keeps failing as the neighborhood grows further into it) — **or** once the gap is already sufficient and stability is holding.
7. `patience` is not a fixed constant: it is derived per point from `noise_ratio = sqrt(normal_error) / local_scale`. Clean, curve-like data gives `noise_ratio ≈ 0` and `patience = 1` (react immediately, same as reacting to the first bad reading); noisy, blob-like data gives a larger `noise_ratio`, and `patience` grows towards a small fixed cap, giving the small-sample tangent estimate room to settle before a bad reading is trusted.
8. A `scale_factor` safety cap (`sigma_i > scale_factor * local_scale`) stops growth outright if the radius has grown far beyond the point's own local scale, regardless of gap or stability — a backstop against runaway growth.

### Why this matters downstream

`sphere_radii`, `densities`, and the anchors selected from them (Step 5) are only as good as the neighborhoods they come from. A neighborhood that silently crossed into a different branch produces an oversized sphere, which produces an anchor that covers more than one branch, which produces spurious intersections and mixed-branch stitching. Reacting to *instability*, not just to a *insufficient gap*, and keeping the best neighborhood seen (not the largest) both exist specifically to keep that error from propagating downstream into the atlas and the stitching step.

---

## 5. Density Estimation

For each local sphere, the algorithm estimates density using a Gaussian-like local weighting:

```text
rho_i = sum exp(-d_j² / (2 sigma_i²))
density_i = rho_i / radius_i^manifold_dim
```

Intuitively:

* points with many close neighbors get higher density,
* large-radius sparse neighborhoods get lower density.

This density is later used both for anchor selection and for stitching.

---

## 6. Step 4: Sorting Points by Density

After computing densities, the algorithm sorts all points from highest to lowest density.

The goal is to construct the atlas starting from the most reliable local regions.

High-density regions are preferred because their tangent estimates are usually more stable.

---

## 7. Step 5: Iterative Atlas Construction

The algorithm now selects anchors.

An anchor is selected if:

* it is not already covered,
* it has high density,
* and its sphere has an acceptable overlap with the current atlas.

The overlap is controlled by:

```text
o_min <= overlap_ratio <= o_max
```

This prevents two problematic cases:

1. **Too little overlap**
   The new sphere would be disconnected from the existing atlas.

2. **Too much overlap**
   The new sphere would be redundant.

If no valid connected candidate is found, the code allows a new seed. This is useful for disconnected components or possible bifurcation-like structures.

After selecting an anchor, all points inside its sphere are marked as covered.

The loop continues until all points are covered or no valid candidate remains.

---

## 8. Step 5b: Orphans Pass

After atlas construction, there may still be a few uncovered points.

These are called orphans.

For each orphan, the algorithm finds the nearest existing anchor and grows that anchor radius just enough to include the orphan.

This is a final cleanup step to avoid leaving isolated uncovered points.

---

## 9. Step 6: Computing SVD for Final Anchor Spheres

Implemented in `compute_anchor_svd` (`src/lomanle.F90`). After the final anchor set and radii are known, the algorithm recomputes the local tangent information only for the selected anchors -- "recomputes" in the algorithmic sense (Steps 1-3 already computed something similar per candidate neighborhood; this step redoes it once, precisely, over each anchor's final sphere). 

This matters because some radii may have changed during atlas construction or during the orphan pass.

For each anchor:

1. Collect all points inside the final sphere.
2. Recompute the local center.
3. Recompute the covariance matrix.
4. Recompute the tangent direction.
5. Store the anchor center and tangent scale.
6. Recompute the anchor's own `normal_error`: the mean squared distance of every point in the *full, final* sphere from this tangent subspace — a broader, more reliable fit-quality measure than the `normal_error` computed during the narrower Steps 1–3 neighborhood search.

These final anchor manifolds are the ones used during stitching.

This anchor-level `normal_error` is not just a diagnostic: it is exactly the quantity Step 10's inverse-variance weighting (`weight_k = 1 / normal_error_k`) uses to decide how much each covering anchor should pull a point's stitched position. The recomputed anchor center (`anchor_centers`) is also returned as an explicit output of `lomanle_pass`/`lomanle_compute` (it used to be a purely internal, local variable) — it was exposed specifically because Step 11's backbone construction (§16) needs a stable anchor position to build the anchor-to-anchor graph and to accumulate distance along a branch.

---

## 10. Step 6.5: Anchor Mapping

The code creates mappings between:

```text
anchor index → original point index
original point index → anchor index
```

This is needed because not every point is an anchor.

For example:

```text
anchor 1 may correspond to original point 17
anchor 2 may correspond to original point 42
```

The mapping allows the algorithm to work with a compact anchor list while still preserving the original point identities.

---

## 11. Step 7: Membership Matrix

The membership matrix is one of the central structures for stitching.

```fortran
point_in_anchor(point, anchor)
```

It stores whether a point lies inside a given anchor sphere.

Example:

```text
              Anchor 1   Anchor 2   Anchor 3
Point 1          T          F          F
Point 2          T          T          F
Point 3          F          T          T
Point 4          F          F          T
```

This tells us:

* Point 1 belongs only to Anchor 1.
* Point 2 is in the intersection of Anchor 1 and Anchor 2.
* Point 3 is in the intersection of Anchor 2 and Anchor 3.

Points belonging to more than one anchor are the main stitching candidates.

The code also computes:

```fortran
anchor_count(i)
```

which stores how many anchor spheres contain point `i`.

---

## 12. Step 8: Building Pairwise Anchor Intersections

The algorithm then identifies pairs of anchors that overlap.

An edge is created between two anchors if they share at least one point:

```text
edge = Anchor A intersects Anchor B
```

Formally, two anchors are connected if:

```fortran
any(point_in_anchor(:, A) .and. point_in_anchor(:, B))
```

This creates a graph of anchor intersections.

Example:

```text
Anchor 1 ----- Anchor 2 ----- Anchor 3
```

means:

* Anchor 1 overlaps with Anchor 2,
* Anchor 2 overlaps with Anchor 3.

This graph is the basis for detecting larger multi-sphere intersection regions.

---

## 13. Why CSR Is Used

CSR means Compressed Sparse Row.

In this implementation, CSR is used to store sparse relationships efficiently.

The code builds two CSR structures:

### Point → Edges

For each point, we can quickly ask:

```text
Which anchor-intersection edges does this point belong to?
```

Example:

```text
Point 10 → Edge 1, Edge 4
Point 11 → Edge 1
Point 12 → Edge 4, Edge 5
```

### Edge → Points

For each edge, we can quickly ask:

```text
Which points are shared by the two anchors of this edge?
```

Example:

```text
Edge 1 → Point 10, Point 11
Edge 4 → Point 10, Point 12
```

The reason for using CSR is efficiency.

Without CSR, the BFS would repeatedly scan many irrelevant points or edges. With CSR, each point only visits its own edges, and each edge only expands to its own points.

The intended cost becomes closer to:

```text
O(number of edges + number of intersection points)
```

instead of repeatedly scanning:

```text
O(number of points × number of edges)
```

---

## 14. Step 9: BFS for Intersection Regions

BFS means Breadth First Search.

Here, BFS is used to grow connected intersection regions.

The issue describes this as:

```text
Grow or union intersections that overlap.
```

The BFS implementation does exactly this.

It starts from a point that belongs to at least one intersection edge. Then:

1. Visit the point.
2. Find all intersection edges containing that point.
3. For each edge, find all points in that edge.
4. Add those points to the queue.
5. Continue until no connected points or edges remain.

All points reached in the same BFS traversal receive the same label.

This produces connected intersection components:

```text
label 1 → one unified intersection region
label 2 → another unified intersection region
label 3 → another one
```

This is how the implementation handles multi-sphere intersections.

Instead of treating each pairwise overlap independently, connected overlaps are grouped into one unified region.

---

## 15. Step 10: Density-Weighted Stitching

After identifying the local intersections, the algorithm updates point positions.

There are three cases.

---

### Case 1: Point Belongs to Multiple Anchors

If:

```fortran
anchor_count(i) >= 2
```

then the point lies in an overlap region.

The algorithm projects the point onto each local manifold whose sphere contains it.

For each covering anchor:

1. Compute the difference from the anchor center:

```text
p_diff = point - anchor_center
```

2. Project onto the local tangent direction:

```text
proj_val = dot(p_diff, tangent_basis)
```

3. Reconstruct the projected point:

```text
projection = anchor_center + proj_val * tangent_basis
```

This gives several possible projected positions for the same original point, one from each local anchor manifold.

Then the final stitched position is computed as a weighted average of these projections.

---

### Original Weighting Idea

Following the issue, the first approach used density-dependent weights:

```text
stitched_point = Σ density_k * projection_k / Σ density_k
```

This follows the idea that denser local manifolds should contribute more strongly to the final position.

However, in practice, this sometimes produced several visible lines instead of collapsing the overlap into one clean local trend.

This happened because anchors with moderately similar densities still contributed similar weights, so the final stitched structure could preserve multiple competing local projections.

An intermediate version tried making the densest sphere dominate more sharply by raising density to a high power (`density_k^4`), which helped but was still not quite right: a small density advantage between two nearby anchors doesn't necessarily mean one of them actually fits the local geometry better -- density alone says nothing about the quality of the tangent-plane fit.

---

### Current Weighting: Inverse-Variance by Fit Quality

The current code (`lomanle.F90`, Step 10) weights each covering anchor by the inverse of its **own reconstruction error**, not by density:

```text
weight_k = 1 / normal_error_k
```

where `normal_error_k` is the anchor's own mean squared residual (Step 6: how far the points in its sphere sit off its recomputed tangent plane). The stitched position is:

```text
stitched_point =
    Σ (projection_k / normal_error_k)
    ---------------------------------
        Σ (1 / normal_error_k)
```

The idea: trust the anchor whose tangent plane genuinely fits its neighborhood, regardless of how many nearby points it happens to have. A dense anchor whose tangent plane is a poor fit (e.g. because its sphere leans into a curve or a nearby second branch) should not dominate just because it is dense; a sparser but cleanly-fit anchor should. This also collapses ordinary overlap points onto a single central trend while still letting genuine bifurcations, where two anchors fit about equally well, land roughly between both branches instead of being forced onto one.

As a side effect, Step 10 also keeps track of which two anchors give the **highest** and **second-highest** weight for each point (`primary_anchor_ids`, `secondary_anchor_ids`) -- not just to blend the position, but because that same pair is exactly what the backbone/edge construction (Step 11, see below) uses to build the topology graph between anchors.

---

### Case 2: Point Belongs to One Anchor

If:

```fortran
anchor_count(i) == 1
```

the point is projected only onto the local manifold of that anchor.

This collapses non-intersection points onto their local central manifold.

---

### Case 3: Point Belongs to No Anchor

If:

```fortran
anchor_count(i) == 0
```

the point is kept in its original position.

This should be rare after the orphan pass.

---

## 16. Step 11: Backbone / Skeleton Edge Construction

Stitching (Step 10) updates every point's position, but on its own it does not say which points should be drawn as connected on the final skeleton. This step answers exactly that question. It is implemented in `build_skeleton_edges_alloc`, which allocates the working buffers (anchor mapping, MST candidate lists, member chains, branch-adjacency CSR) and delegates the actual branch walk to the pure `build_skeleton_edges` (both in `src/lomanle.F90`). `lomanle_compute` calls `build_skeleton_edges_alloc` automatically, once for the iteration-1 snapshot and once for the converged result -- it is not something the test driver or the plotting code has to do. The rest of this section refers to the machinery collectively as "`build_skeleton_edges`" for readability.

### 16.1 Why not just connect nearby points?

The first implementation of this step (now removed; it lived in `test_aux/test_lomanle.f90` as `compute_backbone_edges`) had each point greedily pick its own "successor": look for the nearest tangent-aligned not-yet-visited neighbor, with several fallback stages (a relaxed version inside the same intersection label, then a distance-only bridge for any point still without an edge, then an unbounded directional cone as a last resort) -- plus an explicit walk of the whole graph before accepting each candidate edge, just to check it would not close a loop.

This worked on very simple cases, but fights the shape of the actual problem:

* **Euclidean closeness is not the same as belonging to the same branch.** Two points on nearby-but-different branches (parallel curves, or the two arms of a bifurcation) can be closer to each other than to their own branch's continuation. A rule based on distance and local tangent alignment alone cannot tell the difference.
* **A "next point" model gives every point at most one outgoing edge**, but a genuine branch point structurally needs three or more. Representing that needs extra bookkeeping bolted on afterwards (a "convergence point" special case, inferred by counting incoming edges) instead of falling naturally out of the data structure.

### 16.2 The core idea: the topology lives on the anchors, not on individual points

Step 10 already records, for every point sitting in an overlap, exactly which two anchors it blends between (`primary_anchor_ids`, `secondary_anchor_ids` -- see the weighting section above). Whenever a point uses two anchors as its primary and secondary, those two anchors are, by construction, part of the same connected piece of the manifold. That is a topological fact already sitting in memory; `build_skeleton_edges` reads it out and turns it into a graph over **anchors**, instead of re-deriving connectivity point by point.

### 16.3 Tier 1: candidate anchor-anchor edges from Step 10's own pairing

For every point with both `primary_anchor_ids(i) > 0` and `secondary_anchor_ids(i) > 0`, the pair (primary anchor, secondary anchor) becomes a candidate edge between those two anchors, weighted by the Euclidean distance between their Step 6 centers. This reuses information Step 10 already computed; no new scan over all points is needed.

### 16.4 Tier 2: geometric refinement for thin overlaps

Tier 1 alone can miss a real connection: if two neighboring anchors overlap only slightly, every point in that thin overlap zone might be dominated by a third anchor and never actually rank that specific pair as its own top-2 -- so no point ever "votes" for the edge, even though the anchors are genuinely adjacent. This showed up in practice as small, visible gaps in an otherwise-correct backbone.

Tier 2 catches this with a purely geometric, closed-form test: two anchors are candidates if their spheres (centered at their Step 6 centroids, radius `sphere_radii`) overlap, i.e. `distance(center_i, center_j) <= radius_i + radius_j` -- the same "shares at least one point" criterion Step 8 already uses for the pairwise-intersection edge list, just as a distance check instead of a scan over every point. By the triangle inequality, this is a strict superset of tier 1 (any point-based pair from tier 1 is automatically also a sphere-overlap pair), so tier 2 never removes a tier-1 edge; it only adds edges tier 1 missed.

### 16.5 Minimum spanning tree over the candidate edges

Both tiers feed the same Kruskal / union-find pass, continuing from wherever tier 1 left off: candidate edges are sorted by distance and added greedily, skipping any edge that would connect two anchors already in the same component. Two things fall out of this for free:

* **The backbone is acyclic by construction** -- no separate loop-detection code (the expensive part of the old point-level approach) is ever needed downstream.
* **Tier 2 only ever runs on pairs tier 1 left in separate components**, so it can only close a genuine gap; it never overrides or competes with a connection tier 1 already made.

### 16.6 Anchor roles fall out of the tree structure

Once the MST is built, each anchor's degree in it classifies it structurally, not heuristically:

```text
degree <= 1  -> endpoint
degree == 2  -> pass-through (mid-branch)
degree >= 3  -> branch / junction
```

This directly replaces the old point-level `gap_kind` classification, which *inferred* "convergence point" by counting how many points pointed at a given point, and "true endpoint" only after an unbounded directional search failed. Here it is simply read off the tree.

### 16.7 From anchor tree to point-level edges: three attempts, three bugs found and fixed

Turning the anchor-level tree into an actual polyline through the stitched points took three iterations; each of the first two produced a visibly wrong picture that pointed directly at the next fix.

**Attempt 1 — per-anchor local chain + single bridge point.** Every point was threaded into the chain of *both* its primary and secondary anchor (matching how Step 10 blends its position), sorted within each anchor by projection onto that anchor's own tangent axis, and consecutive anchors in the tree were joined by one edge between their mutually-closest members.

*What went wrong:* a point in the overlap of two anchors got ordered independently in **two** chains, each using a different local axis that need not agree. Visually this produced short, duplicated, overlapping fragments right at every chart boundary — confirmed by counting edges in the output CSV: for a graph that should have been a single tree, `n_edges` came out far above `n_nodes - 1`.

*Fix:* thread every point into exactly **one** chain — that of its primary anchor only. `secondary_anchor_ids` is still used for tier-1 candidates (16.3) and, further down, to smooth the coordinate at chart boundaries (attempt 2 below), but it no longer causes a point to be sorted twice.

**Attempt 2 — one chain per anchor, still ordered locally.** With duplication fixed, a new symptom appeared: a visible zig-zag/staircase texture, especially where several short-radius anchors sit close together along a curve. Each anchor's own tangent line is estimated from a small, finite, noisy neighborhood, so two adjacent anchors' lines can disagree slightly in orientation — and because *all* of one anchor's points were drawn before *any* of the next anchor's, points that were genuinely interleaved in true position near the shared boundary could not be, which is what produced the zig-zag.

*Fix:* stop giving every anchor its own local coordinate that resets to zero at the next chart. Instead, walk each maximal run of pass-through anchors between two "special" anchors (endpoints/junctions) as one **branch**, and give every point along that branch a single, monotonically increasing coordinate: each anchor's own tangent axis is oriented to point "forward" along the branch, and a running offset (plain center-to-center distance) is added on top so anchor *k*'s local coordinate continues exactly where anchor *k-1*'s left off. Sorting the *whole branch* by this one shared coordinate — instead of "all of anchor A, then all of anchor B" — lets points near a boundary interleave in their true order even when the two anchors' local axes disagree. For a point that also blends into a secondary anchor on the same branch, its coordinate is the *average* of both anchors' values rather than the primary anchor's alone — the same smoothing role Step 10 already gives that point's position, just applied to its ordering too.

**Attempt 3 — branches naturally re-threading their own endpoints.** A branch walk starts and ends at a "special" anchor; a junction, by definition, is where three or more branches meet. The first branch-based implementation pooled *every* anchor along its path — including both special ends — into that branch's own sorted chain. A junction's member points therefore got pooled and re-sorted once *per incident branch*: the same edge-count check as before (`n_edges` vs `n_nodes - 1`) showed dozens of redundant edges concentrated exactly at junctions and endpoints.

*Fix:* thread every "special" anchor's own member cluster into its own simple local chain exactly **once** (it is small, so it does not need a branch-wide coordinate), and have each branch pool *only* the interior (pass-through) anchors' points. Each end of a branch's interior chain is then bridged into its special anchor's nearest member, so a hub's cluster is referenced by every incident branch but never re-threaded by more than one.

### 16.8 What the final graph looks like

After all three fixes, on a `bifurcation_2way` test case the final-stage graph is a single connected component with `n_edges == n_nodes - 1` (a true tree, no residual cycles) -- a handful of near-zero-length edges can still appear where two independently-stitched points collapse onto almost the same position at a hub; these are harmless and invisible at plotting scale. Rendered, it shows one continuous line with a genuine Y-shaped bifurcation, matching the dataset's known structure.

### 16.9 Output format

`lomanle_compute` returns, for both the iteration-1 snapshot and the converged result: `edge_from`/`edge_to` (point-index pairs, sized by a caller-provided `max_edges`), `n_edges`, and `anchor_role` (0 = not an anchor; 1/2/3 as in 16.6). `test_lomanle.f90` writes these directly to a companion `..._lomanle_edges.csv` file (`stage,edge_id,x,y,(z,)xend,yend,(zend)`) next to the usual point-level CSV. No connectivity logic runs outside of `lomanle.F90` any more.

### 16.10 Visualization changes

`r/plot_lomanle_spheres.R` no longer builds or orders any edges itself — the old `build_edge_df`, which joined points using `next_iter1`/`next_final` columns, is gone along with those columns. It reads the companion edges CSV directly and draws it with `geom_segment(..., lineend = "round")`: the earlier default ("butt") caps left a visible notch at every join between two independently-drawn segments, which made an already-correct graph look like separate little pieces stuck together; rounding the caps removes that purely cosmetic artifact. The 3D report gained two matching panels ("Backbone Skeleton", iteration 1 and final), built the same way with `plotly` line traces.

### 16.11 Ambient dimension (`dim`) is no longer hardcoded to 2D in the test driver

`build_skeleton_edges` itself was always written generically in `dim` (unlike the old `compute_backbone_edges`, which hardcoded indices assuming `dim == 2`). `test_aux/test_lomanle.f90` can now be run with `dim = 3`: it reads all `dim` coordinate columns, and writes the extra `z` / `sk_z` / `sk_z_f` / `v{j}_z` columns the R script already expected for its 3D report.

### 16.12 The backbone is inherently 1D — reconstructing a surface/mesh for `manifold_dim > 1` is still open

`dim` (ambient dimension, handled above) and `manifold_dim` (intrinsic dimension of what we're trying to recover) are two different knobs, and they get very different treatment today:

* **Step 10's stitched point positions already respect `manifold_dim` fully.** The projection loops over `base_idx = 1, manifold_dim` (not just the first tangent direction), so with `manifold_dim = 2` on 3D input, each point is genuinely projected onto the local 2-D tangent plane at its anchor(s) and blended across overlaps the same way as in the 1-D case. So if the goal is a denoised point cloud lying on the surface, that's already what we get out of `skeleton_coords`/`work_coords`.
* **Step 11's backbone (`build_skeleton_edges`) is 1D by construction, regardless of `manifold_dim`.** The candidate-edge tiers (16.3-16.4), the MST (16.5), and the branch-threading (16.7) all produce a single ordering coordinate per anchor via `tangent_bases(:, 1, k)` -- only the first tangent direction -- because the whole design in §16.2-16.7 is "thread points into a chain along a branch." A chain is inherently one-dimensional: there is no notion of a second ordering axis, no face/triangle structure, and no code path that would let three or more points at the same "branch position" form a 2-D patch instead of colliding onto a single polyline vertex.

So today: for `manifold_dim = 2` (or higher), we can already see the surface in the *stitched points themselves*, but asking `lomanle_edges.csv` for a mesh or an outline of that surface will still just draw a line/tree threading through the point cloud -- the same kind of output as for a 1-D curve, just embedded in the surface. Getting an actual mesh (faces connecting anchors in two directions, not just a spanning tree) is unimplemented and would need new design work: e.g. connecting anchors along *two* independent tangent directions instead of one, or triangulating within/across overlapping charts, rather than reusing the MST/branch machinery as-is.

---

## 17. What Is Currently Implemented from the Issue

The current implementation covers the main structure proposed in the issue:

* build local spheres,
* detect sphere overlaps,
* represent pairwise overlaps as edges,
* identify multi-sphere intersection regions,
* use CSR to make traversal efficient,
* use BFS to union overlapping intersections,
* project points onto local manifolds,
* update points using an inverse-variance-weighted stitching rule,
* iterate until convergence,
* build a topological anchor graph (tiers 1-2 + MST) and turn it into point-level backbone edges, with structural endpoint/pass-through/branch roles.


# Open Questions

## 1. How should we connect edges after stitching? — Resolved

This used to say the current code updated point coordinates but did not construct the final graph connectivity, and that we had no definitive connectivity rule for intersection regions. That is no longer the case: **see Step 11 (§16) above** for the full implementation (`build_skeleton_edges`, `src/lomanle.F90`), returning an explicit `edge_from`/`edge_to` list plus a structural `anchor_role` (endpoint / pass-through / branch) for every anchor.


## 2. Possible Strategy for Edge Reconstruction — Superseded by what was actually built — Resolved

This section used to sketch a candidate strategy (project + sort within a single anchor; nearest-neighbor plus a shared-anchor/same-label filter inside intersection regions). That specific strategy was **not** what ended up being implemented, and for a concrete reason worth recording: nearest-neighbor-in-stitched-space, even filtered by a shared anchor or label, is still fundamentally a distance rule, and distance is exactly what fails to distinguish "same branch" from "nearby-but-different branch" (see §16.1). The approach that actually shipped instead builds the topology on the **anchor graph** (§16.2-16.6: Step 10's own primary/secondary pairing, a geometric fallback tier, and a minimum spanning tree), and only afterwards turns that anchor-level tree into point-level edges (§16.7) -- which is also why it took three iterations to get the point-level threading right rather than getting it directly from a single nearest-neighbor pass.


## 3. Bifurcations — Validated on synthetic 2way/3way datasets, in 2D and 3D

This used to say only that "preliminary experiments look promising, but more bifurcation datasets are needed." Since then:

* `bifurcation_2way` and `bifurcation_3way` synthetic datasets (2D, and 3D via `results/data/3d/`, at `noise_low`/`noise_medium`/`noise_high`) have been run end to end.
* The backbone construction (§16) was iterated on precisely *because* these datasets exposed real problems that a simple line or a non-bifurcating curve would not have: the duplicate-fragment bug (attempt 1), the zig-zag bug (attempt 2), and the hub-duplication bug (attempt 3) were all found by looking at `bifurcation_2way` output, not hypothesized in the abstract.
* On the final implementation, `anchor_role` correctly marks the branch point as a junction (degree >= 3 in the MST), and the rendered backbone shows one continuous line with a genuine Y-shaped split, matching the dataset's known 2-branch structure.
* Remaining validation gap: T-shaped structures, more than one branch point in the same dataset, and deliberately uneven-density bifurcations have not specifically been tested yet -- the fixes above should generalize (nothing in §16.2-16.7 assumes exactly two branches), but that is an assumption, not something confirmed on those specific shapes.


## 4. Surface/mesh reconstruction for `manifold_dim > 1` — Open, not started

See §16.12 for the detailed explanation. In short: Step 10's stitched point positions already respect `manifold_dim` correctly (so a `manifold_dim = 2` run on 3D data does give back a point cloud lying on the local surface), but Step 11's backbone construction (§16) only ever produces a 1-D chain/tree of point-to-point edges, no matter what `manifold_dim` is -- it was designed for curve/branch extraction, not for surface topology. There is currently no mesh, triangulation, or two-direction connectivity between anchors. If we want an actual reconstructed surface (not just the surface's points, but a mesh/outline of it), that needs new design work -- for example connecting anchors along a second tangent direction, or triangulating within/across overlapping charts -- rather than reusing the MST/branch-threading machinery from §16.3-16.7 as-is.

## 5. Choosing k

The issue mentions the problem of choosing a good value or range for `k`.

This has not been implemented yet.

Currently, the algorithm uses an adaptive local neighborhood process, but we have not implemented a systematic elbow method or score-based selection for `k`.

The issue suggests using a score based on:

* residuals,
* number of overlapping spheres,
* and possibly how much noise is learned.

This remains open.

A possible next step is to run the algorithm across a range of `k_min` values and evaluate:

```text
score(k) =
    projection residual
    + overlap penalty
    + instability penalty
```

Then choose a value near the elbow of the curve, where increasing `k` no longer improves the skeleton meaningfully but starts increasing overlap or smoothing too much.


## 6. Summary of Remaining Work

1. ~~Define how to reconstruct final edges after stitching.~~ Done -- see §16.
2. ~~Decide how intersection labels should be used for edge reconstruction.~~ Superseded -- the anchor-graph/MST approach in §16 does not use the BFS intersection labels for connectivity at all (only the plain point-to-anchor pairing and, as a fallback, sphere overlap); the BFS labels remain useful as a diagnostic/visualization grouping (see the "Intersection Regions" plots), just not as the connectivity rule.
3. Validate bifurcation behavior on more controlled datasets -- partially done (§ Open Questions 3): 2way/3way, 2D and 3D, at three noise levels. Still open: T-shapes, multiple simultaneous branch points, deliberately uneven-density bifurcations.
4. Reconstruct an actual surface/mesh for `manifold_dim > 1`, not just the stitched points -- still open, untouched (§16.12, § Open Questions 4).
5. Implement a systematic method for choosing `k`. Still open, untouched.
6. Define and test a score function for parameter selection. Still open, untouched.

The current implementation covers the full pipeline end to end: local atlas construction, overlap detection, CSR/BFS-based intersection grouping, inverse-variance-weighted point updates, and topological backbone/edge construction with structural endpoint/pass-through/branch classification (§16).

---

## 18. `run_lomanle_tests.sh`: Automated Pipeline Script

This bash script automates the compilation, execution, and visualization of the **LoManLe** (Local Manifold Learning) algorithm across multiple datasets. It handles the batch processing of CSV files and dynamically generates reports using R.

#### 1. Features

* **Automatic Compilation**: Links the Fortran source with LAPACK, BLAS, and LOESS libraries.
* **Batch Processing**: Use the `all` keyword to process every CSV in the data directory.
* **Dynamic Parameterization**: Allows testing of multiple $k$ (nearest neighbors) values in a single run.
* **Integrated Visualization**: Automatically calls R scripts to generate PDF reports after each execution.

#### 2. Usage

```bash
./run_lomanle_tests.sh <input_file|all> <k_min_list> [manifold_dim] [g_threshold_list] [o_max_list] [o_min_list] [stability_threshold_list] [scale_factor_list] [max_iterations_list] [relative_conv_tol_list] [dim]
```

##### Arguments

Every argument except `input_file` and `k_min_list` is optional and falls back to the default shown below. Every list-valued argument (anything ending in `_list`) accepts a comma-separated list (e.g. `10,20,30`); the script then runs the full grid (nested loop) over every list-valued argument, not just `k_min`.

| Argument | Position | Type | Description | Default |
| --- | --- | --- | --- | --- |
| `input_file` | 1 | String | Path to a `.csv` file, or the keyword `all` to process every non-output `.csv` under `results/data/2d/`. | Mandatory |
| `k_min_list` | 2 | List | Minimum neighborhood size to start adaptive growth from (§4). | Mandatory |
| `manifold_dim` | 3 | Integer | Intrinsic dimension to fit ($1$ = curve, $2$ = surface). | `1` |
| `g_threshold_list` | 4 | List | Spectral-gap threshold that stops neighborhood growth (§4). | `1.0` |
| `o_max_list` | 5 | List | Maximum allowed anchor-sphere overlap ratio (§7). | `0.30` |
| `o_min_list` | 6 | List | Minimum required anchor-sphere overlap ratio (§7). | `0.05` |
| `stability_threshold_list` | 7 | List | Minimum tangent-basis stability to keep growing a neighborhood (§4). | `0.90` |
| `scale_factor_list` | 8 | List | Cap on sphere radius as a multiple of the point's local scale (§4). | `2.5` |
| `max_iterations_list` | 9 | List | Outer convergence-loop iteration cap (§2). | `50` |
| `relative_conv_tol_list` | 10 | List | Convergence tolerance, as a fraction of the median nearest-neighbor distance (§2). | `0.01` |
| `dim` | 11 | Integer | Ambient dimension of the input coordinates: `2` for datasets under `results/data/2d/`, `3` for `results/data/3d/`. | `2` |

#### 3. Workflow Logic

1. **Compilation**: The script runs `./build.sh`, then compiles `test_aux/test_lomanle.f90` against the resulting `.mod`/library, linked with LAPACK, BLAS, and the LOESS libraries.
2. **Execution Loop**: For every combination in the cartesian product of all list-valued arguments (`k_min` × `g_threshold` × `o_max` × `o_min` × `stability` × `scale_factor` × `max_iterations` × `relative_conv_tol`), the script runs the Fortran binary once.
3. **File Management**: It renames the generic `lomanle_output.csv` / `lomanle_edges.csv` to unique names incorporating every parameter of that run (e.g. `circular_arc_noise_high_k30_g3.0_omax0.3_st0.90_sf2.5_mi50_ct0.01_lomanle.csv` and its `..._edges.csv` companion).
4. **Plotting**: Executes `r/plot_lomanle_spheres.R` on that pair of CSVs -- see "How the plots are generated" right below.

#### 4. How the plots are generated

A single R script, `r/plot_lomanle_spheres.R`, handles both 2D and 3D. It is called once per run, automatically, from step 4 of the workflow above -- there is no separate plotting command to run by hand. It always reads two files: the point-level CSV (`..._lomanle.csv`) and its companion backbone-edges CSV (`..._lomanle_edges.csv`, written directly by `lomanle_compute`/`build_skeleton_edges`, see §16.9).

* **2D report (always generated):** `plot_2d_report()` renders a multi-page PDF to `results/plots/<dataset>_2d.pdf` using `ggplot2`/`ggforce`/`viridis` -- one page per diagnostic: all adaptive spheres, point density, the global tangent field (PC1), the selected atlas anchors, the anchors' SVD segments, intersection regions (§11-§14), the stitched backbone at iteration 1 and at convergence (§15), the backbone edges themselves (§16) faceted by stage, and (if the `patchwork` package is installed) a side-by-side iteration-1-vs-converged-vs-displacement convergence panel. For a 3D dataset this report still gets generated, using only the `x`/`y` columns (a 2D projection) -- it is not skipped.
* **3D report (only when the CSV has a `z` column):** `plot_3d_report()` renders an interactive, self-contained HTML file to `results/plots/<dataset>_3d.html` using `plotly`, with six `scatter3d` panels covering the same ground in three dimensions (original vs. final stitched, intersection regions, convergence displacement, anchor tangent field, and the backbone skeleton at iteration 1 and at convergence).

So a 2D input dataset (e.g. anything under `results/data/2d/`) produces only the `_2d.pdf`; a 3D input dataset (`results/data/3d/`) produces both the `_2d.pdf` (as an x/y projection) and the `_3d.html` (the full 3D view). Both land in the same flat `results/plots/` directory -- there is no `results/plots/2d`/`results/plots/3d` split; the filename suffix is what tells them apart.

#### 5. Example Command

To process a specific dataset with a single parameter combination:

```bash
./run_lomanle_tests.sh results/data/2d/circular_arc_noise_high.csv 30 1 3.0 0.3 0.10 0.90 2.5 50 0.01 2
```
