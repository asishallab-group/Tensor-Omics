# Shape Truthful Clustering (STC) — Current Algorithm Draft

> **Status:** Current agreed design (research prototype).

## Motivation

STC is an exploratory clustering algorithm inspired by statistical physics and renormalization. It grows ensembles from dense seeds while monitoring the stability of macroscopic observables. Growth stops once the ensemble's effective description changes beyond a configurable threshold.

---

## Definitions

- Ambient vectors:
  \[
  \mathcal{V}=\{\mathbf{v}_1,\ldots,\mathbf{v}_N\}.
  \]

- Ensemble after iteration \(t\):
  \[
  \mathcal{E}_t\subseteq\mathcal{V}.
  \]

- Ambient density label:
  \[
  \rho_i=\sum_{j=1}^{N}\mathbf{1}\!\left(d(\mathbf{v}_i,\mathbf{v}_j)\le r\right),
  \]
  where \(r\) is a fixed neighborhood radius.

---

## Algorithm

### 1. Compute ambient density labels

Compute \(\rho_i\) for every ambient vector.

### 2. Seed selection

Choose high-density seed vectors (currently the top 15%).

Each seed initializes one ensemble.

### 3. Candidate generation

Determine all neighboring vectors at the ensemble surface.

A candidate is accepted into the candidate set iff

\[
\left|
\rho_i-\operatorname{median}(\rho_{\mathcal E_t})
\right|
\le
\alpha_{\mathrm{MAD}}
\cdot
\mathrm{MAD}_{\mathrm{ambient}},
\]

with

\[
\mathrm{MAD}_{\mathrm{ambient}}
=
\operatorname{median}
\left(
\left|
\rho_i-
\operatorname{median}(\rho)
\right|
\right).
\]

### 4. Atomic growth

Add **all** candidate vectors simultaneously

\[
\mathcal E_{t+1}
=
\mathcal E_t
\cup
\mathcal C_t.
\]

No sequential ordering of candidates is performed.

### 5. Compute observables

Current observable vector:

- arithmetic mean density
- harmonic mean density
- density heterogeneity
  \[
  H_\rho=
  \frac{\rho_{\mathrm{arith}}}
       {\rho_{\mathrm{harm}}}
  \]
- ensemble size
- number of candidates

Observable history is stored as a real matrix

\[
O(:,t).
\]

### 6. Accept / reject

Reject growth if

\[
\left|
\log_2
\frac{H_{\rho,t+1}}
     {H_{\rho,t}}
\right|
\ge
\alpha_{\mathrm{accept}},
\]

(default \(\alpha_{\mathrm{accept}}=1\), corresponding to a two-fold change),

or if the maximum number of iterations is reached.

If rejected, return the previous ensemble

\[
\mathcal E_t.
\]

Otherwise continue growing.

---

## API

```text
seed_ensembles()

    ↓

grow_ensemble()

    ↓

compute_ensemble_observable()

    ↓

accept_ensemble()
```

`accept_ensemble()` receives only the observable trajectory and user parameters.

---

## Output

Return

- logical membership matrix
- observable trajectories
- optional merged ensembles (post-processing)

---

# Future directions

Planned observables include

- entropy
- energy
- angular coherence
- manifold descriptors
- semantic descriptors

Each iteration yields an observable vector

\[
\Omega_t
=
(\rho,
H_\rho,
S,
E,
\ldots).
\]

The trajectory

\[
\Omega_1,\ldots,\Omega_T
\]

forms a point cloud in observable space.

Future STC versions may analyse this trajectory using

- covariance matrices,
- bootstrap confidence regions,
- Mahalanobis distances,
- Hotelling's \(T^2\),
- manifold learning,

to discover scale-invariant and discriminative observables.
