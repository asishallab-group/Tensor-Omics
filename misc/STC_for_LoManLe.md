# Shape Truthful Clustering for Local Manifold Learning (LoManLe)

## Overview

This document describes the current STC formulation for constructing LoManLe atlas charts.

Given a point cloud $$X=\{\mathbf{x}_1,\ldots,\mathbf{x}_n\}\subset\mathbb{R}^D,$$ STC grows local ensembles around density seeds until geometric observables become scale invariant.

## Observable

For growth iteration $t$ define

$$O_t=(U_t,d_t,G_t),$$

where $U_t$ is the local tangent subspace, $d_t$ the inferred intrinsic dimension and $G_t$ the spectral-gap statistic.

## Difference operator

$$\mathcal{D}(O_t,O_{t+1})=(\delta_U,\delta_d,\delta_G).$$

### Tangent-space change

Given orthonormal tangent bases $$U_t,U_{t+1}\in\mathbb{R}^{D\times d},$$ compute

$$M=U_t^TU_{t+1}.$$

Compute the SVD $$M=P\Sigma Q^T.$$ If $\sigma_i$ are the singular values, then the principal angles are

$$\theta_i=\arccos(\sigma_i).$$

Define

$$\delta_U=\max_i\theta_i.$$

### Dimension change

$$\delta_d=|d_{t+1}-d_t|.$$

### Spectral-gap change

For $$G_t=\lambda_d/(\lambda_{d+1}+\varepsilon),$$ define

$$\delta_G=\left|\log\frac{G_{t+1}}{G_t}\right|.$$

## Accept criterion

Accept the smallest ensemble satisfying

$$\delta_U<\tau_U,\qquad \delta_d=0,\qquad \delta_G<\tau_G.$$

## Pseudocode

```text
for each seed
    repeat
        grow ensemble by fixed radius
        compute local covariance/SVD
        compute O_t=(U_t,d_t,G_t)
        if previous observable exists:
            Δ=D(O_t,O_{t+1})
            if Δ.U<τU and Δ.d==0 and Δ.G<τG:
                accept
                break
```

Each accepted ensemble becomes one atlas chart for LoManLe.