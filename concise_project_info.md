Tensor Omics is a high-performance framework for explainable, geometry-based analysis of multimodal omics and related high-dimensional datasets. Instead of relying on black-box models, it treats expression profiles, clinical measures, or socioeconomic indicators as vectors in semantically meaningful spaces (e.g. tissues, disease stages, conditions). By measuring distances, angles, projections, and trajectories in these spaces, Tensor Omics enables direct comparison of activity across genes, paralogs, sexes, species, or patient groups. This geometric approach makes complex multivariate patterns interpretable and reproducible while remaining robust to sparsity and noise.

Designed for distributed high-performance computing, Tensor Omics is implemented in Fortran and C with OpenMP parallelisation, SIMD optimisation, and Fortran Coarrays, making the algorithms embarrassingly parallel and suitable for federated datasets where privacy and efficiency are critical. Scientific use cases include: detecting disease biomarkers and subtype-specific trajectories in medical data; quantifying divergence and neofunctionalization of gene duplicates in plant and animal transcriptomes; and reconstructing global gender-equality trajectories from socioeconomic indicators. Across these domains, Tensor Omics provides a unified, geometry-driven methodology for discovering explanatory patterns in heterogeneous, high-dimensional data.

## Key Features

* **Geometry-based analysis**: distances, angles, projections, and trajectory shifts are used as primary primitives.
* **Explainable outputs**: results are interpretable in terms of vector geometry rather than opaque model coefficients.
* **Multi-modal integration**: unifies transcriptomics, proteomics, metabolomics, clinical, or socioeconomic data within one framework.
* **Parallel and federated**: implemented in Fortran/C with OpenMP, SIMD, and coarrays for efficient large-scale computation on distributed datasets.
* **Robust to sparsity and noise**: percentile-based empirical thresholds and local geometric measures enable stability.
* **Broad applications**: demonstrated on medical biomarker discovery, gene duplication outcomes, developmental trajectories, and socioeconomic indicators.

## How the API is shaped

Every procedure is published under one name, the same in Fortran, C, Python and R. Its documentation, found on these pages, is the documentation in every language.

**Two tiers, where there is something to take over.** `foo` is the entry point to reach for first: it validates its arguments, allocates the work arrays and prepares what it can derive itself (a sorted permutation, a workspace size, a threshold), then computes. `foo_expert`, where it exists, validates and computes with what you hand it. It allocates nothing and prepares nothing, so you can reuse your own buffers or supply a different sort order. Python and R get the expert tier only where it offers more than buffers, because they allocate the work arrays themselves either way.

**The calling language allocates the outputs.** No procedure allocates an array and hands it back. In Fortran and C you pass the output arrays in. In Python and R the binding allocates them as NumPy arrays or R vectors before the call, and the library writes into that memory directly. Every result is therefore an ordinary object of your language: its garbage collector frees it, and there is nothing to release by hand. The only memory the library allocates itself is scratch space inside `foo`, and it is released before `foo` returns.

This is deliberate. An all-in-one entry point that allocates its outputs inside the library and returns them would save the bindings a few lines. But memory the library allocated must also be freed by the library, so every language would need its own finaliser calling back into it, or a copy of each result into a native array. The cleanup would still be needed, only harder to get right and in more places.


*Tensor Omics shows that geometry, when treated not as preprocessing but as the central instrument of analysis, can open entirely new ways to read complex biological and social data — simple, transparent, and surprisingly powerful.*
