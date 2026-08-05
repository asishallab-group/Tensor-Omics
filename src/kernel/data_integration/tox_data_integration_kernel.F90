!> # Jensen-Shannon-Divergence (JSD) Compatibility Test (JCT)
!|
!| In multi-study omics analyses, it is often unclear whether biological replicates originating from different studies can be safely treated as sampling the same biological condition.
!| Even when studies nominally target the same tissue and condition, differences in sample handling, sequencing technologies, preprocessing pipelines, or cohort,
!| composition can introduce batch effects that are not easily detectable from mean expression levels alone.
!|
!| This ambiguity has direct consequences for downstream analyses in Tensor Omics. Integrating incompatible replicate sets can:
!|
!|  - distort expression spaces,
!|  - affect distance-based analyses,
!|  - bias machine learning models,
!|
!| while unnecessarily separating compatible datasets reduces statistical power.
!|
!| To address this, we introduce a Jensen-Shannon-Divergence based compatibility test (JSD-Comp-Test)
!| that empirically evaluates whether two sets of biological replicates exhibit comparable replicate-level variability.
!| Rather than comparing mean expression values, the method focuses on the distribution of signed residuals (replicate deviations from the gene-wise mean),
!| conditioned on mean expression levels to account for heteroscedasticity, which is a well-known property of omics data.
!|
!| The family is split over four kernel modules -- preprocessing, JSD calculation, per-family
!| JSD calculation and the permutation test -- which this module gathers, so a caller reaches
!| the whole pipeline through it.
module tox_data_integration_kernel
    use tox_data_integration_preprocessing_kernel
    use tox_data_integration_jsd_kernel
    use tox_data_integration_per_family_kernel
    use tox_data_integration_stats_kernel
end module tox_data_integration_kernel
