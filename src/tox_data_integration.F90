!> In multi-study omics analyses, it is often unclear whether biological replicates originating from different studies can be safely treated as sampling the same biological condition.
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
!| To address this, we introduce a Jensen–Shannon-Divergence based compatibility test (JSD-Comp-Test)
!| that empirically evaluates whether two sets of biological replicates exhibit comparable replicate-level variability.
!| Rather than comparing mean expression values, the method focuses on the distribution of signed residuals (replicate deviations from the gene-wise mean),
!| conditioned on mean expression levels to account for heteroscedasticity, which is a well-known property of omics data.
!|
!| The goal of this issue is to define, implement, and validate this compatibility test as a diagnostic tool that can be applied prior to data integration.
!| The test is intended to support principled decisions
!| on whether replicate sets from different studies should be merged or treated as distinct conditions within Tensor Omics workflows.
module tox_data_integration
    use safeguard
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use tox_data_integration_preprocessing, only: compute_gene_means_impl => compute_gene_means,&
        compute_residuals_impl => compute_residuals,&
        pool_means_alloc_impl => pool_means_alloc,&
        pool_means_impl => pool_means,&
        calc_neighborhood_size_impl => calc_neighborhood_size,&
        construct_neighborhoods_alloc_impl => construct_neighborhoods_alloc,&
        construct_neighborhoods_impl => construct_neighborhoods
    use tox_data_integration_jsd, only: build_residual_histograms_impl => build_residual_histograms,&
        compute_divergence_per_reference_point_impl => compute_divergence_per_reference_point,&
        compute_weighted_global_divergence_impl => compute_weighted_global_divergence
    implicit none

    interface compute_gene_means
        module procedure compute_gene_means_impl
    end interface compute_gene_means

    interface compute_residuals
        module procedure compute_residuals_impl
    end interface compute_residuals

    interface pool_means_alloc
        module procedure pool_means_alloc_impl
    end interface pool_means_alloc

    interface pool_means
        module procedure pool_means_impl
    end interface pool_means

    interface calc_neighborhood_size
        module procedure calc_neighborhood_size_impl
    end interface calc_neighborhood_size

    interface construct_neighborhoods_alloc
        module procedure construct_neighborhoods_alloc_impl
    end interface construct_neighborhoods_alloc

    interface construct_neighborhoods
        module procedure construct_neighborhoods_impl
    end interface construct_neighborhoods

    interface build_residual_histograms
        module procedure build_residual_histograms_impl
    end interface build_residual_histograms

    interface compute_divergence_per_reference_point
        module procedure compute_divergence_per_reference_point_impl
    end interface compute_divergence_per_reference_point

    interface compute_weighted_global_divergence
        module procedure compute_weighted_global_divergence_impl
    end interface compute_weighted_global_divergence
end module tox_data_integration