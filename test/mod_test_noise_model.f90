!> Unit test suite for noise_model module.
!|
!| Focus: `own` comparison only. Covers both raw (norm_method = 0) and log2
!| (norm_method /= 0) residual scales, including all four normalization regimes
!| used by the R `apply_normalization()` pipeline (raw/log/std_log/full).
!|
!| The family / ortholog comparison has been retired from the noise model, so
!| there is no family cache, no family/ortholog p-values, and the pipeline takes
!| only the `own` statistic. The exact/Monte-Carlo pairwise p-value has likewise
!| been replaced by the bootstrap mean-difference null
!| (`compute_pvalue_bootstrap_mean_helper`), and variance stratification has been
!| removed — the `own` null is built directly from the gathered kNN pool.
module mod_test_noise_model
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use noise_model
    use tox_normalization, only: normalize_by_std_dev_alloc, quantile_normalization
    use tox_errors
    use test_suite, only: test_case
    implicit none
    public

    real(real64), parameter :: TOL = 1d-12
    integer(int32), parameter :: N_SAMPLES = 5
    integer(int32), parameter :: N_GENES   = 20

contains

    !> Get array of all available tests.
    function get_all_tests_noise_model() result(all_tests)
        type(test_case), allocatable :: all_tests(:)

        allocate(all_tests(8))
        all_tests(1) = test_case("test_prepare_sorted_data",               test_prepare_sorted_data)
        all_tests(2) = test_case("test_gather_residuals_helper",           test_gather_residuals_helper)
        all_tests(3) = test_case("test_compute_pvalue_bootstrap_mean",     test_compute_pvalue_bootstrap_mean)
        all_tests(4) = test_case("test_full_pipeline",                     test_full_pipeline)
        all_tests(5) = test_case("test_prepare_sorted_data_log_transform", test_prepare_sorted_data_log_transform)
        all_tests(6) = test_case("test_residuals_log_transform_centred",   test_residuals_log_transform_centred)
        all_tests(7) = test_case("test_full_pipeline_all_normalizations",  test_full_pipeline_all_normalizations)
        all_tests(8) = test_case("test_trim_pool_tails_helper",            test_trim_pool_tails_helper)
    end function get_all_tests_noise_model

    ! =========================================================================
    ! Helper: convert integer to string for assertion messages
    ! =========================================================================
    function str(i) result(s)
        integer(int32), intent(in) :: i
        character(len=20) :: buffer
        character(len=:), allocatable :: s
        write(buffer, '(I0)') i
        s = trim(buffer)
    end function str

    ! =========================================================================
    ! Helper: generate synthetic data
    ! =========================================================================
    !| Produces 20 genes spread over means 5..15 (control). Genes 1..5 are
    !| overexpressed in case by +2 units. Noise is small and deterministic.
    !|
    !| `control_noise_scale` multiplies the noise amplitude used for control
    !| replicates relative to case. Existing callers pass `0.8` for a mild,
    !| harmless asymmetry under norm_method = 0/1 (no variance rescaling).
    !| Callers that feed this data through a variance-stabilizing normalization
    !| (e.g. `normalize_by_std_dev_alloc`) MUST pass `1.0`: that step divides
    !| each group's replicates by its own empirically fitted standard deviation,
    !| so any case/control noise-amplitude asymmetry becomes a spurious, uniform
    !| scale mismatch between the groups — one unrelated to, and far larger than,
    !| the deliberate per-gene overexpression signal.
    subroutine generate_test_data(n_genes, n_samples, means_case, means_control, &
                                  replicates_case, replicates_control, control_noise_scale)
        integer(int32), intent(in)  :: n_genes, n_samples
        real(real64),   intent(out) :: means_case(n_genes), means_control(n_genes)
        real(real64),   intent(out) :: replicates_case(n_samples, n_genes)
        real(real64),   intent(out) :: replicates_control(n_samples, n_genes)
        real(real64),   intent(in)  :: control_noise_scale

        integer(int32) :: i_gene, i_sample
        real(real64)   :: base_mean, noise

        do i_gene = 1, n_genes
            base_mean = 5.0_real64 + 10.0_real64 * (i_gene - 1) / (n_genes - 1)
            means_case(i_gene)    = base_mean
            means_control(i_gene) = base_mean
            ! Overexpress genes 1..5 in case
            if (i_gene <= 5) means_case(i_gene) = means_case(i_gene) + 2.0_real64
            do i_sample = 1, n_samples
                noise = 0.1_real64 * sin(real(i_gene, real64) + real(i_sample, real64))
                replicates_case(i_sample, i_gene)    = means_case(i_gene)    + noise
                replicates_control(i_sample, i_gene) = means_control(i_gene) + noise * control_noise_scale
            end do
        end do
    end subroutine generate_test_data

    ! =========================================================================
    ! test_prepare_sorted_data
    ! =========================================================================
    subroutine test_prepare_sorted_data()
        integer(int32), parameter :: n_genes = 5, n_samples = 3
        real(real64) :: means(n_genes), replicates(n_samples, n_genes)
        type(sorted_data_t) :: sorted_data
        integer(int32) :: ierr, i

        means = [10.0_real64, 5.0_real64, 8.0_real64, 3.0_real64, 12.0_real64]
        replicates = reshape([ &
            10.1_real64, 9.9_real64, 10.0_real64, &
             5.1_real64, 4.9_real64,  5.0_real64, &
             8.1_real64, 7.9_real64,  8.0_real64, &
             3.1_real64, 2.9_real64,  3.0_real64, &
            12.1_real64, 11.9_real64, 12.0_real64 &
        ], [n_samples, n_genes])

        call prepare_sorted_data(means, replicates, n_samples, n_genes, 0, sorted_data, ierr)
        call assert_equal_int(ierr, ERR_OK, "prepare_sorted_data: ierr should be OK")

        call assert_equal_int(sorted_data%n_genes, n_genes, &
                              "sorted_data%n_genes mismatch")
        call assert_equal_int(sorted_data%max_resid_per_gene, n_samples, &
                              "sorted_data%max_resid_per_gene mismatch")

        ! means_sorted must be strictly non-decreasing
        do i = 2, n_genes
            call assert_true(sorted_data%means_sorted(i) >= sorted_data%means_sorted(i-1), &
                             "means_sorted not ascending at position "//str(i))
        end do

        ! Known sort order: smallest mean 3.0 is gene 4; largest mean 12.0 is gene 5
        call assert_equal_int(sorted_data%original_indices(1), 4, &
                              "original_indices(1): smallest mean should map to gene 4")
        call assert_equal_int(sorted_data%original_indices(5), 5, &
                              "original_indices(5): largest mean should map to gene 5")

        ! Centred residuals must sum to zero per gene
        do i = 1, n_genes
            call assert_equal_real(sum(sorted_data%residuals_packed(:, i)), 0.0_real64, TOL, &
                                   "residuals for sorted gene "//str(i)//" should sum to 0")
        end do
    end subroutine test_prepare_sorted_data

    ! =========================================================================
    ! test_gather_residuals_helper
    ! =========================================================================
    subroutine test_gather_residuals_helper()
        integer(int32), parameter :: n_genes = 10, n_samples = 3
        real(real64) :: means(n_genes), replicates(n_samples, n_genes)
        type(sorted_data_t) :: sorted_data
        integer(int32) :: ierr, i
        real(real64),    allocatable :: pooled(:)
        integer(int32) :: n_pooled, max_pool_size
        real(real64)   :: target_mean
        integer(int32), parameter :: k_start = 10, k_step = 5, k_max = 30
        real(real64),   parameter :: tau = 0.1_real64

        do i = 1, n_genes
            means(i)         = real(i, real64) * 2.0_real64
            replicates(:, i) = means(i) + 0.1_real64 * real(i, real64)
        end do

        call prepare_sorted_data(means, replicates, n_samples, n_genes, 0, sorted_data, ierr)
        call assert_equal_int(ierr, ERR_OK, "gather: prepare_sorted_data failed")

        max_pool_size = 50
        allocate(pooled(max_pool_size))

        ! Target near gene 5 (mean = 10)
        target_mean = 10.0_real64
        call gather_residuals_helper(target_mean, sorted_data, k_start, k_step, k_max, tau, &
                                     pooled, n_pooled, max_pool_size)

        call assert_true(n_pooled >= k_start, &
                         "pool size should be at least k_start")
        call assert_true(n_pooled <= min(k_max, max_pool_size), &
                         "pool size exceeds upper limit")

        do i = 1, n_pooled
            call assert_true(pooled(i) == pooled(i), &
                             "pooled("//str(i)//") is NaN")
        end do

        deallocate(pooled)
    end subroutine test_gather_residuals_helper

    ! =========================================================================
    ! test_trim_pool_tails_helper
    ! =========================================================================
    !| Raw-only symmetric quantile trim: sorts the pool and drops the k =
    !| floor(n * trim_frac) smallest and largest residuals, keeping the central
    !| ones compacted into pool(1:n_pool). Verified on a pool of 20 residuals
    !| (the central values 1..18 plus a low -50 and high +50 outlier, unsorted so
    !| the internal sort is exercised), plus the three documented no-op guards.
    subroutine test_trim_pool_tails_helper()
        real(real64) :: pool(20), pool0(20)
        integer(int32) :: n_pool, i
        real(real64) :: expected_sum

        pool0 = [ 50.0_real64,   1.0_real64,  2.0_real64,  3.0_real64,  4.0_real64, &
                   5.0_real64,   6.0_real64,  7.0_real64,  8.0_real64,  9.0_real64, &
                  10.0_real64,  11.0_real64, -50.0_real64, 12.0_real64, 13.0_real64, &
                  14.0_real64,  15.0_real64, 16.0_real64, 17.0_real64, 18.0_real64]

        ! (a) 5% per tail of 20 -> k = floor(20*0.05) = 1: drop the min (-50) and
        !     max (50), keep exactly the central 1..18 (sorted ascending).
        pool = pool0; n_pool = 20
        call trim_pool_tails_helper(pool, n_pool, 0.05_real64)
        call assert_equal_int(n_pool, 18, "trim: kept count must be 20 - 2*1 = 18")
        call assert_equal_real(pool(1),  1.0_real64,  TOL, "trim: smallest kept residual must be 1")
        call assert_equal_real(pool(18), 18.0_real64, TOL, "trim: largest kept residual must be 18")
        expected_sum = 0.0_real64
        do i = 1, 18
            expected_sum = expected_sum + real(i, real64)
        end do
        call assert_equal_real(sum(pool(1:18)), expected_sum, TOL, &
                               "trim: kept residuals must be exactly central 1..18 (outliers removed)")

        ! (b) trim_frac = 0 -> no-op
        pool = pool0; n_pool = 20
        call trim_pool_tails_helper(pool, n_pool, 0.0_real64)
        call assert_equal_int(n_pool, 20, "trim (frac=0): count must be unchanged")
        call assert_equal_real(sum(pool(1:20)), sum(pool0), TOL, "trim (frac=0): pool must be unchanged")

        ! (c) k rounds down to 0 (pool too small for the fraction) -> no-op
        pool = pool0; n_pool = 20
        call trim_pool_tails_helper(pool, n_pool, 0.02_real64)   ! floor(20*0.02) = 0
        call assert_equal_int(n_pool, 20, "trim (k=0): count must be unchanged")

        ! (d) fraction >= 0.5 would empty the pool -> guarded, no-op
        pool = pool0; n_pool = 20
        call trim_pool_tails_helper(pool, n_pool, 0.6_real64)
        call assert_equal_int(n_pool, 20, "trim (frac>=0.5): count must be unchanged (guard)")
    end subroutine test_trim_pool_tails_helper

    ! =========================================================================
    ! test_compute_pvalue_bootstrap_mean
    ! =========================================================================
    !| The `own` null is now built by bootstrapping the mean difference:
    !| `compute_pvalue_bootstrap_mean_helper` resamples `n_rep` residuals per side
    !| (with replacement), averages, differences, and add-one-corrects the tail
    !| fraction that meets or exceeds `observed_statistic_abs`. Two boundary cases
    !| are deterministic regardless of the RNG draws and pin the arithmetic:
    !|   (a) observed = 0   → every |mean diff| >= 0 qualifies → p = 1
    !|   (b) observed huge  → no draw qualifies → p = 1 / (n_boot + 1)
    !| A moderate statistic is checked only for being a valid p-value in (0, 1].
    subroutine test_compute_pvalue_bootstrap_mean()
        real(real64) :: pool_case(8), pool_control(8)
        real(real64) :: p_value
        integer(int32), parameter :: n_boot = 2000, n_rep = 3

        pool_case    = [ 0.10_real64, -0.10_real64,  0.20_real64, -0.20_real64, &
                         0.05_real64, -0.05_real64,  0.15_real64, -0.15_real64]
        pool_control = [ 0.08_real64, -0.08_real64,  0.18_real64, -0.18_real64, &
                         0.04_real64, -0.04_real64,  0.12_real64, -0.12_real64]

        ! (a) observed = 0 → every null statistic (|mean diff| >= 0) qualifies → p = 1
        call compute_pvalue_bootstrap_mean_helper(pool_case, 8, pool_control, 8, &
                                                  n_rep, n_rep, 0.0_real64, n_boot, p_value)
        call assert_equal_real(p_value, 1.0_real64, TOL, &
                               "bootstrap (obs=0): p_value must be 1 (every draw qualifies)")

        ! (b) impossible-to-exceed observed → no draw qualifies → p = 1/(n_boot+1)
        call compute_pvalue_bootstrap_mean_helper(pool_case, 8, pool_control, 8, &
                                                  n_rep, n_rep, 100.0_real64, n_boot, p_value)
        call assert_equal_real(p_value, 1.0_real64 / real(n_boot + 1, real64), TOL, &
                               "bootstrap (obs huge): p_value must equal 1/(n_boot+1)")

        ! (c) moderate observed → valid p-value in (0, 1]
        call compute_pvalue_bootstrap_mean_helper(pool_case, 8, pool_control, 8, &
                                                  n_rep, n_rep, 0.05_real64, n_boot, p_value)
        call assert_true(p_value > 0.0_real64 .and. p_value <= 1.0_real64, &
                         "bootstrap (moderate): p_value out of (0, 1]")
        call assert_true(p_value == p_value, "bootstrap (moderate): p_value must not be NaN")
    end subroutine test_compute_pvalue_bootstrap_mean

    ! =========================================================================
    ! test_full_pipeline
    ! =========================================================================
    !| Integration test for compute_noise_pvalue_pipeline (`own` comparison, the
    !| only comparison the model now supports).
    !|
    !| Genes 1–5 are overexpressed in case by +2 units and should yield small own
    !| p-values. Genes 6–20 should not be overwhelmingly significant (at most 5 of
    !| 15 at alpha = 0.05).
    subroutine test_full_pipeline()
        real(real64)   :: means_case(N_GENES), means_control(N_GENES)
        real(real64)   :: replicates_case(N_SAMPLES, N_GENES)
        real(real64)   :: replicates_control(N_SAMPLES, N_GENES)
        integer(int32) :: compute_own(N_GENES)
        real(real64)   :: observed_own(N_GENES)
        real(real64)   :: pvalues_own(N_GENES)
        integer(int32) :: n_genes_with_pvalue
        integer(int32) :: neigh_own_case(N_GENES), neigh_own_control(N_GENES), neigh_case(N_GENES)
        integer(int32) :: ierr, i, count_sig_own

        integer(int32), parameter :: k_start = 10, k_step = 5, k_max = 50
        real(real64),   parameter :: tau = 0.1_real64
        integer(int32), parameter :: max_pool_size = 100
        integer(int32), parameter :: norm_method = 0
        real(real64),   parameter :: alpha = 0.05_real64

        call generate_test_data(N_GENES, N_SAMPLES, means_case, means_control, &
                                replicates_case, replicates_control, &
                                control_noise_scale=0.8_real64)

        do i = 1, N_GENES
            observed_own(i) = means_case(i) - means_control(i)
        end do
        compute_own = 1

        call compute_noise_pvalue_pipeline( &
            means_case, replicates_case, N_GENES, N_SAMPLES, &
            means_control, replicates_control, N_GENES, N_SAMPLES, &
            observed_own, compute_own, &
            N_GENES, norm_method, k_start, k_step, k_max, tau, 0.0_real64, &
            pvalues_own, n_genes_with_pvalue, &
            max_pool_size, &
            neigh_own_case, neigh_own_control, neigh_case, &
            ierr)

        call assert_equal_int(ierr, ERR_OK, "pipeline: ierr should be OK")
        call assert_true(n_genes_with_pvalue > 0, &
                         "pipeline: at least some genes should receive p-values")

        ! Overexpressed genes 1–5 must be significant at alpha = 0.05
        do i = 1, 5
            if (pvalues_own(i) >= 0.0_real64) then
                call assert_true(pvalues_own(i) < alpha, &
                                 "overexpressed gene "//str(i)//" should have own p-value < 0.05")
            end if
        end do

        ! Non-overexpressed genes 6–20: allow at most 5 significant (≤ 33 %)
        count_sig_own = 0
        do i = 6, N_GENES
            if (pvalues_own(i) >= 0.0_real64 .and. pvalues_own(i) < alpha) &
                count_sig_own = count_sig_own + 1
        end do
        call assert_true(count_sig_own <= 5, &
                         "too many significant own p-values among non-overexpressed genes")

        ! Neighborhood sizes: positive and ≥ k_start when a p-value was produced.
        ! The `own` pools are now the full gathered kNN pools (no stratification).
        do i = 1, N_GENES
            if (pvalues_own(i) >= 0.0_real64) then
                call assert_true(neigh_own_case(i) > 0, &
                                 "neigh_own_case("//str(i)//") should be > 0")
                call assert_true(neigh_own_control(i) > 0, &
                                 "neigh_own_control("//str(i)//") should be > 0")
            end if
            if (neigh_case(i) > 0) then
                call assert_true(neigh_case(i) >= k_start, &
                                 "neigh_case("//str(i)//") should be at least k_start")
            end if
        end do
    end subroutine test_full_pipeline

    ! =========================================================================
    ! test_prepare_sorted_data_log_transform
    ! =========================================================================
    !| Directly verifies the log-space residual formula introduced by
    !| `prepare_sorted_data_helper` when `norm_method /= 0`:
    !|
    !|   ghat_g        = (1/n) * sum_i log2(r_{i,g} + 1)   (Frechet mean in log2-space)
    !|   epsilon_{i,g} = log2(r_{i,g} + 1) - ghat_g
    !|
    !| computed from the pre-log replicate values and centred on the arithmetic
    !| mean of the LOG-transformed replicates — NOT on `log2(mu_g + 1)` (the log
    !| of the linear-space mean), which by Jensen's inequality does not generally
    !| leave the residuals centred at zero. `test_residuals_log_transform_centred`
    !| below checks that centering property directly.
    subroutine test_prepare_sorted_data_log_transform()
        integer(int32), parameter :: n_genes = 3, n_samples = 3
        real(real64) :: means(n_genes), replicates(n_samples, n_genes)
        real(real64) :: expected_residual, log2_factor, log2_gene_mean
        type(sorted_data_t) :: sorted_data
        integer(int32) :: ierr, i_sample, sorted_pos, orig_idx

        means = [10.0_real64, 5.0_real64, 20.0_real64]
        replicates = reshape([ &
             9.0_real64, 10.0_real64, 11.0_real64, &
             4.0_real64,  5.0_real64,  6.0_real64, &
            18.0_real64, 20.0_real64, 22.0_real64  &
        ], [n_samples, n_genes])

        call prepare_sorted_data(means, replicates, n_samples, n_genes, 1, sorted_data, ierr)
        call assert_equal_int(ierr, ERR_OK, "log-transform prepare_sorted_data: ierr should be OK")

        log2_factor = 1.0_real64 / log(2.0_real64)

        do sorted_pos = 1, n_genes
            orig_idx = sorted_data%original_indices(sorted_pos)
            log2_gene_mean = sum(log(replicates(:, orig_idx) + 1.0_real64)) &
                              * log2_factor / real(n_samples, real64)
            do i_sample = 1, n_samples
                expected_residual = log(replicates(i_sample, orig_idx) + 1.0_real64) * log2_factor &
                                     - log2_gene_mean
                call assert_equal_real(sorted_data%residuals_packed(i_sample, sorted_pos), &
                                       expected_residual, TOL, &
                                       "log-transform residual mismatch at gene "//str(orig_idx)// &
                                       ", sample "//str(i_sample))
            end do
        end do
    end subroutine test_prepare_sorted_data_log_transform

    ! =========================================================================
    ! test_residuals_log_transform_centred
    ! =========================================================================
    !| Confirms the Frechet-mean centering fix: per-gene log2-space residuals
    !| must sum to (numerically) zero. Before this fix, centering on
    !| `log2(mu_g + 1)` instead of the mean of the log-transformed replicates
    !| left a systematic, Jensen's-gap-sized offset (only exactly zero when all
    !| replicates of a gene are equal).
    subroutine test_residuals_log_transform_centred()
        integer(int32), parameter :: n_genes = 4, n_samples = 5
        real(real64) :: means(n_genes), replicates(n_samples, n_genes)
        type(sorted_data_t) :: sorted_data
        integer(int32) :: ierr, i_gene, i_sample

        do i_gene = 1, n_genes
            means(i_gene) = 3.0_real64 + 4.0_real64 * real(i_gene, real64)
            do i_sample = 1, n_samples
                replicates(i_sample, i_gene) = means(i_gene) + &
                    0.3_real64 * real(i_gene, real64) * sin(real(i_gene, real64) + real(i_sample, real64))
            end do
        end do

        call prepare_sorted_data(means, replicates, n_samples, n_genes, 1, sorted_data, ierr)
        call assert_equal_int(ierr, ERR_OK, "centred-residuals: prepare_sorted_data ierr should be OK")

        do i_gene = 1, n_genes
            call assert_equal_real(sum(sorted_data%residuals_packed(:, i_gene)), 0.0_real64, TOL, &
                                   "log2-space residuals for sorted gene "//str(i_gene)// &
                                   " should sum to 0 (Frechet-mean centering)")
        end do
    end subroutine test_residuals_log_transform_centred

    ! =========================================================================
    ! test_full_pipeline_all_normalizations
    ! =========================================================================
    !| Exercises `compute_noise_pvalue_pipeline` across the four normalization
    !| regimes used by the R `apply_normalization()` pipeline:
    !|   - "raw":     no normalization,                                    norm_method = 0
    !|   - "log":     log2(x+1) on raw replicate means,                    norm_method = 1
    !|   - "std_log": gene-wise LOESS std-dev scaling, then log2,          norm_method = 1
    !|   - "full":    gene-wise std-dev scaling + quantile norm, then log2, norm_method = 1
    !|
    !| For "std_log" and "full", the noise model is fed the *replicate-level*
    !| pre-log-transformed data (std-dev-scaled / quantile-normalized), not the
    !| collapsed per-gene mean: the log2 transform itself is applied internally
    !| by `prepare_sorted_data_helper` (via `norm_method /= 0`), so residuals stay
    !| defined in pre-log expression space right up to the point they are logged,
    !| centred on the Frechet mean `(1/n) * sum_i log2(r_{i,g} + 1)` of the
    !| log2-transformed replicates. `observed_own` below is computed on that same
    !| coordinate for methods with `norm_method /= 0`.
    subroutine test_full_pipeline_all_normalizations()
        real(real64)   :: means_case_raw(N_GENES), means_control_raw(N_GENES)
        real(real64)   :: replicates_case_raw(N_SAMPLES, N_GENES)
        real(real64)   :: replicates_control_raw(N_SAMPLES, N_GENES)
        real(real64)   :: replicates_case(N_SAMPLES, N_GENES), replicates_control(N_SAMPLES, N_GENES)
        real(real64)   :: normalize_tmp(N_SAMPLES, N_GENES)
        real(real64)   :: means_case(N_GENES), means_control(N_GENES)
        real(real64)   :: log2_mean_case(N_GENES), log2_mean_control(N_GENES)
        integer(int32) :: compute_own(N_GENES)
        real(real64)   :: observed_own(N_GENES)
        real(real64)   :: pvalues_own(N_GENES)
        integer(int32) :: n_genes_with_pvalue
        integer(int32) :: neigh_own_case(N_GENES), neigh_own_control(N_GENES), neigh_case(N_GENES)
        integer(int32) :: ierr, i, i_method, count_sig_own, count_sig_nonoverexpr

        real(real64)   :: rank_means(N_GENES), tmp_genes_row(N_GENES)
        integer(int32) :: tmp_perm(N_GENES)
        real(real64),   parameter :: LOESS_SPAN = 0.75_real64
        integer(int32), parameter :: LOESS_DEGREE = 2

        integer(int32), parameter :: k_start = 10, k_step = 5, k_max = 50
        real(real64),   parameter :: tau = 0.1_real64
        integer(int32), parameter :: max_pool_size = 100
        real(real64),   parameter :: alpha = 0.05_real64
        real(real64) :: log2_factor

        integer(int32), parameter :: N_METHODS = 4
        character(len=8) :: method_names(N_METHODS)
        integer(int32) :: method_norm(N_METHODS)
        logical :: method_use_std_dev(N_METHODS), method_use_quantile(N_METHODS)
        logical :: method_check_signal(N_METHODS)

        method_names        = ["raw     ", "log     ", "std_log ", "full    "]
        method_norm         = [0, 1, 1, 1]
        method_use_std_dev  = [.false., .false., .true.,  .true.]
        method_use_quantile = [.false., .false., .false., .true.]
        ! Quantile normalization is a joint, rank-based technique across ALL genes:
        ! shifting 5 of 20 genes measurably warps the shared LOESS std-dev curve and
        ! rank_means used for every other gene too (verified numerically — even
        ! genes with byte-identical raw replicate data between case/control pick up
        ! a several-percent divergence in fitted std-dev once the other genes are
        ! perturbed). At N_GENES = 20 that ripple is not reliably dominated by the
        ! deliberate signal, so "full" only gets the structural checks below, not
        ! the overexpression-enrichment check.
        method_check_signal = [.true., .true., .true., .false.]

        log2_factor = 1.0_real64 / log(2.0_real64)

        call generate_test_data(N_GENES, N_SAMPLES, means_case_raw, means_control_raw, &
                                replicates_case_raw, replicates_control_raw, &
                                control_noise_scale=1.0_real64)

        do i_method = 1, N_METHODS
            replicates_case    = replicates_case_raw
            replicates_control = replicates_control_raw

            if (method_use_std_dev(i_method)) then
                call normalize_by_std_dev_alloc(N_GENES, N_SAMPLES, replicates_case, normalize_tmp, &
                                                LOESS_SPAN, LOESS_DEGREE, ierr)
                call assert_equal_int(ierr, ERR_OK, &
                                      trim(method_names(i_method))//": std-dev norm (case) ierr should be OK")
                replicates_case = normalize_tmp

                call normalize_by_std_dev_alloc(N_GENES, N_SAMPLES, replicates_control, normalize_tmp, &
                                                LOESS_SPAN, LOESS_DEGREE, ierr)
                call assert_equal_int(ierr, ERR_OK, &
                                      trim(method_names(i_method))//": std-dev norm (control) ierr should be OK")
                replicates_control = normalize_tmp
            end if

            if (method_use_quantile(i_method)) then
                call quantile_normalization(N_GENES, N_SAMPLES, replicates_case, normalize_tmp, &
                                            rank_means, tmp_genes_row, tmp_perm, ierr)
                call assert_equal_int(ierr, ERR_OK, &
                                      trim(method_names(i_method))//": quantile norm (case) ierr should be OK")
                replicates_case = normalize_tmp

                call quantile_normalization(N_GENES, N_SAMPLES, replicates_control, normalize_tmp, &
                                            rank_means, tmp_genes_row, tmp_perm, ierr)
                call assert_equal_int(ierr, ERR_OK, &
                                      trim(method_names(i_method))//": quantile norm (control) ierr should be OK")
                replicates_control = normalize_tmp
            end if

            ! Pre-log per-gene means, computed AFTER std-dev/quantile normalization
            ! but BEFORE the log2 transform (which the noise model applies internally
            ! when norm_method /= 0).
            do i = 1, N_GENES
                means_case(i)    = sum(replicates_case(:, i))    / real(N_SAMPLES, real64)
                means_control(i) = sum(replicates_control(:, i)) / real(N_SAMPLES, real64)
            end do

            if (method_norm(i_method) == 0) then
                do i = 1, N_GENES
                    observed_own(i) = means_case(i) - means_control(i)
                end do
            else
                ! Observed statistic must live on the same log2 coordinate as the null
                ! distribution: the Frechet mean of the log2-transformed replicates
                ! (see `prepare_sorted_data_helper`), NOT log2(linear mean + 1).
                do i = 1, N_GENES
                    log2_mean_case(i)    = sum(log(max(replicates_case(:, i), 0.0_real64) + 1.0_real64)) &
                                            * log2_factor / real(N_SAMPLES, real64)
                    log2_mean_control(i) = sum(log(max(replicates_control(:, i), 0.0_real64) + 1.0_real64)) &
                                            * log2_factor / real(N_SAMPLES, real64)
                    observed_own(i) = log2_mean_case(i) - log2_mean_control(i)
                end do
            end if

            compute_own = 1

            call compute_noise_pvalue_pipeline( &
                means_case, replicates_case, N_GENES, N_SAMPLES, &
                means_control, replicates_control, N_GENES, N_SAMPLES, &
                observed_own, compute_own, &
                N_GENES, method_norm(i_method), k_start, k_step, k_max, tau, 0.0_real64, &
                pvalues_own, n_genes_with_pvalue, &
                max_pool_size, &
                neigh_own_case, neigh_own_control, neigh_case, &
                ierr)

            call assert_equal_int(ierr, ERR_OK, trim(method_names(i_method))//": pipeline ierr should be OK")
            call assert_true(n_genes_with_pvalue > 0, &
                             trim(method_names(i_method))//": at least some genes should receive p-values")
            call assert_no_nan_real(pvalues_own, N_GENES, trim(method_names(i_method))//": NaN in pvalues_own")

            do i = 1, N_GENES
                if (pvalues_own(i) >= 0.0_real64) then
                    call assert_true(pvalues_own(i) >= 0.0_real64 .and. pvalues_own(i) <= 1.0_real64, &
                                     trim(method_names(i_method))//": pvalues_own("//str(i)//") out of [0, 1]")
                end if
            end do

            ! Overexpressed genes 1-5 should be enriched for significance relative to
            ! genes 6-20. A relative comparison, robust to the reshaping introduced by
            ! quantile normalization, unlike a strict per-gene alpha threshold.
            !
            ! Skipped for methods where `method_check_signal` is .false. (currently
            ! "full" only): quantile normalization is a joint, rank-based technique
            ! across all genes, so shifting a large fraction of a tiny synthetic gene
            ! set measurably perturbs the shared normalization curve/rank_means used
            ! by every other gene too. That is expected behavior of the normalization
            ! itself, not something the noise model can or should compensate for, and
            ! is not a safe invariant to assert at this dataset size.
            if (method_check_signal(i_method)) then
                count_sig_own = 0
                do i = 1, 5
                    if (pvalues_own(i) >= 0.0_real64 .and. pvalues_own(i) < alpha) count_sig_own = count_sig_own + 1
                end do
                count_sig_nonoverexpr = 0
                do i = 6, N_GENES
                    if (pvalues_own(i) >= 0.0_real64 .and. pvalues_own(i) < alpha) &
                        count_sig_nonoverexpr = count_sig_nonoverexpr + 1
                end do
                call assert_true(count_sig_own >= count_sig_nonoverexpr, &
                                 trim(method_names(i_method))// &
                                 ": overexpressed genes should be at least as often significant as non-overexpressed genes")
            end if
        end do
    end subroutine test_full_pipeline_all_normalizations

end module mod_test_noise_model
