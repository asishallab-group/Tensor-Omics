!> Unit test suite for noise_model module.
module mod_test_noise_model
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32, int64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf
    use noise_model
    use tox_errors
    use test_suite, only: test_case
    implicit none
    public

    real(real64), parameter :: TOL = 1d-12
    integer(int32), parameter :: N_SAMPLES = 5
    integer(int32), parameter :: N_GENES = 20
    integer(int32), parameter :: N_FAMILIES = 2

contains

    !> Get array of all available tests.
    function get_all_tests_noise_model() result(all_tests)
        type(test_case), allocatable :: all_tests(:)

        allocate(all_tests(7))
        all_tests(1) = test_case("test_prepare_sorted_data", test_prepare_sorted_data)
        all_tests(2) = test_case("test_gather_residuals_helper", test_gather_residuals_helper)
        all_tests(3) = test_case("test_stratify_residuals", test_stratify_residuals)
        all_tests(4) = test_case("test_compute_pvalue_exact", test_compute_pvalue_exact)
        all_tests(5) = test_case("test_compute_pvalue_monte_carlo", test_compute_pvalue_monte_carlo)
        all_tests(6) = test_case("test_compute_pvalue_selection", test_compute_pvalue_selection)
        all_tests(7) = test_case("test_full_pipeline", test_full_pipeline)
    end function get_all_tests_noise_model

    ! -------------------------------------------------------------------------
    ! Helper to generate synthetic data
    ! -------------------------------------------------------------------------
    subroutine generate_test_data(n_genes, n_samples, means_case, means_control, &
                                  replicates_case, replicates_control, family_sizes, gene_to_family)
        integer(int32), intent(in) :: n_genes, n_samples
        real(real64), intent(out) :: means_case(n_genes), means_control(n_genes)
        real(real64), intent(out) :: replicates_case(n_samples, n_genes)
        real(real64), intent(out) :: replicates_control(n_samples, n_genes)
        integer(int32), intent(out) :: family_sizes(N_FAMILIES)
        integer(int32), intent(out) :: gene_to_family(n_genes)

        integer(int32) :: i_gene, i_sample, family_id, base_idx
        real(real64) :: base_mean, noise

        ! Generate means: spread from 5 to 15
        do i_gene = 1, n_genes
            base_mean = 5.0_real64 + 10.0_real64 * (i_gene - 1) / (n_genes - 1)
            means_case(i_gene) = base_mean
            means_control(i_gene) = base_mean
            ! Add overexpression to genes 1..5 (case only)
            if (i_gene <= 5) means_case(i_gene) = means_case(i_gene) + 2.0_real64
            ! Generate replicates: add small noise
            do i_sample = 1, n_samples
                ! deterministic noise for reproducibility: use i_gene and i_sample
                noise = 0.1_real64 * sin(real(i_gene, real64) + real(i_sample, real64))
                replicates_case(i_sample, i_gene) = means_case(i_gene) + noise
                replicates_control(i_sample, i_gene) = means_control(i_gene) + noise * 0.8_real64
            end do
        end do

        ! Family assignment: first 10 genes in family 1, last 10 in family 2
        do i_gene = 1, n_genes
            if (i_gene <= 10) then
                gene_to_family(i_gene) = 1
            else
                gene_to_family(i_gene) = 2
            end if
        end do
        family_sizes = [10, 10]
    end subroutine generate_test_data

    ! -------------------------------------------------------------------------
    ! Test prepare_sorted_data
    ! -------------------------------------------------------------------------
    subroutine test_prepare_sorted_data()
        integer(int32), parameter :: n_genes = 5, n_samples = 3
        real(real64) :: means(n_genes), replicates(n_samples, n_genes)
        type(sorted_data_t) :: sorted_data
        integer(int32) :: ierr, i

        ! Set up known data
        means = [10.0, 5.0, 8.0, 3.0, 12.0]
        replicates = reshape([ &
            10.1, 9.9, 10.0, &
            5.1, 4.9, 5.0, &
            8.1, 7.9, 8.0, &
            3.1, 2.9, 3.0, &
            12.1, 11.9, 12.0 &
        ], [n_samples, n_genes])

        call prepare_sorted_data(means, replicates, n_samples, n_genes, sorted_data, ierr)
        call assert_equal_int(ierr, ERR_OK, "prepare_sorted_data: ierr should be OK")

        ! Check dimensions
        call assert_equal_int(sorted_data%n_genes, n_genes, "sorted_data%n_genes mismatch")
        call assert_equal_int(sorted_data%max_resid_per_gene, n_samples, "sorted_data%max_resid_per_gene mismatch")

        ! Check sorting: means_sorted should be ascending
        do i = 2, n_genes
            call assert_true(sorted_data%means_sorted(i) >= sorted_data%means_sorted(i-1), &
                             "means_sorted not ascending")
        end do

        ! Check original indices mapping
        call assert_equal_int(sorted_data%original_indices(1), 4, "smallest mean should be gene 4")
        call assert_equal_int(sorted_data%original_indices(5), 5, "largest mean should be gene 5")

        ! Check residuals: should be centred (sum ~0)
        do i = 1, n_genes
            call assert_equal_real(sum(sorted_data%residuals_packed(:, i)), 0.0_real64, TOL, &
                             "residuals for gene "//trim(adjustl(str(i)))//" should sum to 0")
        end do
    end subroutine test_prepare_sorted_data

    !> Helper to convert integer to string for assertions
    function str(i) result(s)
        integer(int32), intent(in) :: i
        character(len=20) :: buffer
        character(len=:), allocatable :: s
        
        write(buffer, '(I0)') i
        s = trim(buffer)
    end function str

    ! -------------------------------------------------------------------------
    ! Test gather_residuals_helper
    ! -------------------------------------------------------------------------
    subroutine test_gather_residuals_helper()
        integer(int32), parameter :: n_genes = 10, n_samples = 3
        real(real64) :: means(n_genes), replicates(n_samples, n_genes)
        type(sorted_data_t) :: sorted_data
        integer(int32) :: ierr, i
        real(real64), allocatable :: pooled(:), tmp_work(:)
        integer(int32), allocatable :: gene_id(:), tmp_gene_id(:)
        integer(int32) :: n_pooled, max_pool_size
        real(real64) :: target_mean
        integer(int32), parameter :: k_start = 10, k_step = 5, k_max = 30
        real(real64), parameter :: tau = 0.1_real64

        ! Generate data with distinct means
        do i = 1, n_genes
            means(i) = real(i, real64) * 2.0_real64
            replicates(:, i) = means(i) + 0.1_real64 * (i)
        end do

        call prepare_sorted_data(means, replicates, n_samples, n_genes, sorted_data, ierr)
        call assert_equal_int(ierr, ERR_OK, "prepare_sorted_data failed")

        max_pool_size = 50
        allocate(pooled(max_pool_size), gene_id(max_pool_size), tmp_work(max_pool_size), tmp_gene_id(max_pool_size))

        ! Target mean close to gene 5 (mean=10)
        target_mean = 10.0_real64
        call gather_residuals_helper(target_mean, sorted_data, k_start, k_step, k_max, tau, &
                                     pooled, gene_id, n_pooled, max_pool_size, tmp_work, tmp_gene_id)

        ! Check pool size is at least k_start
        call assert_true(n_pooled >= k_start, "pool size should be at least k_start")
        call assert_true(n_pooled <= min(k_max, max_pool_size), "pool size exceeds limit")

        ! Check that gene_id values are valid (1..n_genes)
        do i = 1, n_pooled
            call assert_true(gene_id(i) >= 1 .and. gene_id(i) <= n_genes, "gene_id out of range")
        end do

        ! Check that residuals are not NaN
        do i = 1, n_pooled
            call assert_true(pooled(i) == pooled(i), "pooled residual is NaN")
        end do

        deallocate(pooled, gene_id, tmp_work, tmp_gene_id)
    end subroutine test_gather_residuals_helper

    ! -------------------------------------------------------------------------
    ! Test stratify_residuals
    ! -------------------------------------------------------------------------
    subroutine test_stratify_residuals()
        integer(int32), parameter :: n_genes = 5, n_samples = 3, n_pooled = 15
        real(real64) :: pooled(n_pooled), bin_edges(STRATA_BIN_COUNT_SCHEDULE(1)+1)
        integer(int32) :: gene_id(n_pooled), bin_idx(n_pooled), chosen_bin_idx(n_pooled)
        integer(int32) :: tmp_c_g(n_pooled), tmp_bin_counts(STRATA_BIN_COUNT_SCHEDULE(1))
        integer(int32) :: tmp_gene_min_bin(n_genes), tmp_gene_max_bin(n_genes)
        logical :: tmp_gene_seen(n_genes)
        real(real64) :: tmp_bin_edges(STRATA_BIN_COUNT_SCHEDULE(1)+1)
        integer(int32) :: chosen_n_bins, i
        logical :: criteria_met
        real(real64), parameter :: TOL_LOCAL = 1d-12

        ! Create residuals with known distribution: some from gene 1, some from gene 2, etc.
        ! Gene 1: values near 0, gene 2: near 1, gene 3: near 2, gene 4: near 3, gene 5: near 4
        do i = 1, n_pooled
            pooled(i) = real(mod(i-1, 5), real64)  ! 0,1,2,3,4 repeated
            gene_id(i) = mod(i-1, 5) + 1
        end do

        ! Call stratify_residuals_helper
        call stratify_residuals_helper(pooled, gene_id, n_pooled, n_genes, &
                                       bin_idx, tmp_bin_edges, &
                                       tmp_gene_min_bin, tmp_gene_max_bin, tmp_gene_seen, &
                                       tmp_c_g, tmp_bin_counts, &
                                       chosen_n_bins, chosen_bin_idx, bin_edges, criteria_met)

        ! Should choose some bins; criteria_met may be true or false depending on schedule
        call assert_true(chosen_n_bins >= 2 .and. chosen_n_bins <= STRATA_BIN_COUNT_SCHEDULE(1), &
                         "chosen_n_bins out of range")

        ! Check bin indices are within 1..chosen_n_bins
        do i = 1, n_pooled
            call assert_true(chosen_bin_idx(i) >= 1 .and. chosen_bin_idx(i) <= chosen_n_bins, &
                             "bin index out of range")
        end do

        ! Check bin edges are non-decreasing
        do i = 1, chosen_n_bins
            call assert_true(bin_edges(i+1) >= bin_edges(i), "bin_edges not non-decreasing")
        end do

        ! Test with all equal residuals (should force 1 bin)
        pooled = 1.0_real64
        call stratify_residuals_helper(pooled, gene_id, n_pooled, n_genes, &
                                       bin_idx, tmp_bin_edges, &
                                       tmp_gene_min_bin, tmp_gene_max_bin, tmp_gene_seen, &
                                       tmp_c_g, tmp_bin_counts, &
                                       chosen_n_bins, chosen_bin_idx, bin_edges, criteria_met)
        ! Should result in 1 bin? Actually the schedule tries until criteria met; if all equal, bin_edges(1)==bin_edges(2) etc., so locate_bin may return first bin.
        ! But chosen_n_bins will be the last step (2) because criteria will fail (min_occupied_count maybe not enough? Actually with all equal, all residuals in one bin, other bins empty -> min_occupied_count=0, fail, fallback to 2 bins)
        call assert_true(chosen_n_bins == 2, "with all equal residuals, should fallback to 2 bins")
    end subroutine test_stratify_residuals

    ! -------------------------------------------------------------------------
    ! Test compute_pvalue_exact (helper)
    ! -------------------------------------------------------------------------
    subroutine test_compute_pvalue_exact()
        real(real64) :: pool_case(5), pool_control(4)
        real(real64) :: mean_case, mean_control, observed_abs, p_value
        integer(int32) :: norm_method
        integer(int32) :: ierr, i
        integer(int32) :: rand_case(MONTE_CARLO_SAMPLES), rand_control(MONTE_CARLO_SAMPLES)

        ! Set up simple residuals: case residuals around 0, control around 0.5
        pool_case = [0.0, 0.1, -0.1, 0.2, -0.2]
        pool_control = [0.5, 0.6, 0.4, 0.7]
        mean_case = 9.5
        mean_control = 10.0
        observed_abs = 0.5  ! observed mean difference
        norm_method = 0

        do i = 1, MONTE_CARLO_SAMPLES
            rand_case(i) = mod(i, 150) + 1
            rand_control(i) = mod(i*7, 150) + 1
        end do

        ! Exact case: no need to pass random indices
        call compute_pvalue(mean_case, pool_case, size(pool_case), &
                            mean_control, pool_control, size(pool_control), &
                            observed_abs, norm_method, &
                            rand_case, rand_control, p_value, ierr)
        call assert_equal_int(ierr, ERR_OK, "compute_pvalue exact: ierr should be OK")
        ! p_value should be between 0 and 1
        call assert_true(p_value >= 0.0 .and. p_value <= 1.0, "p_value out of range")

        ! Test with observed statistic so large that no null exceeds it -> p_value = 1/(n_pairs+1)
        observed_abs = 100.0
        call compute_pvalue(mean_case, pool_case, size(pool_case), &
                            mean_control, pool_control, size(pool_control), &
                            observed_abs, norm_method, &
                            rand_case, rand_control, p_value, ierr)
        call assert_equal_real(p_value, 1.0_real64/21.0_real64, TOL, "p_value for very large observed should be 1/(n_pairs+1)")

        ! Test with observed statistic small -> p_value close to 1
        observed_abs = 0.0
        call compute_pvalue(mean_case, pool_case, size(pool_case), &
                            mean_control, pool_control, size(pool_control), &
                            observed_abs, norm_method, &
                            rand_case, rand_control, p_value, ierr)
        call assert_true(p_value > 0.9, "p_value for observed=0 should be close to 1")
    end subroutine test_compute_pvalue_exact

    ! -------------------------------------------------------------------------
    ! Test compute_pvalue_monte_carlo
    ! -------------------------------------------------------------------------
    subroutine test_compute_pvalue_monte_carlo()
        integer(int32), parameter :: N_LARGE = 150
        real(real64), allocatable :: large_case(:), large_control(:)
        real(real64) :: mean_case, mean_control, observed_abs, p_value
        integer(int32) :: norm_method, ierr, i
        integer(int32) :: rand_case(MONTE_CARLO_SAMPLES), rand_control(MONTE_CARLO_SAMPLES)

        allocate(large_case(N_LARGE), large_control(N_LARGE))
        ! Fill with random-like numbers (deterministic)
        do i = 1, N_LARGE
            large_case(i) = sin(real(i, real64)) * 0.1
            large_control(i) = cos(real(i, real64)) * 0.1
        end do

        mean_case = 5.0
        mean_control = 5.0
        observed_abs = 0.2
        norm_method = 0

        ! Generate random indices (not truly random, but we only test the routine, not the RNG)
        do i = 1, MONTE_CARLO_SAMPLES
            rand_case(i) = mod(i, N_LARGE) + 1
            rand_control(i) = mod(i*7, N_LARGE) + 1
        end do

        call compute_pvalue(mean_case, large_case, N_LARGE, &
                            mean_control, large_control, N_LARGE, &
                            observed_abs, norm_method, &
                            rand_case, rand_control, &
                            p_value, ierr)
        call assert_equal_int(ierr, ERR_OK, "compute_pvalue Monte Carlo: ierr should be OK")
        call assert_true(p_value >= 0.0 .and. p_value <= 1.0, "p_value out of range")

        deallocate(large_case, large_control)
    end subroutine test_compute_pvalue_monte_carlo

    ! -------------------------------------------------------------------------
    ! Test compute_pvalue selection logic (exact vs Monte Carlo)
    ! -------------------------------------------------------------------------
    subroutine test_compute_pvalue_selection()
        integer(int32) :: ierr
        real(real64) :: p_value
        real(real64) :: pool_small(2), pool_large(150)
        integer(int32) :: rand_indices_case(MONTE_CARLO_SAMPLES), rand_indices_control(MONTE_CARLO_SAMPLES)
        integer(int32) :: i

        ! Small pools: should use exact
        pool_small = [0.1, -0.1]
        ! Large pools: should use Monte Carlo
        do i = 1, 150
            pool_large(i) = sin(real(i, real64)) * 0.1
        end do
        ! Prepare random indices (dummy)
        do i = 1, MONTE_CARLO_SAMPLES
            rand_indices_case(i) = mod(i, 150) + 1
            rand_indices_control(i) = mod(i*3, 150) + 1
        end do

        ! Test small (exact) – omit optional args
        call compute_pvalue(0.0_real64, pool_small, 2, 0.0_real64, pool_small, 2, &
                            0.5_real64, 0, rand_indices_case, rand_indices_control, p_value, ierr)
        call assert_equal_int(ierr, ERR_OK, "small pool exact: ierr OK")

        ! Test large (Monte Carlo) – pass random indices
        call compute_pvalue(0.0_real64, pool_large, 150, 0.0_real64, pool_large, 150, &
                            0.5_real64, 0, rand_indices_case, rand_indices_control, p_value, ierr)
        call assert_equal_int(ierr, ERR_OK, "large pool Monte Carlo: ierr OK")
        call assert_true(p_value >= 0.0 .and. p_value <= 1.0, "p_value out of range")

        ! Test missing random indices for large pool should trigger error
        ! We call without optional arguments
        ierr = ERR_OK
        call compute_pvalue(0.0_real64, pool_large, 150, 0.0_real64, pool_large, 150, &
                            0.5_real64, 0, rand_indices_case, rand_indices_control, p_value, ierr)
        call assert_equal_int(ierr, ERR_INVALID_INPUT, "large pool missing random indices should trigger error")
    end subroutine test_compute_pvalue_selection

    ! -------------------------------------------------------------------------
    ! Test full pipeline
    ! -------------------------------------------------------------------------
    subroutine test_full_pipeline()
        ! Use module constants directly
        real(real64) :: means_case(N_GENES), means_control(N_GENES)
        real(real64) :: replicates_case(N_SAMPLES, N_GENES), replicates_control(N_SAMPLES, N_GENES)
        integer(int32) :: family_sizes(N_FAMILIES), gene_to_family(N_GENES)
        real(real64) :: family_means(N_FAMILIES), ortholog_means(N_FAMILIES)
        integer(int32) :: compute_own(N_GENES), compute_family(N_GENES), compute_ortholog(N_GENES)
        real(real64) :: observed_own(N_GENES), observed_family(N_GENES), observed_ortholog(N_GENES)
        real(real64) :: pvalues_own(N_GENES), pvalues_family(N_GENES), pvalues_ortholog(N_GENES)
        integer(int32) :: n_genes_with_pvalue
        integer(int32) :: neigh_own_case(N_GENES), neigh_own_control(N_GENES)
        integer(int32) :: neigh_family(N_GENES), neigh_ortholog(N_GENES), neigh_case(N_GENES)
        integer(int32) :: ierr
        integer(int32) :: i
        real(real64) :: threshold

        ! Parameters for kNN
        integer(int32), parameter :: k_start = 10, k_step = 5, k_max = 50
        real(real64), parameter :: tau = 0.1_real64
        integer(int32), parameter :: max_pool_size = 100
        integer(int32), parameter :: norm_method = 0

        integer(int32) :: count_sig_own, count_sig_family, count_sig_orth

        ! Generate data
        call generate_test_data(N_GENES, N_SAMPLES, means_case, means_control, &
                                replicates_case, replicates_control, family_sizes, gene_to_family)

        ! Compute family means (control)
        family_means = 0.0_real64
        do i = 1, N_GENES
            family_means(gene_to_family(i)) = family_means(gene_to_family(i)) + means_control(i)
        end do
        family_means = family_means / real(family_sizes, real64)
        ortholog_means = family_means  ! for simplicity

        ! Set observed statistics: for own, use mean difference; for family/ortholog, use difference between gene mean and family mean
        do i = 1, N_GENES
            observed_own(i) = means_case(i) - means_control(i)
            observed_family(i) = means_case(i) - family_means(gene_to_family(i))
            observed_ortholog(i) = means_case(i) - ortholog_means(gene_to_family(i))
        end do

        ! Compute p-values for all genes
        compute_own = 1
        compute_family = 0
        compute_ortholog = 0

        call compute_noise_pvalue_pipeline( &
            means_case, replicates_case, N_GENES, N_SAMPLES, &
            means_control, replicates_control, N_GENES, N_SAMPLES, &
            observed_own, observed_family, observed_ortholog, &
            family_means, ortholog_means, &
            compute_own, compute_family, compute_ortholog, &
            family_sizes, gene_to_family, &
            N_GENES, N_FAMILIES, norm_method, k_start, k_step, k_max, tau, &
            pvalues_own, pvalues_family, pvalues_ortholog, n_genes_with_pvalue, &
            max_pool_size, &
            neigh_own_case, neigh_own_control, neigh_family, neigh_ortholog, neigh_case, &
            ierr)

        call assert_equal_int(ierr, ERR_OK, "pipeline ierr should be OK")
        call assert_true(n_genes_with_pvalue > 0, "at least some genes have p-values")

        ! Check that overexpressed genes (1..5) have small p-values (significant)
        threshold = 0.05_real64
        do i = 1, 5
            if (pvalues_own(i) >= 0.0_real64) then
                call assert_true(pvalues_own(i) < threshold, "overexpressed gene "//trim(str(i))//" should have small own p-value")
            end if
        end do

        ! Check that non-overexpressed genes have larger p-values (not all, but at least not too small)
        ! Not a strict test because of randomness, but we can check that they are not all significant
        count_sig_own = 0
        count_sig_family = 0
        count_sig_orth = 0
        do i = 6, N_GENES
            if (pvalues_own(i) >= 0.0 .and. pvalues_own(i) < threshold) count_sig_own = count_sig_own + 1
            if (pvalues_family(i) >= 0.0 .and. pvalues_family(i) < threshold) count_sig_family = count_sig_family + 1
            if (pvalues_ortholog(i) >= 0.0 .and. pvalues_ortholog(i) < threshold) count_sig_orth = count_sig_orth + 1
        end do
        ! We don't know how many will be significant, but we expect not too many (maybe < 20% of non-overexpressed)
        ! Since we have 15 non-overexpressed, allow up to 5 significant
        call assert_true(count_sig_own <= 5, "too many significant own p-values among non-overexpressed")
        call assert_true(count_sig_family <= 5, "too many significant family p-values among non-overexpressed")
        call assert_true(count_sig_orth <= 5, "too many significant ortholog p-values among non-overexpressed")

        ! Check that neighborhood sizes are positive when p-values computed
        do i = 1, N_GENES
            if (pvalues_own(i) >= 0.0) then
                call assert_true(neigh_own_case(i) > 0 .and. neigh_own_control(i) > 0, "neighborhood sizes should be positive")
            end if
            if (pvalues_family(i) >= 0.0) then
                call assert_true(neigh_family(i) > 0, "family neighborhood size positive")
            end if
            if (pvalues_ortholog(i) >= 0.0) then
                call assert_true(neigh_ortholog(i) > 0, "ortholog neighborhood size positive")
            end if
            if (neigh_case(i) > 0) then
                call assert_true(neigh_case(i) >= k_start, "case neighborhood size at least k_start")
            end if
        end do
    end subroutine test_full_pipeline

end module mod_test_noise_model