!> Unit test suite for tox_data_integration routine.
module mod_test_data_integration
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use tox_data_integration
    use tox_errors
    use f42_utils, only: above, below, init_random, shuffle_vector
    implicit none

    ! Abstract interface for all test procedures
    abstract interface
        subroutine test_interface()
        end subroutine test_interface
    end interface

    ! Type to hold test name and procedure pointer
    type :: test_case
        character(len=128) :: name
        procedure(test_interface), pointer, nopass :: test_proc => null()
    end type test_case

    real(real64), parameter :: TOL = 1d-12

contains

    !> Get array of all available tests.
    function get_all_tests() result(all_tests)
        type(test_case) :: all_tests(21)
        all_tests(1) = test_case("test_determine_shared_residual_range", test_determine_shared_residual_range)
        all_tests(2) = test_case("test_build_residual_histograms", test_build_residual_histograms)
        all_tests(3) = test_case("test_compute_divergence_per_reference_point", test_compute_divergence_per_reference_point)
        all_tests(4) = test_case("test_compute_weighted_global_divergence", test_compute_weighted_global_divergence)
        all_tests(5) = test_case("test_shuffle_reference_point_helper", test_shuffle_reference_point_helper)
        all_tests(6) = test_case("test_gjct_permutation_test", test_gjct_permutation_test)

        all_tests(7) = test_case("test_compute_gene_means_basic", test_compute_gene_means_basic)
        all_tests(8) = test_case("test_compute_gene_means_with_nan", test_compute_gene_means_with_nan)
        all_tests(9) = test_case("test_compute_gene_means_all_nan", test_compute_gene_means_all_nan)
        all_tests(10) = test_case("test_compute_gene_means_invalid_input", test_compute_gene_means_invalid_input)
        
        all_tests(11) = test_case("test_compute_residuals_basic", test_compute_residuals_basic)
        all_tests(12) = test_case("test_compute_residuals_with_nan", test_compute_residuals_with_nan)
        all_tests(13) = test_case("test_compute_residuals_all_nan", test_compute_residuals_all_nan)
        all_tests(14) = test_case("test_compute_residuals_invalid_input", test_compute_residuals_invalid_input)
        
        all_tests(15) = test_case("test_pool_means_alloc_basic", test_pool_means_alloc_basic)
        all_tests(16) = test_case("test_pool_means_alloc_with_nan", test_pool_means_alloc_with_nan)
        all_tests(17) = test_case("test_pool_means_alloc_single_study", test_pool_means_alloc_single_study)
        all_tests(18) = test_case("test_pool_means_alloc_invalid_input", test_pool_means_alloc_invalid_input)
        
        all_tests(19) = test_case("test_construct_neighborhoods_basic", test_construct_neighborhoods_basic)
        all_tests(20) = test_case("test_construct_neighborhoods_nan_means", test_construct_neighborhoods_nan_means)

        all_tests(21) = test_case("test_fjct", test_fjct)
        ! all_tests(22) = test_case("test_fjct_compute_contribution_scores", test_fjct_compute_contribution_scores)
    end function get_all_tests

    !> Run all tox_data_integration tests.
    subroutine run_all_tests_tox_data_integration
        type(test_case), allocatable :: all_tests(:)
        integer(int32) :: i

        all_tests = get_all_tests()

        do i = 1, size(all_tests)
            call all_tests(i)%test_proc()
            print "(' ',A,' passed.')", trim(all_tests(i)%name)
        end do
        print *, "All tox_data_integration tests passed successfully."
    end subroutine run_all_tests_tox_data_integration

    !> Run specific tox_data_integration tests by name.
    subroutine run_named_tests_tox_data_integration(test_names)
        character(len=*), intent(in) :: test_names(:)
        type(test_case), allocatable :: all_tests(:)
        integer(int32) :: i, j
        logical :: found

        all_tests = get_all_tests()

        do i = 1, size(test_names)
            found = .false.
            do j = 1, size(all_tests)
                if (trim(test_names(i)) == trim(all_tests(j)%name)) then
                    call all_tests(j)%test_proc()
                    print "(' ',A,' passed.')", trim(test_names(i))
                    found = .true.
                    exit
                end if
            end do
            if (.not. found) then
                print *, "Unknown test: ", trim(test_names(i))
            end if
        end do
    end subroutine run_named_tests_tox_data_integration

    subroutine test_fjct
        integer(int32), parameter :: n_reps_S1 = 3, n_reps_S2 = 4, n_neighbors = 5, n_points = 3, k_families = 2
        integer(int32), parameter :: n_bins = 4, n_genes_S1 = 10, n_genes_S2 = 10
        integer(int32), dimension(k_families), parameter :: included_families = [1, 2]
        integer(int32), dimension(n_genes_S1), parameter :: gene_to_family_S1 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 2]
        integer(int32), dimension(n_genes_S2), parameter :: gene_to_family_S2 = [1, 0, 0, 0, 0, 0, 0, 0, 3, 0]
        integer(int32), dimension(n_neighbors, n_points) :: neighborhood_genes_S1, neighborhood_genes_S2
        real(real64), dimension(n_reps_S1, n_neighbors, n_points) :: neighborhood_residuals_S1
        real(real64), dimension(n_reps_S2, n_neighbors, n_points) :: neighborhood_residuals_S2
        integer(int32), dimension(n_points) :: included_n_reps_S1, included_n_reps_S2, expected_included_n_reps_S1, expected_included_n_reps_S2
        real(real64), dimension(n_points) :: weights, expected_weights, js_divergences, expected_js_divergences
        integer(int32), dimension(k_families) :: total_included_n_reps
        real(real64), dimension(k_families) :: global_js_divergence, support_weights, expected_support_weights, contribution_scores, expected_contribution_scores
        real(real64) :: shared_residual_range
        integer(int32) :: ierr, expected_total_included_n_reps, family_idx

        ! ============================================================
        ! Test 1 — Simple symmetric case, no NaNs
        ! ============================================================
        !
        ! shared_residual_range = 2, M = 4 → bin width w = 1
        ! Bins: [-2,-1), [-1,0), [0,1), [1,2]
        !
        shared_residual_range = 2.0_real64

        neighborhood_residuals_S1 = 0.0_real64
        neighborhood_residuals_S1(:, 1,1) = [-2.0, -0.5, 0.2]
        neighborhood_residuals_S1(:, 1,3) = [2.5, -3.0, 1.2] ! (clamping applies -> [2,-2,1.2])
        neighborhood_genes_S1(:, 1) = [1, 2, 3, 4, 5]
        neighborhood_genes_S1(:, 2) = [2, 3, 4, 5, 6]
        neighborhood_genes_S1(:, 3) = [10, 2, 3, 4, 5]

        neighborhood_residuals_S2 = 0.0_real64
        neighborhood_residuals_S2(:, 1,1) = [-2.0_real64, ieee_value(1.0_real64, ieee_quiet_nan), -0.5_real64, 0.2_real64]
        neighborhood_residuals_S2(:, 3,3) = [0.4_real64, -0.1_real64, 0.0_real64, ieee_value(1.0_real64, ieee_quiet_nan)]
        neighborhood_genes_S2(:, 1) = [1, 2, 3, 4, 5]
        neighborhood_genes_S2(:, 2) = [2, 3, 4, 5, 6]
        neighborhood_genes_S2(:, 3) = [9, 2, 3, 4, 5]

        ! Values fall into bins:
        ! [-2,-1): -2, -1.2 → 2
        ! [-1,0): -0.5 → 1
        ! [0,1): 0.2, 0.9 → 2
        ! [1,2]: 1.7 → 1
        ! 
        ! point 2 — all zeros → all in bin [0,1)
        ! 
        ! point 3 — clamping:
        ! 2.5 → 2
        ! -3 → -2
        ! bins:
        ! [-2,-1): -2 → 1
        ! [-1,0): -0.1 → 1
        ! [0,1): 0.4, 0.0 → 2
        ! [1,2]: 1.2, 2 → 2

        family_idx = 1
        expected_included_n_reps_S1 = [3, 0, 0] ! first neighbor has the genes of fam 1
        ! -> pmf values point 1 = [1/3, 1/3, 1/3]
        expected_included_n_reps_S2 = [3, 0, 0] ! first neighbor has the genes of fam 1
        ! -> pmf values point 1 = [1/3, 1/3, 1/3]
        expected_total_included_n_reps = 6_int32
        expected_js_divergences = 0.0_real64

        expected_weights = [1.0_real64, 0.0_real64, 0.0_real64]
        call fjct_compute_jsd_alloc(&
            family_idx, gene_to_family_S1, gene_to_family_S2, n_genes_S1, n_genes_S2, neighborhood_residuals_S1, neighborhood_residuals_S2, &
            neighborhood_genes_S1, neighborhood_genes_S2, n_reps_S1, n_reps_S2, n_neighbors, n_points, n_bins, shared_residual_range, js_divergences, &
            included_n_reps_S1, included_n_reps_S2, total_included_n_reps(family_idx), global_js_divergence(family_idx), weights, ierr &
        )
        call assert_equal_int(ierr, ERR_OK, "test_fjct: Family 1: ierr should be OK")
        call assert_equal_int(total_included_n_reps(family_idx), expected_total_included_n_reps, "test_fjct: Family 1: total count of included replicates mismatch")
        call assert_equal_array_int(included_n_reps_S1, expected_included_n_reps_S1, n_points, "test_fjct: Family 1: included_n_reps_S1 mismatch")
        call assert_equal_array_int(included_n_reps_S2, expected_included_n_reps_S2, n_points, "test_fjct: Family 1: included_n_reps_S2 mismatch")
        call assert_equal_array_real(js_divergences, expected_js_divergences, n_points, TOL, "test_fjct: Family 1: js_divergences mismatch")
        call assert_equal_array_real(weights, expected_weights, n_points, TOL, "test_fjct: Family 1: weights mismatch")
        call assert_equal_real(global_js_divergence(family_idx), expected_weights(3) * expected_js_divergences(3), TOL, "test_fjct: Family 1: global_js_divergence mismatch")

        family_idx = 2
        expected_included_n_reps_S1 = [0, 0, 3] ! third neighbor has the gene of fam 2
        ! -> pmf values point 3 = [0, 0, 1.0]
        expected_included_n_reps_S2 = [0, 0, 0] ! no neighbor included
        expected_total_included_n_reps = 3_int32
        expected_js_divergences = [0.0_real64, 0.0_real64, 0.5_real64*log(2.0_real64)] ! s1 pmf val is 1.0, s2 is 0 -> mean is 0.5 -> jsd is 0.5*log(2)
        expected_weights = [0.0_real64, 0.0_real64, 1.0_real64]
        call fjct_compute_jsd_alloc(&
            family_idx, gene_to_family_S1, gene_to_family_S2, n_genes_S1, n_genes_S2, neighborhood_residuals_S1, neighborhood_residuals_S2, &
            neighborhood_genes_S1, neighborhood_genes_S2, n_reps_S1, n_reps_S2, n_neighbors, n_points, n_bins, shared_residual_range, js_divergences, &
            included_n_reps_S1, included_n_reps_S2, total_included_n_reps(family_idx), global_js_divergence(family_idx), weights, ierr &
        )
        call assert_equal_int(ierr, ERR_OK, "test_fjct: Family 2: ierr should be OK")
        call assert_equal_int(total_included_n_reps(family_idx), expected_total_included_n_reps, "test_fjct: Family 2: total count of included replicates mismatch")
        call assert_equal_array_int(included_n_reps_S1, expected_included_n_reps_S1, n_points, "test_fjct: Family 2: included_n_reps_S1 mismatch")
        call assert_equal_array_int(included_n_reps_S2, expected_included_n_reps_S2, n_points, "test_fjct: Family 2: included_n_reps_S2 mismatch")
        call assert_equal_array_real(js_divergences, expected_js_divergences, n_points, TOL, "test_fjct: Family 2: js_divergences mismatch")
        call assert_equal_array_real(weights, expected_weights, n_points, TOL, "test_fjct: Family 2: weights mismatch")
        call assert_equal_real(global_js_divergence(family_idx), expected_weights(3) * expected_js_divergences(3), TOL, "test_fjct: Family 2: global_js_divergence mismatch")


        expected_support_weights = [(real(total_included_n_reps(family_idx), kind=real64) / 9.0_real64, family_idx = 1, k_families)]
        expected_contribution_scores = [(expected_support_weights(family_idx) * global_js_divergence(family_idx), family_idx = 1, k_families)]
        call fjct_compute_contribution_scores(global_js_divergence, total_included_n_reps, k_families, support_weights, contribution_scores, ierr)

        call assert_equal_int(ierr, ERR_OK, "test_fjct: Contribution scores: ierr should be OK")
        call assert_equal_array_real(support_weights, expected_support_weights, k_families, TOL, "test_fjct: Contribution scores: support weights mismatch")
        call assert_equal_array_real(contribution_scores, expected_contribution_scores, k_families, TOL, "test_fjct: Contribution scores: support weights mismatch")
    end subroutine test_fjct

    subroutine test_gjct_permutation_test
        integer(int32), parameter :: n_reps_S1 = 4, n_reps_S2 = 3, n_neighbors = 1, n_points = 2, n_permutations = 2, n_bins = 4
        real(real64), dimension((n_reps_S1 + n_reps_S2) * n_neighbors * n_points ), target :: S_12, expected_S_12
        real(real64), dimension(:, :, :), pointer :: S1, S2
        real(real64), dimension(:, :), pointer :: tmp_pool
        real(real64), dimension(:), pointer :: S1_flat, S2_flat
        integer(int32), parameter :: random_seed = 666
        integer(int32) :: ierr, i_permutation
        real(real64) :: p_value, global_jsd_observed
        real(real64), dimension(n_permutations) :: jsd_null
        real(real64), dimension((n_reps_S1 + n_reps_S2) * n_neighbors), target :: tmp_pool_flat
        integer(int32), dimension(n_points, n_bins) :: tmp_counts
        real(real64), dimension(n_points, n_bins) :: tmp_pmf_S1
        real(real64), dimension(n_points, n_bins) :: tmp_pmf_S2
        integer(int32), dimension(n_points) :: tmp_included_n_reps_S1
        integer(int32), dimension(n_points) :: tmp_included_n_reps_S2
        real(real64), dimension(n_points) :: tmp_js_divergences
        real(real64), dimension(n_points) :: tmp_weights

        call init_random(random_seed)

        tmp_pool(1:n_reps_S1+n_reps_S2, 1:n_neighbors) => tmp_pool_flat

        ! ============================================================
        ! Test 1 — Test randomness with seed
        ! ============================================================
        !
        S_12 = [ 1,2,3,4,  5,6,-7,8, 2,-4,6,8,  1,3 ]

        ! Simulate two permutations 
        expected_S_12 = S_12
        S1(1:n_reps_S1, 1:n_neighbors, 1:n_points) => expected_S_12(1:n_reps_S1*n_neighbors*n_points)
        S1_flat(1:n_reps_S1 * n_neighbors * n_points) => expected_S_12(1:n_reps_S1*n_neighbors*n_points)
        S2(1:n_reps_S2, 1:n_neighbors, 1:n_points) => expected_S_12(n_reps_S1*n_neighbors*n_points+1:)
        S2_flat(1:n_reps_S2 * n_neighbors * n_points) => expected_S_12(n_reps_S1*n_neighbors*n_points+1:)
        do i_permutation = 1, n_permutations
            ! First point
            tmp_pool_flat(1:n_reps_S1*n_neighbors) = S1_flat(1:n_reps_S1*n_neighbors)
            tmp_pool_flat(n_reps_S1*n_neighbors+1:) = S2_flat(1:n_reps_S2*n_neighbors)
            call shuffle_vector(tmp_pool_flat)
            S1_flat(1:n_reps_S1*n_neighbors) = tmp_pool_flat(1:n_reps_S1*n_neighbors)
            S2_flat(1:n_reps_S2*n_neighbors) = tmp_pool_flat(n_reps_S1*n_neighbors+1:)

            ! Second point
            tmp_pool_flat(1:n_reps_S1*n_neighbors) = S1_flat(n_reps_S1*n_neighbors+1:)
            tmp_pool_flat(n_reps_S1*n_neighbors+1:) = S2_flat(n_reps_S2*n_neighbors+1:)
            call shuffle_vector(tmp_pool_flat)
            S1_flat(n_reps_S1*n_neighbors+1:) = tmp_pool_flat(1:n_reps_S1*n_neighbors)
            S2_flat(n_reps_S2*n_neighbors+1:) = tmp_pool_flat(n_reps_S1*n_neighbors+1:)
        end do

        ! for ifx: expected_shuffle = [2, 3, 6, 1, 5, 6, 1, -7, -4, 2, 4, 8, 3, 8]
       
        S1(1:n_reps_S1, 1:n_neighbors, 1:n_points) => S_12(1:n_reps_S1*n_neighbors*n_points)
        S2(1:n_reps_S2, 1:n_neighbors, 1:n_points) => S_12(n_reps_S1*n_neighbors*n_points+1:)
        
        global_jsd_observed = 0.0_real64
        call gjct_permutation_test(S1, S2, n_reps_S1, n_reps_S2, n_neighbors, n_points, global_jsd_observed, n_bins, shared_residual_range=10.0_real64, n_permutations=2_int32, jsd_null=jsd_null, p_value=p_value, ierr=ierr, random_seed=random_seed, tmp_pool=tmp_pool, tmp_counts=tmp_counts, tmp_pmf_S1=tmp_pmf_S1, tmp_pmf_S2=tmp_pmf_S2, tmp_included_n_reps_S1=tmp_included_n_reps_S1, tmp_included_n_reps_S2=tmp_included_n_reps_S2, tmp_js_divergences=tmp_js_divergences, tmp_weights=tmp_weights)
        call assert_equal_int(ierr, ERR_OK, "test_gjct_permutation_test: Test 1: unexpected error")

        call assert_equal_array_real(S_12, expected_S_12, size(S_12, kind=int32), 0.0_real64, "test_gjct_permutation_test: Test 1: concatenated S1, S2 does not match the expected permutation")

        ! ============================================================
        ! Test 2 — Test randomness without seed
        ! ============================================================
        !
        S_12 = [ 1,2,3,4,  5,6,-7,8, 2,-4,6,8,  1,3 ]
        S1(1:n_reps_S1, 1:n_neighbors, 1:n_points) => S_12(1:n_reps_S1*n_neighbors*n_points)
        S2(1:n_reps_S2, 1:n_neighbors, 1:n_points) => S_12(n_reps_S1*n_neighbors*n_points+1:)

        call init_random(random_seed)
        call gjct_permutation_test(S1, S2, n_reps_S1, n_reps_S2, n_neighbors, n_points, global_jsd_observed, n_bins, shared_residual_range=10.0_real64, n_permutations=2_int32, jsd_null=jsd_null, p_value=p_value, ierr=ierr, tmp_pool=tmp_pool, tmp_counts=tmp_counts, tmp_pmf_S1=tmp_pmf_S1, tmp_pmf_S2=tmp_pmf_S2, tmp_included_n_reps_S1=tmp_included_n_reps_S1, tmp_included_n_reps_S2=tmp_included_n_reps_S2, tmp_js_divergences=tmp_js_divergences, tmp_weights=tmp_weights)
        call assert_equal_int(ierr, ERR_OK, "test_gjct_permutation_test: Test 2: 1. call, unexpected error")
        call assert_equal_array_real(S_12, expected_S_12, size(S_12, kind=int32), 0.0_real64, "test_gjct_permutation_test: Test 2: 1. call, concatenated S1, S2 does not match the expected permutation")

        call gjct_permutation_test(S1, S2, n_reps_S1, n_reps_S2, n_neighbors, n_points, global_jsd_observed, n_bins, shared_residual_range=10.0_real64, n_permutations=2_int32, jsd_null=jsd_null, p_value=p_value, ierr=ierr, tmp_pool=tmp_pool, tmp_counts=tmp_counts, tmp_pmf_S1=tmp_pmf_S1, tmp_pmf_S2=tmp_pmf_S2, tmp_included_n_reps_S1=tmp_included_n_reps_S1, tmp_included_n_reps_S2=tmp_included_n_reps_S2, tmp_js_divergences=tmp_js_divergences, tmp_weights=tmp_weights)
        call assert_equal_int(ierr, ERR_OK, "test_gjct_permutation_test: Test 2: 2. call, unexpected error")
        call assert_true(any(S_12 /= expected_S_12), "test_gjct_permutation_test: Test 2: 2. call, concatenated S1, S2 should not match the expected permutation")


        ! ============================================================
        ! Test 3 — Test p value
        ! ============================================================
        !
        S_12 = [ 1,2,3,4,  5,6,-7,8, 2,-4,6,8,  1,3 ]
        S1(1:n_reps_S1, 1:n_neighbors, 1:n_points) => S_12(1:n_reps_S1*n_neighbors*n_points)
        S2(1:n_reps_S2, 1:n_neighbors, 1:n_points) => S_12(n_reps_S1*n_neighbors*n_points+1:)

        ! all should be greater or equal -> p_value=1
        global_jsd_observed = 0.0_real64

        call gjct_permutation_test_alloc(S1, S2, n_reps_S1, n_reps_S2, n_neighbors, n_points, global_jsd_observed, n_bins, shared_residual_range=10.0_real64, n_permutations=2_int32, jsd_null=jsd_null, p_value=p_value, ierr=ierr)
        call assert_equal_int(ierr, ERR_OK, "test_gjct_permutation_test: Test 3: unexpected error")

        call assert_equal_real(p_value, 1.0_real64, TOL, "test_gjct_permutation_test: Test 3: for zero observed jsd, p_value should be 1")

        ! no one should be greater or equal -> p_value=1/(n_permutations+1)=1/3
        global_jsd_observed = huge(0.0_real64)
        call gjct_permutation_test_alloc(S1, S2, n_reps_S1, n_reps_S2, n_neighbors, n_points, global_jsd_observed, n_bins, shared_residual_range=10.0_real64, n_permutations=2_int32, jsd_null=jsd_null, p_value=p_value, ierr=ierr)
        call assert_equal_int(ierr, ERR_OK, "test_gjct_permutation_test: Test 3: unexpected error")

        call assert_equal_real(p_value, 1.0_real64 / 3.0_real64, TOL, "test_gjct_permutation_test: Test 3: for max observed jsd, p_value should be 1/3")
    end subroutine test_gjct_permutation_test

    subroutine test_shuffle_reference_point_helper
        integer(int32), parameter :: n_reps_S1 = 4, n_reps_S2 = 3, n_neighbors = 2
        real(real64), dimension((n_reps_S1 + n_reps_S2) * n_neighbors ), target :: S_12, expected_S_12
        real(real64), dimension(:), pointer :: S1, S2
        real(real64), dimension((n_reps_S1 + n_reps_S2) * n_neighbors) :: pool_flat
        integer(int32), parameter :: random_seed = 666

        call init_random(random_seed)

        S_12 = [ 1,2,3,4,  5,6,-7,8, 2,-4,6,8,  1,3 ]
        S1(1:n_reps_S1*n_neighbors) => S_12(1:n_reps_S1*n_neighbors)
        S2(1:n_reps_S2*n_neighbors) => S_12(n_reps_S1*n_neighbors+1:)

        expected_S_12 = S_12
        call shuffle_vector(expected_S_12)
        ! for ifx: expected_S_12 = [-7, 4, 6, 3, 2, 8, 8, 5, -4, 6, 2, 3, 1, 1]

        call init_random(random_seed)

        call shuffle_reference_point_helper(S1, S2, n_reps_S1, n_reps_S2, n_neighbors, pool_flat)

        call assert_equal_array_real(pool_flat, S_12, size(S_12, kind=int32), TOL, "test_shuffle_reference_point_helper: pool should match concatenated S1, S2")
        call assert_equal_array_real(S_12, expected_S_12, size(S_12, kind=int32), TOL, "test_shuffle_reference_point_helper: concatenated S1, S2 does not match the expected permutation")
    end subroutine test_shuffle_reference_point_helper

    subroutine test_determine_shared_residual_range
        integer(int32), parameter :: n_reps_S1 = 4, n_reps_S2 = 3, n_neighbors = 2, n_points = 2
        real(real64), dimension(n_reps_S1, n_neighbors, n_points) :: S1
        real(real64), dimension(n_reps_S2, n_neighbors, n_points) :: S2
        real(real64) :: R
        integer(int32) :: ierr
        real(real64) :: q

        ! ============================================================
        ! Test 1 — Basic correctness with simple values
        ! ============================================================
        !
        S1 = reshape([ &
            1,2,3,4,  5,6,-7,8,  9,10,11,12, 1,1,1,1 ], shape(S1))
        S2 = reshape([ &
            2,-4,6,8,  1,3,5,7,  9,0,1,2 ], shape(S2))

        ! Sorted pool:
        !   [0, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 11, 12]
        !
        ! 95% quantile → 0.95 * (28 - 1) + 1 = 26.65
        ! sorted(26) = 10
        ! sorted(27) = 11
        !
        ! Expected R = 10 + 0.65 * (11-10) = 10.65
        !

        call determine_shared_residual_range_alloc(S1, S2, n_reps_S1, n_reps_S2, n_neighbors, n_points, R, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_determine_shared_residual_range: Test 1: ierr should be OK")
        call assert_equal_real(R, 10.65_real64, TOL, "test_determine_shared_residual_range: Test 1: R should be 10.65")

        ! ============================================================
        ! Test 2 — Custom quantile (50%)
        ! ============================================================
        !
        ! Median of sorted array above = 0.5 * (sorted(12) + sorted(13)) = 5.5
        !
        q = 50.0_real64
        call determine_shared_residual_range_alloc(S1, S2, n_reps_S1, n_reps_S2, n_neighbors, n_points, R, ierr, q)
        call assert_equal_int(ierr, ERR_OK, "test_determine_shared_residual_range: Test 2: ierr should be OK")
        call assert_equal_real(R, 4.0_real64, TOL, "test_determine_shared_residual_range: Test 2: R should be 4.0")

        ! ============================================================
        ! Test 3 — Quantile < 0 → error
        ! ============================================================
        q = below(0.0_real64)
        call determine_shared_residual_range_alloc(S1, S2, n_reps_S1, n_reps_S2, n_neighbors, n_points, R, ierr, q)
        call assert_equal_int(ierr, ERR_INVALID_INPUT, "test_determine_shared_residual_range: Test 3: ierr should be INVALID_INPUT")

        ! ============================================================
        ! Test 4 — Quantile > 100 → error
        ! ============================================================
        q = above(100.0_real64)
        call determine_shared_residual_range_alloc(S1, S2, n_reps_S1, n_reps_S2, n_neighbors, n_points, R, ierr, q)
        call assert_equal_int(ierr, ERR_INVALID_INPUT, "test_determine_shared_residual_range: Test 4: ierr should be INVALID_INPUT")

        ! ============================================================
        ! Test 5 — NaNs must be ignored
        ! ============================================================
        !
        ! Replace some values with NaN; remaining values should determine R.
        !
        S1 = reshape([ &
            -1, 2, 3, 4, &
            5, 6, -7, 8, &
            9, 10, -11, 12, &
            1,1,1,1 ], shape(S1))
        S2 = reshape([ &
            -1, 2, 3, 4, &
            5, 6, -7, 8, &
            9, 10, -11, 12 ], shape(S2))

        S1(1,1,1) = ieee_value(1.0_real64, ieee_quiet_nan)
        S2(3,2,2) = ieee_value(1.0_real64, ieee_quiet_nan)

        ! Pool now excludes two NaNs → 26 values
        ! sorted = [1, 1, 1, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12]
        ! 95% quantile → 0.95*(26-1)+1=24.75
        ! sorted(24) = 11
        ! sorted(25) = 11
        ! -> R = 11
        !
        call determine_shared_residual_range_alloc(S1, S2, n_reps_S1, n_reps_S2, n_neighbors, n_points, R, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_determine_shared_residual_range: Test 5: ierr should be OK")
        call assert_equal_real(R, 11.0_real64, TOL, "test_determine_shared_residual_range: Test 5: R should ignore NaNs")

        ! ============================================================
        ! Test 6 — All zeros
        ! ============================================================
        S1 = 0.0_real64
        S2 = 0.0_real64
        call determine_shared_residual_range_alloc(S1, S2, n_reps_S1, n_reps_S2, n_neighbors, n_points, R, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_determine_shared_residual_range: Test 6: ierr should be OK")
        call assert_equal_real(R, 0.0_real64, TOL, "test_determine_shared_residual_range: Test 6: R should be zero")

        ! ============================================================
        ! Test 7 — Single residual (n_reps_S1=1, n_reps_S2=1, n_neighbors=1, n_points=1)
        ! ============================================================
        S1 = 3.0_real64
        S2 = -4.0_real64
        ! sorted = [3, 4]
        ! rank = 0.95 * (2-1) + 1 = 1.95
        ! R = 3 + (4-3)*0.95 = 3.95
        call determine_shared_residual_range_alloc(S1, S2, 1_int32, 1_int32, 1_int32, 1_int32, R, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_determine_shared_residual_range: Test 6: ierr should be OK")
        call assert_equal_real(R, 3.95_real64, TOL, "test_determine_shared_residual_range: Test 7: R should be 3.95")

    end subroutine test_determine_shared_residual_range

    subroutine test_build_residual_histograms
        integer(int32), parameter :: n_reps = 3, n_neighbors = 2, n_points = 3
        integer(int32), parameter :: n_bins = 4
        real(real64), dimension(n_reps, n_neighbors, n_points) :: E
        real(real64), dimension(n_points, n_bins) :: pmf, expected_pmf
        integer(int32), dimension(n_points, n_bins) :: counts, expected_counts
        integer(int32), dimension(n_points) :: included
        real(real64) :: R
        integer(int32) :: ierr

        ! ============================================================
        ! Test 1 — Simple symmetric case, no NaNs
        ! ============================================================
        !
        ! R = 2, M = 4 → bin width w = 1
        ! Bins: [-2,-1), [-1,0), [0,1), [1,2]
        !
        R = 2.0_real64

        E(:, 1,1) = [-2.0, -0.5, 0.2]
        E(:, 2,1) = [1.7, 0.9, -1.2]
        E(:, :,2) = 0.0_real64
        E(:, 1,3) = [2.5, -3.0, 1.2] ! (clamping applies -> [2,-2,1.2])
        E(:, 2,3) = [0.4, -0.1, 0.0]

        call build_residual_histograms(E, n_reps, n_neighbors, n_points, R, n_bins, &
                                       counts, pmf, included, ierr)

        call assert_equal_int(ierr, ERR_OK, "test_build_residual_histograms: Test 1: ierr should be OK")

        ! point 1
        ! Values fall into bins:
        ! [-2,-1): -2, -1.2 → 2
        ! [-1,0): -0.5 → 1
        ! [0,1): 0.2, 0.9 → 2
        ! [1,2]: 1.7 → 1
        ! 
        ! point 2 — all zeros → all in bin [0,1)
        ! 
        ! point 3 — clamping:
        ! 2.5 → 2
        ! -3 → -2
        ! bins:
        ! [-2,-1): -2 → 1
        ! [-1,0): -0.1 → 1
        ! [0,1): 0.4, 0.0 → 2
        ! [1,2]: 1.2, 2 → 2
        expected_counts = reshape([&
            2, 0, 1,&
            1, 0, 1,&
            2, 6, 2,&
            1, 0, 2&
        ], [n_points, n_bins])
        expected_pmf = reshape([&
            0.3333333333333333_real64,  0.0_real64, 0.16666666666666666_real64,&
            0.16666666666666666_real64, 0.0_real64, 0.16666666666666666_real64,&
            0.3333333333333333_real64,  1.0_real64, 0.3333333333333333_real64,&
            0.16666666666666666_real64, 0.0_real64, 0.3333333333333333_real64&
        ], [n_points, n_bins])

        call assert_equal_array_int(counts, expected_counts, size(counts, kind=int32), "test_build_residual_histograms: Test 1: counts don't match")
        call assert_equal_array_real(pmf, expected_pmf, size(pmf, kind=int32), TOL, "test_build_residual_histograms: Test 1: pmf don't match")
        call assert_equal_int(included(1), 6, "test_build_residual_histograms: Test 1: included row 1")
        call assert_equal_int(included(2), 6, "test_build_residual_histograms: Test 1: included row 2")
        call assert_equal_int(included(3), 6, "test_build_residual_histograms: Test 1: included row 3")

        ! ============================================================
        ! Test 2 — NaNs must be ignored
        ! ============================================================
        E = 0.0_real64
        E(1,1,1) = ieee_value(1.0_real64, ieee_quiet_nan)
        E(3,1,1) = ieee_value(1.0_real64, ieee_quiet_nan)
        E(2,1,2) = ieee_value(1.0_real64, ieee_quiet_nan)
        E(3,2,3) = ieee_value(1.0_real64, ieee_quiet_nan)

        call build_residual_histograms(E, n_reps, n_neighbors, n_points, R, n_bins, &
                                       counts, pmf, included, ierr)

        call assert_equal_int(ierr, ERR_OK, "test_build_residual_histograms: Test 2: ierr should be OK")

        expected_counts = reshape([&
            0, 0, 0,&
            0, 0, 0,&
            4, 5, 5,&
            0, 0, 0&
        ], [n_points, n_bins])
        expected_pmf = reshape([&
            0.0_real64,  0.0_real64, 0.0_real64,&
            0.0_real64,  0.0_real64, 0.0_real64,&
            1.0_real64,  1.0_real64, 1.0_real64,&
            0.0_real64,  0.0_real64, 0.0_real64&
        ], [n_points, n_bins])

        call assert_equal_array_int(counts, expected_counts, size(counts, kind=int32), "test_build_residual_histograms: Test 2: counts don't match")
        call assert_equal_array_real(pmf, expected_pmf, size(pmf, kind=int32), TOL, "test_build_residual_histograms: Test 2: pmf don't match")
        call assert_equal_int(included(1), 4, "test_build_residual_histograms: Test 2: included row 1")
        call assert_equal_int(included(2), 5, "test_build_residual_histograms: Test 2: included row 2")
        call assert_equal_int(included(3), 5, "test_build_residual_histograms: Test 2: included row 3")

        ! ============================================================
        ! Test 3 — All NaN → pmf = 0, counts = 0, included = 0
        ! ============================================================
        E = ieee_value(1.0_real64, ieee_quiet_nan)

        call build_residual_histograms(E, n_reps, n_neighbors, n_points, R, n_bins, &
                                       counts, pmf, included, ierr)

        call assert_equal_int(ierr, ERR_OK, "test_build_residual_histograms: Test 3: ierr should be OK")

        expected_counts = 0.0_real64
        expected_pmf = 0.0_real64
        call assert_equal_array_int(counts, expected_counts, size(counts, kind=int32), "test_build_residual_histograms: Test 3: counts don't match")
        call assert_equal_array_real(pmf, expected_pmf, size(pmf, kind=int32), TOL, "test_build_residual_histograms: Test 3: pmf don't match")
        call assert_equal_int(included(1), 0, "test_build_residual_histograms: Test 2: included row 1")
        call assert_equal_int(included(2), 0, "test_build_residual_histograms: Test 2: included row 2")
        call assert_equal_int(included(3), 0, "test_build_residual_histograms: Test 2: included row 3")

        ! ============================================================
        ! Test 4 — Residuals exactly on boundaries
        ! ============================================================
        !
        ! R = 2, bins = 4, width = 1
        ! Values: -2, -1, 0, 1, 2
        !
        E = reshape([ -2.0, -1.0, 0.0, 1.0, 2.0, 0.0, &
                      -2.0, -1.0, 0.0, 1.0, 2.0, 0.0, &
                      -2.0, -1.0, 0.0, 1.0, 2.0, 0.0 ], shape(E))

        call build_residual_histograms(E, n_reps, n_neighbors, n_points, R, n_bins, &
                                       counts, pmf, included, ierr)

        call assert_equal_int(ierr, ERR_OK, "test_build_residual_histograms: Test 4: ierr should be OK")

        ! Expected binning:
        ! -2 → bin 1
        ! -1 → bin 2
        !  0 → bin 3
        !  1 → bin 4
        !  2 → bin 4 (right boundary included)
        !
        expected_counts = reshape([&
            1, 1, 1,&
            1, 1, 1,&
            2, 2, 2,&
            2, 2, 2&
        ], [n_points, n_bins])
        expected_pmf = reshape([&
            0.16666666666666666_real64, 0.16666666666666666_real64, 0.16666666666666666_real64,&
            0.16666666666666666_real64, 0.16666666666666666_real64, 0.16666666666666666_real64,&
            0.3333333333333333_real64, 0.3333333333333333_real64, 0.3333333333333333_real64,&
            0.3333333333333333_real64, 0.3333333333333333_real64, 0.3333333333333333_real64&
        ], [n_points, n_bins])

        call assert_equal_array_int(counts, expected_counts, size(counts, kind=int32), "test_build_residual_histograms: Test 3: counts don't match")
        call assert_equal_array_real(pmf, expected_pmf, size(pmf, kind=int32), TOL, "test_build_residual_histograms: Test 3: pmf don't match")
        call assert_equal_int(included(1), 6, "test_build_residual_histograms: Test 3: included row 1")
        call assert_equal_int(included(2), 6, "test_build_residual_histograms: Test 3: included row 2")
        call assert_equal_int(included(3), 6, "test_build_residual_histograms: Test 3: included row 3")

    end subroutine test_build_residual_histograms

    subroutine test_compute_divergence_per_reference_point
        integer(int32), parameter :: n_points = 3, n_bins = 4
        real(real64), dimension(n_points, n_bins) :: p, q
        real(real64), dimension(n_points) :: jsd, expected_jsd
        integer(int32) :: ierr
        real(real64) :: tol

        ! ============================================================
        ! Test 1 — Identical PMFs → JSD = 0
        ! ============================================================
        p = reshape([0.1, 0.2, 0.3, 0.4, &
                     0.25,0.25,0.25,0.25, &
                     1.0, 0.0, 0.0, 0.0], [n_points,n_bins])
        q = p

        call compute_divergence_per_reference_point(p, q, n_points, n_bins, jsd, ierr)

        call assert_equal_int(ierr, ERR_OK, "test_compute_divergence_per_reference_point: Test 1: ierr OK")

        expected_jsd = 0.0_real64
        call assert_equal_array_real(jsd, expected_jsd, size(jsd, kind=int32), TOL, "test_compute_divergence_per_reference_point: Test 1: identical PMFs → JSD=0")

        ! ============================================================
        ! Test 2 — Completely disjoint PMFs → JSD = log(2)
        ! ============================================================
        !
        ! p = [1,0,0,0]
        ! q = [0,1,0,0]
        !
        ! mix = [0.5,0.5,0,0]
        !
        ! KL(p||mix) = 1 * log(1/0.5) = log(2)
        ! KL(q||mix) = log(2)
        !
        ! JSD = 0.5*(log2 + log2) = log(2)
        !
        p = 0.0_real64
        q = 0.0_real64
        p(1,1) = 1.0
        q(1,2) = 1.0

        call compute_divergence_per_reference_point(p, q, n_points, n_bins, jsd, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_compute_divergence_per_reference_point: Test 2: ierr OK")

        expected_jsd = 0.0_real64
        expected_jsd(1) = log(2.0_real64)
        call assert_equal_array_real(jsd, expected_jsd, size(jsd, kind=int32), TOL, "test_compute_divergence_per_reference_point: Test 2: disjoint PMFs → JSD=log(2)")
        call assert_equal_real(jsd(2), 0.0_real64, TOL, "test_compute_divergence_per_reference_point: Test 2: rows 2,3 are zero PMFs → JSD=0")
        call assert_equal_real(jsd(3), 0.0_real64, TOL, "test_compute_divergence_per_reference_point: Test 2: rows 2,3 are zero PMFs → JSD=0")

        ! ============================================================
        ! Test 3 — Partially overlapping PMFs (analytic check)
        ! ============================================================
        !
        ! p = [0.5, 0.5, 0, 0]
        ! q = [0.0, 1.0, 0, 0]
        !
        ! mix = [0.25, 0.75, 0, 0]
        !
        ! KL(p||mix) = 0.5*log(0.5/0.25) + 0.5*log(0.5/0.75)
        !            = 0.5*log(2) + 0.5*log(2/3)
        !
        ! KL(q||mix) = 1.0*log(1.0/0.75)
        !
        ! JSD = 0.5*(KL_p + KL_q)
        !
        p = 0.0_real64
        q = 0.0_real64
        p(1,:) = [0.5, 0.5, 0.0, 0.0]
        q(1,:) = [0.0, 1.0, 0.0, 0.0]

        call compute_divergence_per_reference_point(p, q, n_points, n_bins, jsd, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_compute_divergence_per_reference_point: Test 3: ierr OK")

        ! Compute expected value analytically
        expected_jsd = 0.0_real64
        expected_jsd(1) = 0.5_real64 * ( &
            0.5_real64*log(2.0_real64) + 0.5_real64*log(2.0_real64/3.0_real64) &
          + log(1.0_real64/0.75_real64) )
        call assert_equal_array_real(jsd, expected_jsd, size(jsd, kind=int32), TOL, "test_compute_divergence_per_reference_point: Test 3: analytic partial-overlap JSD")

        ! ============================================================
        ! Test 4 — Zero-probability bins handled correctly
        ! ============================================================
        !
        ! p = [1,0,0,0]
        ! q = [0,0,1,0]
        !
        ! mix = [0.5,0,0.5,0]
        !
        ! KL(p||mix) = log(1/0.5) = log(2)
        ! KL(q||mix) = log(1/0.5) = log(2)
        !
        ! JSD = log(2)
        !
        p = 0.0_real64
        q = 0.0_real64
        p(1,1) = 1.0
        q(1,3) = 1.0

        call compute_divergence_per_reference_point(p, q, n_points, n_bins, jsd, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_compute_divergence_per_reference_point: Test 4: ierr OK")

        expected_jsd = 0.0_real64
        expected_jsd(1) = log(2.0_real64)
        call assert_equal_array_real(jsd, expected_jsd, size(jsd, kind=int32), TOL, "test_compute_divergence_per_reference_point: Test 4: zero-probability bins handled correctly")

        ! ============================================================
        ! Test 5 — Multiple neighbors, mixed patterns
        ! ============================================================
        p = reshape([ &
            0.2,0.3,0.5,&
            0.0,1.0,0.0,&
            0.0,0.0,0.25,&
            0.25,0.25,0.25&
        ], [n_points,n_bins])

        q = reshape([ &
            0.2,0.3,0.5,&
            0.0,0.0,1.0,&
            0.0,0.0,0.25,&
            0.25,0.25,0.25&
        ], [n_points,n_bins])

        call compute_divergence_per_reference_point(p, q, n_points, n_bins, jsd, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_compute_divergence_per_reference_point: Test 5: ierr OK")

        expected_jsd = 0.0_real64
        expected_jsd(2) = 0.5 * (1.0_real64 * log(1.0_real64 / 0.5_real64))
        expected_jsd(3) = 0.5 * (1.0_real64 * log(1.0_real64 / 0.5_real64))
        call assert_equal_array_real(jsd, expected_jsd, size(jsd, kind=int32), TOL, "test_compute_divergence_per_reference_point: Test 5: mixed patterns")
    end subroutine test_compute_divergence_per_reference_point

    subroutine test_compute_weighted_global_divergence
        integer(int32), parameter :: n_points = 4
        real(real64), dimension(n_points) :: jsd
        integer(int32), dimension(n_points) :: n1, n2
        real(real64), dimension(n_points) :: w, expected_weights
        real(real64) :: global_jsd
        integer(int32) :: ierr
        real(real64) :: expected, tol

        ! ============================================================
        ! Test 1 — Simple case: equal sample counts → uniform weights
        ! ============================================================
        !
        ! jsd = [0.1, 0.2, 0.3, 0.4]
        ! n1 = [5,5,5,5]
        ! n2 = [5,5,5,5]
        !
        ! n_j = [10,10,10,10]
        ! T = 40
        ! w = [0.25,0.25,0.25,0.25]
        !
        ! global_jsd = 0.25*(0.1+0.2+0.3+0.4) = 0.25
        !
        jsd = [0.1_real64, 0.2_real64, 0.3_real64, 0.4_real64]
        n1  = [5_int32,5_int32,5_int32,5_int32]
        n2  = [5_int32,5_int32,5_int32,5_int32]

        call compute_weighted_global_divergence(jsd, n_points, n1, n2, &
                                                global_jsd, w, ierr)

        expected_weights = 0.25_real64
        call assert_equal_int(ierr, ERR_OK, "test_compute_weighted_global_divergence: Test 1: ierr OK")

        ! Strangely this array comparison signals DENORMAL in gfortran if a test fails, seems to be an optimization bug
        call assert_equal_array_real(w, expected_weights, size(w, kind=int32), TOL, "test_compute_weighted_global_divergence: Test 1: uniform weights")
        call assert_equal_real(global_jsd, 0.25_real64, TOL, "test_compute_weighted_global_divergence: Test 1: global JSD")

        ! ============================================================
        ! Test 2 — Unequal sample counts → weighted average
        ! ============================================================
        !
        ! jsd = [1.0, 2.0, 3.0, 4.0]
        ! n1 = [10, 20, 30, 40]
        ! n2 = [ 0, 10, 10, 10]
        !
        ! n_j = [10, 30, 40, 50]
        ! T = 130
        ! w = [10/130, 30/130, 40/130, 50/130]
        !
        ! global_jsd = sum_j w(j)*jsd(j)
        !
        jsd = [1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64]
        n1  = [10_int32,20_int32,30_int32,40_int32]
        n2  = [ 0_int32,10_int32,10_int32,10_int32]

        call compute_weighted_global_divergence(jsd, n_points, n1, n2, &
                                                global_jsd, w, ierr)

        call assert_equal_int(ierr, ERR_OK, "test_compute_weighted_global_divergence: Test 2: ierr OK")
        
        expected_weights = [10.0_real64, 30.0_real64, 40.0_real64, 50.0_real64] / 130.0_real64
        call assert_equal_array_real(w, expected_weights, size(w, kind=int32), TOL, "test_compute_weighted_global_divergence: Test 2: wweights")

        expected = (1.0_real64*(10.0_real64/130.0_real64) + &
                    2.0_real64*(30.0_real64/130.0_real64) + &
                    3.0_real64*(40.0_real64/130.0_real64) + &
                    4.0_real64*(50.0_real64/130.0_real64))

        call assert_equal_real(global_jsd, expected, TOL, "test_compute_weighted_global_divergence: Test 2: weighted global JSD")

        ! ============================================================
        ! Test 3 — Some neighborhoods have zero samples → weight = 0
        ! ============================================================
        !
        ! jsd = [0.5, 1.0, 2.0, 4.0]
        ! n1 = [0, 10, 0, 5]
        ! n2 = [0,  0, 0, 5]
        !
        ! n_j = [0,10,0,10]
        ! T = 20
        ! w = [0, 10/20, 0, 10/20] = [0,0.5,0,0.5]
        !
        ! global_jsd = 0.5*1.0 + 0.5*4.0 = 2.5
        !
        jsd = [0.5_real64, 1.0_real64, 2.0_real64, 4.0_real64]
        n1  = [0_int32,10_int32,0_int32,5_int32]
        n2  = [0_int32, 0_int32,0_int32,5_int32]

        call compute_weighted_global_divergence(jsd, n_points, n1, n2, &
                                                global_jsd, w, ierr)

        call assert_equal_int(ierr, ERR_OK, "test_compute_weighted_global_divergence: Test 3: ierr OK")

        expected_weights = [0.0_real64,0.5_real64,0.0_real64,0.5_real64]
        call assert_equal_array_real(w, expected_weights, size(w, kind=int32), TOL, "test_compute_weighted_global_divergence: Test 3: weights with zero-sample neighborhoods")
        call assert_equal_real(global_jsd, 2.5_real64, TOL, "test_compute_weighted_global_divergence: Test 3: weighted global JSD")

        ! ============================================================
        ! Test 4 — All neighborhoods have zero samples → weights=0, global JSD=0
        ! ============================================================
        !
        ! jsd = [1,2,3,4]
        ! n1 = [0,0,0,0]
        ! n2 = [0,0,0,0]
        !
        ! n_j = [0,0,0,0]
        ! T = 0 → weights = 0
        ! global_jsd = 0
        !
        jsd = [1.0_real64,2.0_real64,3.0_real64,4.0_real64]
        n1  = 0
        n2  = 0

        call compute_weighted_global_divergence(jsd, n_points, n1, n2, &
                                                global_jsd, w, ierr)

        expected_weights = 0.0_real64
        call assert_equal_int(ierr, ERR_OK, "test_compute_weighted_global_divergence: Test 4: ierr OK")
        call assert_no_nan_real(w, size(w, kind=int32), "test_compute_weighted_global_divergence: Test 4: all weights not NaN")
        call assert_equal_array_real(w, expected_weights, size(w, kind=int32), TOL, "test_compute_weighted_global_divergence: Test 4: all weights zero")
        call assert_equal_real(global_jsd, 0.0_real64, TOL, "test_compute_weighted_global_divergence: Test 4: global JSD zero when no samples")

        ! ============================================================
        ! Test 5 — Mixed jsd values, mixed sample counts
        ! ============================================================
        !
        ! jsd = [0.0, 0.5, 1.0, 2.0]
        ! n1 = [5, 0, 10, 5]
        ! n2 = [5, 5,  0, 5]
        !
        ! n_j = [10,5,10,10]
        ! T = 35
        ! w = [10/35, 5/35, 10/35, 10/35]
        !
        ! global_jsd = sum_j w(j)*jsd(j)
        !
        jsd = [0.0_real64, 0.5_real64, 1.0_real64, 2.0_real64]
        n1  = [5_int32,0_int32,10_int32,5_int32]
        n2  = [5_int32,5_int32, 0_int32,5_int32]

        call compute_weighted_global_divergence(jsd, n_points, n1, n2, &
                                                global_jsd, w, ierr)

        call assert_equal_int(ierr, ERR_OK, "test_compute_weighted_global_divergence: Test 5: ierr OK")

        expected_weights = [10.0_real64, 5.0_real64, 10.0_real64, 10.0_real64] / 35.0_real64

        call assert_equal_array_real(w, expected_weights, size(w, kind=int32), TOL, "test_compute_weighted_global_divergence: Test 5: all weights zero")

        expected = (0.0_real64*(10.0_real64/35.0_real64) + &
                    0.5_real64*( 5.0_real64/35.0_real64) + &
                    1.0_real64*(10.0_real64/35.0_real64) + &
                    2.0_real64*(10.0_real64/35.0_real64))

        call assert_equal_real(global_jsd, expected, TOL, "test_compute_weighted_global_divergence: Test 5: weighted global JSD")

    end subroutine test_compute_weighted_global_divergence

    ! --------------------------------------------------------------------------
    ! Test Cases for compute_gene_means
    ! --------------------------------------------------------------------------

    ! Test case 1: Basic compute_gene_means functionality.
    subroutine test_compute_gene_means_basic()
        integer, parameter :: n_genes = 4, n_reps = 3
        real(real64) :: expr(n_reps, n_genes), means(n_genes)
        real(real64) :: expected_means(n_genes)
        integer(int32) :: ierr
        
        ! Test data
        expr = reshape([1.0, 2.0, 3.0,    &   ! Gene 1: mean = 2.0
                        4.0, 5.0, 6.0,    &   ! Gene 2: mean = 5.0
                        10.0, 20.0, 30.0, &   ! Gene 3: mean = 20.0
                        0.0, 0.0, 0.0],   &   ! Gene 4: mean = 0.0
                       [n_reps, n_genes])
        
        expected_means = [2.0, 5.0, 20.0, 0.0]
        
        call compute_gene_means(n_genes, n_reps, expr, means, ierr)
        
        call assert_equal_int(ierr, ERR_OK, "test_compute_gene_means_basic: should succeed")
        call assert_allclose_array_real(means, expected_means, n_genes, 0.0_real64, &
                                        TOL, "test_compute_gene_means_basic: means")
    end subroutine test_compute_gene_means_basic

    ! Test case 2: compute_gene_means with NaN values.
    subroutine test_compute_gene_means_with_nan()
        integer, parameter :: n_genes = 3, n_reps = 4
        real(real64) :: expr(n_reps, n_genes), means(n_genes)
        integer(int32) :: ierr
        
        expr(:, 1) = [1.0_real64, 2.0_real64, ieee_value(1.0_real64, ieee_quiet_nan), 3.0_real64]  ! mean = (1+2+3)/3 = 2.0
        expr(:, 2) = [ieee_value(1.0_real64, ieee_quiet_nan), 5.0_real64, 7.0_real64, 9.0_real64]  ! mean = (5+7+9)/3 = 7.0
        expr(:, 3) = [10.0, 20.0, 30.0, 40.0]  ! mean = 25.0
        
        call compute_gene_means(n_genes, n_reps, expr, means, ierr)
        
        call assert_equal_int(ierr, ERR_OK, "test_compute_gene_means_with_nan: should succeed")
        call assert_equal_real(means(1), 2.0_real64, TOL, "test_compute_gene_means_with_nan: gene 1 mean")
        call assert_equal_real(means(2), 7.0_real64, TOL, "test_compute_gene_means_with_nan: gene 2 mean")
        call assert_equal_real(means(3), 25.0_real64, TOL, "test_compute_gene_means_with_nan: gene 3 mean")
    end subroutine test_compute_gene_means_with_nan

    ! Test case 3: compute_gene_means with all NaN values for a gene.
    subroutine test_compute_gene_means_all_nan()
        integer, parameter :: n_genes = 2, n_reps = 3
        real(real64) :: expr(n_reps, n_genes), means(n_genes)
        integer(int32) :: ierr
        
        expr(:, 1) = [1.0, 2.0, 3.0]  ! Normal gene
        expr(:, 2) = ieee_value(1.0_real64, ieee_quiet_nan)  ! All NaN gene
        
        call compute_gene_means(n_genes, n_reps, expr, means, ierr)
        
        call assert_equal_int(ierr, ERR_OK, "test_compute_gene_means_all_nan: should succeed")
        call assert_equal_real(means(1), 2.0_real64, TOL, "test_compute_gene_means_all_nan: gene 1 mean")
        call assert_true(ieee_is_nan(means(2)), "test_compute_gene_means_all_nan: gene 2 should be NaN")
    end subroutine test_compute_gene_means_all_nan

    ! Test case 4: compute_gene_means with invalid input.
    subroutine test_compute_gene_means_invalid_input()
        integer, parameter :: n_genes = 0, n_reps = 3, n_genes_neg = -1
        real(real64) :: expr(3, 1), means(1)
        integer(int32) :: ierr
        
        ! Test with zero genes
        call compute_gene_means(n_genes, n_reps, expr, means, ierr)
        call assert_not_equal_int(ierr, ERR_OK, "test_compute_gene_means_invalid_input: zero genes should fail")
        
        ! Test with negative genes
        call compute_gene_means(n_genes_neg, n_reps, expr, means, ierr)
        call assert_not_equal_int(ierr, ERR_OK, "test_compute_gene_means_invalid_input: negative genes should fail")
        
        ! Test with zero replicates
        call compute_gene_means(n_genes, 0, expr, means, ierr)
        call assert_not_equal_int(ierr, ERR_OK, "test_compute_gene_means_invalid_input: zero replicates should fail")
    end subroutine test_compute_gene_means_invalid_input

    ! --------------------------------------------------------------------------
    ! Test Cases for compute_residuals
    ! --------------------------------------------------------------------------

    ! Test case 5: Basic compute_residuals functionality.
    subroutine test_compute_residuals_basic()
        integer, parameter :: n_genes = 4, n_reps = 3
        real(real64) :: expr(n_reps, n_genes), means(n_genes), resid(n_reps, n_genes)
        real(real64) :: expected_resid(n_reps, n_genes)
        integer(int32) :: ierr
        
        expr = reshape([1.0, 2.0, 3.0,    &   ! Gene 1
                        4.0, 5.0, 6.0,    &   ! Gene 2
                        10.0, 20.0, 30.0, &   ! Gene 3
                        0.0, 0.0, 0.0],   &   ! Gene 4
                       [n_reps, n_genes])
        
        means = [2.0, 5.0, 20.0, 0.0]
        expected_resid = reshape([-1.0, 0.0, 1.0,     &   ! Gene 1 residuals
                                  -1.0, 0.0, 1.0,     &   ! Gene 2 residuals
                                  -10.0, 0.0, 10.0,   &   ! Gene 3 residuals
                                  0.0, 0.0, 0.0],     &   ! Gene 4 residuals
                                 [n_reps, n_genes])
        
        call compute_residuals(n_genes, n_reps, expr, means, resid, ierr)
        
        call assert_equal_int(ierr, ERR_OK, "test_compute_residuals_basic: should succeed")
        call assert_allclose_array_real(reshape(resid, [n_reps*n_genes]), &
                                        reshape(expected_resid, [n_reps*n_genes]), &
                                        n_reps*n_genes, 0.0_real64, TOL, &
                                        "test_compute_residuals_basic: residuals")
    end subroutine test_compute_residuals_basic

    ! Test case 6: compute_residuals with NaN values.
    subroutine test_compute_residuals_with_nan()
        integer, parameter :: n_genes = 2, n_reps = 4
        real(real64) :: expr(n_reps, n_genes), means(n_genes), resid(n_reps, n_genes)
        integer(int32) :: ierr
        
        expr(:, 1) = [1.0_real64, 2.0_real64, ieee_value(1.0_real64, ieee_quiet_nan), 3.0_real64]
        expr(:, 2) = [ieee_value(1.0_real64, ieee_quiet_nan), 5.0_real64, 7.0_real64, 9.0_real64]
        means = [2.0_real64, 7.0_real64]
        
        call compute_residuals(n_genes, n_reps, expr, means, resid, ierr)
        
        call assert_equal_int(ierr, ERR_OK, "test_compute_residuals_with_nan: should succeed")
        ! Check specific values
        call assert_equal_real(resid(1, 1), -1.0_real64, TOL, "test_compute_residuals_with_nan: resid(1,1)")
        call assert_equal_real(resid(2, 1), 0.0_real64, TOL, "test_compute_residuals_with_nan: resid(2,1)")
        call assert_true(ieee_is_nan(resid(3, 1)), "test_compute_residuals_with_nan: resid(3,1) should be NaN")
        call assert_equal_real(resid(4, 1), 1.0_real64, TOL, "test_compute_residuals_with_nan: resid(4,1)")
        
        call assert_true(ieee_is_nan(resid(1, 2)), "test_compute_residuals_with_nan: resid(1,2) should be NaN")
        call assert_equal_real(resid(2, 2), -2.0_real64, TOL, "test_compute_residuals_with_nan: resid(2,2)")
        call assert_equal_real(resid(3, 2), 0.0_real64, TOL, "test_compute_residuals_with_nan: resid(3,2)")
        call assert_equal_real(resid(4, 2), 2.0_real64, TOL, "test_compute_residuals_with_nan: resid(4,2)")
    end subroutine test_compute_residuals_with_nan

    ! Test case 7: compute_residuals with all NaN values.
    subroutine test_compute_residuals_all_nan()
        integer, parameter :: n_genes = 2, n_reps = 3
        real(real64) :: expr(n_reps, n_genes), means(n_genes), resid(n_reps, n_genes)
        integer(int32) :: ierr
        
        expr(:, 1) = [1.0, 2.0, 3.0]
        expr(:, 2) = ieee_value(1.0_real64, ieee_quiet_nan)  ! All NaN
        means = [2.0, 0.0]  ! Second mean is irrelevant
        
        call compute_residuals(n_genes, n_reps, expr, means, resid, ierr)
        
        call assert_equal_int(ierr, ERR_OK, "test_compute_residuals_all_nan: should succeed")
        ! All residuals for gene 2 should be NaN
        call assert_true(all(ieee_is_nan(resid(:, 2))), "test_compute_residuals_all_nan: all residuals for NaN gene should be NaN")
    end subroutine test_compute_residuals_all_nan

    ! Test case 8: compute_residuals with invalid input.
    subroutine test_compute_residuals_invalid_input()
        integer, parameter :: n_genes = 0, n_reps = 3
        real(real64) :: expr(3, 1), means(1), resid(3, 1)
        integer(int32) :: ierr
        
        ! Test with zero genes
        call compute_residuals(n_genes, n_reps, expr, means, resid, ierr)
        call assert_not_equal_int(ierr, ERR_OK, "test_compute_residuals_invalid_input: zero genes should fail")

    end subroutine test_compute_residuals_invalid_input

    ! --------------------------------------------------------------------------
    ! Test Cases for pool_means_alloc
    ! --------------------------------------------------------------------------

    ! Test case 9: Basic pool_means_alloc functionality.
    subroutine test_pool_means_alloc_basic()
        integer, parameter :: n_genes_S1 = 5, n_genes_S2 = 5, n_points = 3
        real(real64) :: mean_S1(n_genes_S1), mean_S2(n_genes_S2), x_star(n_points)
        integer(int32) :: N_pool, ierr
        
        mean_S1 = [1.0, 3.0, 5.0, 7.0, 9.0]
        mean_S2 = [2.0, 4.0, 6.0, 8.0, 10.0]
        
        call pool_means_alloc(n_genes_S1, mean_S1, n_genes_S2, mean_S2, n_points, N_pool, x_star, ierr)
        
        call assert_equal_int(ierr, ERR_OK, "test_pool_means_alloc_basic: should succeed")
        call assert_equal_int(N_pool, 10, "test_pool_means_alloc_basic: N_pool should be 10")
        
        ! Check that x_star contains quantiles from pooled data
        ! Pooled data: [1,2,3,4,5,6,7,8,9,10]
        ! For n_points=3, quantiles at positions: 10/4=2.5, 20/4=5.0, 30/4=7.5
        ! Floored: 2, 5, 7 -> values: 2, 5, 7 -> interpolation to 3.25, 5.5 and 7.75
        call assert_equal_real(x_star(1), 3.25_real64, TOL, "test_pool_means_alloc_basic: first quantile")
        call assert_equal_real(x_star(2), 5.5_real64, TOL, "test_pool_means_alloc_basic: second quantile")
        call assert_equal_real(x_star(3), 7.75_real64, TOL, "test_pool_means_alloc_basic: third quantile")
    end subroutine test_pool_means_alloc_basic

    ! Test case 10: pool_means_alloc with NaN values.
    subroutine test_pool_means_alloc_with_nan()
        integer, parameter :: n_genes_S1 = 4, n_genes_S2 = 4, n_points = 2
        real(real64) :: mean_S1(n_genes_S1), mean_S2(n_genes_S2), x_star(n_points)
        integer(int32) :: N_pool, ierr
        
        mean_S1 = [1.0_real64, ieee_value(1.0_real64, ieee_quiet_nan), 3.0_real64, 5.0_real64]
        mean_S2 = [2.0_real64, 4.0_real64, ieee_value(1.0_real64, ieee_quiet_nan), 6.0_real64]
        
        call pool_means_alloc(n_genes_S1, mean_S1, n_genes_S2, mean_S2, n_points, N_pool, x_star, ierr)
        
        call assert_equal_int(ierr, ERR_OK, "test_pool_means_alloc_with_nan: should succeed")
        call assert_equal_int(N_pool, 6, "test_pool_means_alloc_with_nan: N_pool should exclude NaN values")
        
        ! Pooled data (excluding NaN): [1,2,3,4,5,6]
        ! Values: 2.666, 4.3333 -> interpolation
        call assert_equal_real(x_star(1), 2.0_real64 + 2.0_real64/3.0_real64, TOL, "test_pool_means_alloc_with_nan: first quantile")
        call assert_equal_real(x_star(2), 4.0_real64 + 1.0_real64/3.0_real64, TOL, "test_pool_means_alloc_with_nan: second quantile")
    end subroutine test_pool_means_alloc_with_nan

    ! Test case 11: pool_means_alloc with single study.
    subroutine test_pool_means_alloc_single_study()
        integer, parameter :: n_genes_S1 = 5, n_genes_S2 = 1, n_points = 3
        real(real64) :: mean_S1(n_genes_S1), mean_S2(1), x_star(n_points)
        integer(int32) :: N_pool, ierr
        
        mean_S1 = [1.0, 2.0, 3.0, 4.0, 5.0]
        mean_S2 = [0.0]  ! Dummy
        
        call pool_means_alloc(n_genes_S1, mean_S1, n_genes_S2, mean_S2, n_points, N_pool, x_star, ierr)
        
        call assert_equal_int(ierr, ERR_OK, "test_pool_means_alloc_single_study: should succeed")
        ! Does not work with a single study
    end subroutine test_pool_means_alloc_single_study

    ! Test case 12: pool_means_alloc with invalid input.
    subroutine test_pool_means_alloc_invalid_input()
        integer, parameter :: n_genes_S1 = 0, n_genes_S2 = 5, n_points = 3
        real(real64) :: mean_S1(1), mean_S2(5), x_star(n_points)
        integer(int32) :: N_pool, ierr
        
        mean_S2 = [1.0, 2.0, 3.0, 4.0, 5.0]
        
        ! Test with zero genes in S1
        call pool_means_alloc(n_genes_S1, mean_S1, n_genes_S2, mean_S2, n_points, N_pool, x_star, ierr)
        call assert_not_equal_int(ierr, ERR_OK, "test_pool_means_alloc_invalid_input: zero genes in S1 should fail")
        
        ! Test with zero points
        call pool_means_alloc(5, mean_S2, n_genes_S2, mean_S2, 0, N_pool, x_star, ierr)
        call assert_not_equal_int(ierr, ERR_OK, "test_pool_means_alloc_invalid_input: zero points should fail")
    end subroutine test_pool_means_alloc_invalid_input

    subroutine test_construct_neighborhoods_basic()
          integer(int32), parameter :: n_points    = 2
          integer(int32), parameter :: n_genes_S   = 5
          integer(int32), parameter :: n_reps_S    = 3
          integer(int32), parameter :: n_neighbors = 2

          integer(int32) :: ierr
          real(real64) :: x_star(n_points)
          real(real64) :: mean_S(n_genes_S)
          real(real64) :: resid_S(n_reps_S, n_genes_S)
          real(real64) :: tmp_distances(n_genes_S)
          integer(int32) :: tmp_distances_perm(n_genes_S)
          real(real64) :: neighborhood_residuals(n_reps_S, n_neighbors, n_points)
          integer(int32) :: neighborhood_indices(n_neighbors, n_points)

          ! -----------------------------
          ! Inputs
          ! -----------------------------
          x_star = [ 2.0_real64, 10.0_real64 ]
          mean_S = [ 1.0, 2.5, 9.0, 10.5, 20.0 ]

          resid_S = reshape([ &
              1.0,  2.0,  3.0,  4.0,  5.0, &   ! rep 1
              10.0,20.0,30.0,40.0,50.0, &   ! rep 2
              -1.0,-2.0,-3.0,-4.0,-5.0  &    ! rep 3
          ], shape(resid_S))

          ! -----------------------------
          ! Call routine
          ! -----------------------------
          call construct_neighborhoods( &
              n_points, x_star, n_genes_S, mean_S, n_reps_S, resid_S, &
              tmp_distances, tmp_distances_perm, &
              neighborhood_residuals, neighborhood_indices, &
              n_neighbors, ierr )

          call assert_equal_int(ierr, ERR_OK, "ierr must be ERR_OK")

          ! -----------------------------
          ! Expected neighbors
          !
          ! For x_star(1)=2.0:
          !   distances = [1.0, 0.5, 7.0, 8.5, 18.0]
          !   sorted → gene 2, gene 1
          !   tie-breaking: ascending gene index  <-- IMPORTANT
          !
          ! For x_star(2)=10.0:
          !   distances = [9.0, 7.5, 1.0, 0.5, 10.0]
          !   sorted → gene 4, gene 3
          ! -----------------------------

          call assert_equal_array_int( neighborhood_indices(:,1), [2,1], n_neighbors, &
              "test_construct_neighborhoods_basic: Incorrect neighborhood indices for point 1" )

          call assert_equal_array_int( neighborhood_indices(:,2), [4,3], n_neighbors, &
              "test_construct_neighborhoods_basic: Incorrect neighborhood indices for point 2" )

          ! -----------------------------
          ! Expected residuals
          ! -----------------------------
          call assert_equal_array_real( neighborhood_residuals(:,1,1), resid_S(:, 2), n_reps_S, TOL, &
              "test_construct_neighborhoods_basic: Incorrect residuals(:,1,1)" )

          call assert_equal_array_real( neighborhood_residuals(:,2,1), resid_S(:, 1), n_reps_S, TOL, &
              "test_construct_neighborhoods_basic: Incorrect residuals(:,2,1)" )

          call assert_equal_array_real( neighborhood_residuals(:,1,2), resid_S(:, 4), n_reps_S, TOL, &
              "test_construct_neighborhoods_basic: Incorrect residuals(:,1,2)" )

          call assert_equal_array_real( neighborhood_residuals(:,2,2), resid_S(:, 3), n_reps_S, TOL, &
              "test_construct_neighborhoods_basic: Incorrect residuals(:,2,2)" )

          call construct_neighborhoods( &
              0_int32, x_star, n_genes_S, mean_S, n_reps_S, resid_S, &
              tmp_distances, tmp_distances_perm, &
              neighborhood_residuals, neighborhood_indices, &
              n_neighbors, ierr )

          call assert_equal_int(ierr, ERR_EMPTY_INPUT, "ierr must be ERR_EMPTY_INPUT")

          call construct_neighborhoods( &
              n_points, x_star, 0_int32, mean_S, n_reps_S, resid_S, &
              tmp_distances, tmp_distances_perm, &
              neighborhood_residuals, neighborhood_indices, &
              n_neighbors, ierr )

          call assert_equal_int(ierr, ERR_EMPTY_INPUT, "ierr must be ERR_EMPTY_INPUT")

          call construct_neighborhoods( &
              n_points, x_star, n_genes_S, mean_S, 0_int32, resid_S, &
              tmp_distances, tmp_distances_perm, &
              neighborhood_residuals, neighborhood_indices, &
              n_neighbors, ierr )

          call assert_equal_int(ierr, ERR_EMPTY_INPUT, "ierr must be ERR_EMPTY_INPUT")

          call construct_neighborhoods( &
              n_points, x_star, n_genes_S, mean_S, n_reps_S, resid_S, &
              tmp_distances, tmp_distances_perm, &
              neighborhood_residuals, neighborhood_indices, &
              0_int32, ierr )

          call assert_equal_int(ierr, ERR_EMPTY_INPUT, "ierr must be ERR_EMPTY_INPUT")
    end subroutine test_construct_neighborhoods_basic

    subroutine test_construct_neighborhoods_nan_means()
          integer(int32), parameter :: n_points    = 1
          integer(int32), parameter :: n_genes_S   = 4
          integer(int32), parameter :: n_reps_S    = 2
          integer(int32), parameter :: n_neighbors = 2

          integer(int32) :: ierr
          real(real64) :: x_star(n_points)
          real(real64) :: mean_S(n_genes_S)
          real(real64) :: resid_S(n_reps_S, n_genes_S)
          real(real64) :: tmp_distances(n_genes_S)
          integer(int32) :: tmp_distances_perm(n_genes_S)
          real(real64) :: neighborhood_residuals(n_reps_S, n_neighbors, n_points)
          integer(int32) :: neighborhood_indices(n_neighbors, n_points)

          x_star = [ 5.0_real64 ]
          mean_S = [ 4.0_real64, ieee_value(1.0_real64, ieee_quiet_nan), 6.0_real64, ieee_value(1.0_real64, ieee_quiet_nan) ]

          resid_S = reshape([ &
              1.0, 2.0, 3.0, 4.0, &
              10.0,20.0,30.0,40.0 &
          ], shape(resid_S))

          call construct_neighborhoods( &
              n_points, x_star, n_genes_S, mean_S, n_reps_S, resid_S, &
              tmp_distances, tmp_distances_perm, &
              neighborhood_residuals, neighborhood_indices, &
              n_neighbors, ierr )

          call assert_equal_int(ierr, ERR_OK, "ierr must be ERR_OK")

          ! Only genes 1 and 3 are valid (non-NaN)
          call assert_equal_array_int( neighborhood_indices(:,1), [1,3], n_neighbors, &
              "test_construct_neighborhoods_nan_means: NaN mean handling incorrect" )
    end subroutine test_construct_neighborhoods_nan_means

end module mod_test_data_integration
