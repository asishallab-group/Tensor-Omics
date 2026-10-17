!> Unit test suite for the study-label Permutation-Test of the JS-Comp-Test
!| (tox_data_integration_js_comp_test_permutation_impl): one assignment of samples to studies
!| (k-NN and adaptive), the whole test with its observed value, null distribution and p-values,
!| the rerun of the parameter search, and the two sizing routines.
module mod_test_data_integration_js_comp_test_permutation
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: iso_c_binding, only: c_bool
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_nan
    use tox_errors
    use test_suite, only: test_case
    use f42_random_gsl, only: rng_t, create_rng, destroy_rng, random_uniform
    use tox_data_integration_preprocessing, only: compute_gene_means, compute_residuals, pool_means, &
                                                  construct_adaptive_neighborhoods
    use tox_data_integration_js_comp_test_run, only: run_js_comp_test, run_js_comp_test_adaptive
    use tox_data_integration_js_comp_test_fixed_k_search, only: run_js_comp_test_parameter_search
    use tox_data_integration_js_comp_test_adaptive_search, only: run_js_comp_test_adaptive_parameter_search
    use tox_data_integration_js_comp_test_fixed_k_search_impl, only: calc_js_comp_test_candidate_bounds
    use tox_data_integration_js_comp_test_adaptive_search_impl, only: calc_adaptive_js_comp_test_bounds
    use tox_data_integration_js_comp_test_impl, only: MAX_N_BINS, MAX_CANDIDATE_PAIRS, METHOD_JOIN_MIN
    use tox_data_integration_js_comp_test_permutation, only: run_js_comp_test_on_sample_assignment, &
                                                             run_js_comp_test_adaptive_on_sample_assignment, &
                                                             run_js_comp_test_label_permutation_test, &
                                                             run_js_comp_test_adaptive_label_permutation_test
    use tox_data_integration_js_comp_test_permutation_impl, only: calc_js_comp_test_permutation_bounds, &
                                                                  calc_adaptive_js_comp_test_permutation_bounds, &
                                                                  run_js_comp_test_on_sample_assignment_impl, &
                                                                  run_js_comp_test_adaptive_on_sample_assignment_impl, &
                                                                  run_js_comp_test_label_permutation_test_impl, &
                                                                  run_js_comp_test_adaptive_label_permutation_test_impl

    implicit none

    ! The fixed parameters every k-NN test uses, and the adaptive growth knobs.
    integer(int32), parameter :: N_POINTS_KNN = 10_int32, N_NEIGHBORS_KNN = 30_int32
    integer(int32), parameter :: K_START = 30_int32, K_STEP = 10_int32, K_MAX = 90_int32

contains

    !> Get array of all available tests.
    function get_all_tests_data_integration_js_comp_test_permutation() result(all_tests)
        type(test_case), allocatable :: all_tests(:)
        allocate (all_tests(17))

        all_tests(1) = test_case("test_knn_identity_matches_plain_pipeline", test_knn_identity_matches_plain_pipeline)
        all_tests(2) = test_case("test_adaptive_identity_matches_plain_pipeline", &
                                 test_adaptive_identity_matches_plain_pipeline)
        all_tests(3) = test_case("test_knn_swapping_two_studies_swaps_their_values", &
                                 test_knn_swapping_two_studies_swaps_their_values)
        all_tests(4) = test_case("test_knn_within_study_permutation_keeps_values", &
                                 test_knn_within_study_permutation_keeps_values)
        all_tests(5) = test_case("test_adaptive_swap_and_within_study_permutation", &
                                 test_adaptive_swap_and_within_study_permutation)
        all_tests(6) = test_case("test_knn_on_assignment_rejects_invalid_input", test_knn_on_assignment_rejects_invalid_input)
        all_tests(7) = test_case("test_adaptive_on_assignment_rejects_invalid_input", &
                                 test_adaptive_on_assignment_rejects_invalid_input)
        all_tests(8) = test_case("test_knn_label_test_reproducible_rows_match_level_a", &
                                 test_knn_label_test_reproducible_rows_match_level_a)
        all_tests(9) = test_case("test_adaptive_label_test_reproducible_rows_match_level_a", &
                                 test_adaptive_label_test_reproducible_rows_match_level_a)
        all_tests(10) = test_case("test_knn_label_test_null_and_batch", test_knn_label_test_null_and_batch)
        all_tests(11) = test_case("test_adaptive_label_test_null_and_batch", test_adaptive_label_test_null_and_batch)
        all_tests(12) = test_case("test_knn_rerun_parameter_search", test_knn_rerun_parameter_search)
        all_tests(13) = test_case("test_adaptive_rerun_parameter_search", test_adaptive_rerun_parameter_search)
        all_tests(14) = test_case("test_permutation_bounds_by_hand", test_permutation_bounds_by_hand)
        all_tests(15) = test_case("test_label_test_zero_permutations_direct_call", &
                                  test_label_test_zero_permutations_direct_call)
        all_tests(16) = test_case("test_adaptive_label_test_failed_permutations_bookkeeping", &
                                  test_adaptive_label_test_failed_permutations_bookkeeping)
        all_tests(17) = test_case("test_label_test_rejects_invalid_input", test_label_test_rejects_invalid_input)
    end function get_all_tests_data_integration_js_comp_test_permutation

    ! ========================================================================================
    ! Data and reference helpers
    ! ========================================================================================

    !> Test helper: an expression array of `size(n_reps)` studies, gene `g` centred at
    !| `2 + 6*(g-1)/(n_genes-1)` with Gaussian noise of SD `noise_sd(s)` (GSL generator, `seed`);
    !| replicate slots beyond `n_reps(s)` are NaN padding. `dyadic` rounds every value to a
    !| multiple of 1/8, so the gene-mean sums are exact.
    subroutine make_expression(n_genes, n_reps, noise_sd, seed, expression, dyadic)
        integer(int32), intent(in) :: n_genes
            !! Number of genes
        integer(int32), intent(in) :: n_reps(:)
            !! Number of samples of each study
        real(real64), intent(in) :: noise_sd(:)
            !! Noise standard deviation of each study
        integer(int32), intent(in) :: seed
            !! Seed of the generator
        real(real64), allocatable, intent(out) :: expression(:, :, :)
            !! The data, `(maxval(n_reps), n_genes, size(n_reps))`
        logical, intent(in), optional :: dyadic
            !! Round every value to a multiple of 1/8
        type(rng_t) :: rng
        integer(int32) :: ierr, i_study, i_gene, i_rep
        real(real64) :: u1, u2, value
        real(real64), parameter :: TWO_PI = 6.283185307179586_real64

        allocate (expression(maxval(n_reps), n_genes, size(n_reps)))
        expression = ieee_value(1.0_real64, ieee_quiet_nan)
        rng = create_rng(ierr, seed)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "make_expression: rng")
        do i_study = 1, size(n_reps)
            do i_gene = 1, n_genes
                do i_rep = 1, n_reps(i_study)
                    u1 = random_uniform(rng, ierr)
                    u2 = random_uniform(rng, ierr)
                    value = 2.0_real64 + 6.0_real64*real(i_gene - 1, real64)/real(n_genes - 1, real64) + &
                            noise_sd(i_study)*sqrt(-2.0_real64*log(1.0_real64 - u1))*cos(TWO_PI*u2)
                    if (present(dyadic)) then
                        if (dyadic) value = anint(8.0_real64*value)/8.0_real64
                    end if
                    expression(i_rep, i_gene, i_study) = value
                end do
            end do
        end do
        call destroy_rng(rng)
    end subroutine make_expression

    !> Test helper: every study's gene means and residuals, exactly as a caller computes them per study.
    subroutine reference_means_residuals(expression, gene_means, residuals)
        real(real64), intent(in) :: expression(:, :, :)
            !! The data
        real(real64), allocatable, intent(out) :: gene_means(:, :)
            !! Per-gene mean of every study
        real(real64), allocatable, intent(out) :: residuals(:, :, :)
            !! Residuals of every study
        integer(int32) :: n_reps, n_genes, n_studies, i_study, ierr

        n_reps = size(expression, 1, kind=int32)
        n_genes = size(expression, 2, kind=int32)
        n_studies = size(expression, 3, kind=int32)
        allocate (gene_means(n_genes, n_studies), residuals(n_reps, n_genes, n_studies))
        do i_study = 1, n_studies
            call compute_gene_means(n_genes, n_reps, expression(:, :, i_study), gene_means(:, i_study), ierr)
            call assert_equal_int(get_err_code(ierr), ERR_OK, "reference_means_residuals: means")
            call compute_residuals(n_genes, n_reps, expression(:, :, i_study), gene_means(:, i_study), &
                                   residuals(:, :, i_study), ierr)
            call assert_equal_int(get_err_code(ierr), ERR_OK, "reference_means_residuals: residuals")
        end do
    end subroutine reference_means_residuals

    !> Test helper: the plain k-NN pipeline a caller runs today: means/residuals per study, x_star
    !| from the pooled means, then run_js_comp_test without its permutation test.
    subroutine reference_knn(expression, n_points, n_neighbors, global_jsd, ierr)
        real(real64), intent(in) :: expression(:, :, :)
            !! The data
        integer(int32), intent(in) :: n_points
            !! Number of reference points
        integer(int32), intent(in) :: n_neighbors
            !! Number of neighbors
        real(real64), intent(out) :: global_jsd(:)
            !! Global JSD per study
        integer(int32), intent(out) :: ierr
            !! Error code of the last call
        real(real64), allocatable :: gene_means(:, :), residuals(:, :, :), x_star(:), range_low(:), range_high(:)
        real(real64), allocatable :: mean_occ(:), pmfs(:, :, :), mean_pmf(:, :), js(:, :), weights(:, :), p(:), p_obs(:)
        integer(int32), allocatable :: indices(:, :, :), ranges(:, :, :), n_bins(:), n_pooled(:), min_occ(:), max_occ(:)
        integer(int32), allocatable :: sturges(:), fd(:), counts(:, :, :), included(:, :), mean_counts(:, :)
        integer(int32), allocatable :: mean_included(:)
        logical(c_bool), allocatable :: occ_failed(:)
        integer(int32) :: n_studies, n_genes, n_reps, n_pool, n_bins_used

        call reference_means_residuals(expression, gene_means, residuals)
        n_reps = size(expression, 1, kind=int32)
        n_genes = size(expression, 2, kind=int32)
        n_studies = size(expression, 3, kind=int32)
        allocate (x_star(n_points))
        call pool_means(reshape(gene_means, [n_genes*n_studies]), n_genes*n_studies, n_points, n_pool, x_star, ierr)
        if (is_err(ierr)) return
        allocate (indices(n_neighbors, n_points, n_studies), ranges(2, n_points, n_studies), n_bins(n_points), &
                  n_pooled(n_points), min_occ(n_points), max_occ(n_points), sturges(n_points), fd(n_points), &
                  counts(n_points, MAX_N_BINS, n_studies), included(n_points, n_studies), &
                  mean_counts(n_points, MAX_N_BINS), mean_included(n_points), occ_failed(n_points), &
                  range_low(n_points), range_high(n_points), mean_occ(n_points), pmfs(n_points, MAX_N_BINS, n_studies), &
                  mean_pmf(n_points, MAX_N_BINS), js(n_points, n_studies), weights(n_points, n_studies), &
                  p(n_studies), p_obs(n_studies))
        call run_js_comp_test(n_studies, n_genes, n_reps, n_points, n_neighbors, gene_means, residuals, x_star, indices, &
                              ranges, n_bins, range_low, range_high, n_bins_used, occ_failed, n_pooled, min_occ, mean_occ, &
                              max_occ, sturges, fd, pmfs, counts, included, mean_pmf, mean_counts, mean_included, js, &
                              weights, global_jsd, p, p_obs, n_permutations=0_int32, ierr=ierr)
    end subroutine reference_knn

    !> Test helper: the plain adaptive pipeline a caller runs today: means/residuals per study,
    !| construct_adaptive_neighborhoods, then run_js_comp_test_adaptive without its permutation test.
    subroutine reference_adaptive(expression, k_start, k_step, k_max, global_jsd, n_points, ierr)
        real(real64), intent(in) :: expression(:, :, :)
            !! The data
        integer(int32), intent(in) :: k_start
            !! Growth knob
        integer(int32), intent(in) :: k_step
            !! Growth knob
        integer(int32), intent(in) :: k_max
            !! Growth knob
        real(real64), intent(out) :: global_jsd(:)
            !! Global JSD per study
        integer(int32), intent(out) :: n_points
            !! Number of neighborhoods built
        integer(int32), intent(out) :: ierr
            !! Error code of the last call
        real(real64), allocatable :: gene_means(:, :), residuals(:, :, :), x_star(:), dispersion(:), mads(:)
        real(real64), allocatable :: range_low(:), range_high(:), mean_occ(:), pmfs(:, :, :), mean_pmf(:, :)
        real(real64), allocatable :: js(:, :), weights(:, :), p(:), p_obs(:)
        integer(int32), allocatable :: ranges(:, :), n_per_point(:, :), stops(:), n_bins(:), n_pooled(:), min_occ(:)
        integer(int32), allocatable :: max_occ(:), sturges(:), fd(:), counts(:, :, :), included(:, :), mean_counts(:, :)
        integer(int32), allocatable :: mean_included(:), run_n_per_point(:, :)
        logical(c_bool), allocatable :: occ_failed(:)
        integer(int32) :: n_studies, n_genes, n_reps, n_pool, max_nn, status, n_bins_used

        call reference_means_residuals(expression, gene_means, residuals)
        n_reps = size(expression, 1, kind=int32)
        n_genes = size(expression, 2, kind=int32)
        n_studies = size(expression, 3, kind=int32)
        n_pool = n_genes*n_studies
        allocate (x_star(n_pool), ranges(2, n_pool), n_per_point(n_studies, n_pool), stops(n_pool), dispersion(n_pool), &
                  mads(n_pool))
        call construct_adaptive_neighborhoods(n_studies, n_genes, n_reps, gene_means, residuals, k_start, k_step, k_max, &
                                              n_points, x_star, ranges, n_per_point, stops, dispersion, mads, max_nn, &
                                              status, ierr=ierr)
        if (is_err(ierr)) return
        allocate (run_n_per_point(n_studies, n_points), n_bins(n_points), n_pooled(n_points), min_occ(n_points), &
                  max_occ(n_points), sturges(n_points), fd(n_points), counts(n_points, MAX_N_BINS, n_studies), &
                  included(n_points, n_studies), mean_counts(n_points, MAX_N_BINS), mean_included(n_points), &
                  occ_failed(n_points), range_low(n_points), range_high(n_points), mean_occ(n_points), &
                  pmfs(n_points, MAX_N_BINS, n_studies), mean_pmf(n_points, MAX_N_BINS), js(n_points, n_studies), &
                  weights(n_points, n_studies), p(n_studies), p_obs(n_studies))
        call run_js_comp_test_adaptive(n_studies, n_genes, n_reps, n_points, gene_means, residuals, ranges(:, 1:n_points), &
                                       run_n_per_point, n_bins, range_low, range_high, n_bins_used, occ_failed, n_pooled, &
                                       min_occ, mean_occ, max_occ, sturges, fd, pmfs, counts, included, mean_pmf, &
                                       mean_counts, mean_included, js, weights, global_jsd, p, p_obs, &
                                       n_permutations=0_int32, ierr=ierr)
    end subroutine reference_adaptive

    !> Test helper: run_js_comp_test_parameter_search with `n_bootstraps` and otherwise default
    !| options; returns the `(n_points, n_neighbors)` it selects.
    subroutine reference_fixed_k_search(gene_means, residuals, n_bootstraps, n_points, n_neighbors)
        real(real64), intent(in) :: gene_means(:, :)
            !! Per-gene mean of every study
        real(real64), intent(in) :: residuals(:, :, :)
            !! Residuals of every study
        integer(int32), intent(in) :: n_bootstraps
            !! Number of bootstraps
        integer(int32), intent(out) :: n_points
            !! Selected number of reference points
        integer(int32), intent(out) :: n_neighbors
            !! Selected number of neighbors
        integer(int32) :: n_studies, n_genes, n_reps, max_p, max_k, ierr, n_admissible
        integer(int32) :: trace_np(MAX_CANDIDATE_PAIRS), trace_nn(MAX_CANDIDATE_PAIRS)
        real(real64) :: trace_dmed(MAX_CANDIDATE_PAIRS), trace_dmax(MAX_CANDIDATE_PAIRS)
        logical(c_bool) :: plateau
        integer(int32), allocatable :: n_bins(:), t_bins(:, :), t_pooled(:, :), t_min(:, :), t_max(:, :), t_st(:, :)
        integer(int32), allocatable :: t_fd(:, :)
        logical(c_bool), allocatable :: t_failed(:, :)
        real(real64), allocatable :: r_low(:), r_high(:), ci(:, :), t_jsd(:, :), t_lo(:, :), t_up(:, :), t_w(:, :)
        real(real64), allocatable :: t_wr(:, :), t_ov(:, :), t_d(:, :), t_mean(:, :), t_rl(:, :), t_rh(:, :)

        n_genes = size(gene_means, 1, kind=int32)
        n_studies = size(gene_means, 2, kind=int32)
        n_reps = size(residuals, 1, kind=int32)
        call calc_js_comp_test_candidate_bounds(n_genes, max_p, max_k, ierr)
        allocate (n_bins(max_p), r_low(max_p), r_high(max_p), ci(2, n_studies), &
                  t_jsd(n_studies, MAX_CANDIDATE_PAIRS), t_lo(n_studies, MAX_CANDIDATE_PAIRS), &
                  t_up(n_studies, MAX_CANDIDATE_PAIRS), t_w(n_studies, MAX_CANDIDATE_PAIRS), &
                  t_wr(n_studies, MAX_CANDIDATE_PAIRS), t_ov(n_studies, MAX_CANDIDATE_PAIRS), &
                  t_d(n_studies, MAX_CANDIDATE_PAIRS), t_bins(max_p, MAX_CANDIDATE_PAIRS), &
                  t_failed(max_p, MAX_CANDIDATE_PAIRS), t_pooled(max_p, MAX_CANDIDATE_PAIRS), &
                  t_min(max_p, MAX_CANDIDATE_PAIRS), t_mean(max_p, MAX_CANDIDATE_PAIRS), &
                  t_max(max_p, MAX_CANDIDATE_PAIRS), t_st(max_p, MAX_CANDIDATE_PAIRS), t_fd(max_p, MAX_CANDIDATE_PAIRS), &
                  t_rl(max_p, MAX_CANDIDATE_PAIRS), t_rh(max_p, MAX_CANDIDATE_PAIRS))
        call run_js_comp_test_parameter_search(n_studies, n_genes, n_reps, gene_means, residuals, n_bootstraps, &
                                               METHOD_JOIN_MIN, max_p, max_k, n_points, n_neighbors, n_bins, r_low, &
                                               r_high, ci, plateau, n_admissible, trace_np, trace_nn, t_jsd, t_lo, t_up, &
                                               t_w, t_wr, t_ov, t_d, trace_dmed, trace_dmax, t_bins, t_failed, t_pooled, &
                                               t_min, t_mean, t_max, t_st, t_fd, t_rl, t_rh, ierr=ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "reference_fixed_k_search: ierr")
    end subroutine reference_fixed_k_search

    !> Test helper: run_js_comp_test_adaptive_parameter_search with `n_bootstraps` and otherwise
    !| default options and its own capacity; returns the knobs it selects.
    subroutine reference_adaptive_search(gene_means, residuals, n_bootstraps, k_start, k_step, k_max)
        real(real64), intent(in) :: gene_means(:, :)
            !! Per-gene mean of every study
        real(real64), intent(in) :: residuals(:, :, :)
            !! Residuals of every study
        integer(int32), intent(in) :: n_bootstraps
            !! Number of bootstraps
        integer(int32), intent(out) :: k_start
            !! Selected knob
        integer(int32), intent(out) :: k_step
            !! Selected knob
        integer(int32), intent(out) :: k_max
            !! Selected knob
        integer(int32) :: n_studies, n_genes, n_reps, cap, ierr, n_points, n_admissible, n_tried
        integer(int32), dimension(MAX_CANDIDATE_PAIRS) :: t_ks, t_kst, t_km, t_np, c_ks, c_np, c_status
        real(real64) :: trace_dmed(MAX_CANDIDATE_PAIRS), trace_dmax(MAX_CANDIDATE_PAIRS)
        logical(c_bool) :: plateau
        integer(int32), allocatable :: n_bins(:), t_bins(:, :), t_pooled(:, :), t_min(:, :), t_max(:, :), t_st(:, :)
        integer(int32), allocatable :: t_fd(:, :)
        logical(c_bool), allocatable :: t_failed(:, :)
        real(real64), allocatable :: r_low(:), r_high(:), ci(:, :), t_jsd(:, :), t_lo(:, :), t_up(:, :), t_w(:, :)
        real(real64), allocatable :: t_wr(:, :), t_ov(:, :), t_d(:, :), t_mean(:, :), t_rl(:, :), t_rh(:, :)

        n_genes = size(gene_means, 1, kind=int32)
        n_studies = size(gene_means, 2, kind=int32)
        n_reps = size(residuals, 1, kind=int32)
        call calc_adaptive_js_comp_test_bounds(n_genes, n_studies, cap, ierr)
        allocate (n_bins(cap), r_low(cap), r_high(cap), ci(2, n_studies), &
                  t_jsd(n_studies, MAX_CANDIDATE_PAIRS), t_lo(n_studies, MAX_CANDIDATE_PAIRS), &
                  t_up(n_studies, MAX_CANDIDATE_PAIRS), t_w(n_studies, MAX_CANDIDATE_PAIRS), &
                  t_wr(n_studies, MAX_CANDIDATE_PAIRS), t_ov(n_studies, MAX_CANDIDATE_PAIRS), &
                  t_d(n_studies, MAX_CANDIDATE_PAIRS), t_bins(cap, MAX_CANDIDATE_PAIRS), &
                  t_failed(cap, MAX_CANDIDATE_PAIRS), t_pooled(cap, MAX_CANDIDATE_PAIRS), &
                  t_min(cap, MAX_CANDIDATE_PAIRS), t_mean(cap, MAX_CANDIDATE_PAIRS), &
                  t_max(cap, MAX_CANDIDATE_PAIRS), t_st(cap, MAX_CANDIDATE_PAIRS), t_fd(cap, MAX_CANDIDATE_PAIRS), &
                  t_rl(cap, MAX_CANDIDATE_PAIRS), t_rh(cap, MAX_CANDIDATE_PAIRS))
        call run_js_comp_test_adaptive_parameter_search(n_studies, n_genes, n_reps, gene_means, residuals, n_bootstraps, &
                                                        METHOD_JOIN_MIN, cap, k_start, k_step, k_max, n_points, n_bins, &
                                                        r_low, r_high, ci, plateau, n_admissible, t_ks, t_kst, t_km, t_np, &
                                                        t_jsd, t_lo, t_up, t_w, t_wr, t_ov, t_d, trace_dmed, trace_dmax, &
                                                        t_bins, t_failed, t_pooled, t_min, t_mean, t_max, t_st, t_fd, &
                                                        t_rl, t_rh, n_tried, c_ks, c_np, c_status, ierr=ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "reference_adaptive_search: ierr")
    end subroutine reference_adaptive_search

    !> Test helper: the label draws of the Permutation-Test, re-implemented here (unbiased
    !| Fisher-Yates on one GSL generator seeded with `seed`): column `i` of `assignments` is draw `i`.
    subroutine redraw_assignments(seed, n_samples, n_draws, assignments)
        integer(int32), intent(in) :: seed
            !! Seed of the generator
        integer(int32), intent(in) :: n_samples
            !! Number of samples
        integer(int32), intent(in) :: n_draws
            !! Number of draws
        integer(int32), intent(out) :: assignments(n_samples, n_draws)
            !! The drawn permutations, one per column
        type(rng_t) :: rng
        integer(int32) :: ierr, i_draw, i, j, swapped
        integer(int32) :: perm(n_samples)

        rng = create_rng(ierr, seed)
        do i_draw = 1, n_draws
            perm = [(i, i=1, n_samples)]
            do i = n_samples, 2, -1
                j = min(1_int32 + int(random_uniform(rng, ierr)*real(i, real64), kind=int32), i)
                swapped = perm(i)
                perm(i) = perm(j)
                perm(j) = swapped
            end do
            assignments(:, i_draw) = perm
        end do
        call destroy_rng(rng)
    end subroutine redraw_assignments

    !> Test helper: level (a) k-NN through the plain wrapper, fixed parameters.
    subroutine knn_on_assignment(expression, assignment, global_jsd, ierr)
        real(real64), intent(in) :: expression(:, :, :)
            !! The data
        integer(int32), intent(in) :: assignment(:)
            !! The assignment
        real(real64), intent(out) :: global_jsd(:)
            !! Global JSD per study
        integer(int32), intent(out) :: ierr
            !! Error code
        integer(int32) :: used_np, used_nn

        call run_js_comp_test_on_sample_assignment(size(expression, 3, kind=int32), size(expression, 2, kind=int32), &
                                                   size(expression, 1, kind=int32), expression, &
                                                   size(assignment, kind=int32), assignment, N_POINTS_KNN, &
                                                   N_NEIGHBORS_KNN, global_jsd, used_np, used_nn, ierr=ierr)
    end subroutine knn_on_assignment

    !> Test helper: level (a) adaptive through the plain wrapper, fixed knobs.
    subroutine adaptive_on_assignment(expression, assignment, global_jsd, ierr)
        real(real64), intent(in) :: expression(:, :, :)
            !! The data
        integer(int32), intent(in) :: assignment(:)
            !! The assignment
        real(real64), intent(out) :: global_jsd(:)
            !! Global JSD per study
        integer(int32), intent(out) :: ierr
            !! Error code
        integer(int32) :: used_ks, used_kst, used_km, n_points

        call run_js_comp_test_adaptive_on_sample_assignment(size(expression, 3, kind=int32), &
                                                            size(expression, 2, kind=int32), &
                                                            size(expression, 1, kind=int32), expression, &
                                                            size(assignment, kind=int32), assignment, K_START, K_STEP, &
                                                            K_MAX, global_jsd, used_ks, used_kst, used_km, n_points, &
                                                            ierr=ierr)
    end subroutine adaptive_on_assignment

    !> Test helper: the identity assignment of `n` samples.
    pure function identity(n) result(perm)
        integer(int32), intent(in) :: n
            !! Number of samples
        integer(int32) :: perm(n)
        integer(int32) :: i

        perm = [(i, i=1, n)]
    end function identity

    !> Test helper: `.true.` iff `a` and `b` have the same size and are equal element by element; any NaN makes
    !| it `.false.` (unlike `same_values`, so a comparison of two failed results cannot pass).
    pure logical function equal_and_not_nan(a, b)
        real(real64), intent(in) :: a(:)
            !! First array
        real(real64), intent(in) :: b(:)
            !! Second array

        equal_and_not_nan = size(a) == size(b)
        if (.not. equal_and_not_nan) return
        equal_and_not_nan = all(a == b)
    end function equal_and_not_nan

    !> Test helper: `.true.` iff `a` and `b` are bit-identical (both NaN counts as equal).
    pure logical function same_values(a, b)
        real(real64), intent(in) :: a(:)
            !! First array
        real(real64), intent(in) :: b(:)
            !! Second array
        integer :: i

        same_values = size(a) == size(b)
        if (.not. same_values) return
        do i = 1, size(a)
            if (ieee_is_nan(a(i)) .and. ieee_is_nan(b(i))) cycle
            if (a(i) /= b(i)) then
                same_values = .false.
                return
            end if
        end do
    end function same_values

    ! ========================================================================================
    ! Level (a): one assignment
    ! ========================================================================================

    !> The identity assignment reproduces, bit for bit, the plain k-NN pipeline on the observed
    !| data, with unequal sample counts (NaN padding) and gene 7 missing in study 2.
    subroutine test_knn_identity_matches_plain_pipeline()
        character(*), parameter :: name = "test_knn_identity_matches_plain_pipeline: "
        real(real64), allocatable :: expression(:, :, :)
        real(real64) :: reference(3), result(3)
        integer(int32) :: ierr, used_np, used_nn, n_samples

        call make_expression(300_int32, [6_int32, 5_int32, 4_int32], [0.5_real64, 0.5_real64, 0.5_real64], 11_int32, &
                             expression)
        expression(:, 7, 2) = ieee_value(1.0_real64, ieee_quiet_nan)
        n_samples = 15_int32

        call reference_knn(expression, N_POINTS_KNN, N_NEIGHBORS_KNN, reference, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"reference ierr")
        call run_js_comp_test_on_sample_assignment(3_int32, 300_int32, 6_int32, expression, n_samples, &
                                                   identity(n_samples), N_POINTS_KNN, N_NEIGHBORS_KNN, result, used_np, &
                                                   used_nn, ierr=ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"ierr")
        call assert_true(same_values(result, reference), name//"identity equals the plain pipeline exactly")
        call assert_no_nan_real(result, 3_int32, name//"no NaN")
        call assert_equal_int(used_np, N_POINTS_KNN, name//"used_n_points")
        call assert_equal_int(used_nn, N_NEIGHBORS_KNN, name//"used_n_neighbors")
    end subroutine test_knn_identity_matches_plain_pipeline

    !> The adaptive identity assignment reproduces the plain construction + run exactly.
    subroutine test_adaptive_identity_matches_plain_pipeline()
        character(*), parameter :: name = "test_adaptive_identity_matches_plain_pipeline: "
        real(real64), allocatable :: expression(:, :, :)
        real(real64) :: reference(3), result(3)
        integer(int32) :: ierr, n_points_reference, used_ks, used_kst, used_km, n_points, n_samples

        call make_expression(300_int32, [6_int32, 5_int32, 4_int32], [0.5_real64, 0.5_real64, 0.5_real64], 12_int32, &
                             expression)
        expression(:, 7, 2) = ieee_value(1.0_real64, ieee_quiet_nan)
        n_samples = 15_int32

        call reference_adaptive(expression, K_START, K_STEP, K_MAX, reference, n_points_reference, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"reference ierr")
        call run_js_comp_test_adaptive_on_sample_assignment(3_int32, 300_int32, 6_int32, expression, n_samples, &
                                                            identity(n_samples), K_START, K_STEP, K_MAX, result, used_ks, &
                                                            used_kst, used_km, n_points, ierr=ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"ierr")
        call assert_true(same_values(result, reference), name//"identity equals the plain pipeline exactly")
        call assert_no_nan_real(result, 3_int32, name//"no NaN")
        call assert_equal_int(n_points, n_points_reference, name//"n_points")
        call assert_equal_int(used_ks, K_START, name//"used_k_start")
        call assert_equal_int(used_kst, K_STEP, name//"used_k_step")
        call assert_equal_int(used_km, K_MAX, name//"used_k_max")
    end subroutine test_adaptive_identity_matches_plain_pipeline

    !> Two studies of 4 samples each: the assignment that swaps the studies' samples swaps their
    !| values exactly (J1' = J2, J2' = J1).
    subroutine test_knn_swapping_two_studies_swaps_their_values()
        character(*), parameter :: name = "test_knn_swapping_two_studies_swaps_their_values: "
        real(real64), allocatable :: expression(:, :, :)
        real(real64) :: observed(2), swapped(2)
        integer(int32) :: ierr

        call make_expression(200_int32, [4_int32, 4_int32], [0.5_real64, 0.8_real64], 21_int32, expression, dyadic=.true.)
        call knn_on_assignment(expression, identity(8_int32), observed, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"observed ierr")
        call knn_on_assignment(expression, [5_int32, 6_int32, 7_int32, 8_int32, 1_int32, 2_int32, 3_int32, 4_int32], &
                               swapped, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"swapped ierr")
        call assert_true(observed(1) /= observed(2), name//"the two studies differ")
        call assert_true(swapped(1) == observed(2), name//"J1' == J2")
        call assert_true(swapped(2) == observed(1), name//"J2' == J1")
    end subroutine test_knn_swapping_two_studies_swaps_their_values

    !> An assignment that only reorders the samples inside each study leaves every value unchanged
    !| (dyadic data, so the gene means are exact in any order).
    subroutine test_knn_within_study_permutation_keeps_values()
        character(*), parameter :: name = "test_knn_within_study_permutation_keeps_values: "
        real(real64), allocatable :: expression(:, :, :)
        real(real64) :: observed(2), permuted(2)
        integer(int32) :: ierr

        call make_expression(200_int32, [4_int32, 4_int32], [0.5_real64, 0.8_real64], 22_int32, expression, dyadic=.true.)
        call knn_on_assignment(expression, identity(8_int32), observed, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"observed ierr")
        call knn_on_assignment(expression, [3_int32, 1_int32, 4_int32, 2_int32, 8_int32, 7_int32, 6_int32, 5_int32], &
                               permuted, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"permuted ierr")
        call assert_true(same_values(permuted, observed), name//"values unchanged")
    end subroutine test_knn_within_study_permutation_keeps_values

    !> The adaptive counterpart of the two hand-checkable cases: a swap swaps the values, a
    !| within-study reorder keeps them. The swap uses continuous data: the construction sorts all
    !| studies' means pooled, and dyadic means would tie across studies, where the order among tied
    !| entries follows their (swapped) pooled positions.
    subroutine test_adaptive_swap_and_within_study_permutation()
        character(*), parameter :: name = "test_adaptive_swap_and_within_study_permutation: "
        real(real64), allocatable :: expression(:, :, :)
        real(real64) :: observed(2), swapped(2), permuted(2)
        integer(int32) :: ierr

        call make_expression(200_int32, [4_int32, 4_int32], [0.5_real64, 0.8_real64], 23_int32, expression)
        call adaptive_on_assignment(expression, identity(8_int32), observed, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"observed ierr")
        call adaptive_on_assignment(expression, [5_int32, 6_int32, 7_int32, 8_int32, 1_int32, 2_int32, 3_int32, 4_int32], &
                                    swapped, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"swapped ierr")
        call assert_true(observed(1) /= observed(2), name//"the two studies differ")
        call assert_true(swapped(1) == observed(2), name//"J1' == J2")
        call assert_true(swapped(2) == observed(1), name//"J2' == J1")

        call make_expression(200_int32, [4_int32, 4_int32], [0.5_real64, 0.8_real64], 24_int32, expression, dyadic=.true.)
        call adaptive_on_assignment(expression, identity(8_int32), observed, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"dyadic observed ierr")
        call adaptive_on_assignment(expression, [2_int32, 4_int32, 1_int32, 3_int32, 6_int32, 5_int32, 8_int32, 7_int32], &
                                    permuted, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"permuted ierr")
        call assert_true(same_values(permuted, observed), name//"within-study reorder keeps the values")
    end subroutine test_adaptive_swap_and_within_study_permutation

    !> Invalid inputs of the k-NN level (a): through the plain wrapper (runtime errors carry no
    !| position there) and through the implementation (each names its argument).
    subroutine test_knn_on_assignment_rejects_invalid_input()
        character(*), parameter :: name = "test_knn_on_assignment_rejects_invalid_input: "
        real(real64), allocatable :: expression(:, :, :), empty_study(:, :, :)
        real(real64) :: result(3)
        integer(int32) :: ierr, used_np, used_nn, max_p, max_k, top_k

        call make_expression(300_int32, [4_int32, 4_int32, 4_int32], [0.5_real64, 0.5_real64, 0.5_real64], 31_int32, &
                             expression)
        empty_study = expression
        empty_study(:, :, 2) = ieee_value(1.0_real64, ieee_quiet_nan)

        ! Plain wrapper.
        call knn_on_assignment(expression, [1_int32, 1_int32, 3_int32, 4_int32, 5_int32, 6_int32, 7_int32, 8_int32, &
                                            9_int32, 10_int32, 11_int32, 12_int32], result, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, name//"repeated value rejected", arg_pos=0_int32)
        call assert_true(all(ieee_is_nan(result)), name//"result NaN on error")
        call knn_on_assignment(expression, identity(11_int32), result, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, name//"n_samples below the sample count rejected", arg_pos=0_int32)
        call knn_on_assignment(expression, identity(13_int32), result, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, name//"n_samples above the sample count rejected", arg_pos=0_int32)
        call knn_on_assignment(empty_study, identity(8_int32), result, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, name//"study without samples rejected", arg_pos=0_int32)
        call knn_on_assignment(expression, [0_int32, 2_int32, 3_int32, 4_int32, 5_int32, 6_int32, 7_int32, 8_int32, &
                                            9_int32, 10_int32, 11_int32, 12_int32], result, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, name//"value 0 rejected by the wrapper", arg_pos=6_int32)

        ! Implementation: positions, and the rerun capacity checks.
        call calc_js_comp_test_permutation_bounds(300_int32, N_POINTS_KNN, N_NEIGHBORS_KNN, max_n_points_candidate=max_p, &
                                                  max_n_neighbors_candidate=max_k, n_bootstrapping_top_k_jsds=top_k, &
                                                  ierr=ierr)
        call direct_knn_on_assignment(3_int32, 300_int32, 4_int32, expression, 12_int32, &
                                      [2_int32, 2_int32, 3_int32, 4_int32, 5_int32, 6_int32, 7_int32, 8_int32, 9_int32, &
                                       10_int32, 11_int32, 12_int32], N_POINTS_KNN, N_NEIGHBORS_KNN, max_p, max_k, top_k, &
                                      result, used_np, used_nn, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, name//"impl: repeated value", arg_pos=6_int32)
        call direct_knn_on_assignment(3_int32, 300_int32, 4_int32, expression, 11_int32, identity(11_int32), N_POINTS_KNN, &
                                      N_NEIGHBORS_KNN, max_p, max_k, top_k, result, used_np, used_nn, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, name//"impl: n_samples", arg_pos=5_int32)
        call direct_knn_on_assignment(3_int32, 300_int32, 4_int32, empty_study, 8_int32, identity(8_int32), N_POINTS_KNN, &
                                      N_NEIGHBORS_KNN, max_p, max_k, top_k, result, used_np, used_nn, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, name//"impl: empty study", arg_pos=4_int32)
        call direct_knn_on_assignment(3_int32, 300_int32, 4_int32, expression, 12_int32, identity(12_int32), N_POINTS_KNN, &
                                      N_NEIGHBORS_KNN, max_p - 1_int32, max_k, top_k, result, used_np, used_nn, ierr, &
                                      rerun_parameter_search=.true._c_bool)
        call assert_err(ierr, ERR_INVALID_INPUT, name//"impl: point capacity below the search's", arg_pos=9_int32)
        call direct_knn_on_assignment(3_int32, 300_int32, 4_int32, expression, 12_int32, identity(12_int32), 4_int32, &
                                      4_int32, max_p, 3_int32, top_k, result, used_np, used_nn, ierr, &
                                      rerun_parameter_search=.true._c_bool)
        call assert_err(ierr, ERR_INVALID_INPUT, name//"impl: neighbor capacity below the search's", arg_pos=10_int32)
        ! The same small capacities without the rerun are fine.
        call direct_knn_on_assignment(3_int32, 300_int32, 4_int32, expression, 12_int32, identity(12_int32), N_POINTS_KNN, &
                                      N_NEIGHBORS_KNN, N_POINTS_KNN, N_NEIGHBORS_KNN, top_k, result, used_np, used_nn, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"impl: exact capacities without rerun")
    end subroutine test_knn_on_assignment_rejects_invalid_input

    !> Invalid inputs of the adaptive level (a), and a construction with no neighborhood
    !| (`k_start` beyond the data) or more than the capacity.
    subroutine test_adaptive_on_assignment_rejects_invalid_input()
        character(*), parameter :: name = "test_adaptive_on_assignment_rejects_invalid_input: "
        real(real64), allocatable :: expression(:, :, :), empty_study(:, :, :)
        real(real64) :: result(3)
        integer(int32) :: ierr, used_ks, used_kst, used_km, n_points, cap, top_k, search_cap

        call make_expression(300_int32, [4_int32, 4_int32, 4_int32], [0.5_real64, 0.5_real64, 0.5_real64], 32_int32, &
                             expression)
        empty_study = expression
        empty_study(:, :, 3) = ieee_value(1.0_real64, ieee_quiet_nan)
        call calc_adaptive_js_comp_test_permutation_bounds(300_int32, 3_int32, K_START, max_n_points_candidate=cap, &
                                                           n_bootstrapping_top_k_jsds=top_k, ierr=ierr)

        call direct_adaptive_on_assignment(3_int32, 300_int32, 4_int32, expression, 12_int32, &
                                           [1_int32, 2_int32, 3_int32, 4_int32, 5_int32, 6_int32, 7_int32, 8_int32, &
                                            9_int32, 10_int32, 12_int32, 12_int32], K_START, K_STEP, K_MAX, cap, top_k, &
                                           result, used_ks, used_kst, used_km, n_points, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, name//"impl: repeated value", arg_pos=6_int32)
        call direct_adaptive_on_assignment(3_int32, 300_int32, 4_int32, expression, 13_int32, identity(13_int32), K_START, &
                                           K_STEP, K_MAX, cap, top_k, result, used_ks, used_kst, used_km, n_points, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, name//"impl: n_samples", arg_pos=5_int32)
        call direct_adaptive_on_assignment(3_int32, 300_int32, 4_int32, empty_study, 8_int32, identity(8_int32), K_START, &
                                           K_STEP, K_MAX, cap, top_k, result, used_ks, used_kst, used_km, n_points, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, name//"impl: empty study", arg_pos=4_int32)
        call calc_adaptive_js_comp_test_bounds(300_int32, 3_int32, search_cap, ierr)
        call direct_adaptive_on_assignment(3_int32, 300_int32, 4_int32, expression, 12_int32, identity(12_int32), K_START, &
                                           K_STEP, K_MAX, search_cap - 1_int32, top_k, result, used_ks, used_kst, used_km, &
                                           n_points, ierr, rerun_parameter_search=.true._c_bool)
        call assert_err(ierr, ERR_INVALID_INPUT, name//"impl: capacity below the search's", arg_pos=10_int32)
        ! Capacity 1: the construction builds more neighborhoods.
        call direct_adaptive_on_assignment(3_int32, 300_int32, 4_int32, expression, 12_int32, identity(12_int32), K_START, &
                                           K_STEP, K_MAX, 1_int32, top_k, result, used_ks, used_kst, used_km, n_points, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, name//"impl: more neighborhoods than the capacity", arg_pos=0_int32)
        call assert_true(n_points > 1_int32, name//"n_points reports the construction's count")
        call assert_true(all(ieee_is_nan(result)), name//"result NaN on error")

        ! k_start beyond the 900 pooled means: the construction builds nothing (plain wrapper).
        call run_js_comp_test_adaptive_on_sample_assignment(3_int32, 300_int32, 4_int32, expression, 12_int32, &
                                                            identity(12_int32), 1000_int32, 10_int32, 1000_int32, result, &
                                                            used_ks, used_kst, used_km, n_points, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, name//"no neighborhood built", arg_pos=0_int32)
        call assert_equal_int(n_points, 0_int32, name//"n_points 0")
        ! Plain wrapper: a repeated value.
        call adaptive_on_assignment(expression, [1_int32, 2_int32, 3_int32, 4_int32, 5_int32, 6_int32, 7_int32, 8_int32, &
                                                 9_int32, 9_int32, 11_int32, 12_int32], result, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, name//"wrapper: repeated value", arg_pos=0_int32)
    end subroutine test_adaptive_on_assignment_rejects_invalid_input

    ! ========================================================================================
    ! Level (b): the whole test
    ! ========================================================================================

    !> k-NN: a fixed seed reproduces the test; the p-values lie in `[1/(n+1), 1]`; the observed
    !| value is level (a) on the identity; every row of `jsd_null` is level (a) on the assignment
    !| the test drew (re-drawn here with the same Fisher-Yates and seed).
    subroutine test_knn_label_test_reproducible_rows_match_level_a()
        character(*), parameter :: name = "test_knn_label_test_reproducible_rows_match_level_a: "
        integer(int32), parameter :: n = 6_int32, seed = 7_int32
        real(real64), allocatable :: expression(:, :, :)
        real(real64) :: observed(3), p(3), null(n, 3), observed2(3), p2(3), null2(n, 3), row(3)
        integer(int32) :: ierr, n_failed, n_failed2, assignments(15, n), i

        call make_expression(300_int32, [6_int32, 5_int32, 4_int32], [0.5_real64, 0.5_real64, 0.5_real64], 41_int32, &
                             expression)
        call run_js_comp_test_label_permutation_test(3_int32, 300_int32, 6_int32, expression, N_POINTS_KNN, &
                                                     N_NEIGHBORS_KNN, n, observed, p, null, n_failed, random_seed=seed, &
                                                     ierr=ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"ierr")
        call run_js_comp_test_label_permutation_test(3_int32, 300_int32, 6_int32, expression, N_POINTS_KNN, &
                                                     N_NEIGHBORS_KNN, n, observed2, p2, null2, n_failed2, random_seed=seed, &
                                                     ierr=ierr)
        call assert_true(same_values(p, p2) .and. same_values(observed, observed2) .and. &
                         same_values(reshape(null, [3*n]), reshape(null2, [3*n])), name//"reproducible")
        call assert_equal_int(n_failed, 0_int32, name//"no failed permutation")
        do i = 1, 3
            call assert_in_range_real(p(i), 1.0_real64/real(n + 1, real64), 1.0_real64, name//"p in [1/(n+1), 1]")
        end do

        call knn_on_assignment(expression, identity(15_int32), row, ierr)
        call assert_true(same_values(observed, row), name//"observed is level (a) on the identity")
        call redraw_assignments(seed, 15_int32, n, assignments)
        do i = 1, n
            call knn_on_assignment(expression, assignments(:, i), row, ierr)
            call assert_equal_int(get_err_code(ierr), ERR_OK, name//"level (a) ierr")
            call assert_true(same_values(null(i, :), row), name//"row equals level (a) on the drawn assignment")
        end do
        call assert_true(any(assignments(:, 1) /= identity(15_int32)), name//"the draw is not the identity")
    end subroutine test_knn_label_test_reproducible_rows_match_level_a

    !> The adaptive counterpart of the reproducibility and row checks.
    subroutine test_adaptive_label_test_reproducible_rows_match_level_a()
        character(*), parameter :: name = "test_adaptive_label_test_reproducible_rows_match_level_a: "
        integer(int32), parameter :: n = 5_int32, seed = 8_int32
        real(real64), allocatable :: expression(:, :, :)
        real(real64) :: observed(3), p(3), null(n, 3), observed2(3), p2(3), null2(n, 3), row(3)
        integer(int32) :: ierr, n_failed, n_failed2, assignments(15, n), i

        call make_expression(300_int32, [6_int32, 5_int32, 4_int32], [0.5_real64, 0.5_real64, 0.5_real64], 42_int32, &
                             expression)
        call run_js_comp_test_adaptive_label_permutation_test(3_int32, 300_int32, 6_int32, expression, K_START, K_STEP, &
                                                              K_MAX, n, observed, p, null, n_failed, random_seed=seed, &
                                                              ierr=ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"ierr")
        call run_js_comp_test_adaptive_label_permutation_test(3_int32, 300_int32, 6_int32, expression, K_START, K_STEP, &
                                                              K_MAX, n, observed2, p2, null2, n_failed2, random_seed=seed, &
                                                              ierr=ierr)
        call assert_true(same_values(p, p2) .and. same_values(observed, observed2) .and. &
                         same_values(reshape(null, [3*n]), reshape(null2, [3*n])), name//"reproducible")
        call assert_equal_int(n_failed, 0_int32, name//"no failed permutation")
        do i = 1, 3
            call assert_in_range_real(p(i), 1.0_real64/real(n + 1, real64), 1.0_real64, name//"p in [1/(n+1), 1]")
        end do

        call adaptive_on_assignment(expression, identity(15_int32), row, ierr)
        call assert_true(same_values(observed, row), name//"observed is level (a) on the identity")
        call redraw_assignments(seed, 15_int32, n, assignments)
        do i = 1, n
            call adaptive_on_assignment(expression, assignments(:, i), row, ierr)
            call assert_equal_int(get_err_code(ierr), ERR_OK, name//"level (a) ierr")
            call assert_true(same_values(null(i, :), row), name//"row equals level (a) on the drawn assignment")
        end do
    end subroutine test_adaptive_label_test_reproducible_rows_match_level_a

    !> k-NN: on null data (every study drawn from the same distribution) the p-values are mostly
    !| well above the minimum `1/(n+1)`; a study with four times the noise gets exactly `1/(n+1)`.
    subroutine test_knn_label_test_null_and_batch()
        character(*), parameter :: name = "test_knn_label_test_null_and_batch: "
        integer(int32), parameter :: n = 19_int32
        real(real64), allocatable :: expression(:, :, :)
        real(real64) :: observed(3), p(3), null(n, 3), p_min
        integer(int32) :: ierr, n_failed, i_round, n_above

        p_min = 1.0_real64/real(n + 1, real64)
        n_above = 0_int32
        do i_round = 1, 4
            call make_expression(300_int32, [5_int32, 5_int32, 5_int32], [0.5_real64, 0.5_real64, 0.5_real64], &
                                 100_int32 + i_round, expression)
            call run_js_comp_test_label_permutation_test(3_int32, 300_int32, 5_int32, expression, N_POINTS_KNN, &
                                                         N_NEIGHBORS_KNN, n, observed, p, null, n_failed, &
                                                         random_seed=i_round, ierr=ierr)
            call assert_equal_int(get_err_code(ierr), ERR_OK, name//"null ierr")
            n_above = n_above + count(p > p_min + 1.0e-12_real64, kind=int32)
        end do
        call assert_true(n_above >= 9_int32, name//"null p-values mostly above the minimum (got "//n_above//" of 12)")

        call make_expression(300_int32, [5_int32, 5_int32, 5_int32], [0.5_real64, 0.5_real64, 2.0_real64], 109_int32, &
                             expression)
        call run_js_comp_test_label_permutation_test(3_int32, 300_int32, 5_int32, expression, N_POINTS_KNN, &
                                                     N_NEIGHBORS_KNN, n, observed, p, null, n_failed, random_seed=5_int32, &
                                                     ierr=ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"batch ierr")
        call assert_true(p(3) == p_min, name//"batch study at the minimum p-value")
    end subroutine test_knn_label_test_null_and_batch

    !> The adaptive counterpart of the null/batch check.
    subroutine test_adaptive_label_test_null_and_batch()
        character(*), parameter :: name = "test_adaptive_label_test_null_and_batch: "
        integer(int32), parameter :: n = 19_int32
        real(real64), allocatable :: expression(:, :, :)
        real(real64) :: observed(3), p(3), null(n, 3), p_min
        integer(int32) :: ierr, n_failed, i_round, n_above

        p_min = 1.0_real64/real(n + 1, real64)
        n_above = 0_int32
        do i_round = 1, 4
            call make_expression(300_int32, [5_int32, 5_int32, 5_int32], [0.5_real64, 0.5_real64, 0.5_real64], &
                                 200_int32 + i_round, expression)
            call run_js_comp_test_adaptive_label_permutation_test(3_int32, 300_int32, 5_int32, expression, K_START, &
                                                                  K_STEP, K_MAX, n, observed, p, null, n_failed, &
                                                                  random_seed=i_round, ierr=ierr)
            call assert_equal_int(get_err_code(ierr), ERR_OK, name//"null ierr")
            n_above = n_above + count(p > p_min + 1.0e-12_real64, kind=int32)
        end do
        call assert_true(n_above >= 9_int32, name//"null p-values mostly above the minimum (got "//n_above//" of 12)")

        call make_expression(300_int32, [5_int32, 5_int32, 5_int32], [0.5_real64, 0.5_real64, 2.0_real64], 209_int32, &
                             expression)
        call run_js_comp_test_adaptive_label_permutation_test(3_int32, 300_int32, 5_int32, expression, K_START, K_STEP, &
                                                              K_MAX, n, observed, p, null, n_failed, random_seed=5_int32, &
                                                              ierr=ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"batch ierr")
        call assert_true(p(3) == p_min, name//"batch study at the minimum p-value")
    end subroutine test_adaptive_label_test_null_and_batch

    !> k-NN rerun option: on, the result equals the search on the identity data followed by the run
    !| with the parameters it selected; off, the fixed-parameter run (test 1). The whole test with
    !| the rerun option reports the same observed value.
    subroutine test_knn_rerun_parameter_search()
        character(*), parameter :: name = "test_knn_rerun_parameter_search: "
        real(real64), allocatable :: expression(:, :, :), gene_means(:, :), residuals(:, :, :)
        real(real64) :: reference(3), result(3), fixed(3), observed(3), p(3), null(1, 3)
        integer(int32) :: ierr, used_np, used_nn, search_np, search_nn, n_failed

        call make_expression(300_int32, [5_int32, 4_int32, 4_int32], [0.5_real64, 0.5_real64, 0.6_real64], 51_int32, &
                             expression)
        call reference_means_residuals(expression, gene_means, residuals)
        call reference_fixed_k_search(gene_means, residuals, 20_int32, search_np, search_nn)
        call reference_knn(expression, search_np, search_nn, reference, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"reference ierr")

        call run_js_comp_test_on_sample_assignment(3_int32, 300_int32, 5_int32, expression, 13_int32, identity(13_int32), &
                                                   N_POINTS_KNN, N_NEIGHBORS_KNN, result, used_np, used_nn, &
                                                   rerun_parameter_search=.true._c_bool, n_bootstraps=20_int32, ierr=ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"rerun ierr")
        call assert_equal_int(used_np, search_np, name//"used_n_points is the search's")
        call assert_equal_int(used_nn, search_nn, name//"used_n_neighbors is the search's")
        call assert_true(same_values(result, reference), name//"rerun equals search + run")

        call run_js_comp_test_on_sample_assignment(3_int32, 300_int32, 5_int32, expression, 13_int32, identity(13_int32), &
                                                   N_POINTS_KNN, N_NEIGHBORS_KNN, fixed, used_np, used_nn, &
                                                   rerun_parameter_search=.false._c_bool, n_bootstraps=20_int32, ierr=ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"off ierr")
        call reference_knn(expression, N_POINTS_KNN, N_NEIGHBORS_KNN, reference, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"off reference ierr")
        call assert_true(equal_and_not_nan(fixed, reference), name//"off equals the fixed-parameter run")
        call assert_equal_int(used_np, N_POINTS_KNN, name//"off: used_n_points")

        call run_js_comp_test_label_permutation_test(3_int32, 300_int32, 5_int32, expression, N_POINTS_KNN, &
                                                     N_NEIGHBORS_KNN, 1_int32, observed, p, null, n_failed, &
                                                     rerun_parameter_search=.true._c_bool, n_bootstraps=20_int32, ierr=ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"label test with rerun ierr")
        call assert_true(same_values(observed, result), name//"label test with rerun: observed value")
    end subroutine test_knn_rerun_parameter_search

    !> Adaptive rerun option: on, the result equals the adaptive search on the identity data followed
    !| by the construction and run with the knobs it selected; off, the fixed-knob run.
    subroutine test_adaptive_rerun_parameter_search()
        character(*), parameter :: name = "test_adaptive_rerun_parameter_search: "
        real(real64), allocatable :: expression(:, :, :), gene_means(:, :), residuals(:, :, :)
        real(real64) :: reference(3), result(3), observed(3), p(3), null(1, 3)
        integer(int32) :: ierr, used_ks, used_kst, used_km, n_points, ks, kst, km, n_points_reference, n_failed

        call make_expression(300_int32, [5_int32, 4_int32, 4_int32], [0.5_real64, 0.5_real64, 0.6_real64], 52_int32, &
                             expression)
        call reference_means_residuals(expression, gene_means, residuals)
        call reference_adaptive_search(gene_means, residuals, 20_int32, ks, kst, km)
        call reference_adaptive(expression, ks, kst, km, reference, n_points_reference, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"reference ierr")

        call run_js_comp_test_adaptive_on_sample_assignment(3_int32, 300_int32, 5_int32, expression, 13_int32, &
                                                            identity(13_int32), K_START, K_STEP, K_MAX, result, used_ks, &
                                                            used_kst, used_km, n_points, &
                                                            rerun_parameter_search=.true._c_bool, n_bootstraps=20_int32, &
                                                            ierr=ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"rerun ierr")
        call assert_equal_int(used_ks, ks, name//"used_k_start is the search's")
        call assert_equal_int(used_kst, kst, name//"used_k_step is the search's")
        call assert_equal_int(used_km, km, name//"used_k_max is the search's")
        call assert_equal_int(n_points, n_points_reference, name//"n_points")
        call assert_true(same_values(result, reference), name//"rerun equals search + construction + run")

        call run_js_comp_test_adaptive_on_sample_assignment(3_int32, 300_int32, 5_int32, expression, 13_int32, &
                                                            identity(13_int32), K_START, K_STEP, K_MAX, result, used_ks, &
                                                            used_kst, used_km, n_points, &
                                                            rerun_parameter_search=.false._c_bool, ierr=ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"off ierr")
        call reference_adaptive(expression, K_START, K_STEP, K_MAX, reference, n_points_reference, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"off reference ierr")
        call assert_true(equal_and_not_nan(result, reference), name//"off equals the fixed-knob run")
        call assert_equal_int(used_ks, K_START, name//"off: used_k_start")

        call run_js_comp_test_adaptive_on_sample_assignment(3_int32, 300_int32, 5_int32, expression, 13_int32, &
                                                            identity(13_int32), K_START, K_STEP, K_MAX, result, used_ks, &
                                                            used_kst, used_km, n_points, &
                                                            rerun_parameter_search=.true._c_bool, n_bootstraps=20_int32, &
                                                            ierr=ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"second rerun ierr")
        call run_js_comp_test_adaptive_label_permutation_test(3_int32, 300_int32, 5_int32, expression, K_START, K_STEP, &
                                                              K_MAX, 1_int32, observed, p, null, n_failed, &
                                                              rerun_parameter_search=.true._c_bool, &
                                                              n_bootstraps=20_int32, ierr=ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"label test with rerun ierr")
        call assert_true(same_values(observed, result), name//"label test with rerun: observed value")
    end subroutine test_adaptive_rerun_parameter_search

    !> The sizing routines, by hand for G = 300 (S = 3): the fixed-k grid's bounds are 300 points
    !| (`clamp(ceiling(4*sqrt(300)) = 70, 300, 1500)`) and `floor(300/(0.25*300)) = 4` neighbors;
    !| the adaptive search's capacity is `ceiling(2*900/11) + 1 = 165` (smallest `k_start` 11).
    subroutine test_permutation_bounds_by_hand()
        character(*), parameter :: name = "test_permutation_bounds_by_hand: "
        integer(int32) :: max_p, max_k, top_k, ierr

        call calc_js_comp_test_permutation_bounds(300_int32, 10_int32, 30_int32, max_n_points_candidate=max_p, &
                                                  max_n_neighbors_candidate=max_k, n_bootstrapping_top_k_jsds=top_k, &
                                                  ierr=ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"knn ierr")
        call assert_equal_int(max_p, 300_int32, name//"max(10, 300)")
        call assert_equal_int(max_k, 30_int32, name//"max(30, 4)")
        call assert_equal_int(top_k, 25_int32, name//"default heap: floor(2.5% of 1000)")
        call calc_js_comp_test_permutation_bounds(300_int32, 400_int32, 2_int32, 100_int32, 5.0_real64, max_p, max_k, &
                                                  top_k, ierr)
        call assert_equal_int(max_p, 400_int32, name//"max(400, 300)")
        call assert_equal_int(max_k, 4_int32, name//"max(2, 4)")
        call assert_equal_int(top_k, 5_int32, name//"floor(5% of 100)")
        call calc_js_comp_test_permutation_bounds(0_int32, 10_int32, 1_int32, max_n_points_candidate=max_p, &
                                                  max_n_neighbors_candidate=max_k, n_bootstrapping_top_k_jsds=top_k, &
                                                  ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, name//"G = 0", arg_pos=1_int32)
        call calc_js_comp_test_permutation_bounds(300_int32, 10_int32, 301_int32, max_n_points_candidate=max_p, &
                                                  max_n_neighbors_candidate=max_k, n_bootstrapping_top_k_jsds=top_k, &
                                                  ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, name//"n_neighbors > G", arg_pos=3_int32)
        call calc_js_comp_test_permutation_bounds(300_int32, 10_int32, 30_int32, 0_int32, max_n_points_candidate=max_p, &
                                                  max_n_neighbors_candidate=max_k, n_bootstrapping_top_k_jsds=top_k, &
                                                  ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, name//"n_bootstraps 0", arg_pos=4_int32)
        call calc_js_comp_test_permutation_bounds(300_int32, 10_int32, 30_int32, 100_int32, 101.0_real64, max_p, max_k, &
                                                  top_k, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, name//"level 101", arg_pos=5_int32)

        call calc_adaptive_js_comp_test_permutation_bounds(300_int32, 3_int32, 30_int32, max_n_points_candidate=max_p, &
                                                           n_bootstrapping_top_k_jsds=top_k, ierr=ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"adaptive ierr")
        call assert_equal_int(max_p, 165_int32, name//"max(165, ceiling(2*900/30) + 1 = 61)")
        call assert_equal_int(top_k, 25_int32, name//"adaptive default heap")
        call calc_adaptive_js_comp_test_permutation_bounds(300_int32, 3_int32, 5_int32, 200_int32, 2.5_real64, max_p, &
                                                           top_k, ierr)
        call assert_equal_int(max_p, 361_int32, name//"max(165, ceiling(2*900/5) + 1 = 361)")
        call assert_equal_int(top_k, 5_int32, name//"floor(2.5% of 200)")
        call calc_adaptive_js_comp_test_permutation_bounds(300_int32, 3_int32, 1_int32, max_n_points_candidate=max_p, &
                                                           n_bootstrapping_top_k_jsds=top_k, ierr=ierr)
        call assert_equal_int(max_p, 900_int32, name//"capped at N")
        call calc_adaptive_js_comp_test_permutation_bounds(300_int32, 3_int32, 0_int32, max_n_points_candidate=max_p, &
                                                           n_bootstrapping_top_k_jsds=top_k, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, name//"k_start 0", arg_pos=3_int32)
        call calc_adaptive_js_comp_test_permutation_bounds(300_int32, 0_int32, 30_int32, max_n_points_candidate=max_p, &
                                                           n_bootstrapping_top_k_jsds=top_k, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, name//"n_studies 0", arg_pos=2_int32)
    end subroutine test_permutation_bounds_by_hand

    !> A direct implementation call with `n_permutations = 0` (both methods): the observed value is
    !| computed, the p-values are NaN, nothing failed.
    subroutine test_label_test_zero_permutations_direct_call()
        character(*), parameter :: name = "test_label_test_zero_permutations_direct_call: "
        real(real64), allocatable :: expression(:, :, :), null(:, :)
        real(real64) :: observed(3), p(3), reference(3)
        integer(int32) :: ierr, n_failed, max_p, max_k, top_k, n_points

        call make_expression(300_int32, [4_int32, 4_int32, 4_int32], [0.5_real64, 0.5_real64, 0.5_real64], 61_int32, &
                             expression)
        allocate (null(0, 3))
        call calc_js_comp_test_permutation_bounds(300_int32, N_POINTS_KNN, N_NEIGHBORS_KNN, max_n_points_candidate=max_p, &
                                                  max_n_neighbors_candidate=max_k, n_bootstrapping_top_k_jsds=top_k, &
                                                  ierr=ierr)
        call direct_knn_label_test(3_int32, 300_int32, 4_int32, expression, N_POINTS_KNN, N_NEIGHBORS_KNN, max_p, max_k, &
                                   top_k, 0_int32, observed, p, null, n_failed, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"knn ierr")
        call assert_true(all(ieee_is_nan(p)), name//"knn p-values NaN")
        call assert_equal_int(n_failed, 0_int32, name//"knn n_failed")
        call reference_knn(expression, N_POINTS_KNN, N_NEIGHBORS_KNN, reference, ierr)
        call assert_true(same_values(observed, reference), name//"knn observed")

        call calc_adaptive_js_comp_test_permutation_bounds(300_int32, 3_int32, K_START, max_n_points_candidate=max_p, &
                                                           n_bootstrapping_top_k_jsds=top_k, ierr=ierr)
        call direct_adaptive_label_test(3_int32, 300_int32, 4_int32, expression, K_START, K_STEP, K_MAX, max_p, top_k, &
                                        0_int32, observed, p, null, n_failed, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"adaptive ierr")
        call assert_true(all(ieee_is_nan(p)), name//"adaptive p-values NaN")
        call reference_adaptive(expression, K_START, K_STEP, K_MAX, reference, n_points, ierr)
        call assert_true(same_values(observed, reference), name//"adaptive observed")
    end subroutine test_label_test_zero_permutations_direct_call

    !> Failed permutations, a natural case: every study has 2 of its 5 samples shifted by +20, so
    !| on the observed data every study's gene means are shifted alike. A shuffle that gives the
    !| studies unequal numbers of shifted samples separates their means, some adaptive neighborhood
    !| then holds no gene of some study, and the run rejects it with invalid input. Such a
    !| permutation's `jsd_null` row is NaN, it is counted in `n_failed_permutations`, and the
    !| p-values count only the valid rows; when every permutation fails, the p-values are NaN.
    subroutine test_adaptive_label_test_failed_permutations_bookkeeping()
        character(*), parameter :: name = "test_adaptive_label_test_failed_permutations_bookkeeping: "
        integer(int32), parameter :: n = 19_int32
        real(real64), allocatable :: expression(:, :, :)
        real(real64) :: observed(3), p(3), null(n, 3), expected, row(3)
        integer(int32) :: ierr, n_failed, i_study, n_valid, n_at_least, i, assignments(15, n), i_seed
        logical :: is_failed(n), found
        real(real64) :: null_one(1, 3)

        call make_expression(300_int32, [5_int32, 5_int32, 5_int32], [0.5_real64, 0.5_real64, 0.5_real64], 71_int32, &
                             expression)
        expression(1:2, :, :) = expression(1:2, :, :) + 20.0_real64
        call run_js_comp_test_adaptive_label_permutation_test(3_int32, 300_int32, 5_int32, expression, K_START, K_STEP, &
                                                              K_MAX, n, observed, p, null, n_failed, random_seed=3_int32, &
                                                              ierr=ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"ierr")
        call assert_no_nan_real(observed, 3_int32, name//"observed data run")
        do i = 1, n
            is_failed(i) = all(ieee_is_nan(null(i, :)))
            call assert_true(is_failed(i) .or. .not. any(ieee_is_nan(null(i, :))), name//"a row is all NaN or none")
        end do
        call assert_equal_int(n_failed, count(is_failed, kind=int32), name//"n_failed counts the NaN rows")
        call assert_true(n_failed > 0_int32 .and. n_failed < n, name//"some but not all permutations failed (got "// &
                         n_failed//")")
        ! A failed row is exactly a level-(a) call that fails with invalid input.
        call redraw_assignments(3_int32, 15_int32, n, assignments)
        do i = 1, n
            call adaptive_on_assignment(expression, assignments(:, i), row, ierr)
            if (is_failed(i)) then
                call assert_err(ierr, ERR_INVALID_INPUT, name//"failed row: level (a) rejects the assignment")
            else
                call assert_true(same_values(null(i, :), row), name//"valid row equals level (a)")
            end if
        end do
        n_valid = count(.not. is_failed, kind=int32)
        do i_study = 1, 3
            n_at_least = count((.not. is_failed) .and. null(:, i_study) >= observed(i_study), kind=int32)
            expected = real(n_at_least + 1, real64)/real(n_valid + 1, real64)
            call assert_true(p(i_study) == expected, name//"p over the valid permutations only")
        end do
        ! Every permutation failed (n_valid = 0): one permutation, and the first seed whose single
        ! draw fails on this data. The p-values are then NaN, the observed values are not.
        found = .false.
        do i_seed = 1, 40
            call run_js_comp_test_adaptive_label_permutation_test(3_int32, 300_int32, 5_int32, expression, K_START, &
                                                                  K_STEP, K_MAX, 1_int32, observed, p, null_one, &
                                                                  n_failed, random_seed=i_seed, ierr=ierr)
            call assert_equal_int(get_err_code(ierr), ERR_OK, name//"single permutation ierr")
            if (n_failed == 1_int32) then
                found = .true.
                exit
            end if
        end do
        call assert_true(found, name//"some seed makes the single permutation fail")
        call assert_true(all(ieee_is_nan(p)), name//"all permutations failed: p NaN")
        call assert_true(all(ieee_is_nan(null_one)), name//"all permutations failed: the null row NaN")
        call assert_no_nan_real(observed, 3_int32, name//"all permutations failed: observed kept")
    end subroutine test_adaptive_label_test_failed_permutations_bookkeeping

    !> Level (b) rejects what does not depend on the assignment with its own positions (a study
    !| without samples, a rerun capacity below the search's), and returns an error on the observed
    !| data at once with NaN outputs.
    subroutine test_label_test_rejects_invalid_input()
        character(*), parameter :: name = "test_label_test_rejects_invalid_input: "
        real(real64), allocatable :: expression(:, :, :), empty_study(:, :, :)
        real(real64) :: observed(3), p(3), null(2, 3)
        integer(int32) :: ierr, n_failed, max_p, max_k, top_k, search_cap

        call make_expression(300_int32, [4_int32, 4_int32, 4_int32], [0.5_real64, 0.5_real64, 0.5_real64], 81_int32, &
                             expression)
        empty_study = expression
        empty_study(:, :, 1) = ieee_value(1.0_real64, ieee_quiet_nan)
        call calc_js_comp_test_permutation_bounds(300_int32, N_POINTS_KNN, N_NEIGHBORS_KNN, max_n_points_candidate=max_p, &
                                                  max_n_neighbors_candidate=max_k, n_bootstrapping_top_k_jsds=top_k, &
                                                  ierr=ierr)
        call direct_knn_label_test(3_int32, 300_int32, 4_int32, empty_study, N_POINTS_KNN, N_NEIGHBORS_KNN, max_p, max_k, &
                                   top_k, 2_int32, observed, p, null, n_failed, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, name//"knn: empty study", arg_pos=4_int32)
        call assert_true(all(ieee_is_nan(p)) .and. all(ieee_is_nan(null)), name//"knn: NaN outputs")
        call direct_knn_label_test(3_int32, 300_int32, 4_int32, expression, N_POINTS_KNN, N_NEIGHBORS_KNN, max_p - 1_int32, &
                                   max_k, top_k, 2_int32, observed, p, null, n_failed, ierr, &
                                   rerun_parameter_search=.true._c_bool)
        call assert_err(ierr, ERR_INVALID_INPUT, name//"knn: point capacity", arg_pos=7_int32)
        call direct_knn_label_test(3_int32, 300_int32, 4_int32, expression, 4_int32, 4_int32, max_p, 3_int32, top_k, &
                                   2_int32, observed, p, null, n_failed, ierr, rerun_parameter_search=.true._c_bool)
        call assert_err(ierr, ERR_INVALID_INPUT, name//"knn: neighbor capacity", arg_pos=8_int32)

        call calc_adaptive_js_comp_test_permutation_bounds(300_int32, 3_int32, K_START, max_n_points_candidate=max_p, &
                                                           n_bootstrapping_top_k_jsds=top_k, ierr=ierr)
        call direct_adaptive_label_test(3_int32, 300_int32, 4_int32, empty_study, K_START, K_STEP, K_MAX, max_p, top_k, &
                                        2_int32, observed, p, null, n_failed, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, name//"adaptive: empty study", arg_pos=4_int32)
        call calc_adaptive_js_comp_test_bounds(300_int32, 3_int32, search_cap, ierr)
        call direct_adaptive_label_test(3_int32, 300_int32, 4_int32, expression, K_START, K_STEP, K_MAX, &
                                        search_cap - 1_int32, top_k, 2_int32, observed, p, null, n_failed, ierr, &
                                        rerun_parameter_search=.true._c_bool)
        call assert_err(ierr, ERR_INVALID_INPUT, name//"adaptive: capacity", arg_pos=8_int32)

        ! The observed data fail (k_start beyond the data): returned at once, every output NaN.
        call run_js_comp_test_adaptive_label_permutation_test(3_int32, 300_int32, 4_int32, expression, 1000_int32, &
                                                              10_int32, 1000_int32, 2_int32, observed, p, null, n_failed, &
                                                              ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, name//"observed data fail", arg_pos=0_int32)
        call assert_true(all(ieee_is_nan(observed)) .and. all(ieee_is_nan(p)) .and. all(ieee_is_nan(null)), &
                         name//"observed failure: NaN outputs")
        call assert_equal_int(n_failed, 0_int32, name//"observed failure: nothing counted")
    end subroutine test_label_test_rejects_invalid_input

    ! ========================================================================================
    ! Direct implementation calls (generated from the implementations' dummy lists)
    ! ========================================================================================

    !> Test helper: calls run_js_comp_test_on_sample_assignment_impl directly, every work array allocated here.
    subroutine direct_knn_on_assignment(n_studies, max_n_genes_all_studies, max_n_reps_all_studies, expression, &
                                        n_samples, sample_assignment, n_points, n_neighbors, &
                                        max_n_points_candidate, max_n_neighbors_candidate, &
                                        n_bootstrapping_top_k_jsds, global_js_divergence, used_n_points, &
                                        used_n_neighbors, ierr, rerun_parameter_search, n_bootstraps, random_seed)
        integer(int32), intent(in) :: n_studies
            !! Forwarded to run_js_comp_test_on_sample_assignment_impl
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Forwarded to run_js_comp_test_on_sample_assignment_impl
        integer(int32), intent(in) :: max_n_reps_all_studies
            !! Forwarded to run_js_comp_test_on_sample_assignment_impl
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(in) :: expression
            !! Forwarded to run_js_comp_test_on_sample_assignment_impl
        integer(int32), intent(in) :: n_samples
            !! Forwarded to run_js_comp_test_on_sample_assignment_impl
        integer(int32), dimension(n_samples), intent(in) :: sample_assignment
            !! Forwarded to run_js_comp_test_on_sample_assignment_impl
        integer(int32), intent(in) :: n_points
            !! Forwarded to run_js_comp_test_on_sample_assignment_impl
        integer(int32), intent(in) :: n_neighbors
            !! Forwarded to run_js_comp_test_on_sample_assignment_impl
        integer(int32), intent(in) :: max_n_points_candidate
            !! Forwarded to run_js_comp_test_on_sample_assignment_impl
        integer(int32), intent(in) :: max_n_neighbors_candidate
            !! Forwarded to run_js_comp_test_on_sample_assignment_impl
        integer(int32), intent(in) :: n_bootstrapping_top_k_jsds
            !! Forwarded to run_js_comp_test_on_sample_assignment_impl
        real(real64), dimension(n_studies), intent(out) :: global_js_divergence
            !! Forwarded to run_js_comp_test_on_sample_assignment_impl
        integer(int32), intent(out) :: used_n_points
            !! Forwarded to run_js_comp_test_on_sample_assignment_impl
        integer(int32), intent(out) :: used_n_neighbors
            !! Forwarded to run_js_comp_test_on_sample_assignment_impl
        integer(int32), intent(out) :: ierr
            !! Error code of the implementation
        logical(c_bool), intent(in), optional :: rerun_parameter_search
            !! Forwarded to run_js_comp_test_on_sample_assignment_impl
        integer(int32), intent(in), optional :: n_bootstraps
            !! Forwarded to run_js_comp_test_on_sample_assignment_impl
        integer(int32), intent(in), optional :: random_seed
            !! Forwarded to run_js_comp_test_on_sample_assignment_impl
        integer(int32), allocatable :: tmp_sample_slots(:)
        logical(c_bool), allocatable :: tmp_sample_mark(:)
        real(real64), allocatable :: tmp_permuted_expression(:, :, :)
        real(real64), allocatable :: tmp_gene_means(:, :)
        real(real64), allocatable :: tmp_residuals(:, :, :)
        integer(int32), allocatable :: tmp_gene_means_perms(:, :)
        integer(int32), allocatable :: tmp_gene_means_perm_all(:)
        real(real64), allocatable :: tmp_x_star(:)
        integer(int32), allocatable :: tmp_neighborhood_indices(:, :, :)
        integer(int32), allocatable :: tmp_neighborhood_range(:, :, :)
        real(real64), allocatable :: tmp_neighborhood_residuals_gathered(:, :, :)
        integer(int32), allocatable :: tmp_point_pool(:)
        integer(int32), allocatable :: tmp_point_draw(:)
        integer(int32), allocatable :: tmp_n_bins_per_point(:)
        real(real64), allocatable :: tmp_shared_residual_range_low(:)
        real(real64), allocatable :: tmp_shared_residual_range_high(:)
        logical(c_bool), allocatable :: tmp_occupancy_failed(:)
        integer(int32), allocatable :: tmp_n_pooled_residuals(:)
        integer(int32), allocatable :: tmp_min_bin_occupancy(:)
        real(real64), allocatable :: tmp_mean_bin_occupancy(:)
        integer(int32), allocatable :: tmp_max_bin_occupancy(:)
        integer(int32), allocatable :: tmp_sturges_bins(:)
        integer(int32), allocatable :: tmp_fd_bins(:)
        real(real64), allocatable :: tmp_pmfs(:, :, :)
        integer(int32), allocatable :: tmp_counts(:, :, :)
        integer(int32), allocatable :: tmp_included_n_reps(:, :)
        real(real64), allocatable :: tmp_mean_pmf(:, :)
        integer(int32), allocatable :: tmp_mean_pmf_counts(:, :)
        integer(int32), allocatable :: tmp_mean_pmf_included_n_reps(:)
        real(real64), allocatable :: tmp_js_divergences(:, :)
        real(real64), allocatable :: tmp_weights(:, :)
        real(real64), allocatable :: tmp_global_js_divergence(:)
        real(real64), allocatable :: tmp_pooled_residuals(:)
        integer(int32), allocatable :: tmp_pooled_residuals_perm(:)
        integer(int32), allocatable :: tmp_residual_to_study_perm(:)
        integer(int32), allocatable :: tmp_bin_counts_search(:)
        integer(int32), allocatable :: tmp_study_block_ends(:)
        integer(int32), allocatable :: tmp_study_block_sizes(:)
        integer(int32), allocatable :: tmp_point_n_neighbors(:)
        real(real64), allocatable :: tmp_p_values(:)
        real(real64), allocatable :: tmp_p_values_observed_consensus(:)
        real(real64), allocatable :: tmp_permutation_mean_pmf(:, :)
        integer(int32), allocatable :: tmp_permutation_mean_pmf_counts(:, :)
        integer(int32), allocatable :: tmp_permutation_counts(:, :)
        real(real64), allocatable :: tmp_permutation_pmfs(:, :, :)
        real(real64), allocatable :: tmp_permutation_js_divergences(:, :)
        real(real64), allocatable :: tmp_permutation_weights(:, :)
        integer(int32), allocatable :: tmp_search_n_bins_per_point(:)
        real(real64), allocatable :: tmp_search_shared_residual_range_low(:)
        real(real64), allocatable :: tmp_search_shared_residual_range_high(:)
        real(real64), allocatable :: tmp_search_confidence_interval(:, :)
        integer(int32), allocatable :: tmp_search_trace_n_points(:)
        integer(int32), allocatable :: tmp_search_trace_n_neighbors(:)
        real(real64), allocatable :: tmp_search_trace_global_js_divergence(:, :)
        real(real64), allocatable :: tmp_search_trace_ci_lower(:, :)
        real(real64), allocatable :: tmp_search_trace_ci_upper(:, :)
        real(real64), allocatable :: tmp_search_trace_ci_width(:, :)
        real(real64), allocatable :: tmp_search_trace_ci_width_relative(:, :)
        real(real64), allocatable :: tmp_search_trace_ci_overlap(:, :)
        real(real64), allocatable :: tmp_search_trace_delta(:, :)
        real(real64), allocatable :: tmp_search_trace_delta_median(:)
        real(real64), allocatable :: tmp_search_trace_delta_max(:)
        integer(int32), allocatable :: tmp_search_trace_selected_n_bins(:, :)
        logical(c_bool), allocatable :: tmp_search_trace_occupancy_failed(:, :)
        integer(int32), allocatable :: tmp_search_trace_n_pooled_residuals(:, :)
        integer(int32), allocatable :: tmp_search_trace_min_bin_occupancy(:, :)
        real(real64), allocatable :: tmp_search_trace_mean_bin_occupancy(:, :)
        integer(int32), allocatable :: tmp_search_trace_max_bin_occupancy(:, :)
        integer(int32), allocatable :: tmp_search_trace_sturges_bins(:, :)
        integer(int32), allocatable :: tmp_search_trace_fd_bins(:, :)
        real(real64), allocatable :: tmp_search_trace_shared_residual_range_low(:, :)
        real(real64), allocatable :: tmp_search_trace_shared_residual_range_high(:, :)
        real(real64), allocatable :: tmp_candidate_confidence_interval(:, :)
        real(real64), allocatable :: tmp_bootstrapping_top_k_jsds(:, :, :)
        real(real64), allocatable :: tmp_prev_global_js_divergence(:)
        real(real64), allocatable :: tmp_best_uncertainty_confidence_interval(:, :)
        integer(int32), allocatable :: tmp_best_n_bins_per_point(:)
        integer(int32), allocatable :: tmp_best_uncertainty_n_bins_per_point(:)
        real(real64), allocatable :: tmp_best_shared_residual_range_low(:)
        real(real64), allocatable :: tmp_best_shared_residual_range_high(:)
        real(real64), allocatable :: tmp_best_uncertainty_shared_residual_range_low(:)
        real(real64), allocatable :: tmp_best_uncertainty_shared_residual_range_high(:)
        allocate (tmp_sample_slots(max_n_reps_all_studies*n_studies))
        allocate (tmp_sample_mark(n_samples))
        allocate (tmp_permuted_expression(max_n_reps_all_studies, max_n_genes_all_studies, n_studies))
        allocate (tmp_gene_means(max_n_genes_all_studies, n_studies))
        allocate (tmp_residuals(max_n_reps_all_studies, max_n_genes_all_studies, n_studies))
        allocate (tmp_gene_means_perms(max_n_genes_all_studies, n_studies))
        allocate (tmp_gene_means_perm_all(max_n_genes_all_studies*n_studies))
        allocate (tmp_x_star(max_n_points_candidate))
        allocate (tmp_neighborhood_indices(max_n_neighbors_candidate, max_n_points_candidate, n_studies))
        allocate (tmp_neighborhood_range(2, max_n_points_candidate, n_studies))
        allocate (tmp_neighborhood_residuals_gathered(max_n_reps_all_studies, max_n_neighbors_candidate, max_n_points_candidate))
        allocate (tmp_point_pool(MAX_N_BINS))
        allocate (tmp_point_draw(MAX_N_BINS))
        allocate (tmp_n_bins_per_point(max_n_points_candidate))
        allocate (tmp_shared_residual_range_low(max_n_points_candidate))
        allocate (tmp_shared_residual_range_high(max_n_points_candidate))
        allocate (tmp_occupancy_failed(max_n_points_candidate))
        allocate (tmp_n_pooled_residuals(max_n_points_candidate))
        allocate (tmp_min_bin_occupancy(max_n_points_candidate))
        allocate (tmp_mean_bin_occupancy(max_n_points_candidate))
        allocate (tmp_max_bin_occupancy(max_n_points_candidate))
        allocate (tmp_sturges_bins(max_n_points_candidate))
        allocate (tmp_fd_bins(max_n_points_candidate))
        allocate (tmp_pmfs(max_n_points_candidate, MAX_N_BINS, n_studies))
        allocate (tmp_counts(max_n_points_candidate, MAX_N_BINS, n_studies))
        allocate (tmp_included_n_reps(max_n_points_candidate, n_studies))
        allocate (tmp_mean_pmf(max_n_points_candidate, MAX_N_BINS))
        allocate (tmp_mean_pmf_counts(max_n_points_candidate, MAX_N_BINS))
        allocate (tmp_mean_pmf_included_n_reps(max_n_points_candidate))
        allocate (tmp_js_divergences(max_n_points_candidate, n_studies))
        allocate (tmp_weights(max_n_points_candidate, n_studies))
        allocate (tmp_global_js_divergence(n_studies))
        allocate (tmp_pooled_residuals(max_n_reps_all_studies*max_n_neighbors_candidate*n_studies))
        allocate (tmp_pooled_residuals_perm(max_n_reps_all_studies*max_n_neighbors_candidate*n_studies))
        allocate (tmp_residual_to_study_perm(max_n_reps_all_studies*max_n_neighbors_candidate*n_studies))
        allocate (tmp_bin_counts_search(MAX_N_BINS))
        allocate (tmp_study_block_ends(n_studies))
        allocate (tmp_study_block_sizes(n_studies))
        allocate (tmp_point_n_neighbors(n_studies))
        allocate (tmp_p_values(n_studies))
        allocate (tmp_p_values_observed_consensus(n_studies))
        allocate (tmp_permutation_mean_pmf(max_n_points_candidate, MAX_N_BINS))
        allocate (tmp_permutation_mean_pmf_counts(max_n_points_candidate, MAX_N_BINS))
        allocate (tmp_permutation_counts(max_n_points_candidate, MAX_N_BINS))
        allocate (tmp_permutation_pmfs(max_n_points_candidate, MAX_N_BINS, n_studies))
        allocate (tmp_permutation_js_divergences(max_n_points_candidate, n_studies))
        allocate (tmp_permutation_weights(max_n_points_candidate, n_studies))
        allocate (tmp_search_n_bins_per_point(max_n_points_candidate))
        allocate (tmp_search_shared_residual_range_low(max_n_points_candidate))
        allocate (tmp_search_shared_residual_range_high(max_n_points_candidate))
        allocate (tmp_search_confidence_interval(2, n_studies))
        allocate (tmp_search_trace_n_points(MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_n_neighbors(MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_global_js_divergence(n_studies, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_ci_lower(n_studies, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_ci_upper(n_studies, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_ci_width(n_studies, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_ci_width_relative(n_studies, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_ci_overlap(n_studies, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_delta(n_studies, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_delta_median(MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_delta_max(MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_selected_n_bins(max_n_points_candidate, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_occupancy_failed(max_n_points_candidate, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_n_pooled_residuals(max_n_points_candidate, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_min_bin_occupancy(max_n_points_candidate, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_mean_bin_occupancy(max_n_points_candidate, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_max_bin_occupancy(max_n_points_candidate, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_sturges_bins(max_n_points_candidate, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_fd_bins(max_n_points_candidate, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_shared_residual_range_low(max_n_points_candidate, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_shared_residual_range_high(max_n_points_candidate, MAX_CANDIDATE_PAIRS))
        allocate (tmp_candidate_confidence_interval(2, n_studies))
        allocate (tmp_bootstrapping_top_k_jsds(n_bootstrapping_top_k_jsds, 2, n_studies))
        allocate (tmp_prev_global_js_divergence(n_studies))
        allocate (tmp_best_uncertainty_confidence_interval(2, n_studies))
        allocate (tmp_best_n_bins_per_point(max_n_points_candidate))
        allocate (tmp_best_uncertainty_n_bins_per_point(max_n_points_candidate))
        allocate (tmp_best_shared_residual_range_low(max_n_points_candidate))
        allocate (tmp_best_shared_residual_range_high(max_n_points_candidate))
        allocate (tmp_best_uncertainty_shared_residual_range_low(max_n_points_candidate))
        allocate (tmp_best_uncertainty_shared_residual_range_high(max_n_points_candidate))
        call run_js_comp_test_on_sample_assignment_impl(n_studies=n_studies, &
            max_n_genes_all_studies=max_n_genes_all_studies, max_n_reps_all_studies=max_n_reps_all_studies, &
            expression=expression, n_samples=n_samples, sample_assignment=sample_assignment, n_points=n_points, &
            n_neighbors=n_neighbors, max_n_points_candidate=max_n_points_candidate, &
            max_n_neighbors_candidate=max_n_neighbors_candidate, &
            n_bootstrapping_top_k_jsds=n_bootstrapping_top_k_jsds, global_js_divergence=global_js_divergence, &
            used_n_points=used_n_points, used_n_neighbors=used_n_neighbors, tmp_sample_slots=tmp_sample_slots, &
            tmp_sample_mark=tmp_sample_mark, tmp_permuted_expression=tmp_permuted_expression, &
            tmp_gene_means=tmp_gene_means, tmp_residuals=tmp_residuals, tmp_gene_means_perms=tmp_gene_means_perms, &
            tmp_gene_means_perm_all=tmp_gene_means_perm_all, tmp_x_star=tmp_x_star, &
            tmp_neighborhood_indices=tmp_neighborhood_indices, tmp_neighborhood_range=tmp_neighborhood_range, &
            tmp_neighborhood_residuals_gathered=tmp_neighborhood_residuals_gathered, &
            tmp_point_pool=tmp_point_pool, tmp_point_draw=tmp_point_draw, &
            tmp_n_bins_per_point=tmp_n_bins_per_point, &
            tmp_shared_residual_range_low=tmp_shared_residual_range_low, &
            tmp_shared_residual_range_high=tmp_shared_residual_range_high, &
            tmp_occupancy_failed=tmp_occupancy_failed, tmp_n_pooled_residuals=tmp_n_pooled_residuals, &
            tmp_min_bin_occupancy=tmp_min_bin_occupancy, tmp_mean_bin_occupancy=tmp_mean_bin_occupancy, &
            tmp_max_bin_occupancy=tmp_max_bin_occupancy, tmp_sturges_bins=tmp_sturges_bins, &
            tmp_fd_bins=tmp_fd_bins, tmp_pmfs=tmp_pmfs, tmp_counts=tmp_counts, &
            tmp_included_n_reps=tmp_included_n_reps, tmp_mean_pmf=tmp_mean_pmf, &
            tmp_mean_pmf_counts=tmp_mean_pmf_counts, tmp_mean_pmf_included_n_reps=tmp_mean_pmf_included_n_reps, &
            tmp_js_divergences=tmp_js_divergences, tmp_weights=tmp_weights, &
            tmp_global_js_divergence=tmp_global_js_divergence, tmp_pooled_residuals=tmp_pooled_residuals, &
            tmp_pooled_residuals_perm=tmp_pooled_residuals_perm, &
            tmp_residual_to_study_perm=tmp_residual_to_study_perm, tmp_bin_counts_search=tmp_bin_counts_search, &
            tmp_study_block_ends=tmp_study_block_ends, tmp_study_block_sizes=tmp_study_block_sizes, &
            tmp_point_n_neighbors=tmp_point_n_neighbors, tmp_p_values=tmp_p_values, &
            tmp_p_values_observed_consensus=tmp_p_values_observed_consensus, &
            tmp_permutation_mean_pmf=tmp_permutation_mean_pmf, &
            tmp_permutation_mean_pmf_counts=tmp_permutation_mean_pmf_counts, &
            tmp_permutation_counts=tmp_permutation_counts, tmp_permutation_pmfs=tmp_permutation_pmfs, &
            tmp_permutation_js_divergences=tmp_permutation_js_divergences, &
            tmp_permutation_weights=tmp_permutation_weights, &
            tmp_search_n_bins_per_point=tmp_search_n_bins_per_point, &
            tmp_search_shared_residual_range_low=tmp_search_shared_residual_range_low, &
            tmp_search_shared_residual_range_high=tmp_search_shared_residual_range_high, &
            tmp_search_confidence_interval=tmp_search_confidence_interval, &
            tmp_search_trace_n_points=tmp_search_trace_n_points, &
            tmp_search_trace_n_neighbors=tmp_search_trace_n_neighbors, &
            tmp_search_trace_global_js_divergence=tmp_search_trace_global_js_divergence, &
            tmp_search_trace_ci_lower=tmp_search_trace_ci_lower, &
            tmp_search_trace_ci_upper=tmp_search_trace_ci_upper, &
            tmp_search_trace_ci_width=tmp_search_trace_ci_width, &
            tmp_search_trace_ci_width_relative=tmp_search_trace_ci_width_relative, &
            tmp_search_trace_ci_overlap=tmp_search_trace_ci_overlap, &
            tmp_search_trace_delta=tmp_search_trace_delta, &
            tmp_search_trace_delta_median=tmp_search_trace_delta_median, &
            tmp_search_trace_delta_max=tmp_search_trace_delta_max, &
            tmp_search_trace_selected_n_bins=tmp_search_trace_selected_n_bins, &
            tmp_search_trace_occupancy_failed=tmp_search_trace_occupancy_failed, &
            tmp_search_trace_n_pooled_residuals=tmp_search_trace_n_pooled_residuals, &
            tmp_search_trace_min_bin_occupancy=tmp_search_trace_min_bin_occupancy, &
            tmp_search_trace_mean_bin_occupancy=tmp_search_trace_mean_bin_occupancy, &
            tmp_search_trace_max_bin_occupancy=tmp_search_trace_max_bin_occupancy, &
            tmp_search_trace_sturges_bins=tmp_search_trace_sturges_bins, &
            tmp_search_trace_fd_bins=tmp_search_trace_fd_bins, &
            tmp_search_trace_shared_residual_range_low=tmp_search_trace_shared_residual_range_low, &
            tmp_search_trace_shared_residual_range_high=tmp_search_trace_shared_residual_range_high, &
            tmp_candidate_confidence_interval=tmp_candidate_confidence_interval, &
            tmp_bootstrapping_top_k_jsds=tmp_bootstrapping_top_k_jsds, &
            tmp_prev_global_js_divergence=tmp_prev_global_js_divergence, &
            tmp_best_uncertainty_confidence_interval=tmp_best_uncertainty_confidence_interval, &
            tmp_best_n_bins_per_point=tmp_best_n_bins_per_point, &
            tmp_best_uncertainty_n_bins_per_point=tmp_best_uncertainty_n_bins_per_point, &
            tmp_best_shared_residual_range_low=tmp_best_shared_residual_range_low, &
            tmp_best_shared_residual_range_high=tmp_best_shared_residual_range_high, &
            tmp_best_uncertainty_shared_residual_range_low=tmp_best_uncertainty_shared_residual_range_low, &
            tmp_best_uncertainty_shared_residual_range_high=tmp_best_uncertainty_shared_residual_range_high, &
            rerun_parameter_search=rerun_parameter_search, n_bootstraps=n_bootstraps, random_seed=random_seed, &
            ierr=ierr)
    end subroutine direct_knn_on_assignment

    !> Test helper: calls run_js_comp_test_adaptive_on_sample_assignment_impl directly, every work array allocated here.
    subroutine direct_adaptive_on_assignment(n_studies, max_n_genes_all_studies, max_n_reps_all_studies, &
                                             expression, n_samples, sample_assignment, k_start, k_step, k_max, &
                                             max_n_points_candidate, n_bootstrapping_top_k_jsds, &
                                             global_js_divergence, used_k_start, used_k_step, used_k_max, &
                                             n_points, ierr, rerun_parameter_search, n_bootstraps, random_seed)
        integer(int32), intent(in) :: n_studies
            !! Forwarded to run_js_comp_test_adaptive_on_sample_assignment_impl
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Forwarded to run_js_comp_test_adaptive_on_sample_assignment_impl
        integer(int32), intent(in) :: max_n_reps_all_studies
            !! Forwarded to run_js_comp_test_adaptive_on_sample_assignment_impl
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(in) :: expression
            !! Forwarded to run_js_comp_test_adaptive_on_sample_assignment_impl
        integer(int32), intent(in) :: n_samples
            !! Forwarded to run_js_comp_test_adaptive_on_sample_assignment_impl
        integer(int32), dimension(n_samples), intent(in) :: sample_assignment
            !! Forwarded to run_js_comp_test_adaptive_on_sample_assignment_impl
        integer(int32), intent(in) :: k_start
            !! Forwarded to run_js_comp_test_adaptive_on_sample_assignment_impl
        integer(int32), intent(in) :: k_step
            !! Forwarded to run_js_comp_test_adaptive_on_sample_assignment_impl
        integer(int32), intent(in) :: k_max
            !! Forwarded to run_js_comp_test_adaptive_on_sample_assignment_impl
        integer(int32), intent(in) :: max_n_points_candidate
            !! Forwarded to run_js_comp_test_adaptive_on_sample_assignment_impl
        integer(int32), intent(in) :: n_bootstrapping_top_k_jsds
            !! Forwarded to run_js_comp_test_adaptive_on_sample_assignment_impl
        real(real64), dimension(n_studies), intent(out) :: global_js_divergence
            !! Forwarded to run_js_comp_test_adaptive_on_sample_assignment_impl
        integer(int32), intent(out) :: used_k_start
            !! Forwarded to run_js_comp_test_adaptive_on_sample_assignment_impl
        integer(int32), intent(out) :: used_k_step
            !! Forwarded to run_js_comp_test_adaptive_on_sample_assignment_impl
        integer(int32), intent(out) :: used_k_max
            !! Forwarded to run_js_comp_test_adaptive_on_sample_assignment_impl
        integer(int32), intent(out) :: n_points
            !! Forwarded to run_js_comp_test_adaptive_on_sample_assignment_impl
        integer(int32), intent(out) :: ierr
            !! Error code of the implementation
        logical(c_bool), intent(in), optional :: rerun_parameter_search
            !! Forwarded to run_js_comp_test_adaptive_on_sample_assignment_impl
        integer(int32), intent(in), optional :: n_bootstraps
            !! Forwarded to run_js_comp_test_adaptive_on_sample_assignment_impl
        integer(int32), intent(in), optional :: random_seed
            !! Forwarded to run_js_comp_test_adaptive_on_sample_assignment_impl
        integer(int32), allocatable :: tmp_sample_slots(:)
        logical(c_bool), allocatable :: tmp_sample_mark(:)
        real(real64), allocatable :: tmp_permuted_expression(:, :, :)
        real(real64), allocatable :: tmp_gene_means(:, :)
        real(real64), allocatable :: tmp_residuals(:, :, :)
        integer(int32), allocatable :: tmp_gene_means_perm_all(:)
        real(real64), allocatable :: tmp_x_star(:)
        integer(int32), allocatable :: tmp_pooled_neighborhood_range(:, :)
        integer(int32), allocatable :: tmp_n_neighbors_per_point(:, :)
        integer(int32), allocatable :: tmp_stop_reason(:)
        real(real64), allocatable :: tmp_neighborhood_dispersion(:)
        real(real64), allocatable :: tmp_neighborhood_mad(:)
        integer(int32), allocatable :: tmp_point_neighborhood_indices(:, :)
        real(real64), allocatable :: tmp_neighbor_residuals(:, :)
        integer(int32), allocatable :: tmp_point_pool(:)
        integer(int32), allocatable :: tmp_point_draw(:)
        integer(int32), allocatable :: tmp_n_bins_per_point(:)
        real(real64), allocatable :: tmp_shared_residual_range_low(:)
        real(real64), allocatable :: tmp_shared_residual_range_high(:)
        logical(c_bool), allocatable :: tmp_occupancy_failed(:)
        integer(int32), allocatable :: tmp_n_pooled_residuals(:)
        integer(int32), allocatable :: tmp_min_bin_occupancy(:)
        real(real64), allocatable :: tmp_mean_bin_occupancy(:)
        integer(int32), allocatable :: tmp_max_bin_occupancy(:)
        integer(int32), allocatable :: tmp_sturges_bins(:)
        integer(int32), allocatable :: tmp_fd_bins(:)
        real(real64), allocatable :: tmp_pmfs(:, :, :)
        integer(int32), allocatable :: tmp_counts(:, :, :)
        integer(int32), allocatable :: tmp_included_n_reps(:, :)
        real(real64), allocatable :: tmp_mean_pmf(:, :)
        integer(int32), allocatable :: tmp_mean_pmf_counts(:, :)
        integer(int32), allocatable :: tmp_mean_pmf_included_n_reps(:)
        real(real64), allocatable :: tmp_js_divergences(:, :)
        real(real64), allocatable :: tmp_weights(:, :)
        real(real64), allocatable :: tmp_global_js_divergence(:)
        real(real64), allocatable :: tmp_pooled_residuals(:)
        integer(int32), allocatable :: tmp_pooled_residuals_perm(:)
        integer(int32), allocatable :: tmp_residual_to_study_perm(:)
        integer(int32), allocatable :: tmp_bin_counts_search(:)
        integer(int32), allocatable :: tmp_study_block_ends(:)
        integer(int32), allocatable :: tmp_study_block_sizes(:)
        real(real64), allocatable :: tmp_p_values(:)
        real(real64), allocatable :: tmp_p_values_observed_consensus(:)
        real(real64), allocatable :: tmp_permutation_mean_pmf(:, :)
        integer(int32), allocatable :: tmp_permutation_mean_pmf_counts(:, :)
        integer(int32), allocatable :: tmp_permutation_counts(:, :)
        real(real64), allocatable :: tmp_permutation_pmfs(:, :, :)
        real(real64), allocatable :: tmp_permutation_js_divergences(:, :)
        real(real64), allocatable :: tmp_permutation_weights(:, :)
        integer(int32), allocatable :: tmp_search_n_bins_per_point(:)
        real(real64), allocatable :: tmp_search_shared_residual_range_low(:)
        real(real64), allocatable :: tmp_search_shared_residual_range_high(:)
        real(real64), allocatable :: tmp_search_confidence_interval(:, :)
        integer(int32), allocatable :: tmp_search_trace_k_start(:)
        integer(int32), allocatable :: tmp_search_trace_k_step(:)
        integer(int32), allocatable :: tmp_search_trace_k_max(:)
        integer(int32), allocatable :: tmp_search_trace_n_points(:)
        real(real64), allocatable :: tmp_search_trace_global_js_divergence(:, :)
        real(real64), allocatable :: tmp_search_trace_ci_lower(:, :)
        real(real64), allocatable :: tmp_search_trace_ci_upper(:, :)
        real(real64), allocatable :: tmp_search_trace_ci_width(:, :)
        real(real64), allocatable :: tmp_search_trace_ci_width_relative(:, :)
        real(real64), allocatable :: tmp_search_trace_ci_overlap(:, :)
        real(real64), allocatable :: tmp_search_trace_delta(:, :)
        real(real64), allocatable :: tmp_search_trace_delta_median(:)
        real(real64), allocatable :: tmp_search_trace_delta_max(:)
        integer(int32), allocatable :: tmp_search_trace_selected_n_bins(:, :)
        logical(c_bool), allocatable :: tmp_search_trace_occupancy_failed(:, :)
        integer(int32), allocatable :: tmp_search_trace_n_pooled_residuals(:, :)
        integer(int32), allocatable :: tmp_search_trace_min_bin_occupancy(:, :)
        real(real64), allocatable :: tmp_search_trace_mean_bin_occupancy(:, :)
        integer(int32), allocatable :: tmp_search_trace_max_bin_occupancy(:, :)
        integer(int32), allocatable :: tmp_search_trace_sturges_bins(:, :)
        integer(int32), allocatable :: tmp_search_trace_fd_bins(:, :)
        real(real64), allocatable :: tmp_search_trace_shared_residual_range_low(:, :)
        real(real64), allocatable :: tmp_search_trace_shared_residual_range_high(:, :)
        integer(int32), allocatable :: tmp_search_candidate_k_start(:)
        integer(int32), allocatable :: tmp_search_candidate_n_points(:)
        integer(int32), allocatable :: tmp_search_candidate_status(:)
        real(real64), allocatable :: tmp_candidate_confidence_interval(:, :)
        real(real64), allocatable :: tmp_bootstrapping_top_k_jsds(:, :, :)
        real(real64), allocatable :: tmp_prev_global_js_divergence(:)
        integer(int32), allocatable :: tmp_delta_perm(:)
        real(real64), allocatable :: tmp_best_uncertainty_confidence_interval(:, :)
        integer(int32), allocatable :: tmp_best_n_bins_per_point(:)
        integer(int32), allocatable :: tmp_best_uncertainty_n_bins_per_point(:)
        real(real64), allocatable :: tmp_best_shared_residual_range_low(:)
        real(real64), allocatable :: tmp_best_shared_residual_range_high(:)
        real(real64), allocatable :: tmp_best_uncertainty_shared_residual_range_low(:)
        real(real64), allocatable :: tmp_best_uncertainty_shared_residual_range_high(:)
        allocate (tmp_sample_slots(max_n_reps_all_studies*n_studies))
        allocate (tmp_sample_mark(n_samples))
        allocate (tmp_permuted_expression(max_n_reps_all_studies, max_n_genes_all_studies, n_studies))
        allocate (tmp_gene_means(max_n_genes_all_studies, n_studies))
        allocate (tmp_residuals(max_n_reps_all_studies, max_n_genes_all_studies, n_studies))
        allocate (tmp_gene_means_perm_all(max_n_genes_all_studies*n_studies))
        allocate (tmp_x_star(max_n_genes_all_studies*n_studies))
        allocate (tmp_pooled_neighborhood_range(2, max_n_genes_all_studies*n_studies))
        allocate (tmp_n_neighbors_per_point(n_studies, max_n_genes_all_studies*n_studies))
        allocate (tmp_stop_reason(max_n_genes_all_studies*n_studies))
        allocate (tmp_neighborhood_dispersion(max_n_genes_all_studies*n_studies))
        allocate (tmp_neighborhood_mad(max_n_genes_all_studies*n_studies))
        allocate (tmp_point_neighborhood_indices(max_n_genes_all_studies, n_studies))
        allocate (tmp_neighbor_residuals(max_n_reps_all_studies, max_n_genes_all_studies))
        allocate (tmp_point_pool(MAX_N_BINS))
        allocate (tmp_point_draw(MAX_N_BINS))
        allocate (tmp_n_bins_per_point(max_n_points_candidate))
        allocate (tmp_shared_residual_range_low(max_n_points_candidate))
        allocate (tmp_shared_residual_range_high(max_n_points_candidate))
        allocate (tmp_occupancy_failed(max_n_points_candidate))
        allocate (tmp_n_pooled_residuals(max_n_points_candidate))
        allocate (tmp_min_bin_occupancy(max_n_points_candidate))
        allocate (tmp_mean_bin_occupancy(max_n_points_candidate))
        allocate (tmp_max_bin_occupancy(max_n_points_candidate))
        allocate (tmp_sturges_bins(max_n_points_candidate))
        allocate (tmp_fd_bins(max_n_points_candidate))
        allocate (tmp_pmfs(max_n_points_candidate, MAX_N_BINS, n_studies))
        allocate (tmp_counts(max_n_points_candidate, MAX_N_BINS, n_studies))
        allocate (tmp_included_n_reps(max_n_points_candidate, n_studies))
        allocate (tmp_mean_pmf(max_n_points_candidate, MAX_N_BINS))
        allocate (tmp_mean_pmf_counts(max_n_points_candidate, MAX_N_BINS))
        allocate (tmp_mean_pmf_included_n_reps(max_n_points_candidate))
        allocate (tmp_js_divergences(max_n_points_candidate, n_studies))
        allocate (tmp_weights(max_n_points_candidate, n_studies))
        allocate (tmp_global_js_divergence(n_studies))
        allocate (tmp_pooled_residuals(max_n_reps_all_studies*max_n_genes_all_studies*n_studies))
        allocate (tmp_pooled_residuals_perm(max_n_reps_all_studies*max_n_genes_all_studies*n_studies))
        allocate (tmp_residual_to_study_perm(max_n_reps_all_studies*max_n_genes_all_studies*n_studies))
        allocate (tmp_bin_counts_search(MAX_N_BINS))
        allocate (tmp_study_block_ends(n_studies))
        allocate (tmp_study_block_sizes(n_studies))
        allocate (tmp_p_values(n_studies))
        allocate (tmp_p_values_observed_consensus(n_studies))
        allocate (tmp_permutation_mean_pmf(max_n_points_candidate, MAX_N_BINS))
        allocate (tmp_permutation_mean_pmf_counts(max_n_points_candidate, MAX_N_BINS))
        allocate (tmp_permutation_counts(max_n_points_candidate, MAX_N_BINS))
        allocate (tmp_permutation_pmfs(max_n_points_candidate, MAX_N_BINS, n_studies))
        allocate (tmp_permutation_js_divergences(max_n_points_candidate, n_studies))
        allocate (tmp_permutation_weights(max_n_points_candidate, n_studies))
        allocate (tmp_search_n_bins_per_point(max_n_points_candidate))
        allocate (tmp_search_shared_residual_range_low(max_n_points_candidate))
        allocate (tmp_search_shared_residual_range_high(max_n_points_candidate))
        allocate (tmp_search_confidence_interval(2, n_studies))
        allocate (tmp_search_trace_k_start(MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_k_step(MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_k_max(MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_n_points(MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_global_js_divergence(n_studies, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_ci_lower(n_studies, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_ci_upper(n_studies, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_ci_width(n_studies, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_ci_width_relative(n_studies, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_ci_overlap(n_studies, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_delta(n_studies, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_delta_median(MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_delta_max(MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_selected_n_bins(max_n_points_candidate, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_occupancy_failed(max_n_points_candidate, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_n_pooled_residuals(max_n_points_candidate, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_min_bin_occupancy(max_n_points_candidate, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_mean_bin_occupancy(max_n_points_candidate, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_max_bin_occupancy(max_n_points_candidate, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_sturges_bins(max_n_points_candidate, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_fd_bins(max_n_points_candidate, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_shared_residual_range_low(max_n_points_candidate, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_shared_residual_range_high(max_n_points_candidate, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_candidate_k_start(MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_candidate_n_points(MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_candidate_status(MAX_CANDIDATE_PAIRS))
        allocate (tmp_candidate_confidence_interval(2, n_studies))
        allocate (tmp_bootstrapping_top_k_jsds(n_bootstrapping_top_k_jsds, 2, n_studies))
        allocate (tmp_prev_global_js_divergence(n_studies))
        allocate (tmp_delta_perm(n_studies))
        allocate (tmp_best_uncertainty_confidence_interval(2, n_studies))
        allocate (tmp_best_n_bins_per_point(max_n_points_candidate))
        allocate (tmp_best_uncertainty_n_bins_per_point(max_n_points_candidate))
        allocate (tmp_best_shared_residual_range_low(max_n_points_candidate))
        allocate (tmp_best_shared_residual_range_high(max_n_points_candidate))
        allocate (tmp_best_uncertainty_shared_residual_range_low(max_n_points_candidate))
        allocate (tmp_best_uncertainty_shared_residual_range_high(max_n_points_candidate))
        call run_js_comp_test_adaptive_on_sample_assignment_impl(n_studies=n_studies, &
            max_n_genes_all_studies=max_n_genes_all_studies, max_n_reps_all_studies=max_n_reps_all_studies, &
            expression=expression, n_samples=n_samples, sample_assignment=sample_assignment, k_start=k_start, &
            k_step=k_step, k_max=k_max, max_n_points_candidate=max_n_points_candidate, &
            n_bootstrapping_top_k_jsds=n_bootstrapping_top_k_jsds, global_js_divergence=global_js_divergence, &
            used_k_start=used_k_start, used_k_step=used_k_step, used_k_max=used_k_max, n_points=n_points, &
            tmp_sample_slots=tmp_sample_slots, tmp_sample_mark=tmp_sample_mark, &
            tmp_permuted_expression=tmp_permuted_expression, tmp_gene_means=tmp_gene_means, &
            tmp_residuals=tmp_residuals, tmp_gene_means_perm_all=tmp_gene_means_perm_all, tmp_x_star=tmp_x_star, &
            tmp_pooled_neighborhood_range=tmp_pooled_neighborhood_range, &
            tmp_n_neighbors_per_point=tmp_n_neighbors_per_point, tmp_stop_reason=tmp_stop_reason, &
            tmp_neighborhood_dispersion=tmp_neighborhood_dispersion, tmp_neighborhood_mad=tmp_neighborhood_mad, &
            tmp_point_neighborhood_indices=tmp_point_neighborhood_indices, &
            tmp_neighbor_residuals=tmp_neighbor_residuals, tmp_point_pool=tmp_point_pool, &
            tmp_point_draw=tmp_point_draw, tmp_n_bins_per_point=tmp_n_bins_per_point, &
            tmp_shared_residual_range_low=tmp_shared_residual_range_low, &
            tmp_shared_residual_range_high=tmp_shared_residual_range_high, &
            tmp_occupancy_failed=tmp_occupancy_failed, tmp_n_pooled_residuals=tmp_n_pooled_residuals, &
            tmp_min_bin_occupancy=tmp_min_bin_occupancy, tmp_mean_bin_occupancy=tmp_mean_bin_occupancy, &
            tmp_max_bin_occupancy=tmp_max_bin_occupancy, tmp_sturges_bins=tmp_sturges_bins, &
            tmp_fd_bins=tmp_fd_bins, tmp_pmfs=tmp_pmfs, tmp_counts=tmp_counts, &
            tmp_included_n_reps=tmp_included_n_reps, tmp_mean_pmf=tmp_mean_pmf, &
            tmp_mean_pmf_counts=tmp_mean_pmf_counts, tmp_mean_pmf_included_n_reps=tmp_mean_pmf_included_n_reps, &
            tmp_js_divergences=tmp_js_divergences, tmp_weights=tmp_weights, &
            tmp_global_js_divergence=tmp_global_js_divergence, tmp_pooled_residuals=tmp_pooled_residuals, &
            tmp_pooled_residuals_perm=tmp_pooled_residuals_perm, &
            tmp_residual_to_study_perm=tmp_residual_to_study_perm, tmp_bin_counts_search=tmp_bin_counts_search, &
            tmp_study_block_ends=tmp_study_block_ends, tmp_study_block_sizes=tmp_study_block_sizes, &
            tmp_p_values=tmp_p_values, tmp_p_values_observed_consensus=tmp_p_values_observed_consensus, &
            tmp_permutation_mean_pmf=tmp_permutation_mean_pmf, &
            tmp_permutation_mean_pmf_counts=tmp_permutation_mean_pmf_counts, &
            tmp_permutation_counts=tmp_permutation_counts, tmp_permutation_pmfs=tmp_permutation_pmfs, &
            tmp_permutation_js_divergences=tmp_permutation_js_divergences, &
            tmp_permutation_weights=tmp_permutation_weights, &
            tmp_search_n_bins_per_point=tmp_search_n_bins_per_point, &
            tmp_search_shared_residual_range_low=tmp_search_shared_residual_range_low, &
            tmp_search_shared_residual_range_high=tmp_search_shared_residual_range_high, &
            tmp_search_confidence_interval=tmp_search_confidence_interval, &
            tmp_search_trace_k_start=tmp_search_trace_k_start, tmp_search_trace_k_step=tmp_search_trace_k_step, &
            tmp_search_trace_k_max=tmp_search_trace_k_max, tmp_search_trace_n_points=tmp_search_trace_n_points, &
            tmp_search_trace_global_js_divergence=tmp_search_trace_global_js_divergence, &
            tmp_search_trace_ci_lower=tmp_search_trace_ci_lower, &
            tmp_search_trace_ci_upper=tmp_search_trace_ci_upper, &
            tmp_search_trace_ci_width=tmp_search_trace_ci_width, &
            tmp_search_trace_ci_width_relative=tmp_search_trace_ci_width_relative, &
            tmp_search_trace_ci_overlap=tmp_search_trace_ci_overlap, &
            tmp_search_trace_delta=tmp_search_trace_delta, &
            tmp_search_trace_delta_median=tmp_search_trace_delta_median, &
            tmp_search_trace_delta_max=tmp_search_trace_delta_max, &
            tmp_search_trace_selected_n_bins=tmp_search_trace_selected_n_bins, &
            tmp_search_trace_occupancy_failed=tmp_search_trace_occupancy_failed, &
            tmp_search_trace_n_pooled_residuals=tmp_search_trace_n_pooled_residuals, &
            tmp_search_trace_min_bin_occupancy=tmp_search_trace_min_bin_occupancy, &
            tmp_search_trace_mean_bin_occupancy=tmp_search_trace_mean_bin_occupancy, &
            tmp_search_trace_max_bin_occupancy=tmp_search_trace_max_bin_occupancy, &
            tmp_search_trace_sturges_bins=tmp_search_trace_sturges_bins, &
            tmp_search_trace_fd_bins=tmp_search_trace_fd_bins, &
            tmp_search_trace_shared_residual_range_low=tmp_search_trace_shared_residual_range_low, &
            tmp_search_trace_shared_residual_range_high=tmp_search_trace_shared_residual_range_high, &
            tmp_search_candidate_k_start=tmp_search_candidate_k_start, &
            tmp_search_candidate_n_points=tmp_search_candidate_n_points, &
            tmp_search_candidate_status=tmp_search_candidate_status, &
            tmp_candidate_confidence_interval=tmp_candidate_confidence_interval, &
            tmp_bootstrapping_top_k_jsds=tmp_bootstrapping_top_k_jsds, &
            tmp_prev_global_js_divergence=tmp_prev_global_js_divergence, tmp_delta_perm=tmp_delta_perm, &
            tmp_best_uncertainty_confidence_interval=tmp_best_uncertainty_confidence_interval, &
            tmp_best_n_bins_per_point=tmp_best_n_bins_per_point, &
            tmp_best_uncertainty_n_bins_per_point=tmp_best_uncertainty_n_bins_per_point, &
            tmp_best_shared_residual_range_low=tmp_best_shared_residual_range_low, &
            tmp_best_shared_residual_range_high=tmp_best_shared_residual_range_high, &
            tmp_best_uncertainty_shared_residual_range_low=tmp_best_uncertainty_shared_residual_range_low, &
            tmp_best_uncertainty_shared_residual_range_high=tmp_best_uncertainty_shared_residual_range_high, &
            rerun_parameter_search=rerun_parameter_search, n_bootstraps=n_bootstraps, random_seed=random_seed, &
            ierr=ierr)
    end subroutine direct_adaptive_on_assignment

    !> Test helper: calls run_js_comp_test_label_permutation_test_impl directly, every work array allocated here.
    subroutine direct_knn_label_test(n_studies, max_n_genes_all_studies, max_n_reps_all_studies, expression, &
                                     n_points, n_neighbors, max_n_points_candidate, max_n_neighbors_candidate, &
                                     n_bootstrapping_top_k_jsds, n_permutations, global_js_divergence_observed, &
                                     p_values, jsd_null, n_failed_permutations, ierr, rerun_parameter_search, &
                                     n_bootstraps, random_seed)
        integer(int32), intent(in) :: n_studies
            !! Forwarded to run_js_comp_test_label_permutation_test_impl
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Forwarded to run_js_comp_test_label_permutation_test_impl
        integer(int32), intent(in) :: max_n_reps_all_studies
            !! Forwarded to run_js_comp_test_label_permutation_test_impl
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(in) :: expression
            !! Forwarded to run_js_comp_test_label_permutation_test_impl
        integer(int32), intent(in) :: n_points
            !! Forwarded to run_js_comp_test_label_permutation_test_impl
        integer(int32), intent(in) :: n_neighbors
            !! Forwarded to run_js_comp_test_label_permutation_test_impl
        integer(int32), intent(in) :: max_n_points_candidate
            !! Forwarded to run_js_comp_test_label_permutation_test_impl
        integer(int32), intent(in) :: max_n_neighbors_candidate
            !! Forwarded to run_js_comp_test_label_permutation_test_impl
        integer(int32), intent(in) :: n_bootstrapping_top_k_jsds
            !! Forwarded to run_js_comp_test_label_permutation_test_impl
        integer(int32), intent(in) :: n_permutations
            !! Forwarded to run_js_comp_test_label_permutation_test_impl
        real(real64), dimension(n_studies), intent(out) :: global_js_divergence_observed
            !! Forwarded to run_js_comp_test_label_permutation_test_impl
        real(real64), dimension(n_studies), intent(out) :: p_values
            !! Forwarded to run_js_comp_test_label_permutation_test_impl
        real(real64), dimension(n_permutations, n_studies), intent(out) :: jsd_null
            !! Forwarded to run_js_comp_test_label_permutation_test_impl
        integer(int32), intent(out) :: n_failed_permutations
            !! Forwarded to run_js_comp_test_label_permutation_test_impl
        integer(int32), intent(out) :: ierr
            !! Error code of the implementation
        logical(c_bool), intent(in), optional :: rerun_parameter_search
            !! Forwarded to run_js_comp_test_label_permutation_test_impl
        integer(int32), intent(in), optional :: n_bootstraps
            !! Forwarded to run_js_comp_test_label_permutation_test_impl
        integer(int32), intent(in), optional :: random_seed
            !! Forwarded to run_js_comp_test_label_permutation_test_impl
        integer(int32), allocatable :: tmp_sample_slots(:)
        integer(int32), allocatable :: tmp_sample_assignment(:)
        real(real64), allocatable :: tmp_permuted_global_js_divergence(:)
        logical(c_bool), allocatable :: tmp_sample_mark(:)
        real(real64), allocatable :: tmp_permuted_expression(:, :, :)
        real(real64), allocatable :: tmp_gene_means(:, :)
        real(real64), allocatable :: tmp_residuals(:, :, :)
        integer(int32), allocatable :: tmp_gene_means_perms(:, :)
        integer(int32), allocatable :: tmp_gene_means_perm_all(:)
        real(real64), allocatable :: tmp_x_star(:)
        integer(int32), allocatable :: tmp_neighborhood_indices(:, :, :)
        integer(int32), allocatable :: tmp_neighborhood_range(:, :, :)
        real(real64), allocatable :: tmp_neighborhood_residuals_gathered(:, :, :)
        integer(int32), allocatable :: tmp_point_pool(:)
        integer(int32), allocatable :: tmp_point_draw(:)
        integer(int32), allocatable :: tmp_n_bins_per_point(:)
        real(real64), allocatable :: tmp_shared_residual_range_low(:)
        real(real64), allocatable :: tmp_shared_residual_range_high(:)
        logical(c_bool), allocatable :: tmp_occupancy_failed(:)
        integer(int32), allocatable :: tmp_n_pooled_residuals(:)
        integer(int32), allocatable :: tmp_min_bin_occupancy(:)
        real(real64), allocatable :: tmp_mean_bin_occupancy(:)
        integer(int32), allocatable :: tmp_max_bin_occupancy(:)
        integer(int32), allocatable :: tmp_sturges_bins(:)
        integer(int32), allocatable :: tmp_fd_bins(:)
        real(real64), allocatable :: tmp_pmfs(:, :, :)
        integer(int32), allocatable :: tmp_counts(:, :, :)
        integer(int32), allocatable :: tmp_included_n_reps(:, :)
        real(real64), allocatable :: tmp_mean_pmf(:, :)
        integer(int32), allocatable :: tmp_mean_pmf_counts(:, :)
        integer(int32), allocatable :: tmp_mean_pmf_included_n_reps(:)
        real(real64), allocatable :: tmp_js_divergences(:, :)
        real(real64), allocatable :: tmp_weights(:, :)
        real(real64), allocatable :: tmp_global_js_divergence(:)
        real(real64), allocatable :: tmp_pooled_residuals(:)
        integer(int32), allocatable :: tmp_pooled_residuals_perm(:)
        integer(int32), allocatable :: tmp_residual_to_study_perm(:)
        integer(int32), allocatable :: tmp_bin_counts_search(:)
        integer(int32), allocatable :: tmp_study_block_ends(:)
        integer(int32), allocatable :: tmp_study_block_sizes(:)
        integer(int32), allocatable :: tmp_point_n_neighbors(:)
        real(real64), allocatable :: tmp_p_values(:)
        real(real64), allocatable :: tmp_p_values_observed_consensus(:)
        real(real64), allocatable :: tmp_permutation_mean_pmf(:, :)
        integer(int32), allocatable :: tmp_permutation_mean_pmf_counts(:, :)
        integer(int32), allocatable :: tmp_permutation_counts(:, :)
        real(real64), allocatable :: tmp_permutation_pmfs(:, :, :)
        real(real64), allocatable :: tmp_permutation_js_divergences(:, :)
        real(real64), allocatable :: tmp_permutation_weights(:, :)
        integer(int32), allocatable :: tmp_search_n_bins_per_point(:)
        real(real64), allocatable :: tmp_search_shared_residual_range_low(:)
        real(real64), allocatable :: tmp_search_shared_residual_range_high(:)
        real(real64), allocatable :: tmp_search_confidence_interval(:, :)
        integer(int32), allocatable :: tmp_search_trace_n_points(:)
        integer(int32), allocatable :: tmp_search_trace_n_neighbors(:)
        real(real64), allocatable :: tmp_search_trace_global_js_divergence(:, :)
        real(real64), allocatable :: tmp_search_trace_ci_lower(:, :)
        real(real64), allocatable :: tmp_search_trace_ci_upper(:, :)
        real(real64), allocatable :: tmp_search_trace_ci_width(:, :)
        real(real64), allocatable :: tmp_search_trace_ci_width_relative(:, :)
        real(real64), allocatable :: tmp_search_trace_ci_overlap(:, :)
        real(real64), allocatable :: tmp_search_trace_delta(:, :)
        real(real64), allocatable :: tmp_search_trace_delta_median(:)
        real(real64), allocatable :: tmp_search_trace_delta_max(:)
        integer(int32), allocatable :: tmp_search_trace_selected_n_bins(:, :)
        logical(c_bool), allocatable :: tmp_search_trace_occupancy_failed(:, :)
        integer(int32), allocatable :: tmp_search_trace_n_pooled_residuals(:, :)
        integer(int32), allocatable :: tmp_search_trace_min_bin_occupancy(:, :)
        real(real64), allocatable :: tmp_search_trace_mean_bin_occupancy(:, :)
        integer(int32), allocatable :: tmp_search_trace_max_bin_occupancy(:, :)
        integer(int32), allocatable :: tmp_search_trace_sturges_bins(:, :)
        integer(int32), allocatable :: tmp_search_trace_fd_bins(:, :)
        real(real64), allocatable :: tmp_search_trace_shared_residual_range_low(:, :)
        real(real64), allocatable :: tmp_search_trace_shared_residual_range_high(:, :)
        real(real64), allocatable :: tmp_candidate_confidence_interval(:, :)
        real(real64), allocatable :: tmp_bootstrapping_top_k_jsds(:, :, :)
        real(real64), allocatable :: tmp_prev_global_js_divergence(:)
        real(real64), allocatable :: tmp_best_uncertainty_confidence_interval(:, :)
        integer(int32), allocatable :: tmp_best_n_bins_per_point(:)
        integer(int32), allocatable :: tmp_best_uncertainty_n_bins_per_point(:)
        real(real64), allocatable :: tmp_best_shared_residual_range_low(:)
        real(real64), allocatable :: tmp_best_shared_residual_range_high(:)
        real(real64), allocatable :: tmp_best_uncertainty_shared_residual_range_low(:)
        real(real64), allocatable :: tmp_best_uncertainty_shared_residual_range_high(:)
        allocate (tmp_sample_slots(max_n_reps_all_studies*n_studies))
        allocate (tmp_sample_assignment(max_n_reps_all_studies*n_studies))
        allocate (tmp_permuted_global_js_divergence(n_studies))
        allocate (tmp_sample_mark(max_n_reps_all_studies*n_studies))
        allocate (tmp_permuted_expression(max_n_reps_all_studies, max_n_genes_all_studies, n_studies))
        allocate (tmp_gene_means(max_n_genes_all_studies, n_studies))
        allocate (tmp_residuals(max_n_reps_all_studies, max_n_genes_all_studies, n_studies))
        allocate (tmp_gene_means_perms(max_n_genes_all_studies, n_studies))
        allocate (tmp_gene_means_perm_all(max_n_genes_all_studies*n_studies))
        allocate (tmp_x_star(max_n_points_candidate))
        allocate (tmp_neighborhood_indices(max_n_neighbors_candidate, max_n_points_candidate, n_studies))
        allocate (tmp_neighborhood_range(2, max_n_points_candidate, n_studies))
        allocate (tmp_neighborhood_residuals_gathered(max_n_reps_all_studies, max_n_neighbors_candidate, max_n_points_candidate))
        allocate (tmp_point_pool(MAX_N_BINS))
        allocate (tmp_point_draw(MAX_N_BINS))
        allocate (tmp_n_bins_per_point(max_n_points_candidate))
        allocate (tmp_shared_residual_range_low(max_n_points_candidate))
        allocate (tmp_shared_residual_range_high(max_n_points_candidate))
        allocate (tmp_occupancy_failed(max_n_points_candidate))
        allocate (tmp_n_pooled_residuals(max_n_points_candidate))
        allocate (tmp_min_bin_occupancy(max_n_points_candidate))
        allocate (tmp_mean_bin_occupancy(max_n_points_candidate))
        allocate (tmp_max_bin_occupancy(max_n_points_candidate))
        allocate (tmp_sturges_bins(max_n_points_candidate))
        allocate (tmp_fd_bins(max_n_points_candidate))
        allocate (tmp_pmfs(max_n_points_candidate, MAX_N_BINS, n_studies))
        allocate (tmp_counts(max_n_points_candidate, MAX_N_BINS, n_studies))
        allocate (tmp_included_n_reps(max_n_points_candidate, n_studies))
        allocate (tmp_mean_pmf(max_n_points_candidate, MAX_N_BINS))
        allocate (tmp_mean_pmf_counts(max_n_points_candidate, MAX_N_BINS))
        allocate (tmp_mean_pmf_included_n_reps(max_n_points_candidate))
        allocate (tmp_js_divergences(max_n_points_candidate, n_studies))
        allocate (tmp_weights(max_n_points_candidate, n_studies))
        allocate (tmp_global_js_divergence(n_studies))
        allocate (tmp_pooled_residuals(max_n_reps_all_studies*max_n_neighbors_candidate*n_studies))
        allocate (tmp_pooled_residuals_perm(max_n_reps_all_studies*max_n_neighbors_candidate*n_studies))
        allocate (tmp_residual_to_study_perm(max_n_reps_all_studies*max_n_neighbors_candidate*n_studies))
        allocate (tmp_bin_counts_search(MAX_N_BINS))
        allocate (tmp_study_block_ends(n_studies))
        allocate (tmp_study_block_sizes(n_studies))
        allocate (tmp_point_n_neighbors(n_studies))
        allocate (tmp_p_values(n_studies))
        allocate (tmp_p_values_observed_consensus(n_studies))
        allocate (tmp_permutation_mean_pmf(max_n_points_candidate, MAX_N_BINS))
        allocate (tmp_permutation_mean_pmf_counts(max_n_points_candidate, MAX_N_BINS))
        allocate (tmp_permutation_counts(max_n_points_candidate, MAX_N_BINS))
        allocate (tmp_permutation_pmfs(max_n_points_candidate, MAX_N_BINS, n_studies))
        allocate (tmp_permutation_js_divergences(max_n_points_candidate, n_studies))
        allocate (tmp_permutation_weights(max_n_points_candidate, n_studies))
        allocate (tmp_search_n_bins_per_point(max_n_points_candidate))
        allocate (tmp_search_shared_residual_range_low(max_n_points_candidate))
        allocate (tmp_search_shared_residual_range_high(max_n_points_candidate))
        allocate (tmp_search_confidence_interval(2, n_studies))
        allocate (tmp_search_trace_n_points(MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_n_neighbors(MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_global_js_divergence(n_studies, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_ci_lower(n_studies, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_ci_upper(n_studies, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_ci_width(n_studies, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_ci_width_relative(n_studies, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_ci_overlap(n_studies, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_delta(n_studies, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_delta_median(MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_delta_max(MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_selected_n_bins(max_n_points_candidate, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_occupancy_failed(max_n_points_candidate, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_n_pooled_residuals(max_n_points_candidate, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_min_bin_occupancy(max_n_points_candidate, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_mean_bin_occupancy(max_n_points_candidate, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_max_bin_occupancy(max_n_points_candidate, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_sturges_bins(max_n_points_candidate, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_fd_bins(max_n_points_candidate, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_shared_residual_range_low(max_n_points_candidate, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_shared_residual_range_high(max_n_points_candidate, MAX_CANDIDATE_PAIRS))
        allocate (tmp_candidate_confidence_interval(2, n_studies))
        allocate (tmp_bootstrapping_top_k_jsds(n_bootstrapping_top_k_jsds, 2, n_studies))
        allocate (tmp_prev_global_js_divergence(n_studies))
        allocate (tmp_best_uncertainty_confidence_interval(2, n_studies))
        allocate (tmp_best_n_bins_per_point(max_n_points_candidate))
        allocate (tmp_best_uncertainty_n_bins_per_point(max_n_points_candidate))
        allocate (tmp_best_shared_residual_range_low(max_n_points_candidate))
        allocate (tmp_best_shared_residual_range_high(max_n_points_candidate))
        allocate (tmp_best_uncertainty_shared_residual_range_low(max_n_points_candidate))
        allocate (tmp_best_uncertainty_shared_residual_range_high(max_n_points_candidate))
        call run_js_comp_test_label_permutation_test_impl(n_studies=n_studies, &
            max_n_genes_all_studies=max_n_genes_all_studies, max_n_reps_all_studies=max_n_reps_all_studies, &
            expression=expression, n_points=n_points, n_neighbors=n_neighbors, &
            max_n_points_candidate=max_n_points_candidate, max_n_neighbors_candidate=max_n_neighbors_candidate, &
            n_bootstrapping_top_k_jsds=n_bootstrapping_top_k_jsds, n_permutations=n_permutations, &
            global_js_divergence_observed=global_js_divergence_observed, p_values=p_values, jsd_null=jsd_null, &
            n_failed_permutations=n_failed_permutations, tmp_sample_slots=tmp_sample_slots, &
            tmp_sample_assignment=tmp_sample_assignment, &
            tmp_permuted_global_js_divergence=tmp_permuted_global_js_divergence, tmp_sample_mark=tmp_sample_mark, &
            tmp_permuted_expression=tmp_permuted_expression, tmp_gene_means=tmp_gene_means, &
            tmp_residuals=tmp_residuals, tmp_gene_means_perms=tmp_gene_means_perms, &
            tmp_gene_means_perm_all=tmp_gene_means_perm_all, tmp_x_star=tmp_x_star, &
            tmp_neighborhood_indices=tmp_neighborhood_indices, tmp_neighborhood_range=tmp_neighborhood_range, &
            tmp_neighborhood_residuals_gathered=tmp_neighborhood_residuals_gathered, &
            tmp_point_pool=tmp_point_pool, tmp_point_draw=tmp_point_draw, &
            tmp_n_bins_per_point=tmp_n_bins_per_point, &
            tmp_shared_residual_range_low=tmp_shared_residual_range_low, &
            tmp_shared_residual_range_high=tmp_shared_residual_range_high, &
            tmp_occupancy_failed=tmp_occupancy_failed, tmp_n_pooled_residuals=tmp_n_pooled_residuals, &
            tmp_min_bin_occupancy=tmp_min_bin_occupancy, tmp_mean_bin_occupancy=tmp_mean_bin_occupancy, &
            tmp_max_bin_occupancy=tmp_max_bin_occupancy, tmp_sturges_bins=tmp_sturges_bins, &
            tmp_fd_bins=tmp_fd_bins, tmp_pmfs=tmp_pmfs, tmp_counts=tmp_counts, &
            tmp_included_n_reps=tmp_included_n_reps, tmp_mean_pmf=tmp_mean_pmf, &
            tmp_mean_pmf_counts=tmp_mean_pmf_counts, tmp_mean_pmf_included_n_reps=tmp_mean_pmf_included_n_reps, &
            tmp_js_divergences=tmp_js_divergences, tmp_weights=tmp_weights, &
            tmp_global_js_divergence=tmp_global_js_divergence, tmp_pooled_residuals=tmp_pooled_residuals, &
            tmp_pooled_residuals_perm=tmp_pooled_residuals_perm, &
            tmp_residual_to_study_perm=tmp_residual_to_study_perm, tmp_bin_counts_search=tmp_bin_counts_search, &
            tmp_study_block_ends=tmp_study_block_ends, tmp_study_block_sizes=tmp_study_block_sizes, &
            tmp_point_n_neighbors=tmp_point_n_neighbors, tmp_p_values=tmp_p_values, &
            tmp_p_values_observed_consensus=tmp_p_values_observed_consensus, &
            tmp_permutation_mean_pmf=tmp_permutation_mean_pmf, &
            tmp_permutation_mean_pmf_counts=tmp_permutation_mean_pmf_counts, &
            tmp_permutation_counts=tmp_permutation_counts, tmp_permutation_pmfs=tmp_permutation_pmfs, &
            tmp_permutation_js_divergences=tmp_permutation_js_divergences, &
            tmp_permutation_weights=tmp_permutation_weights, &
            tmp_search_n_bins_per_point=tmp_search_n_bins_per_point, &
            tmp_search_shared_residual_range_low=tmp_search_shared_residual_range_low, &
            tmp_search_shared_residual_range_high=tmp_search_shared_residual_range_high, &
            tmp_search_confidence_interval=tmp_search_confidence_interval, &
            tmp_search_trace_n_points=tmp_search_trace_n_points, &
            tmp_search_trace_n_neighbors=tmp_search_trace_n_neighbors, &
            tmp_search_trace_global_js_divergence=tmp_search_trace_global_js_divergence, &
            tmp_search_trace_ci_lower=tmp_search_trace_ci_lower, &
            tmp_search_trace_ci_upper=tmp_search_trace_ci_upper, &
            tmp_search_trace_ci_width=tmp_search_trace_ci_width, &
            tmp_search_trace_ci_width_relative=tmp_search_trace_ci_width_relative, &
            tmp_search_trace_ci_overlap=tmp_search_trace_ci_overlap, &
            tmp_search_trace_delta=tmp_search_trace_delta, &
            tmp_search_trace_delta_median=tmp_search_trace_delta_median, &
            tmp_search_trace_delta_max=tmp_search_trace_delta_max, &
            tmp_search_trace_selected_n_bins=tmp_search_trace_selected_n_bins, &
            tmp_search_trace_occupancy_failed=tmp_search_trace_occupancy_failed, &
            tmp_search_trace_n_pooled_residuals=tmp_search_trace_n_pooled_residuals, &
            tmp_search_trace_min_bin_occupancy=tmp_search_trace_min_bin_occupancy, &
            tmp_search_trace_mean_bin_occupancy=tmp_search_trace_mean_bin_occupancy, &
            tmp_search_trace_max_bin_occupancy=tmp_search_trace_max_bin_occupancy, &
            tmp_search_trace_sturges_bins=tmp_search_trace_sturges_bins, &
            tmp_search_trace_fd_bins=tmp_search_trace_fd_bins, &
            tmp_search_trace_shared_residual_range_low=tmp_search_trace_shared_residual_range_low, &
            tmp_search_trace_shared_residual_range_high=tmp_search_trace_shared_residual_range_high, &
            tmp_candidate_confidence_interval=tmp_candidate_confidence_interval, &
            tmp_bootstrapping_top_k_jsds=tmp_bootstrapping_top_k_jsds, &
            tmp_prev_global_js_divergence=tmp_prev_global_js_divergence, &
            tmp_best_uncertainty_confidence_interval=tmp_best_uncertainty_confidence_interval, &
            tmp_best_n_bins_per_point=tmp_best_n_bins_per_point, &
            tmp_best_uncertainty_n_bins_per_point=tmp_best_uncertainty_n_bins_per_point, &
            tmp_best_shared_residual_range_low=tmp_best_shared_residual_range_low, &
            tmp_best_shared_residual_range_high=tmp_best_shared_residual_range_high, &
            tmp_best_uncertainty_shared_residual_range_low=tmp_best_uncertainty_shared_residual_range_low, &
            tmp_best_uncertainty_shared_residual_range_high=tmp_best_uncertainty_shared_residual_range_high, &
            rerun_parameter_search=rerun_parameter_search, n_bootstraps=n_bootstraps, random_seed=random_seed, &
            ierr=ierr)
    end subroutine direct_knn_label_test

    !> Test helper: calls run_js_comp_test_adaptive_label_permutation_test_impl directly, every work array allocated here.
    subroutine direct_adaptive_label_test(n_studies, max_n_genes_all_studies, max_n_reps_all_studies, expression, &
                                          k_start, k_step, k_max, max_n_points_candidate, &
                                          n_bootstrapping_top_k_jsds, n_permutations, &
                                          global_js_divergence_observed, p_values, jsd_null, &
                                          n_failed_permutations, ierr, rerun_parameter_search, n_bootstraps, &
                                          random_seed)
        integer(int32), intent(in) :: n_studies
            !! Forwarded to run_js_comp_test_adaptive_label_permutation_test_impl
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Forwarded to run_js_comp_test_adaptive_label_permutation_test_impl
        integer(int32), intent(in) :: max_n_reps_all_studies
            !! Forwarded to run_js_comp_test_adaptive_label_permutation_test_impl
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(in) :: expression
            !! Forwarded to run_js_comp_test_adaptive_label_permutation_test_impl
        integer(int32), intent(in) :: k_start
            !! Forwarded to run_js_comp_test_adaptive_label_permutation_test_impl
        integer(int32), intent(in) :: k_step
            !! Forwarded to run_js_comp_test_adaptive_label_permutation_test_impl
        integer(int32), intent(in) :: k_max
            !! Forwarded to run_js_comp_test_adaptive_label_permutation_test_impl
        integer(int32), intent(in) :: max_n_points_candidate
            !! Forwarded to run_js_comp_test_adaptive_label_permutation_test_impl
        integer(int32), intent(in) :: n_bootstrapping_top_k_jsds
            !! Forwarded to run_js_comp_test_adaptive_label_permutation_test_impl
        integer(int32), intent(in) :: n_permutations
            !! Forwarded to run_js_comp_test_adaptive_label_permutation_test_impl
        real(real64), dimension(n_studies), intent(out) :: global_js_divergence_observed
            !! Forwarded to run_js_comp_test_adaptive_label_permutation_test_impl
        real(real64), dimension(n_studies), intent(out) :: p_values
            !! Forwarded to run_js_comp_test_adaptive_label_permutation_test_impl
        real(real64), dimension(n_permutations, n_studies), intent(out) :: jsd_null
            !! Forwarded to run_js_comp_test_adaptive_label_permutation_test_impl
        integer(int32), intent(out) :: n_failed_permutations
            !! Forwarded to run_js_comp_test_adaptive_label_permutation_test_impl
        integer(int32), intent(out) :: ierr
            !! Error code of the implementation
        logical(c_bool), intent(in), optional :: rerun_parameter_search
            !! Forwarded to run_js_comp_test_adaptive_label_permutation_test_impl
        integer(int32), intent(in), optional :: n_bootstraps
            !! Forwarded to run_js_comp_test_adaptive_label_permutation_test_impl
        integer(int32), intent(in), optional :: random_seed
            !! Forwarded to run_js_comp_test_adaptive_label_permutation_test_impl
        integer(int32), allocatable :: tmp_sample_slots(:)
        integer(int32), allocatable :: tmp_sample_assignment(:)
        real(real64), allocatable :: tmp_permuted_global_js_divergence(:)
        logical(c_bool), allocatable :: tmp_sample_mark(:)
        real(real64), allocatable :: tmp_permuted_expression(:, :, :)
        real(real64), allocatable :: tmp_gene_means(:, :)
        real(real64), allocatable :: tmp_residuals(:, :, :)
        integer(int32), allocatable :: tmp_gene_means_perm_all(:)
        real(real64), allocatable :: tmp_x_star(:)
        integer(int32), allocatable :: tmp_pooled_neighborhood_range(:, :)
        integer(int32), allocatable :: tmp_n_neighbors_per_point(:, :)
        integer(int32), allocatable :: tmp_stop_reason(:)
        real(real64), allocatable :: tmp_neighborhood_dispersion(:)
        real(real64), allocatable :: tmp_neighborhood_mad(:)
        integer(int32), allocatable :: tmp_point_neighborhood_indices(:, :)
        real(real64), allocatable :: tmp_neighbor_residuals(:, :)
        integer(int32), allocatable :: tmp_point_pool(:)
        integer(int32), allocatable :: tmp_point_draw(:)
        integer(int32), allocatable :: tmp_n_bins_per_point(:)
        real(real64), allocatable :: tmp_shared_residual_range_low(:)
        real(real64), allocatable :: tmp_shared_residual_range_high(:)
        logical(c_bool), allocatable :: tmp_occupancy_failed(:)
        integer(int32), allocatable :: tmp_n_pooled_residuals(:)
        integer(int32), allocatable :: tmp_min_bin_occupancy(:)
        real(real64), allocatable :: tmp_mean_bin_occupancy(:)
        integer(int32), allocatable :: tmp_max_bin_occupancy(:)
        integer(int32), allocatable :: tmp_sturges_bins(:)
        integer(int32), allocatable :: tmp_fd_bins(:)
        real(real64), allocatable :: tmp_pmfs(:, :, :)
        integer(int32), allocatable :: tmp_counts(:, :, :)
        integer(int32), allocatable :: tmp_included_n_reps(:, :)
        real(real64), allocatable :: tmp_mean_pmf(:, :)
        integer(int32), allocatable :: tmp_mean_pmf_counts(:, :)
        integer(int32), allocatable :: tmp_mean_pmf_included_n_reps(:)
        real(real64), allocatable :: tmp_js_divergences(:, :)
        real(real64), allocatable :: tmp_weights(:, :)
        real(real64), allocatable :: tmp_global_js_divergence(:)
        real(real64), allocatable :: tmp_pooled_residuals(:)
        integer(int32), allocatable :: tmp_pooled_residuals_perm(:)
        integer(int32), allocatable :: tmp_residual_to_study_perm(:)
        integer(int32), allocatable :: tmp_bin_counts_search(:)
        integer(int32), allocatable :: tmp_study_block_ends(:)
        integer(int32), allocatable :: tmp_study_block_sizes(:)
        real(real64), allocatable :: tmp_p_values(:)
        real(real64), allocatable :: tmp_p_values_observed_consensus(:)
        real(real64), allocatable :: tmp_permutation_mean_pmf(:, :)
        integer(int32), allocatable :: tmp_permutation_mean_pmf_counts(:, :)
        integer(int32), allocatable :: tmp_permutation_counts(:, :)
        real(real64), allocatable :: tmp_permutation_pmfs(:, :, :)
        real(real64), allocatable :: tmp_permutation_js_divergences(:, :)
        real(real64), allocatable :: tmp_permutation_weights(:, :)
        integer(int32), allocatable :: tmp_search_n_bins_per_point(:)
        real(real64), allocatable :: tmp_search_shared_residual_range_low(:)
        real(real64), allocatable :: tmp_search_shared_residual_range_high(:)
        real(real64), allocatable :: tmp_search_confidence_interval(:, :)
        integer(int32), allocatable :: tmp_search_trace_k_start(:)
        integer(int32), allocatable :: tmp_search_trace_k_step(:)
        integer(int32), allocatable :: tmp_search_trace_k_max(:)
        integer(int32), allocatable :: tmp_search_trace_n_points(:)
        real(real64), allocatable :: tmp_search_trace_global_js_divergence(:, :)
        real(real64), allocatable :: tmp_search_trace_ci_lower(:, :)
        real(real64), allocatable :: tmp_search_trace_ci_upper(:, :)
        real(real64), allocatable :: tmp_search_trace_ci_width(:, :)
        real(real64), allocatable :: tmp_search_trace_ci_width_relative(:, :)
        real(real64), allocatable :: tmp_search_trace_ci_overlap(:, :)
        real(real64), allocatable :: tmp_search_trace_delta(:, :)
        real(real64), allocatable :: tmp_search_trace_delta_median(:)
        real(real64), allocatable :: tmp_search_trace_delta_max(:)
        integer(int32), allocatable :: tmp_search_trace_selected_n_bins(:, :)
        logical(c_bool), allocatable :: tmp_search_trace_occupancy_failed(:, :)
        integer(int32), allocatable :: tmp_search_trace_n_pooled_residuals(:, :)
        integer(int32), allocatable :: tmp_search_trace_min_bin_occupancy(:, :)
        real(real64), allocatable :: tmp_search_trace_mean_bin_occupancy(:, :)
        integer(int32), allocatable :: tmp_search_trace_max_bin_occupancy(:, :)
        integer(int32), allocatable :: tmp_search_trace_sturges_bins(:, :)
        integer(int32), allocatable :: tmp_search_trace_fd_bins(:, :)
        real(real64), allocatable :: tmp_search_trace_shared_residual_range_low(:, :)
        real(real64), allocatable :: tmp_search_trace_shared_residual_range_high(:, :)
        integer(int32), allocatable :: tmp_search_candidate_k_start(:)
        integer(int32), allocatable :: tmp_search_candidate_n_points(:)
        integer(int32), allocatable :: tmp_search_candidate_status(:)
        real(real64), allocatable :: tmp_candidate_confidence_interval(:, :)
        real(real64), allocatable :: tmp_bootstrapping_top_k_jsds(:, :, :)
        real(real64), allocatable :: tmp_prev_global_js_divergence(:)
        integer(int32), allocatable :: tmp_delta_perm(:)
        real(real64), allocatable :: tmp_best_uncertainty_confidence_interval(:, :)
        integer(int32), allocatable :: tmp_best_n_bins_per_point(:)
        integer(int32), allocatable :: tmp_best_uncertainty_n_bins_per_point(:)
        real(real64), allocatable :: tmp_best_shared_residual_range_low(:)
        real(real64), allocatable :: tmp_best_shared_residual_range_high(:)
        real(real64), allocatable :: tmp_best_uncertainty_shared_residual_range_low(:)
        real(real64), allocatable :: tmp_best_uncertainty_shared_residual_range_high(:)
        allocate (tmp_sample_slots(max_n_reps_all_studies*n_studies))
        allocate (tmp_sample_assignment(max_n_reps_all_studies*n_studies))
        allocate (tmp_permuted_global_js_divergence(n_studies))
        allocate (tmp_sample_mark(max_n_reps_all_studies*n_studies))
        allocate (tmp_permuted_expression(max_n_reps_all_studies, max_n_genes_all_studies, n_studies))
        allocate (tmp_gene_means(max_n_genes_all_studies, n_studies))
        allocate (tmp_residuals(max_n_reps_all_studies, max_n_genes_all_studies, n_studies))
        allocate (tmp_gene_means_perm_all(max_n_genes_all_studies*n_studies))
        allocate (tmp_x_star(max_n_genes_all_studies*n_studies))
        allocate (tmp_pooled_neighborhood_range(2, max_n_genes_all_studies*n_studies))
        allocate (tmp_n_neighbors_per_point(n_studies, max_n_genes_all_studies*n_studies))
        allocate (tmp_stop_reason(max_n_genes_all_studies*n_studies))
        allocate (tmp_neighborhood_dispersion(max_n_genes_all_studies*n_studies))
        allocate (tmp_neighborhood_mad(max_n_genes_all_studies*n_studies))
        allocate (tmp_point_neighborhood_indices(max_n_genes_all_studies, n_studies))
        allocate (tmp_neighbor_residuals(max_n_reps_all_studies, max_n_genes_all_studies))
        allocate (tmp_point_pool(MAX_N_BINS))
        allocate (tmp_point_draw(MAX_N_BINS))
        allocate (tmp_n_bins_per_point(max_n_points_candidate))
        allocate (tmp_shared_residual_range_low(max_n_points_candidate))
        allocate (tmp_shared_residual_range_high(max_n_points_candidate))
        allocate (tmp_occupancy_failed(max_n_points_candidate))
        allocate (tmp_n_pooled_residuals(max_n_points_candidate))
        allocate (tmp_min_bin_occupancy(max_n_points_candidate))
        allocate (tmp_mean_bin_occupancy(max_n_points_candidate))
        allocate (tmp_max_bin_occupancy(max_n_points_candidate))
        allocate (tmp_sturges_bins(max_n_points_candidate))
        allocate (tmp_fd_bins(max_n_points_candidate))
        allocate (tmp_pmfs(max_n_points_candidate, MAX_N_BINS, n_studies))
        allocate (tmp_counts(max_n_points_candidate, MAX_N_BINS, n_studies))
        allocate (tmp_included_n_reps(max_n_points_candidate, n_studies))
        allocate (tmp_mean_pmf(max_n_points_candidate, MAX_N_BINS))
        allocate (tmp_mean_pmf_counts(max_n_points_candidate, MAX_N_BINS))
        allocate (tmp_mean_pmf_included_n_reps(max_n_points_candidate))
        allocate (tmp_js_divergences(max_n_points_candidate, n_studies))
        allocate (tmp_weights(max_n_points_candidate, n_studies))
        allocate (tmp_global_js_divergence(n_studies))
        allocate (tmp_pooled_residuals(max_n_reps_all_studies*max_n_genes_all_studies*n_studies))
        allocate (tmp_pooled_residuals_perm(max_n_reps_all_studies*max_n_genes_all_studies*n_studies))
        allocate (tmp_residual_to_study_perm(max_n_reps_all_studies*max_n_genes_all_studies*n_studies))
        allocate (tmp_bin_counts_search(MAX_N_BINS))
        allocate (tmp_study_block_ends(n_studies))
        allocate (tmp_study_block_sizes(n_studies))
        allocate (tmp_p_values(n_studies))
        allocate (tmp_p_values_observed_consensus(n_studies))
        allocate (tmp_permutation_mean_pmf(max_n_points_candidate, MAX_N_BINS))
        allocate (tmp_permutation_mean_pmf_counts(max_n_points_candidate, MAX_N_BINS))
        allocate (tmp_permutation_counts(max_n_points_candidate, MAX_N_BINS))
        allocate (tmp_permutation_pmfs(max_n_points_candidate, MAX_N_BINS, n_studies))
        allocate (tmp_permutation_js_divergences(max_n_points_candidate, n_studies))
        allocate (tmp_permutation_weights(max_n_points_candidate, n_studies))
        allocate (tmp_search_n_bins_per_point(max_n_points_candidate))
        allocate (tmp_search_shared_residual_range_low(max_n_points_candidate))
        allocate (tmp_search_shared_residual_range_high(max_n_points_candidate))
        allocate (tmp_search_confidence_interval(2, n_studies))
        allocate (tmp_search_trace_k_start(MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_k_step(MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_k_max(MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_n_points(MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_global_js_divergence(n_studies, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_ci_lower(n_studies, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_ci_upper(n_studies, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_ci_width(n_studies, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_ci_width_relative(n_studies, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_ci_overlap(n_studies, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_delta(n_studies, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_delta_median(MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_delta_max(MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_selected_n_bins(max_n_points_candidate, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_occupancy_failed(max_n_points_candidate, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_n_pooled_residuals(max_n_points_candidate, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_min_bin_occupancy(max_n_points_candidate, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_mean_bin_occupancy(max_n_points_candidate, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_max_bin_occupancy(max_n_points_candidate, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_sturges_bins(max_n_points_candidate, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_fd_bins(max_n_points_candidate, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_shared_residual_range_low(max_n_points_candidate, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_shared_residual_range_high(max_n_points_candidate, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_candidate_k_start(MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_candidate_n_points(MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_candidate_status(MAX_CANDIDATE_PAIRS))
        allocate (tmp_candidate_confidence_interval(2, n_studies))
        allocate (tmp_bootstrapping_top_k_jsds(n_bootstrapping_top_k_jsds, 2, n_studies))
        allocate (tmp_prev_global_js_divergence(n_studies))
        allocate (tmp_delta_perm(n_studies))
        allocate (tmp_best_uncertainty_confidence_interval(2, n_studies))
        allocate (tmp_best_n_bins_per_point(max_n_points_candidate))
        allocate (tmp_best_uncertainty_n_bins_per_point(max_n_points_candidate))
        allocate (tmp_best_shared_residual_range_low(max_n_points_candidate))
        allocate (tmp_best_shared_residual_range_high(max_n_points_candidate))
        allocate (tmp_best_uncertainty_shared_residual_range_low(max_n_points_candidate))
        allocate (tmp_best_uncertainty_shared_residual_range_high(max_n_points_candidate))
        call run_js_comp_test_adaptive_label_permutation_test_impl(n_studies=n_studies, &
            max_n_genes_all_studies=max_n_genes_all_studies, max_n_reps_all_studies=max_n_reps_all_studies, &
            expression=expression, k_start=k_start, k_step=k_step, k_max=k_max, &
            max_n_points_candidate=max_n_points_candidate, n_bootstrapping_top_k_jsds=n_bootstrapping_top_k_jsds, &
            n_permutations=n_permutations, global_js_divergence_observed=global_js_divergence_observed, &
            p_values=p_values, jsd_null=jsd_null, n_failed_permutations=n_failed_permutations, &
            tmp_sample_slots=tmp_sample_slots, tmp_sample_assignment=tmp_sample_assignment, &
            tmp_permuted_global_js_divergence=tmp_permuted_global_js_divergence, tmp_sample_mark=tmp_sample_mark, &
            tmp_permuted_expression=tmp_permuted_expression, tmp_gene_means=tmp_gene_means, &
            tmp_residuals=tmp_residuals, tmp_gene_means_perm_all=tmp_gene_means_perm_all, tmp_x_star=tmp_x_star, &
            tmp_pooled_neighborhood_range=tmp_pooled_neighborhood_range, &
            tmp_n_neighbors_per_point=tmp_n_neighbors_per_point, tmp_stop_reason=tmp_stop_reason, &
            tmp_neighborhood_dispersion=tmp_neighborhood_dispersion, tmp_neighborhood_mad=tmp_neighborhood_mad, &
            tmp_point_neighborhood_indices=tmp_point_neighborhood_indices, &
            tmp_neighbor_residuals=tmp_neighbor_residuals, tmp_point_pool=tmp_point_pool, &
            tmp_point_draw=tmp_point_draw, tmp_n_bins_per_point=tmp_n_bins_per_point, &
            tmp_shared_residual_range_low=tmp_shared_residual_range_low, &
            tmp_shared_residual_range_high=tmp_shared_residual_range_high, &
            tmp_occupancy_failed=tmp_occupancy_failed, tmp_n_pooled_residuals=tmp_n_pooled_residuals, &
            tmp_min_bin_occupancy=tmp_min_bin_occupancy, tmp_mean_bin_occupancy=tmp_mean_bin_occupancy, &
            tmp_max_bin_occupancy=tmp_max_bin_occupancy, tmp_sturges_bins=tmp_sturges_bins, &
            tmp_fd_bins=tmp_fd_bins, tmp_pmfs=tmp_pmfs, tmp_counts=tmp_counts, &
            tmp_included_n_reps=tmp_included_n_reps, tmp_mean_pmf=tmp_mean_pmf, &
            tmp_mean_pmf_counts=tmp_mean_pmf_counts, tmp_mean_pmf_included_n_reps=tmp_mean_pmf_included_n_reps, &
            tmp_js_divergences=tmp_js_divergences, tmp_weights=tmp_weights, &
            tmp_global_js_divergence=tmp_global_js_divergence, tmp_pooled_residuals=tmp_pooled_residuals, &
            tmp_pooled_residuals_perm=tmp_pooled_residuals_perm, &
            tmp_residual_to_study_perm=tmp_residual_to_study_perm, tmp_bin_counts_search=tmp_bin_counts_search, &
            tmp_study_block_ends=tmp_study_block_ends, tmp_study_block_sizes=tmp_study_block_sizes, &
            tmp_p_values=tmp_p_values, tmp_p_values_observed_consensus=tmp_p_values_observed_consensus, &
            tmp_permutation_mean_pmf=tmp_permutation_mean_pmf, &
            tmp_permutation_mean_pmf_counts=tmp_permutation_mean_pmf_counts, &
            tmp_permutation_counts=tmp_permutation_counts, tmp_permutation_pmfs=tmp_permutation_pmfs, &
            tmp_permutation_js_divergences=tmp_permutation_js_divergences, &
            tmp_permutation_weights=tmp_permutation_weights, &
            tmp_search_n_bins_per_point=tmp_search_n_bins_per_point, &
            tmp_search_shared_residual_range_low=tmp_search_shared_residual_range_low, &
            tmp_search_shared_residual_range_high=tmp_search_shared_residual_range_high, &
            tmp_search_confidence_interval=tmp_search_confidence_interval, &
            tmp_search_trace_k_start=tmp_search_trace_k_start, tmp_search_trace_k_step=tmp_search_trace_k_step, &
            tmp_search_trace_k_max=tmp_search_trace_k_max, tmp_search_trace_n_points=tmp_search_trace_n_points, &
            tmp_search_trace_global_js_divergence=tmp_search_trace_global_js_divergence, &
            tmp_search_trace_ci_lower=tmp_search_trace_ci_lower, &
            tmp_search_trace_ci_upper=tmp_search_trace_ci_upper, &
            tmp_search_trace_ci_width=tmp_search_trace_ci_width, &
            tmp_search_trace_ci_width_relative=tmp_search_trace_ci_width_relative, &
            tmp_search_trace_ci_overlap=tmp_search_trace_ci_overlap, &
            tmp_search_trace_delta=tmp_search_trace_delta, &
            tmp_search_trace_delta_median=tmp_search_trace_delta_median, &
            tmp_search_trace_delta_max=tmp_search_trace_delta_max, &
            tmp_search_trace_selected_n_bins=tmp_search_trace_selected_n_bins, &
            tmp_search_trace_occupancy_failed=tmp_search_trace_occupancy_failed, &
            tmp_search_trace_n_pooled_residuals=tmp_search_trace_n_pooled_residuals, &
            tmp_search_trace_min_bin_occupancy=tmp_search_trace_min_bin_occupancy, &
            tmp_search_trace_mean_bin_occupancy=tmp_search_trace_mean_bin_occupancy, &
            tmp_search_trace_max_bin_occupancy=tmp_search_trace_max_bin_occupancy, &
            tmp_search_trace_sturges_bins=tmp_search_trace_sturges_bins, &
            tmp_search_trace_fd_bins=tmp_search_trace_fd_bins, &
            tmp_search_trace_shared_residual_range_low=tmp_search_trace_shared_residual_range_low, &
            tmp_search_trace_shared_residual_range_high=tmp_search_trace_shared_residual_range_high, &
            tmp_search_candidate_k_start=tmp_search_candidate_k_start, &
            tmp_search_candidate_n_points=tmp_search_candidate_n_points, &
            tmp_search_candidate_status=tmp_search_candidate_status, &
            tmp_candidate_confidence_interval=tmp_candidate_confidence_interval, &
            tmp_bootstrapping_top_k_jsds=tmp_bootstrapping_top_k_jsds, &
            tmp_prev_global_js_divergence=tmp_prev_global_js_divergence, tmp_delta_perm=tmp_delta_perm, &
            tmp_best_uncertainty_confidence_interval=tmp_best_uncertainty_confidence_interval, &
            tmp_best_n_bins_per_point=tmp_best_n_bins_per_point, &
            tmp_best_uncertainty_n_bins_per_point=tmp_best_uncertainty_n_bins_per_point, &
            tmp_best_shared_residual_range_low=tmp_best_shared_residual_range_low, &
            tmp_best_shared_residual_range_high=tmp_best_shared_residual_range_high, &
            tmp_best_uncertainty_shared_residual_range_low=tmp_best_uncertainty_shared_residual_range_low, &
            tmp_best_uncertainty_shared_residual_range_high=tmp_best_uncertainty_shared_residual_range_high, &
            rerun_parameter_search=rerun_parameter_search, n_bootstraps=n_bootstraps, random_seed=random_seed, &
            ierr=ierr)
    end subroutine direct_adaptive_label_test

end module mod_test_data_integration_js_comp_test_permutation
