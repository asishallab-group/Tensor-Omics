!> Shared fixtures, reference helpers and tests of the study-label Permutation-Test of the JS-Comp-Test: the data
!| and reference helpers every suite uses, the forwarding tests of options that both the k-NN and the adaptive
!| variant take, the sizing routines and the two-study case. The variant suites are
!| `mod_test_data_integration_js_comp_test_permutation_fixed_k` and `mod_test_data_integration_js_comp_test_permutation_adaptive`.
module mod_test_data_integration_js_comp_test_permutation_common
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: iso_c_binding, only: c_bool
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_nan
    use tox_errors
    use test_suite, only: test_case
    use f42_random_gsl, only: rng_t, create_rng, destroy_rng, random_uniform, random_shuffle
    use f42_sort_impl, only: init_perm
    use tox_data_integration_preprocessing, only: compute_gene_means, compute_residuals, pool_means, &
                                                  construct_adaptive_neighborhoods
    use tox_data_integration_js_comp_test_run, only: run_js_comp_test, run_js_comp_test_adaptive
    use tox_data_integration_js_comp_test_fixed_k_search, only: run_js_comp_test_parameter_search, &
                                                                calc_js_comp_test_candidate_bounds
    use tox_data_integration_js_comp_test_adaptive_search, only: run_js_comp_test_adaptive_parameter_search, &
                                                                 calc_adaptive_js_comp_test_bounds
    use tox_data_integration_js_comp_test_impl, only: MAX_N_BINS, MAX_CANDIDATE_PAIRS, METHOD_JOIN_MIN
    use tox_data_integration_jsd_impl, only: MODE_DIVERGENCE_SYMMETRIC, MODE_DIVERGENCE_PAIRWISE_K2
    use tox_data_integration_js_comp_test_binning_impl, only: MODE_BIN_COUNT_OCCUPANCY
    use mod_test_data_integration_js_comp_test, only: assert_equal_array_real_nan_exact
    use tox_data_integration_js_comp_test_permutation_fixed_k, only: run_js_comp_test_on_sample_assignment, &
                                                                     run_js_comp_test_permutation
    use tox_data_integration_js_comp_test_permutation_adaptive, only: run_js_comp_test_adaptive_on_sample_assignment, &
                                                                      run_js_comp_test_adaptive_permutation
    use tox_data_integration_js_comp_test_permutation_bounds, only: calc_js_comp_test_permutation_bounds, &
                                                                    calc_adaptive_js_comp_test_permutation_bounds
    use tox_data_integration_js_comp_test_permutation_common_impl, only: label_permutation_p_values


    implicit none

    ! The fixed parameters every k-NN test uses, and the adaptive growth knobs.
    integer(int32), parameter :: N_POINTS_KNN = 10_int32, N_NEIGHBORS_KNN = 30_int32
    integer(int32), parameter :: K_START = 30_int32, K_STEP = 10_int32, K_MAX = 90_int32
    ! The null/batch tests: 4 null data sets of 3 studies give 12 p-values over n = 19 permutations.
    ! Under the null each p-value is uniform on {1/20, ..., 20/20}, so it sits at the minimum 1/20
    ! with probability 1/20, and the number at the minimum is about Binomial(12, 0.05): at most 3
    ! of the 12 with probability 0.998 (if the p-values were independent). So at least 9 above the
    ! minimum is a loose bar that a working test passes, while a test whose null rows always fall
    ! below the observed value (every p at the minimum) fails it.
    integer(int32), parameter :: N_NULL_P_VALUES = 12_int32, MIN_NULL_ABOVE_MINIMUM = 9_int32
    ! The k-NN forwarding fixture (make_knn_grid_fixture): G = 6400 genes give the default
    ! large-neighborhood grid `n_points_high = ceiling(4*sqrt(6400)) = 320`, then 256, 204.8 -> 205,
    ! 163.84 < 200 stops, plus the floor 200: four candidates (320, 80), (256, 100), (205, ...),
    ! (200, 128), and the small grid has the same `n_points` with half the `n_neighbors`. Data seed
    ! 51 is pinned: on it the k-NN search with effect-size thresholds 1e-300 selects (205, ...) with
    ! `delta_abs_tolerance = 1` and (200, 128) with 0 (found by running seeds 51, 52, 53).
    integer(int32), parameter :: KNN_GRID_GENES = 6400_int32, KNN_GRID_SEED = 51_int32

    !> The options of the plain run that no routine of the suite names elsewhere: an unallocated component is an
    !| absent argument (the routine's default), an allocated one is forwarded.
    type :: option_set_t
        integer(int32), allocatable :: bin_count_mode, m_min, m_max, min_residuals_per_bin
        real(real64), allocatable :: gamma_occupancy, lower_residual_range_quantile, upper_residual_range_quantile
        real(real64), allocatable :: tau, mad_distance_factor
        integer(int32), allocatable :: max_pooled_residuals, min_study_neighbors
    end type option_set_t

contains

    !> Get array of all available tests.
    function get_all_tests_js_comp_test_permutation_common() result(all_tests)
        type(test_case), allocatable :: all_tests(:)
        allocate (all_tests(7))
        all_tests(1) = test_case("test_rerun_forwards_delta_abs_tolerance", test_rerun_forwards_delta_abs_tolerance)
        all_tests(2) = test_case("test_forwards_divergence_mode", test_forwards_divergence_mode)
        all_tests(3) = test_case("test_rerun_forwards_divergence_mode", test_rerun_forwards_divergence_mode)
        all_tests(4) = test_case("test_permutation_bounds_by_hand", test_permutation_bounds_by_hand)
        all_tests(5) = test_case("test_two_studies_default_divergence_null_and_batch", &
                                 test_two_studies_default_divergence_null_and_batch)
        all_tests(6) = test_case("test_label_permutation_p_values_by_hand", test_label_permutation_p_values_by_hand)
        all_tests(7) = test_case("test_forwards_binning_options", test_forwards_binning_options)
    end function get_all_tests_js_comp_test_permutation_common

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
        logical(c_bool), intent(in), optional :: dyadic
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
    !| from the pooled means, then run_js_comp_test (`divergence_mode` forwarded, absent: its default).
    subroutine reference_knn(expression, n_points, n_neighbors, global_jsd, ierr, divergence_mode, options)
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
        integer(int32), intent(in), optional :: divergence_mode
            !! Forwarded to run_js_comp_test
        type(option_set_t), intent(in), optional :: options
            !! The binning options forwarded to run_js_comp_test (absent: all defaults)
        type(option_set_t) :: active
        real(real64), allocatable :: gene_means(:, :), residuals(:, :, :), x_star(:), range_low(:), range_high(:)
        real(real64), allocatable :: mean_occ(:), pmfs(:, :, :), mean_pmf(:, :), js(:, :), weights(:, :)
        integer(int32), allocatable :: indices(:, :, :), ranges(:, :, :), n_bins(:), n_pooled(:), min_occ(:), max_occ(:)
        integer(int32), allocatable :: sturges(:), fd(:), counts(:, :, :), included(:, :), mean_counts(:, :)
        integer(int32), allocatable :: mean_included(:)
        logical(c_bool), allocatable :: occ_failed(:)
        integer(int32) :: n_studies, n_genes, n_reps, n_pool, n_bins_used

        if (present(options)) active = options
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
                  mean_pmf(n_points, MAX_N_BINS), js(n_points, n_studies), weights(n_points, n_studies))
        call run_js_comp_test(n_studies, n_genes, n_reps, n_points, n_neighbors, gene_means, residuals, x_star, indices, &
                              ranges, n_bins, range_low, range_high, n_bins_used, occ_failed, n_pooled, min_occ, mean_occ, &
                              max_occ, sturges, fd, pmfs, counts, included, mean_pmf, mean_counts, mean_included, js, &
                              weights, global_jsd, min_residuals_per_bin=active%min_residuals_per_bin, m_min=active%m_min, &
                              m_max=active%m_max, gamma_occupancy=active%gamma_occupancy, &
                              lower_residual_range_quantile=active%lower_residual_range_quantile, &
                              upper_residual_range_quantile=active%upper_residual_range_quantile, &
                              bin_count_mode=active%bin_count_mode, divergence_mode=divergence_mode, ierr=ierr)
    end subroutine reference_knn

    !> Test helper: the plain adaptive pipeline a caller runs today: means/residuals per study,
    !| construct_adaptive_neighborhoods, then run_js_comp_test_adaptive (`divergence_mode` forwarded,
    !| absent: its default).
    subroutine reference_adaptive(expression, k_start, k_step, k_max, global_jsd, n_points, ierr, divergence_mode, options)
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
        integer(int32), intent(in), optional :: divergence_mode
            !! Forwarded to run_js_comp_test_adaptive
        type(option_set_t), intent(in), optional :: options
            !! The construction and binning options forwarded (absent: all defaults)
        type(option_set_t) :: active
        real(real64), allocatable :: gene_means(:, :), residuals(:, :, :), x_star(:), dispersion(:), mads(:)
        real(real64), allocatable :: range_low(:), range_high(:), mean_occ(:), pmfs(:, :, :), mean_pmf(:, :)
        real(real64), allocatable :: js(:, :), weights(:, :)
        integer(int32), allocatable :: ranges(:, :), n_per_point(:, :), stops(:), n_bins(:), n_pooled(:), min_occ(:)
        integer(int32), allocatable :: max_occ(:), sturges(:), fd(:), counts(:, :, :), included(:, :), mean_counts(:, :)
        integer(int32), allocatable :: mean_included(:), run_n_per_point(:, :)
        logical(c_bool), allocatable :: occ_failed(:)
        integer(int32) :: n_studies, n_genes, n_reps, n_pool, max_nn, status, n_bins_used

        if (present(options)) active = options
        call reference_means_residuals(expression, gene_means, residuals)
        n_reps = size(expression, 1, kind=int32)
        n_genes = size(expression, 2, kind=int32)
        n_studies = size(expression, 3, kind=int32)
        n_pool = n_genes*n_studies
        allocate (x_star(n_pool), ranges(2, n_pool), n_per_point(n_studies, n_pool), stops(n_pool), dispersion(n_pool), &
                  mads(n_pool))
        call construct_adaptive_neighborhoods(n_studies, n_genes, n_reps, gene_means, residuals, k_start, k_step, k_max, &
                                              n_points, x_star, ranges, n_per_point, stops, dispersion, mads, max_nn, &
                                              status, tau=active%tau, mad_distance_factor=active%mad_distance_factor, &
                                              max_pooled_residuals=active%max_pooled_residuals, &
                                              min_study_neighbors=active%min_study_neighbors, ierr=ierr)
        if (is_err(ierr)) return
        allocate (run_n_per_point(n_studies, n_points), n_bins(n_points), n_pooled(n_points), min_occ(n_points), &
                  max_occ(n_points), sturges(n_points), fd(n_points), counts(n_points, MAX_N_BINS, n_studies), &
                  included(n_points, n_studies), mean_counts(n_points, MAX_N_BINS), mean_included(n_points), &
                  occ_failed(n_points), range_low(n_points), range_high(n_points), mean_occ(n_points), &
                  pmfs(n_points, MAX_N_BINS, n_studies), mean_pmf(n_points, MAX_N_BINS), js(n_points, n_studies), &
                  weights(n_points, n_studies))
        call run_js_comp_test_adaptive(n_studies, n_genes, n_reps, n_points, gene_means, residuals, ranges(:, 1:n_points), &
                                       run_n_per_point, n_bins, range_low, range_high, n_bins_used, occ_failed, n_pooled, &
                                       min_occ, mean_occ, max_occ, sturges, fd, pmfs, counts, included, mean_pmf, &
                                       mean_counts, mean_included, js, weights, global_jsd, &
                                       min_residuals_per_bin=active%min_residuals_per_bin, m_min=active%m_min, &
                                       m_max=active%m_max, gamma_occupancy=active%gamma_occupancy, &
                                       lower_residual_range_quantile=active%lower_residual_range_quantile, &
                                       upper_residual_range_quantile=active%upper_residual_range_quantile, &
                                       bin_count_mode=active%bin_count_mode, divergence_mode=divergence_mode, ierr=ierr)
    end subroutine reference_adaptive

    !> Test helper: run_js_comp_test_parameter_search with `n_bootstraps`, the given optionals and
    !| otherwise default options; returns the `(n_points, n_neighbors)` it selects and, if asked,
    !| the number of admissible candidates and their observed global JSDs.
    subroutine reference_fixed_k_search(gene_means, residuals, n_bootstraps, n_points, n_neighbors, &
                                        delta_abs_tolerance, delta_median_threshold, delta_max_threshold, &
                                        plateau_mode, divergence_mode, grid_points_factor, grid_min_points, &
                                        neighborhood_size_mode, delta_min_consecutive_transitions, &
                                        n_admissible_out, trace_jsd_out)
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
        real(real64), intent(in), optional :: delta_abs_tolerance
            !! Forwarded to the search as given
        real(real64), intent(in), optional :: delta_median_threshold
            !! Forwarded to the search as given
        real(real64), intent(in), optional :: delta_max_threshold
            !! Forwarded to the search as given
        integer(int32), intent(in), optional :: plateau_mode
            !! Forwarded to the search as given
        integer(int32), intent(in), optional :: divergence_mode
            !! Forwarded to the search as given
        real(real64), intent(in), optional :: grid_points_factor
            !! Forwarded to the search and its candidate bounds as given
        integer(int32), intent(in), optional :: grid_min_points
            !! Forwarded to the search and its candidate bounds as given
        integer(int32), intent(in), optional :: neighborhood_size_mode
            !! Forwarded to the search as given
        integer(int32), intent(in), optional :: delta_min_consecutive_transitions
            !! Forwarded to the search as given
        integer(int32), intent(out), optional :: n_admissible_out
            !! Number of admissible candidates the search evaluated
        real(real64), intent(out), optional :: trace_jsd_out(:, :)
            !! `(n_studies, MAX_CANDIDATE_PAIRS)`: column `t` is admissible candidate `t`'s observed
            !! global JSD per study (`1:n_admissible_out` meaningful)
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
        call calc_js_comp_test_candidate_bounds(n_genes, max_p, max_k, grid_points_factor=grid_points_factor, &
                                                grid_min_points=grid_min_points, ierr=ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "reference_fixed_k_search: bounds ierr")
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
                                               t_min, t_mean, t_max, t_st, t_fd, t_rl, t_rh, &
                                               delta_median_threshold=delta_median_threshold, &
                                               delta_max_threshold=delta_max_threshold, &
                                               delta_abs_tolerance=delta_abs_tolerance, plateau_mode=plateau_mode, &
                                               divergence_mode=divergence_mode, grid_points_factor=grid_points_factor, &
                                               grid_min_points=grid_min_points, &
                                               neighborhood_size_mode=neighborhood_size_mode, &
                                               delta_min_consecutive_transitions=delta_min_consecutive_transitions, &
                                               ierr=ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "reference_fixed_k_search: ierr")
        if (present(n_admissible_out)) n_admissible_out = n_admissible
        if (present(trace_jsd_out)) trace_jsd_out = t_jsd
    end subroutine reference_fixed_k_search

    !> Test helper: run_js_comp_test_adaptive_parameter_search with `n_bootstraps` and otherwise
    !| default options and its own capacity; returns the knobs it selects and, if asked, the
    !| number of admissible candidates and their observed global JSDs.
    subroutine reference_adaptive_search(gene_means, residuals, n_bootstraps, k_start, k_step, k_max, &
                                         delta_abs_tolerance, delta_median_threshold, delta_max_threshold, &
                                         min_neighbor_overlap, plateau_mode, divergence_mode, &
                                         delta_min_consecutive_transitions, n_admissible_out, trace_jsd_out)
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
        real(real64), intent(in), optional :: delta_abs_tolerance
            !! Forwarded to the search as given
        real(real64), intent(in), optional :: delta_median_threshold
            !! Forwarded to the search as given
        real(real64), intent(in), optional :: delta_max_threshold
            !! Forwarded to the search as given
        real(real64), intent(in), optional :: min_neighbor_overlap
            !! Forwarded to the search as given
        integer(int32), intent(in), optional :: plateau_mode
            !! Forwarded to the search as given
        integer(int32), intent(in), optional :: divergence_mode
            !! Forwarded to the search as given
        integer(int32), intent(in), optional :: delta_min_consecutive_transitions
            !! Forwarded to the search as given
        integer(int32), intent(out), optional :: n_admissible_out
            !! Number of admissible candidates the search evaluated
        real(real64), intent(out), optional :: trace_jsd_out(:, :)
            !! `(n_studies, MAX_CANDIDATE_PAIRS)`: column `t` is admissible candidate `t`'s observed
            !! global JSD per study (`1:n_admissible_out` meaningful)
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
                                                        t_rl, t_rh, n_tried, c_ks, c_np, c_status, &
                                                        delta_median_threshold=delta_median_threshold, &
                                                        delta_max_threshold=delta_max_threshold, &
                                                        delta_abs_tolerance=delta_abs_tolerance, &
                                                        min_neighbor_overlap=min_neighbor_overlap, plateau_mode=plateau_mode, &
                                                        divergence_mode=divergence_mode, &
                                                        delta_min_consecutive_transitions=delta_min_consecutive_transitions, &
                                                        ierr=ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "reference_adaptive_search: ierr")
        if (present(n_admissible_out)) n_admissible_out = n_admissible
        if (present(trace_jsd_out)) trace_jsd_out = t_jsd
    end subroutine reference_adaptive_search

    !> Test helper: replays the label draws of the Permutation-Test through the same public calls the
    !| test makes (one GSL generator seeded with `seed`; each draw `init_perm`, then `random_shuffle`):
    !| column `i` of `assignments` is draw `i`. Not an independent oracle of the drawing rule (the
    !| shuffle itself is tested in mod_test_random_gsl); it checks the plumbing: draw `i` is used for
    !| row `i`, one generator serves all draws in order, and the seed reaches it. Each column is
    !| checked to be a permutation of `1..n_samples`.
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
        integer(int32) :: ierr, i_draw

        rng = create_rng(ierr, seed)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "redraw_assignments: rng")
        do i_draw = 1, n_draws
            call init_perm(assignments(:, i_draw))
            call random_shuffle(rng, n_samples, assignments(:, i_draw))
            call assert_permutation(assignments(:, i_draw), n_samples, "redraw_assignments: draw is a permutation")
        end do
        call destroy_rng(rng)
    end subroutine redraw_assignments

    !> Test helper: the identity assignment of `n` samples (init_perm as a function, for use in calls).
    pure function identity(n) result(perm)
        integer(int32), intent(in) :: n
            !! Number of samples
        integer(int32) :: perm(n)

        call init_perm(perm)
    end function identity

    !> Test helper: `.true.` iff `a` and `b` have the same size and are equal element by element; any NaN makes
    !| it `.false.` (unlike assert_equal_array_real_nan_exact, so a comparison of two failed results cannot pass).
    pure logical function equal_and_not_nan(a, b)
        real(real64), intent(in) :: a(:)
            !! First array
        real(real64), intent(in) :: b(:)
            !! Second array

        equal_and_not_nan = size(a) == size(b)
        if (.not. equal_and_not_nan) return
        equal_and_not_nan = all(a == b)
    end function equal_and_not_nan

    !> Test helper: the k-NN forwarding fixture, three studies of 5, 4 and 4 samples and
    !| KNN_GRID_GENES genes (data seed KNN_GRID_SEED, see there), whose default k-NN grid has four
    !| candidates.
    subroutine make_knn_grid_fixture(expression)
        real(real64), allocatable, intent(out) :: expression(:, :, :)
            !! The data

        call make_expression(KNN_GRID_GENES, [5_int32, 4_int32, 4_int32], [0.5_real64, 0.5_real64, 0.6_real64], &
                             KNN_GRID_SEED, expression)
    end subroutine make_knn_grid_fixture

    !> Test helper: the largest absolute change across studies of the observed global JSD from
    !| admissible candidate `i_transition` to the next.
    function largest_change(jsd, i_transition) result(change)
        real(real64), intent(in) :: jsd(:, :)
            !! `(n_studies, :)`: observed global JSD per study and admissible candidate
        integer(int32), intent(in) :: i_transition
            !! The transition, from candidate `i_transition` to `i_transition + 1`
        real(real64) :: change
            !! The largest absolute change
        integer(int32) :: i_study

        change = 0.0_real64
        do i_study = 1, size(jsd, 1, kind=int32)
            change = max(change, abs(jsd(i_study, i_transition + 1_int32) - jsd(i_study, i_transition)))
        end do
    end function largest_change

    !> Test helper for test_rerun_forwards_divergence_mode: from the observed global JSDs of the
    !| same admissible candidates under the symmetric and the pairwise divergence (two studies),
    !| the absolute tolerance that the symmetric change at transition 1 (candidate 1 to 2) and the
    !| smallest pairwise change at a later transition stay below, and the pairwise change at
    !| transition 1 exceeds. A transition's change is the largest absolute JSD change across the
    !| studies. With `below` the larger of the first two changes and `above` the pairwise change at
    !| transition 1, the tolerance is their geometric mean and `margin = above/below`, so each lies a
    !| factor `sqrt(margin)` from the tolerance; `margin` is 0 with fewer than 3 admissible
    !| candidates.
    subroutine divergence_stop_tolerance(n_admissible, jsd_symmetric, jsd_pairwise, tolerance, margin)
        integer(int32), intent(in) :: n_admissible
            !! Number of admissible candidates (the same in both modes)
        real(real64), intent(in) :: jsd_symmetric(:, :)
            !! `(n_studies, :)`: observed global JSD per study and admissible candidate, symmetric
        real(real64), intent(in) :: jsd_pairwise(:, :)
            !! `(n_studies, :)`: the same with the pairwise divergence
        real(real64), intent(out) :: tolerance
            !! The tolerance described above
        real(real64), intent(out) :: margin
            !! Pairwise change at transition 1 over the larger of the other two changes
        real(real64) :: below, above, smallest_later
        integer(int32) :: i_transition

        tolerance = 0.0_real64
        margin = 0.0_real64
        if (n_admissible < 3_int32) return
        smallest_later = huge(1.0_real64)
        do i_transition = 2, n_admissible - 1_int32
            smallest_later = min(smallest_later, largest_change(jsd_pairwise, i_transition))
        end do
        below = max(largest_change(jsd_symmetric, 1_int32), smallest_later)
        above = largest_change(jsd_pairwise, 1_int32)
        if (below <= 0.0_real64) return
        tolerance = sqrt(below*above)
        margin = above/below
    end subroutine divergence_stop_tolerance

    !> Test helper: the checks of test_two_studies_default_divergence_null_and_batch on one label test.
    subroutine assert_two_study_label_test(name, observed, reference, p, null, n_failed, p_min, batch)
        character(*), intent(in) :: name
            !! Message prefix
        real(real64), intent(in) :: observed(2)
            !! The test's observed global JSDs
        real(real64), intent(in) :: reference(2)
            !! The plain pipeline's pairwise global JSDs
        real(real64), intent(in) :: p(2)
            !! The test's p-values
        real(real64), intent(in) :: null(:, :)
            !! The test's null rows, one column per study
        integer(int32), intent(in) :: n_failed
            !! The test's failed permutations
        real(real64), intent(in) :: p_min
            !! `1/(n+1)`
        logical(c_bool), intent(in) :: batch
            !! `.true.`: a batch between the two studies, both p-values must be the minimum
        integer(int32) :: i_study

        call assert_equal_int(n_failed, 0_int32, name//"no failed permutation")
        call assert_true(equal_and_not_nan(observed, reference), name//"observed is the pairwise value")
        call assert_true(observed(1) == observed(2), name//"observed: one value shared by both studies")
        call assert_true(equal_and_not_nan(null(:, 1), null(:, 2)), name//"null rows: one value shared by both studies")
        call assert_true(p(1) == p(2), name//"the two p-values are equal")
        do i_study = 1, 2
            call assert_equal_real(p(i_study), real(count(null(:, i_study) >= observed(i_study)) + 1, real64)/ &
                                   real(size(null, 1) + 1, real64), 1.0e-15_real64, name//"p recounted from the null rows")
        end do
        ! Exact comparison: p and p_min are both (count + 1)/(n + 1) quotients computed alike.
        if (batch) call assert_true(p(1) == p_min .and. p(2) == p_min, name//"batch: both p-values at the minimum")
    end subroutine assert_two_study_label_test

    !> With the rerun option, all four routines forward `delta_abs_tolerance` to the parameter
    !| search. Effect-size thresholds of 1e-300 make every non-zero change disqualify (the adaptive
    !| part also sets `min_neighbor_overlap = 0`, so that several candidates are admissible), so a
    !| huge tolerance (1, every JSD lies in [0, 1]) is what lets a search find a plateau.
    !| Adaptive: data seed ADAPTIVE_TOLERANCE_SEED (52, a pinned seed: the first of 52..71 on which
    !| the adaptive search called directly selects different knobs with tolerance 1 than with 0,
    !| asserted here); for both tolerances assignment run with rerun reports exactly the knobs and values
    !| of search + construction + run at the same tolerance, and whole test the same observed values.
    !| k-NN: the multi-candidate fixture of make_knn_grid_fixture, where the search selects another
    !| `(n_points, n_neighbors)` with tolerance 1 (a plateau at candidate 3) than with 0 (no
    !| plateau, the smallest-uncertainty candidate), asserted; for both tolerances assignment run and
    !| (b) with rerun equal the k-NN search + run at that tolerance.
    subroutine test_rerun_forwards_delta_abs_tolerance()
        character(*), parameter :: name = "test_rerun_forwards_delta_abs_tolerance: "
        integer(int32), parameter :: ADAPTIVE_TOLERANCE_SEED = 52_int32
        real(real64), parameter :: tiny_threshold = 1.0e-300_real64, tolerances(2) = [0.0_real64, 1.0_real64]
        real(real64), allocatable :: expression(:, :, :), gene_means(:, :), residuals(:, :, :)
        real(real64) :: reference(3), result(3), observed(3), p(3), null(1, 3)
        integer(int32) :: ierr, used_ks, used_kst, used_km, n_points, n_points_reference, n_failed, i_tol
        integer(int32) :: ks(2), kst(2), km(2), search_np(2), search_nn(2), used_np, used_nn, n_samples

        call make_expression(300_int32, [5_int32, 4_int32, 4_int32], [0.5_real64, 0.5_real64, 0.6_real64], &
                             ADAPTIVE_TOLERANCE_SEED, expression)
        call reference_means_residuals(expression, gene_means, residuals)
        do i_tol = 1, 2
            call reference_adaptive_search(gene_means, residuals, 20_int32, ks(i_tol), kst(i_tol), km(i_tol), &
                                           delta_abs_tolerance=tolerances(i_tol), &
                                           delta_median_threshold=tiny_threshold, delta_max_threshold=tiny_threshold, &
                                           min_neighbor_overlap=0.0_real64)
        end do
        call assert_true(ks(1) /= ks(2) .or. kst(1) /= kst(2) .or. km(1) /= km(2), &
                         name//"the tolerance changes the adaptive search's selection on this data set")
        do i_tol = 1, 2
            call reference_adaptive(expression, ks(i_tol), kst(i_tol), km(i_tol), reference, n_points_reference, ierr)
            call assert_equal_int(get_err_code(ierr), ERR_OK, name//"adaptive reference ierr")
            call run_js_comp_test_adaptive_on_sample_assignment(3_int32, 300_int32, 5_int32, expression, 13_int32, &
                                                                identity(13_int32), K_START, K_STEP, K_MAX, result, &
                                                                used_ks, used_kst, used_km, n_points, &
                                                                rerun_parameter_search=.true._c_bool, &
                                                                n_bootstraps=20_int32, &
                                                                delta_median_threshold=tiny_threshold, &
                                                                delta_max_threshold=tiny_threshold, &
                                                                delta_abs_tolerance=tolerances(i_tol), &
                                                                min_neighbor_overlap=0.0_real64, ierr=ierr)
            call assert_equal_int(get_err_code(ierr), ERR_OK, name//"adaptive assignment run ierr")
            call assert_equal_int(used_ks, ks(i_tol), name//"adaptive: used_k_start is the search's")
            call assert_equal_int(used_kst, kst(i_tol), name//"adaptive: used_k_step is the search's")
            call assert_equal_int(used_km, km(i_tol), name//"adaptive: used_k_max is the search's")
            call assert_true(equal_and_not_nan(result, reference), name//"adaptive: equals search + run")
            call run_js_comp_test_adaptive_permutation(3_int32, 300_int32, 5_int32, expression, K_START, &
                                                       K_STEP, K_MAX, 1_int32, observed, p, null, n_failed, &
                                                       rerun_parameter_search=.true._c_bool, &
                                                       n_bootstraps=20_int32, &
                                                       delta_median_threshold=tiny_threshold, &
                                                       delta_max_threshold=tiny_threshold, &
                                                       delta_abs_tolerance=tolerances(i_tol), &
                                                       min_neighbor_overlap=0.0_real64, ierr=ierr)
            call assert_equal_int(get_err_code(ierr), ERR_OK, name//"adaptive whole test ierr")
            call assert_true(equal_and_not_nan(observed, result), name//"adaptive whole test forwards the tolerance")
        end do

        call make_knn_grid_fixture(expression)
        n_samples = count(.not. ieee_is_nan(expression(:, 1, :)), kind=int32)
        call reference_means_residuals(expression, gene_means, residuals)
        do i_tol = 1, 2
            call reference_fixed_k_search(gene_means, residuals, 20_int32, search_np(i_tol), search_nn(i_tol), &
                                          delta_abs_tolerance=tolerances(i_tol), &
                                          delta_median_threshold=tiny_threshold, delta_max_threshold=tiny_threshold)
        end do
        call assert_true(search_np(1) /= search_np(2) .or. search_nn(1) /= search_nn(2), &
                         name//"the tolerance changes the k-NN search's selection on this data set")
        do i_tol = 1, 2
            call reference_knn(expression, search_np(i_tol), search_nn(i_tol), reference, ierr)
            call assert_equal_int(get_err_code(ierr), ERR_OK, name//"k-NN reference ierr")
            call run_js_comp_test_on_sample_assignment(3_int32, size(expression, 2, kind=int32), &
                                                       size(expression, 1, kind=int32), expression, n_samples, &
                                                       identity(n_samples), N_POINTS_KNN, N_NEIGHBORS_KNN, result, &
                                                       used_np, used_nn, rerun_parameter_search=.true._c_bool, &
                                                       n_bootstraps=20_int32, delta_median_threshold=tiny_threshold, &
                                                       delta_max_threshold=tiny_threshold, &
                                                       delta_abs_tolerance=tolerances(i_tol), ierr=ierr)
            call assert_equal_int(get_err_code(ierr), ERR_OK, name//"k-NN assignment run ierr")
            call assert_equal_array_int([used_np, used_nn], [search_np(i_tol), search_nn(i_tol)], 2_int32, &
                                        name//"k-NN: used pair is the search's")
            call assert_true(equal_and_not_nan(result, reference), name//"k-NN: equals search + run")
            call run_js_comp_test_permutation(3_int32, size(expression, 2, kind=int32), size(expression, 1, kind=int32), &
                                              expression, N_POINTS_KNN, N_NEIGHBORS_KNN, 1_int32, observed, p, null, &
                                              n_failed, rerun_parameter_search=.true._c_bool, n_bootstraps=20_int32, &
                                              delta_median_threshold=tiny_threshold, &
                                              delta_max_threshold=tiny_threshold, &
                                              delta_abs_tolerance=tolerances(i_tol), ierr=ierr)
            call assert_equal_int(get_err_code(ierr), ERR_OK, name//"k-NN whole test ierr")
            call assert_true(equal_and_not_nan(observed, result), name//"k-NN whole test forwards the tolerance")
        end do
    end subroutine test_rerun_forwards_delta_abs_tolerance

    !> All four routines forward `divergence_mode` to the run (without the rerun option; the rerun
    !| search is test_rerun_forwards_divergence_mode). Two studies (300 genes, 5 and 4 samples, other
    !| noise levels): the symmetric divergence gives each study its JSD against the consensus (two
    !| different values), the pairwise one (the default) the one JSD of study 1 against study 2,
    !| shared by both. For each mode, assignment run on the identity equals the plain run (k-NN:
    !| run_js_comp_test; adaptive: construction + run_js_comp_test_adaptive) with the same mode,
    !| and whole test reports the same observed values; the two modes' values differ (so the test
    !| fails if the mode is not forwarded), and the default equals the pairwise mode.
    subroutine test_forwards_divergence_mode()
        character(*), parameter :: name = "test_forwards_divergence_mode: "
        integer(int32), parameter :: modes(2) = [MODE_DIVERGENCE_SYMMETRIC, MODE_DIVERGENCE_PAIRWISE_K2]
        integer(int32), parameter :: N_SAMPLES_K2 = 9_int32
        real(real64), allocatable :: expression(:, :, :)
        real(real64) :: knn_values(2, 2), adaptive_values(2, 2), reference(2), observed(2), p(2), null(1, 2)
        real(real64) :: default_values(2)
        integer(int32) :: ierr, i_mode, used_np, used_nn, used_ks, used_kst, used_km, n_points, n_failed

        call make_expression(300_int32, [5_int32, 4_int32], [0.5_real64, 0.8_real64], 91_int32, expression)
        do i_mode = 1, 2
            call run_js_comp_test_on_sample_assignment(2_int32, 300_int32, 5_int32, expression, N_SAMPLES_K2, &
                                                       identity(N_SAMPLES_K2), N_POINTS_KNN, N_NEIGHBORS_KNN, &
                                                       knn_values(:, i_mode), used_np, used_nn, &
                                                       divergence_mode=modes(i_mode), ierr=ierr)
            call assert_equal_int(get_err_code(ierr), ERR_OK, name//"k-NN assignment run ierr")
            call reference_knn(expression, N_POINTS_KNN, N_NEIGHBORS_KNN, reference, ierr, divergence_mode=modes(i_mode))
            call assert_equal_int(get_err_code(ierr), ERR_OK, name//"k-NN reference ierr")
            call assert_true(equal_and_not_nan(knn_values(:, i_mode), reference), name//"k-NN assignment run equals the run")
            call run_js_comp_test_permutation(2_int32, 300_int32, 5_int32, expression, N_POINTS_KNN, N_NEIGHBORS_KNN, &
                                              1_int32, observed, p, null, n_failed, divergence_mode=modes(i_mode), &
                                              ierr=ierr)
            call assert_equal_int(get_err_code(ierr), ERR_OK, name//"k-NN whole test ierr")
            call assert_true(equal_and_not_nan(observed, knn_values(:, i_mode)), name//"k-NN whole test forwards the mode")

            call run_js_comp_test_adaptive_on_sample_assignment(2_int32, 300_int32, 5_int32, expression, N_SAMPLES_K2, &
                                                                identity(N_SAMPLES_K2), K_START, K_STEP, K_MAX, &
                                                                adaptive_values(:, i_mode), used_ks, used_kst, used_km, &
                                                                n_points, divergence_mode=modes(i_mode), ierr=ierr)
            call assert_equal_int(get_err_code(ierr), ERR_OK, name//"adaptive assignment run ierr")
            call reference_adaptive(expression, K_START, K_STEP, K_MAX, reference, n_points, ierr, &
                                    divergence_mode=modes(i_mode))
            call assert_equal_int(get_err_code(ierr), ERR_OK, name//"adaptive reference ierr")
            call assert_true(equal_and_not_nan(adaptive_values(:, i_mode), reference), &
                             name//"adaptive assignment run equals construction + run")
            call run_js_comp_test_adaptive_permutation(2_int32, 300_int32, 5_int32, expression, K_START, K_STEP, K_MAX, &
                                                       1_int32, observed, p, null, n_failed, divergence_mode=modes(i_mode), &
                                                       ierr=ierr)
            call assert_equal_int(get_err_code(ierr), ERR_OK, name//"adaptive whole test ierr")
            call assert_true(equal_and_not_nan(observed, adaptive_values(:, i_mode)), &
                             name//"adaptive whole test forwards the mode")
        end do
        call assert_true(knn_values(1, 1) /= knn_values(2, 1), name//"k-NN symmetric: the two studies differ")
        call assert_true(knn_values(1, 2) == knn_values(2, 2), name//"k-NN pairwise: one shared value")
        call assert_true(knn_values(1, 1) /= knn_values(1, 2), name//"k-NN: the modes give different values")
        call assert_true(adaptive_values(1, 1) /= adaptive_values(2, 1), name//"adaptive symmetric: the two studies differ")
        call assert_true(adaptive_values(1, 2) == adaptive_values(2, 2), name//"adaptive pairwise: one shared value")
        call assert_true(adaptive_values(1, 1) /= adaptive_values(1, 2), name//"adaptive: the modes give different values")

        call run_js_comp_test_on_sample_assignment(2_int32, 300_int32, 5_int32, expression, N_SAMPLES_K2, &
                                                   identity(N_SAMPLES_K2), N_POINTS_KNN, N_NEIGHBORS_KNN, &
                                                   default_values, used_np, used_nn, ierr=ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"k-NN default ierr")
        call assert_true(equal_and_not_nan(default_values, knn_values(:, 2)), name//"k-NN: the default is pairwise")
        call run_js_comp_test_adaptive_on_sample_assignment(2_int32, 300_int32, 5_int32, expression, N_SAMPLES_K2, &
                                                            identity(N_SAMPLES_K2), K_START, K_STEP, K_MAX, &
                                                            default_values, used_ks, used_kst, used_km, n_points, &
                                                            ierr=ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"adaptive default ierr")
        call assert_true(equal_and_not_nan(default_values, adaptive_values(:, 2)), name//"adaptive: the default is pairwise")
    end subroutine test_forwards_divergence_mode

    !> With the rerun option, all four routines forward `divergence_mode` to the rerun parameter
    !| search, not only to the run after it. Two studies (1000 genes, 5 and 4 samples, noise SD 0.5
    !| and 0.8), since with three or more studies both modes coincide. k-NN: data seed
    !| KNN_DIVERGENCE_SEED and `grid_points_factor = 10`, whose grid is (317, 12), (254, 15), (203,
    !| 19), (200, 20); adaptive: data seed ADAPTIVE_DIVERGENCE_SEED, the default knobs and
    !| `min_neighbor_overlap = 0`, so that at least three candidates are admissible.
    !|
    !| The two modes are made to select different candidates by a rule that reads only the observed
    !| JSDs, never the bootstrap: relative thresholds 1e-300 and a single qualifying transition
    !| (`delta_min_consecutive_transitions = 1`), so a search stops at the first transition between
    !| admissible candidates where every study's JSD changes by at most `delta_abs_tolerance`, and
    !| selects the candidate after it. With two studies the symmetric JSD (each study against the
    !| consensus) is roughly a quarter of the pairwise one, and so are its changes. From the observed
    !| JSDs of searches that cannot stop (tolerance 0), divergence_stop_tolerance picks a tolerance
    !| that the symmetric change at transition 1 stays below and the pairwise one exceeds, while a
    !| later pairwise transition stays below it. The test requires these changes to lie at least a
    !| factor DIVERGENCE_STOP_MARGIN apart (each at least its square root from the tolerance), so
    !| rounding differences between compilers or optimisation levels cannot flip a decision; the
    !| seeds are pinned only because they give such a margin (k-NN 3.2, adaptive 3.6 with gfortran).
    !| The symmetric search then selects admissible candidate 2 and the pairwise one a later
    !| candidate (asserted on the direct searches). For each mode, assignment run with rerun reports
    !| exactly the knobs the search selects directly with that mode, and its values equal the run at
    !| them with that mode; whole test with rerun reports the same observed values. If the forwarding
    !| into the rerun search were dropped, the symmetric case would search with the default
    !| (pairwise) divergence, select the later candidate and fail.
    subroutine test_rerun_forwards_divergence_mode()
        character(*), parameter :: name = "test_rerun_forwards_divergence_mode: "
        integer(int32), parameter :: modes(2) = [MODE_DIVERGENCE_SYMMETRIC, MODE_DIVERGENCE_PAIRWISE_K2]
        integer(int32), parameter :: KNN_DIVERGENCE_SEED = 8_int32, ADAPTIVE_DIVERGENCE_SEED = 2_int32
        integer(int32), parameter :: N_GENES_K2 = 1000_int32, N_SAMPLES_K2 = 9_int32, ONE_TRANSITION = 1_int32
        real(real64), parameter :: DIVERGENCE_GRID_FACTOR = 10.0_real64, tiny_threshold = 1.0e-300_real64
        real(real64), parameter :: DIVERGENCE_STOP_MARGIN = 2.0_real64
        real(real64), allocatable :: expression(:, :, :), gene_means(:, :), residuals(:, :, :)
        real(real64) :: reference(2), result(2), observed(2), p(2), null(1, 2), tolerance, margin
        real(real64) :: trace_jsd(2, MAX_CANDIDATE_PAIRS, 2)
        integer(int32) :: ierr, i_mode, used_np, used_nn, used_ks, used_kst, used_km, n_points, n_points_reference
        integer(int32) :: n_failed, search_np(2), search_nn(2), ks(2), kst(2), km(2), n_admissible(2)

        call make_expression(N_GENES_K2, [5_int32, 4_int32], [0.5_real64, 0.8_real64], KNN_DIVERGENCE_SEED, expression)
        call reference_means_residuals(expression, gene_means, residuals)
        do i_mode = 1, 2
            call reference_fixed_k_search(gene_means, residuals, 20_int32, search_np(i_mode), search_nn(i_mode), &
                                          delta_abs_tolerance=0.0_real64, delta_median_threshold=tiny_threshold, &
                                          delta_max_threshold=tiny_threshold, divergence_mode=modes(i_mode), &
                                          grid_points_factor=DIVERGENCE_GRID_FACTOR, &
                                          delta_min_consecutive_transitions=ONE_TRANSITION, &
                                          n_admissible_out=n_admissible(i_mode), trace_jsd_out=trace_jsd(:, :, i_mode))
        end do
        call assert_equal_int(n_admissible(1), n_admissible(2), name//"k-NN: the same admissible candidates")
        call divergence_stop_tolerance(n_admissible(1), trace_jsd(:, :, 1), trace_jsd(:, :, 2), tolerance, margin)
        call assert_true(margin >= DIVERGENCE_STOP_MARGIN, name//"k-NN: the stop decisions have a wide margin")
        do i_mode = 1, 2
            call reference_fixed_k_search(gene_means, residuals, 20_int32, search_np(i_mode), search_nn(i_mode), &
                                          delta_abs_tolerance=tolerance, delta_median_threshold=tiny_threshold, &
                                          delta_max_threshold=tiny_threshold, divergence_mode=modes(i_mode), &
                                          grid_points_factor=DIVERGENCE_GRID_FACTOR, &
                                          delta_min_consecutive_transitions=ONE_TRANSITION)
        end do
        call assert_true(search_np(1) /= search_np(2) .or. search_nn(1) /= search_nn(2), &
                         name//"the mode changes the k-NN search's selection")
        do i_mode = 1, 2
            call run_js_comp_test_on_sample_assignment(2_int32, N_GENES_K2, 5_int32, expression, N_SAMPLES_K2, &
                                                       identity(N_SAMPLES_K2), N_POINTS_KNN, N_NEIGHBORS_KNN, result, &
                                                       used_np, used_nn, rerun_parameter_search=.true._c_bool, &
                                                       n_bootstraps=20_int32, delta_median_threshold=tiny_threshold, &
                                                       delta_max_threshold=tiny_threshold, &
                                                       delta_min_consecutive_transitions=ONE_TRANSITION, &
                                                       delta_abs_tolerance=tolerance, divergence_mode=modes(i_mode), &
                                                       grid_points_factor=DIVERGENCE_GRID_FACTOR, ierr=ierr)
            call assert_equal_int(get_err_code(ierr), ERR_OK, name//"k-NN assignment run ierr")
            call assert_equal_array_int([used_np, used_nn], [search_np(i_mode), search_nn(i_mode)], 2_int32, &
                                        name//"k-NN: used pair is the search's with the same mode")
            call reference_knn(expression, search_np(i_mode), search_nn(i_mode), reference, ierr, &
                               divergence_mode=modes(i_mode))
            call assert_equal_int(get_err_code(ierr), ERR_OK, name//"k-NN reference ierr")
            call assert_true(equal_and_not_nan(result, reference), name//"k-NN: equals search + run with the same mode")
            call run_js_comp_test_permutation(2_int32, N_GENES_K2, 5_int32, expression, N_POINTS_KNN, N_NEIGHBORS_KNN, &
                                              1_int32, observed, p, null, n_failed, rerun_parameter_search=.true._c_bool, &
                                              n_bootstraps=20_int32, delta_median_threshold=tiny_threshold, &
                                              delta_max_threshold=tiny_threshold, &
                                              delta_min_consecutive_transitions=ONE_TRANSITION, &
                                              delta_abs_tolerance=tolerance, divergence_mode=modes(i_mode), &
                                              grid_points_factor=DIVERGENCE_GRID_FACTOR, ierr=ierr)
            call assert_equal_int(get_err_code(ierr), ERR_OK, name//"k-NN whole test ierr")
            call assert_true(equal_and_not_nan(observed, result), name//"k-NN whole test forwards the mode")
        end do

        call make_expression(N_GENES_K2, [5_int32, 4_int32], [0.5_real64, 0.8_real64], ADAPTIVE_DIVERGENCE_SEED, expression)
        call reference_means_residuals(expression, gene_means, residuals)
        do i_mode = 1, 2
            call reference_adaptive_search(gene_means, residuals, 20_int32, ks(i_mode), kst(i_mode), km(i_mode), &
                                           delta_abs_tolerance=0.0_real64, delta_median_threshold=tiny_threshold, &
                                           delta_max_threshold=tiny_threshold, min_neighbor_overlap=0.0_real64, &
                                           divergence_mode=modes(i_mode), &
                                           delta_min_consecutive_transitions=ONE_TRANSITION, &
                                           n_admissible_out=n_admissible(i_mode), trace_jsd_out=trace_jsd(:, :, i_mode))
        end do
        call assert_equal_int(n_admissible(1), n_admissible(2), name//"adaptive: the same admissible candidates")
        call divergence_stop_tolerance(n_admissible(1), trace_jsd(:, :, 1), trace_jsd(:, :, 2), tolerance, margin)
        call assert_true(margin >= DIVERGENCE_STOP_MARGIN, name//"adaptive: the stop decisions have a wide margin")
        do i_mode = 1, 2
            call reference_adaptive_search(gene_means, residuals, 20_int32, ks(i_mode), kst(i_mode), km(i_mode), &
                                           delta_abs_tolerance=tolerance, delta_median_threshold=tiny_threshold, &
                                           delta_max_threshold=tiny_threshold, min_neighbor_overlap=0.0_real64, &
                                           divergence_mode=modes(i_mode), &
                                           delta_min_consecutive_transitions=ONE_TRANSITION)
        end do
        call assert_true(ks(1) /= ks(2) .or. kst(1) /= kst(2) .or. km(1) /= km(2), &
                         name//"the mode changes the adaptive search's selection")
        do i_mode = 1, 2
            call run_js_comp_test_adaptive_on_sample_assignment(2_int32, N_GENES_K2, 5_int32, expression, N_SAMPLES_K2, &
                                                                identity(N_SAMPLES_K2), K_START, K_STEP, K_MAX, result, &
                                                                used_ks, used_kst, used_km, n_points, &
                                                                rerun_parameter_search=.true._c_bool, &
                                                                n_bootstraps=20_int32, &
                                                                delta_median_threshold=tiny_threshold, &
                                                                delta_max_threshold=tiny_threshold, &
                                                                delta_min_consecutive_transitions=ONE_TRANSITION, &
                                                                delta_abs_tolerance=tolerance, &
                                                                min_neighbor_overlap=0.0_real64, &
                                                                divergence_mode=modes(i_mode), ierr=ierr)
            call assert_equal_int(get_err_code(ierr), ERR_OK, name//"adaptive assignment run ierr")
            call assert_equal_array_int([used_ks, used_kst, used_km], [ks(i_mode), kst(i_mode), km(i_mode)], 3_int32, &
                                        name//"adaptive: used knobs are the search's with the same mode")
            call reference_adaptive(expression, ks(i_mode), kst(i_mode), km(i_mode), reference, n_points_reference, ierr, &
                                    divergence_mode=modes(i_mode))
            call assert_equal_int(get_err_code(ierr), ERR_OK, name//"adaptive reference ierr")
            call assert_equal_int(n_points, n_points_reference, name//"adaptive: n_points")
            call assert_true(equal_and_not_nan(result, reference), name//"adaptive: equals search + run with the same mode")
            call run_js_comp_test_adaptive_permutation(2_int32, N_GENES_K2, 5_int32, expression, K_START, K_STEP, K_MAX, &
                                                       1_int32, observed, p, null, n_failed, &
                                                       rerun_parameter_search=.true._c_bool, n_bootstraps=20_int32, &
                                                       delta_median_threshold=tiny_threshold, &
                                                       delta_max_threshold=tiny_threshold, &
                                                       delta_min_consecutive_transitions=ONE_TRANSITION, &
                                                       delta_abs_tolerance=tolerance, min_neighbor_overlap=0.0_real64, &
                                                       divergence_mode=modes(i_mode), ierr=ierr)
            call assert_equal_int(get_err_code(ierr), ERR_OK, name//"adaptive whole test ierr")
            call assert_true(equal_and_not_nan(observed, result), name//"adaptive whole test forwards the mode")
        end do
    end subroutine test_rerun_forwards_divergence_mode

    !> The sizing routines, by hand for G = 300 (S = 3). Default grid: `ceiling(4*sqrt(300)) = 70`
    !| clamps up to the floor 200, so the fixed-k search's bounds are 200 points and
    !| `floor(300/(0.25*200)) = 6` neighbors. With `grid_points_factor = 17.32`, `grid_min_points =
    !| 240`: `ceiling(17.32*sqrt(300)) = ceiling(299.99) = 300` points and `floor(300/(0.25*240)) = 5`
    !| neighbors; with `grid_min_points = 1500` (MAX_POINTS): 1500 points (clamped up) and
    !| `max(1, floor(300/(0.25*1500))) = max(1, 0) = 1` neighbor. The adaptive search's capacity is
    !| `ceiling(2*900/11) + 1 = 165` (smallest `k_start` 11).
    subroutine test_permutation_bounds_by_hand()
        character(*), parameter :: name = "test_permutation_bounds_by_hand: "
        real(real64), parameter :: LARGER_POINTS_FACTOR = 17.32_real64
        integer(int32), parameter :: LARGER_MIN_POINTS = 240_int32, MAX_POINTS_FLOOR = 1500_int32
        integer(int32) :: max_p, max_k, top_k, ierr

        call calc_js_comp_test_permutation_bounds(300_int32, 10_int32, 30_int32, max_n_points_candidate=max_p, &
                                                  max_n_neighbors_candidate=max_k, n_bootstrapping_top_k_jsds=top_k, &
                                                  ierr=ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"knn ierr")
        call assert_equal_int(max_p, 200_int32, name//"max(10, 200)")
        call assert_equal_int(max_k, 30_int32, name//"max(30, 6)")
        call assert_equal_int(top_k, 25_int32, name//"default heap: floor(2.5% of 1000)")
        call calc_js_comp_test_permutation_bounds(300_int32, 400_int32, 2_int32, 100_int32, 5.0_real64, &
                                                  max_n_points_candidate=max_p, max_n_neighbors_candidate=max_k, &
                                                  n_bootstrapping_top_k_jsds=top_k, ierr=ierr)
        call assert_equal_int(max_p, 400_int32, name//"max(400, 200)")
        call assert_equal_int(max_k, 6_int32, name//"max(2, 6)")
        call assert_equal_int(top_k, 5_int32, name//"floor(5% of 100)")
        call calc_js_comp_test_permutation_bounds(300_int32, 10_int32, 2_int32, grid_points_factor=LARGER_POINTS_FACTOR, &
                                                  grid_min_points=LARGER_MIN_POINTS, max_n_points_candidate=max_p, &
                                                  max_n_neighbors_candidate=max_k, n_bootstrapping_top_k_jsds=top_k, &
                                                  ierr=ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"grid options ierr")
        call assert_equal_int(max_p, 300_int32, name//"grid options: max(10, 300)")
        call assert_equal_int(max_k, 5_int32, name//"grid options: max(2, 5)")
        call calc_js_comp_test_permutation_bounds(300_int32, 10_int32, 2_int32, grid_min_points=MAX_POINTS_FLOOR, &
                                                  max_n_points_candidate=max_p, max_n_neighbors_candidate=max_k, &
                                                  n_bootstrapping_top_k_jsds=top_k, ierr=ierr)
        call assert_equal_int(max_p, 1500_int32, name//"floor 1500: max(10, 1500)")
        call assert_equal_int(max_k, 2_int32, name//"floor 1500: max(2, 1)")
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
        call calc_js_comp_test_permutation_bounds(300_int32, 10_int32, 30_int32, 100_int32, 101.0_real64, &
                                                  max_n_points_candidate=max_p, max_n_neighbors_candidate=max_k, &
                                                  n_bootstrapping_top_k_jsds=top_k, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, name//"level 101", arg_pos=5_int32)
        call calc_js_comp_test_permutation_bounds(300_int32, 10_int32, 30_int32, grid_points_factor=0.0_real64, &
                                                  max_n_points_candidate=max_p, max_n_neighbors_candidate=max_k, &
                                                  n_bootstrapping_top_k_jsds=top_k, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, name//"grid_points_factor 0", arg_pos=6_int32)
        call calc_js_comp_test_permutation_bounds(300_int32, 10_int32, 30_int32, grid_min_points=0_int32, &
                                                  max_n_points_candidate=max_p, max_n_neighbors_candidate=max_k, &
                                                  n_bootstrapping_top_k_jsds=top_k, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, name//"grid_min_points 0", arg_pos=7_int32)
        call calc_js_comp_test_permutation_bounds(300_int32, 10_int32, 30_int32, grid_min_points=MAX_POINTS_FLOOR + 1_int32, &
                                                  max_n_points_candidate=max_p, max_n_neighbors_candidate=max_k, &
                                                  n_bootstrapping_top_k_jsds=top_k, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, name//"grid_min_points 1501", arg_pos=7_int32)

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

    !> Two studies with the default divergence (pairwise, one JSD of study 1 against study 2 shared
    !| by both), k-NN and adaptive, 9 permutations. On null data (both studies noise SD 0.5): the
    !| observed values are the pairwise value of the plain pipeline run with
    !| `MODE_DIVERGENCE_PAIRWISE_K2` passed explicitly, equal for both studies; every null row holds
    !| one shared value (`jsd_null(:, 1) == jsd_null(:, 2)` bit for bit), so the two p-values are equal,
    !| and they lie in `[1/(n+1), 1]`. With a batch (study 2 noise SD 2.0) both p-values are the
    !| minimum `1/(n+1)`.
    subroutine test_two_studies_default_divergence_null_and_batch()
        character(*), parameter :: name = "test_two_studies_default_divergence_null_and_batch: "
        integer(int32), parameter :: n = 9_int32
        real(real64), allocatable :: expression(:, :, :)
        real(real64) :: observed(2), p(2), null(n, 2), reference(2), p_min
        integer(int32) :: ierr, n_failed, n_points, i_case
        character(len=8), parameter :: case_names(2) = ["null    ", "batch   "]
        real(real64), parameter :: noise(2, 2) = reshape([0.5_real64, 0.5_real64, 0.5_real64, 2.0_real64], [2, 2])

        p_min = 1.0_real64/real(n + 1, real64)
        do i_case = 1, 2
            call make_expression(300_int32, [5_int32, 5_int32], noise(:, i_case), 300_int32 + i_case, expression)

            call run_js_comp_test_permutation(2_int32, 300_int32, 5_int32, expression, N_POINTS_KNN, N_NEIGHBORS_KNN, n, &
                                              observed, p, null, n_failed, random_seed=i_case, ierr=ierr)
            call assert_equal_int(get_err_code(ierr), ERR_OK, name//"k-NN "//trim(case_names(i_case))//": ierr")
            call reference_knn(expression, N_POINTS_KNN, N_NEIGHBORS_KNN, reference, ierr, &
                               divergence_mode=MODE_DIVERGENCE_PAIRWISE_K2)
            call assert_equal_int(get_err_code(ierr), ERR_OK, name//"k-NN reference ierr")
            call assert_two_study_label_test(name//"k-NN "//trim(case_names(i_case))//": ", observed, reference, p, null, &
                                             n_failed, p_min, logical(i_case == 2, c_bool))

            call run_js_comp_test_adaptive_permutation(2_int32, 300_int32, 5_int32, expression, K_START, K_STEP, K_MAX, n, &
                                                       observed, p, null, n_failed, random_seed=i_case, ierr=ierr)
            call assert_equal_int(get_err_code(ierr), ERR_OK, name//"adaptive "//trim(case_names(i_case))//": ierr")
            call reference_adaptive(expression, K_START, K_STEP, K_MAX, reference, n_points, ierr, &
                                    divergence_mode=MODE_DIVERGENCE_PAIRWISE_K2)
            call assert_equal_int(get_err_code(ierr), ERR_OK, name//"adaptive reference ierr")
            call assert_two_study_label_test(name//"adaptive "//trim(case_names(i_case))//": ", observed, reference, p, &
                                             null, n_failed, p_min, logical(i_case == 2, c_bool))
        end do
    end subroutine test_two_studies_default_divergence_null_and_batch

    !> The p-value routine by hand: `(count + 1)/(n_valid + 1)` with `count` the valid null values that are at least
    !| the observed one (a tie counts), NaN rows left out, and NaN when no row is valid.
    subroutine test_label_permutation_p_values_by_hand()
        character(*), parameter :: name = "test_label_permutation_p_values_by_hand: "
        real(real64) :: observed(3), null(5, 3), p(3), nan
        real(real64) :: null_low(4, 1), null_high(4, 1), observed_one(1), p_one(1)

        nan = ieee_value(1.0_real64, ieee_quiet_nan)
        observed = 0.5_real64
        ! Study 1: 0.5 (a tie), 0.6, 0.7 reach the observed value, 0.4 and 0.3 do not: (3 + 1)/(5 + 1).
        null(:, 1) = [0.5_real64, 0.6_real64, 0.4_real64, 0.3_real64, 0.7_real64]
        ! Study 2: two NaN rows are left out; 0.6 and 0.5 reach it, 0.1 does not: (2 + 1)/(3 + 1).
        null(:, 2) = [nan, 0.6_real64, nan, 0.1_real64, 0.5_real64]
        ! Study 3: no valid row.
        null(:, 3) = nan
        call label_permutation_p_values(5_int32, 3_int32, observed, null, p)
        call assert_equal_real(p(1), 4.0_real64/6.0_real64, 1.0e-15_real64, name//"study 1: a tie counts")
        call assert_equal_real(p(2), 3.0_real64/4.0_real64, 1.0e-15_real64, name//"study 2: NaN rows left out")
        call assert_true(ieee_is_nan(p(3)), name//"study 3: no valid row gives NaN")

        ! The extremes: no null value reaches the observed one -> 1/(n + 1); all do -> 1.
        observed_one = 0.5_real64
        null_low(:, 1) = [0.1_real64, 0.2_real64, 0.3_real64, 0.4_real64]
        null_high(:, 1) = [0.6_real64, 0.7_real64, 0.8_real64, 0.9_real64]
        call label_permutation_p_values(4_int32, 1_int32, observed_one, null_low, p_one)
        call assert_equal_real(p_one(1), 1.0_real64/5.0_real64, 1.0e-15_real64, name//"none reaches: the minimum 1/(n+1)")
        call label_permutation_p_values(4_int32, 1_int32, observed_one, null_high, p_one)
        call assert_equal_real(p_one(1), 1.0_real64, 1.0e-15_real64, name//"all reach: p = 1")
    end subroutine test_label_permutation_p_values_by_hand

    !> Test helper: the binning options of case `i_case` of test_forwards_binning_options. Case 1 switches the bin count
    !| to the occupancy search, every other case sets one more option on top of it, so each one can change the result.
    subroutine binning_case(i_case, options, label)
        integer(int32), intent(in) :: i_case
            !! Case number, 1 to N_BINNING_CASES
        type(option_set_t), intent(out) :: options
            !! The case's options
        character(len=40), intent(out) :: label
            !! The changed option's name

        options%bin_count_mode = MODE_BIN_COUNT_OCCUPANCY
        label = "bin_count_mode"
        select case (i_case)
        case (2)
            options%m_min = 40_int32
            label = "m_min"
        case (3)
            options%m_max = 4_int32
            label = "m_max"
        case (4)
            options%min_residuals_per_bin = 3_int32
            label = "min_residuals_per_bin"
        case (5)
            options%gamma_occupancy = 2.0_real64
            label = "gamma_occupancy"
        case (6)
            options%lower_residual_range_quantile = 0.2_real64
            label = "lower_residual_range_quantile"
        case (7)
            options%upper_residual_range_quantile = 0.8_real64
            label = "upper_residual_range_quantile"
        end select
    end subroutine binning_case

    !> The binning options reach the run through all four routines. For each option, assignment run on the identity and whole test
    !| equal the plain run (run_js_comp_test, construction + run_js_comp_test_adaptive) with that option, and that value
    !| differs from the plain run without the option, so the test fails if the option is dropped or goes to another
    !| argument. Case 1 is the occupancy bin count against the default (Stone); the others are one option on top of it.
    subroutine test_forwards_binning_options()
        character(*), parameter :: name = "test_forwards_binning_options: "
        integer(int32), parameter :: N_SAMPLES = 15_int32
        real(real64), allocatable :: expression(:, :, :)
        real(real64) :: expected_knn(3), baseline_knn(3), expected_adaptive(3), baseline_adaptive(3)
        real(real64) :: result(3), observed(3), p(3), null(1, 3)
        type(option_set_t) :: options, occupancy_only, no_options
        character(len=40) :: label
        integer(int32) :: ierr, i_case, used_np, used_nn, used_ks, used_kst, used_km, n_points, n_failed

        call make_expression(300_int32, [6_int32, 5_int32, 4_int32], [0.5_real64, 0.5_real64, 0.5_real64], 41_int32, &
                             expression)
        call binning_case(1_int32, occupancy_only, label)
        do i_case = 1, 7
            call binning_case(i_case, options, label)
            ! The value the option must change: the default for case 1, the occupancy search alone otherwise.
            if (i_case == 1_int32) then
                call reference_knn(expression, N_POINTS_KNN, N_NEIGHBORS_KNN, baseline_knn, ierr, options=no_options)
                call reference_adaptive(expression, K_START, K_STEP, K_MAX, baseline_adaptive, n_points, ierr, &
                                        options=no_options)
            else
                call reference_knn(expression, N_POINTS_KNN, N_NEIGHBORS_KNN, baseline_knn, ierr, options=occupancy_only)
                call reference_adaptive(expression, K_START, K_STEP, K_MAX, baseline_adaptive, n_points, ierr, &
                                        options=occupancy_only)
            end if
            call reference_knn(expression, N_POINTS_KNN, N_NEIGHBORS_KNN, expected_knn, ierr, options=options)
            call assert_equal_int(get_err_code(ierr), ERR_OK, name//trim(label)//": k-NN reference ierr")
            call reference_adaptive(expression, K_START, K_STEP, K_MAX, expected_adaptive, n_points, ierr, options=options)
            call assert_equal_int(get_err_code(ierr), ERR_OK, name//trim(label)//": adaptive reference ierr")
            call assert_true(any(expected_knn /= baseline_knn), name//trim(label)//": k-NN, the option changes the value")
            call assert_true(any(expected_adaptive /= baseline_adaptive), &
                             name//trim(label)//": adaptive, the option changes the value")

            call run_js_comp_test_on_sample_assignment(3_int32, 300_int32, 6_int32, expression, N_SAMPLES, &
                                                       identity(N_SAMPLES), N_POINTS_KNN, N_NEIGHBORS_KNN, result, &
                                                       used_np, used_nn, min_residuals_per_bin=options%min_residuals_per_bin, &
                                                       m_min=options%m_min, m_max=options%m_max, &
                                                       gamma_occupancy=options%gamma_occupancy, &
                                                       lower_residual_range_quantile=options%lower_residual_range_quantile, &
                                                       upper_residual_range_quantile=options%upper_residual_range_quantile, &
                                                       bin_count_mode=options%bin_count_mode, ierr=ierr)
            call assert_equal_int(get_err_code(ierr), ERR_OK, name//trim(label)//": k-NN assignment run ierr")
            call assert_true(equal_and_not_nan(result, expected_knn), name//trim(label)//": k-NN assignment run")
            call run_js_comp_test_permutation(3_int32, 300_int32, 6_int32, expression, N_POINTS_KNN, N_NEIGHBORS_KNN, &
                                              1_int32, observed, p, null, n_failed, &
                                              min_residuals_per_bin=options%min_residuals_per_bin, m_min=options%m_min, &
                                              m_max=options%m_max, gamma_occupancy=options%gamma_occupancy, &
                                              lower_residual_range_quantile=options%lower_residual_range_quantile, &
                                              upper_residual_range_quantile=options%upper_residual_range_quantile, &
                                              bin_count_mode=options%bin_count_mode, ierr=ierr)
            call assert_equal_int(get_err_code(ierr), ERR_OK, name//trim(label)//": k-NN whole test ierr")
            call assert_true(equal_and_not_nan(observed, expected_knn), name//trim(label)//": k-NN whole test")

            call run_js_comp_test_adaptive_on_sample_assignment(3_int32, 300_int32, 6_int32, expression, N_SAMPLES, &
                                                                identity(N_SAMPLES), K_START, K_STEP, K_MAX, result, &
                                                                used_ks, used_kst, used_km, n_points, &
                                                                min_residuals_per_bin=options%min_residuals_per_bin, &
                                                                m_min=options%m_min, m_max=options%m_max, &
                                                                gamma_occupancy=options%gamma_occupancy, &
                                                                lower_residual_range_quantile=options%lower_residual_range_quantile, &
                                                                upper_residual_range_quantile=options%upper_residual_range_quantile, &
                                                                bin_count_mode=options%bin_count_mode, ierr=ierr)
            call assert_equal_int(get_err_code(ierr), ERR_OK, name//trim(label)//": adaptive assignment run ierr")
            call assert_true(equal_and_not_nan(result, expected_adaptive), name//trim(label)//": adaptive assignment run")
            call run_js_comp_test_adaptive_permutation(3_int32, 300_int32, 6_int32, expression, K_START, K_STEP, K_MAX, &
                                                       1_int32, observed, p, null, n_failed, &
                                                       min_residuals_per_bin=options%min_residuals_per_bin, &
                                                       m_min=options%m_min, m_max=options%m_max, &
                                                       gamma_occupancy=options%gamma_occupancy, &
                                                       lower_residual_range_quantile=options%lower_residual_range_quantile, &
                                                       upper_residual_range_quantile=options%upper_residual_range_quantile, &
                                                       bin_count_mode=options%bin_count_mode, ierr=ierr)
            call assert_equal_int(get_err_code(ierr), ERR_OK, name//trim(label)//": adaptive whole test ierr")
            call assert_true(equal_and_not_nan(observed, expected_adaptive), name//trim(label)//": adaptive whole test")
        end do
    end subroutine test_forwards_binning_options

end module mod_test_data_integration_js_comp_test_permutation_common
