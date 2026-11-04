#include <src/macros.h>
#include <src/tox/data_integration/macros.h>

!> # Study-Label Permutation-Test of the JS-Comp-Test: shared helpers
!|
!| What the k-NN and the adaptive Permutation-Test
!| ([[tox_data_integration_js_comp_test_permutation_fixed_k_impl(module)]], [[tox_data_integration_js_comp_test_permutation_adaptive_impl(module)]]) share: the input checks, the relabelling of the samples and
!| the p-values. The helpers are public for those two modules only and are not published to the wrappers.
module tox_data_integration_js_comp_test_permutation_common_impl
    use f42_safeguard
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use, intrinsic :: iso_c_binding, only: c_bool
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_value, ieee_quiet_nan
    use tox_errors, only: set_ok, set_err_once, is_err, ERR_INVALID_INPUT
    use tox_data_integration_preprocessing_impl, only: compute_gene_means_impl, compute_residuals_impl
    use tox_data_integration_js_comp_test_fixed_k_search_impl, only: calc_js_comp_test_candidate_bounds_impl
    use tox_data_integration_js_comp_test_adaptive_search_impl, only: calc_adaptive_js_comp_test_bounds_impl
    M_IMPLICIT_NONE
    private
    public :: check_permutation_input, assign_samples_and_center, label_permutation_p_values

contains

    !> summary: Find the sample slots of `expression`: every replicate slot that is not all NaN
    !| AUTHOR_LASZLO_LANG
    !| Private. Sample `k` is the `k`-th slot that has at least one non-NaN value, in slot order
    !| (replicate fastest, then study); `sample_slots(k)` is its flat slot index
    !| `(i_study - 1)*max_n_reps_all_studies + i_rep`. Sequential: the sample numbering is the order.
    pure subroutine locate_samples(n_studies, max_n_genes_all_studies, max_n_reps_all_studies, expression, sample_slots, &
                                   n_samples_found, has_empty_study)
        integer(int32), intent(in) :: n_studies
            !! Number of studies
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
        integer(int32), intent(in) :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies
        real(real64), intent(in) :: expression(max_n_reps_all_studies, max_n_genes_all_studies, n_studies)
            !! Expression of every replicate slot of every gene in every study
        integer(int32), intent(out) :: sample_slots(max_n_reps_all_studies*n_studies)
            !! Flat slot index of every sample, in sample order; 0 beyond `n_samples_found`
        integer(int32), intent(out) :: n_samples_found
            !! Number of samples (slots that are not all NaN)
        logical(c_bool), intent(out) :: has_empty_study
            !! `.true.` iff some study has no sample at all

        integer(int32) :: i_study, i_rep, i_gene, n_in_study
        logical(c_bool) :: is_sample

        sample_slots = 0_int32
        n_samples_found = 0_int32
        has_empty_study = .false.
        do i_study = 1, n_studies
            n_in_study = 0_int32
            do i_rep = 1, max_n_reps_all_studies
                is_sample = .false.
                do i_gene = 1, max_n_genes_all_studies
                    if (.not. ieee_is_nan(expression(i_rep, i_gene, i_study))) then
                        is_sample = .true.
                        exit
                    end if
                end do
                if (is_sample) then
                    n_in_study = n_in_study + 1_int32
                    n_samples_found = n_samples_found + 1_int32
                    sample_slots(n_samples_found) = (i_study - 1_int32)*max_n_reps_all_studies + i_rep
                end if
            end do
            if (n_in_study == 0_int32) has_empty_study = .true.
        end do
    end subroutine locate_samples

    !> summary: Check that `sample_assignment` holds every value of `1..n_samples` exactly once
    !| AUTHOR_LASZLO_LANG
    !| Private. Sequential: `sample_mark` records the values seen so far.
    pure subroutine check_sample_assignment(n_samples, sample_assignment, sample_mark, is_permutation)
        integer(int32), intent(in) :: n_samples
            !! Number of samples
        integer(int32), intent(in) :: sample_assignment(n_samples)
            !! The assignment to check
        logical(c_bool), intent(out) :: sample_mark(n_samples)
            !! Working array: `.true.` for every value already seen
        logical(c_bool), intent(out) :: is_permutation
            !! `.true.` iff every value is in `1..n_samples` and none repeats

        integer(int32) :: i_sample, source

        sample_mark = .false.
        is_permutation = .true.
        do i_sample = 1, n_samples
            source = sample_assignment(i_sample)
            if (source < 1_int32 .or. source > n_samples) then
                is_permutation = .false.
                return
            end if
            if (sample_mark(source)) then
                is_permutation = .false.
                return
            end if
            sample_mark(source) = .true.
        end do
    end subroutine check_sample_assignment

    !> summary: Build the relabelled data set and every study's gene means and residuals
    !| AUTHOR_LASZLO_LANG
    !| Internal helper. Destination sample `k` receives all genes of source sample `sample_assignment(k)`;
    !| padding slots stay as they are (all NaN). Then, per study, exactly what a caller does per study
    !| today: [[tox_data_integration_preprocessing_impl(module):compute_gene_means_impl(interface)]] and
    !| [[tox_data_integration_preprocessing_impl(module):compute_residuals_impl(interface)]] on that
    !| study's slice (NaN padding is skipped by both).
    pure subroutine assign_samples_and_center(n_studies, max_n_genes_all_studies, max_n_reps_all_studies, expression, &
                                              n_samples, sample_assignment, sample_slots, permuted_expression, gene_means, &
                                              residuals)
        integer(int32), intent(in) :: n_studies
            !! Number of studies
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
        integer(int32), intent(in) :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies
        real(real64), intent(in) :: expression(max_n_reps_all_studies, max_n_genes_all_studies, n_studies)
            !! Expression of every replicate slot of every gene in every study
        integer(int32), intent(in) :: n_samples
            !! Number of samples
        integer(int32), intent(in) :: sample_assignment(n_samples)
            !! A permutation of `1..n_samples`, already checked
        integer(int32), intent(in) :: sample_slots(n_samples)
            !! Flat slot index of every sample, from locate_samples
        real(real64), intent(out) :: permuted_expression(max_n_reps_all_studies, max_n_genes_all_studies, n_studies)
            !! The relabelled data set
        real(real64), intent(out) :: gene_means(max_n_genes_all_studies, n_studies)
            !! Per-gene mean of every study of the relabelled data set; NaN where a study has no value
        real(real64), intent(out) :: residuals(max_n_reps_all_studies, max_n_genes_all_studies, n_studies)
            !! Signed residuals of the relabelled data set; NaN where it has no value

        integer(int32) :: i_sample, i_gene, i_study, destination, source

        permuted_expression = expression
        ! Independent: the destination slots are distinct, and every source is read from `expression`.
        do concurrent(i_sample=1:n_samples, i_gene=1:max_n_genes_all_studies) local(destination, source) &
            shared(sample_slots, sample_assignment, permuted_expression, expression, max_n_reps_all_studies)
            destination = sample_slots(i_sample) - 1_int32
            source = sample_slots(sample_assignment(i_sample)) - 1_int32
            permuted_expression(mod(destination, max_n_reps_all_studies) + 1_int32, i_gene, &
                                destination/max_n_reps_all_studies + 1_int32) = &
                expression(mod(source, max_n_reps_all_studies) + 1_int32, i_gene, source/max_n_reps_all_studies + 1_int32)
        end do

        do concurrent(i_study=1:n_studies) shared(permuted_expression, gene_means, residuals, max_n_genes_all_studies, &
                                                  max_n_reps_all_studies)
            call compute_gene_means_impl(max_n_genes_all_studies, max_n_reps_all_studies, &
                                         permuted_expression(:, :, i_study), gene_means(:, i_study))
            call compute_residuals_impl(max_n_genes_all_studies, max_n_reps_all_studies, permuted_expression(:, :, i_study), &
                                        gene_means(:, i_study), residuals(:, :, i_study))
        end do
    end subroutine assign_samples_and_center

    !> summary: Check a rerun's capacities against the fixed-k parameter search's own bounds
    !| AUTHOR_LASZLO_LANG
    !| Private. `points_fit`/`neighbors_fit` are `.true.` iff the capacity reaches the bound of
    !| [[tox_data_integration_js_comp_test_fixed_k_search_impl(module):calc_js_comp_test_candidate_bounds_impl(interface)]]
    !| for the grid options `grid_points_factor` and `grid_min_points` (absent: their defaults).
    pure subroutine check_fixed_k_rerun_capacities(max_n_genes_all_studies, max_n_points_candidate, &
                                                   max_n_neighbors_candidate, points_fit, neighbors_fit, &
                                                   grid_points_factor, grid_min_points)
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies, at least 1
        integer(int32), intent(in) :: max_n_points_candidate
            !! Capacity for reference points
        integer(int32), intent(in) :: max_n_neighbors_candidate
            !! Capacity for neighbors per neighborhood
        logical(c_bool), intent(out) :: points_fit
            !! `.true.` iff `max_n_points_candidate` reaches the search's bound
        logical(c_bool), intent(out) :: neighbors_fit
            !! `.true.` iff `max_n_neighbors_candidate` reaches the search's bound
        real(real64), intent(in), optional :: grid_points_factor
            !! The rerun search's `grid_points_factor`
        integer(int32), intent(in), optional :: grid_min_points
            !! The rerun search's `grid_min_points`

        integer(int32) :: search_n_points, search_n_neighbors

        call calc_js_comp_test_candidate_bounds_impl(max_n_genes_all_studies, search_n_points, search_n_neighbors, &
                                                     grid_points_factor, grid_min_points)
        points_fit = max_n_points_candidate >= search_n_points
        neighbors_fit = max_n_neighbors_candidate >= search_n_neighbors
    end subroutine check_fixed_k_rerun_capacities

    !> summary: The input checks every routine of the study-label Permutation-Test shares
    !| AUTHOR_LASZLO_LANG
    !| Internal helper, public for the two Permutation-Test modules and not published. Runtime checks (codegen guide Sec 5.14) of joint conditions no single `DM_` macro can state. Each
    !| failure is invalid input (`ERR_INVALID_INPUT`) at the argument position the caller passes, and returns at once:
    !|
    !| 1. A study without any sample fails at position 4 (`expression`). This step also fills `sample_slots`
    !|    and `n_samples_found`.
    !| 2. With `sample_assignment` (the per-assignment routines): `n_samples` other than the number of samples
    !|    fails at position 5, an assignment that is not a permutation of `1..n_samples` at position 6.
    !| 3. With `rerun`, the capacities against the parameter search's own bounds. Fixed-k (with
    !|    `max_n_neighbors_candidate`): both capacities against
    !|    [[tox_data_integration_js_comp_test_fixed_k_search_impl(module):calc_js_comp_test_candidate_bounds_impl(interface)]]
    !|    for `grid_points_factor` and `grid_min_points`, at `points_arg_pos` and `neighbors_arg_pos`. Adaptive (with
    !|    `search_capacity`): `search_capacity` is set from
    !|    [[tox_data_integration_js_comp_test_adaptive_search_impl(module):calc_adaptive_js_comp_test_bounds_impl(interface)]]
    !|    (its error is returned with its code and no position), and a smaller `max_n_points_candidate` fails at
    !|    `points_arg_pos`.
    pure subroutine check_permutation_input(n_studies, max_n_genes_all_studies, max_n_reps_all_studies, expression, &
                                            rerun, max_n_points_candidate, points_arg_pos, sample_slots, &
                                            n_samples_found, n_samples, sample_assignment, sample_mark, &
                                            max_n_neighbors_candidate, neighbors_arg_pos, search_capacity, &
                                            grid_points_factor, grid_min_points, ierr)
        integer(int32), intent(in) :: n_studies
            !! Number of studies
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
        integer(int32), intent(in) :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies
        real(real64), intent(in) :: expression(max_n_reps_all_studies, max_n_genes_all_studies, n_studies)
            !! Expression of every replicate slot of every gene in every study
        logical(c_bool), intent(in) :: rerun
            !! Whether the parameter search is rerun, so its capacities are checked
        integer(int32), intent(in) :: max_n_points_candidate
            !! Capacity of the per-point work arrays
        integer(int32), intent(in) :: points_arg_pos
            !! Argument position of `max_n_points_candidate` in the caller
        integer(int32), intent(out) :: sample_slots(max_n_reps_all_studies*n_studies)
            !! Flat slot index of every sample, from locate_samples
        integer(int32), intent(out) :: n_samples_found
            !! Number of samples in `expression`
        integer(int32), intent(in), optional :: n_samples
            !! The caller's number of samples; required with `sample_assignment`
        integer(int32), intent(in), optional :: sample_assignment(*)
            !! The assignment to check, of `n_samples` values
        logical(c_bool), intent(out), optional :: sample_mark(*)
            !! Work array of the assignment check, of `n_samples` values; required with `sample_assignment`
        integer(int32), intent(in), optional :: max_n_neighbors_candidate
            !! Capacity of the per-neighbor work arrays (fixed-k)
        integer(int32), intent(in), optional :: neighbors_arg_pos
            !! Argument position of `max_n_neighbors_candidate` in the caller; required with it
        integer(int32), intent(out), optional :: search_capacity
            !! Adaptive: the capacity the adaptive parameter search gets (set only with `rerun`)
        real(real64), intent(in), optional :: grid_points_factor
            !! Fixed-k: the rerun search's `grid_points_factor`
        integer(int32), intent(in), optional :: grid_min_points
            !! Fixed-k: the rerun search's `grid_min_points`
        integer(int32), intent(out) :: ierr
            !! Error code; ERR_INVALID_INPUT at the positions above, or the adaptive bounds' code

        logical(c_bool) :: has_empty_study, is_permutation, points_fit, neighbors_fit

        call set_ok(ierr)
        call locate_samples(n_studies, max_n_genes_all_studies, max_n_reps_all_studies, expression, sample_slots, &
                            n_samples_found, has_empty_study)
        if (has_empty_study) then
            call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=4_int32)
            return
        end if

        if (present(sample_assignment)) then
            if (n_samples_found /= n_samples) then
                call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=5_int32)
                return
            end if
            call check_sample_assignment(n_samples, sample_assignment, sample_mark, is_permutation)
            if (.not. is_permutation) then
                call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=6_int32)
                return
            end if
        end if

        if (.not. rerun) return
        if (present(max_n_neighbors_candidate)) then
            call check_fixed_k_rerun_capacities(max_n_genes_all_studies, max_n_points_candidate, max_n_neighbors_candidate, &
                                                points_fit, neighbors_fit, grid_points_factor, grid_min_points)
            if (.not. points_fit) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=points_arg_pos)
            if (.not. neighbors_fit) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=neighbors_arg_pos)
        else if (present(search_capacity)) then
            ! The search gets the capacity a direct call would, so it rejects the same candidates.
            ! Its joint bounds name no argument (no position).
            call calc_adaptive_js_comp_test_bounds_impl(max_n_genes_all_studies, n_studies, search_capacity, ierr)
            if (is_err(ierr)) return
            if (max_n_points_candidate < search_capacity) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=points_arg_pos)
        end if
    end subroutine check_permutation_input

    !> summary: Turn the permuted values into the study-label Permutation-Test's p-values
    !| AUTHOR_LASZLO_LANG
    !| Internal helper. Per study, `(#{valid permutations with J_perm >= J_obs} + 1) / (n_valid + 1)`, the
    !| Laplace correction (the observed data counts as one data set); NaN rows are failed permutations and do not
    !| count. NaN when no permutation is valid.
    pure subroutine label_permutation_p_values(n_permutations, n_studies, global_js_divergence_observed, jsd_null, &
                                               p_values)
        integer(int32), intent(in) :: n_permutations
            !! Number of permutations (rows of `jsd_null`)
        integer(int32), intent(in) :: n_studies
            !! Number of studies
        real(real64), intent(in) :: global_js_divergence_observed(n_studies)
            !! Observed global JSD per study
        real(real64), intent(in) :: jsd_null(n_permutations, n_studies)
            !! Per-permutation global JSD per study; a NaN row is a failed permutation
        real(real64), intent(out) :: p_values(n_studies)
            !! Empirical p-value per study

        integer(int32) :: i_study, i_permutation, n_valid, n_at_least

        do concurrent(i_study=1:n_studies) local(i_permutation, n_valid, n_at_least) &
            shared(n_permutations, jsd_null, global_js_divergence_observed, p_values)
            n_valid = 0_int32
            n_at_least = 0_int32
            do i_permutation = 1, n_permutations
                if (ieee_is_nan(jsd_null(i_permutation, i_study))) cycle
                n_valid = n_valid + 1_int32
                if (jsd_null(i_permutation, i_study) >= global_js_divergence_observed(i_study)) &
                    n_at_least = n_at_least + 1_int32
            end do
            if (n_valid == 0_int32) then
                p_values(i_study) = M_NAN
            else
                p_values(i_study) = real(n_at_least + 1_int32, real64)/real(n_valid + 1_int32, real64)
            end if
        end do
    end subroutine label_permutation_p_values
end module tox_data_integration_js_comp_test_permutation_common_impl
