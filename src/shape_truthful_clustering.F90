#include <src/macros.h>

!> Shape Truthful Clustering (STC), tangent-space variant.
!| See `misc/mod_STC.md` for the full algorithm definition and
!| `misc/STC_for_LoManLe.md` for how this is intended to integrate with
!| LoManLe (not yet wired in -- see `CLAUDE.md`'s recorded implementation
!| order: this module is built and tested standalone first).
module shape_truthful_clustering

    use, intrinsic :: iso_fortran_env, only: int32, real64
    use tox_errors, only: set_ok, set_err, is_err, ERR_ALLOC_FAIL, ERR_INVALID_INPUT, ERR_INTERNAL, &
                          validate_dimension_size, validate_in_range_int, validate_in_range_real, &
                          validate_all_in_range_real
    use tox_euclidean_distance, only: euclidean_distance
    use f42_utils, only: sort_real_heapsort, calc_percentile
    use kd_tree, only: kd_knn_query, build_kd_index, kd_range_query_list, kd_range_query_mask
    implicit none

    interface
        ! Declared pure to document that dgesdd is a deterministic numerical
        ! routine with no I/O and no global state -- mirroring dsyev's own
        ! purity-as-documentation precedent in src/lomanle.F90. This lets
        ! observable_helper below remain `pure` despite calling an external
        ! LAPACK routine.
        pure subroutine dgesdd(jobz, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, iwork, info)
            import :: int32, real64
            character,      intent(in)    :: jobz
            integer(int32), intent(in)    :: m, n, lda, ldu, ldvt, lwork
            real(real64),   intent(inout) :: a(lda, n)
            real(real64),   intent(out)   :: s(min(m, n))
            real(real64),   intent(out)   :: u(ldu, min(m, n))
            real(real64),   intent(out)   :: vt(ldvt, n)
            real(real64),   intent(out)   :: work(lwork)
            integer(int32), intent(out)   :: iwork(8 * min(m, n))
            integer(int32), intent(out)   :: info
        end subroutine dgesdd

        ! Declared pure for the same reason as dgesdd above. `u`/`vt` are
        ! declared assumed-size (`*`) since this module's only call site
        ! uses JOBU='N', JOBVT='N', for which LAPACK documents them as "not
        ! referenced" -- their exact shape is irrelevant, so assumed-size
        ! sidesteps having to pick one.
        pure subroutine dgesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info)
            import :: int32, real64
            character,      intent(in)    :: jobu, jobvt
            integer(int32), intent(in)    :: m, n, lda, ldu, ldvt, lwork
            real(real64),   intent(inout) :: a(lda, n)
            real(real64),   intent(out)   :: s(min(m, n))
            real(real64),   intent(out)   :: u(ldu, *)
            real(real64),   intent(out)   :: vt(ldvt, *)
            real(real64),   intent(out)   :: work(lwork)
            integer(int32), intent(out)   :: info
        end subroutine dgesvd
    end interface

    private

    public :: normal_error, normal_error_helper
    public :: tangent_scales, tangent_scales_helper
    public :: calculate_density_radius_alloc, calculate_density_radius, calculate_density_radius_helper
    public :: calc_ensemble_growth_radius_alloc, calc_ensemble_growth_radius, calc_ensemble_growth_radius_helper
    public :: density_labels, density_labels_helper
    public :: seeds_alloc, seeds, seeds_helper
    public :: grow_ensemble_alloc, grow_ensemble, grow_ensemble_helper
    public :: observable_alloc, observable, observable_helper
    public :: accept_ensemble_alloc, accept_ensemble, accept_ensemble_helper

contains

    !> Validated Entry Point for SKG `normal_error` (`misc/mod_STC.md`): the
    !| mean squared residual of an ensemble's members off its
    !| `d`-dimensional tangent subspace, obtained directly from the
    !| ensemble covariance's eigenvalues -- no pass over the member vectors
    !| is required, see `misc/mod_STC.md`'s SKG `normal_error`.
    pure subroutine normal_error(d, eigenvalues, n_dimensions, normal_error_value, ierr)

        integer(int32), intent(in) :: n_dimensions
            !! Ambient dimension D
        integer(int32), intent(in) :: d
            !! Intrinsic (tangent) dimension, 0 <= d <= n_dimensions
        real(real64),   intent(in) :: eigenvalues(n_dimensions)
            !! Ensemble covariance eigenvalues, descending: lambda_1 >= ... >= lambda_D >= 0
        real(real64),   intent(out) :: normal_error_value
            !! Mean squared residual off the tangent subspace
        integer(int32), intent(out) :: ierr
            !! Error code: 0 = success

        call set_ok(ierr)
        call validate_dimension_size(n_dimensions, ierr)
        call validate_in_range_int(d, ierr, min=0_int32, max=n_dimensions)
        call validate_all_in_range_real(eigenvalues, n_dimensions, ierr, min=0.0_real64)
        if (is_err(ierr)) return

        call normal_error_helper(d, eigenvalues, n_dimensions, normal_error_value)

    end subroutine normal_error

    !> Core Implementation for SKG `normal_error`.
    pure subroutine normal_error_helper(d, eigenvalues, n_dimensions, normal_error_value)

        integer(int32), intent(in) :: n_dimensions
            !! Ambient dimension D
        integer(int32), intent(in) :: d
            !! Intrinsic (tangent) dimension, 0 <= d <= n_dimensions
        real(real64),   intent(in) :: eigenvalues(n_dimensions)
            !! Ensemble covariance eigenvalues, descending
        real(real64),   intent(out) :: normal_error_value
            !! Mean squared residual off the tangent subspace

        normal_error_value = sum(eigenvalues(d+1:n_dimensions))

    end subroutine normal_error_helper

    !> Validated Entry Point for SKG `tangent_scales` (`misc/mod_STC.md`):
    !| the extent along each tangent direction, obtained directly from the
    !| ensemble covariance's tangent-space eigenvalues.
    pure subroutine tangent_scales(d, eigenvalues, n_dimensions, tangent_scales_value, ierr)

        integer(int32), intent(in) :: n_dimensions
            !! Ambient dimension D
        integer(int32), intent(in) :: d
            !! Intrinsic (tangent) dimension, 0 <= d <= n_dimensions
        real(real64),   intent(in) :: eigenvalues(n_dimensions)
            !! Ensemble covariance eigenvalues, descending: lambda_1 >= ... >= lambda_D >= 0
        real(real64),   intent(out) :: tangent_scales_value(d)
            !! Extent along each of the d tangent directions
        integer(int32), intent(out) :: ierr
            !! Error code: 0 = success

        call set_ok(ierr)
        call validate_dimension_size(n_dimensions, ierr)
        call validate_in_range_int(d, ierr, min=0_int32, max=n_dimensions)
        call validate_all_in_range_real(eigenvalues, n_dimensions, ierr, min=0.0_real64)
        if (is_err(ierr)) return

        call tangent_scales_helper(d, eigenvalues, n_dimensions, tangent_scales_value)

    end subroutine tangent_scales

    !> Core Implementation for SKG `tangent_scales`.
    pure subroutine tangent_scales_helper(d, eigenvalues, n_dimensions, tangent_scales_value)

        integer(int32), intent(in) :: n_dimensions
            !! Ambient dimension D
        integer(int32), intent(in) :: d
            !! Intrinsic (tangent) dimension, 0 <= d <= n_dimensions
        real(real64),   intent(in) :: eigenvalues(n_dimensions)
            !! Ensemble covariance eigenvalues, descending
        real(real64),   intent(out) :: tangent_scales_value(d)
            !! Extent along each of the d tangent directions

        tangent_scales_value = sqrt(eigenvalues(1:d))

    end subroutine tangent_scales_helper

    !> Allocating Wrapper for SKG `calculate_density_radius` (`misc/mod_STC.md`).
    subroutine calculate_density_radius_alloc(vectors, n_dimensions, n_vectors, &
                                              radius, mean_to_other_vecs_dist_quant, ierr)

        integer(int32), intent(in) :: n_dimensions
            !! Ambient dimension D
        integer(int32), intent(in) :: n_vectors
            !! Number of input vectors N
        real(real64),   intent(in) :: vectors(n_dimensions, n_vectors)
            !! Input data matrix
        real(real64),   intent(out) :: radius
            !! Resulting density search radius
        real(real64),   intent(in), optional :: mean_to_other_vecs_dist_quant
            !! Percentile (0.0 to 1.0) of mean-to-other-vector distances used
            !! as the radius; defaults to 0.15 (the 15th percentile)
        integer(int32), intent(out) :: ierr
            !! Error code: 0 = success

        real(real64),   allocatable :: tmp_mean_vec(:)
        real(real64),   allocatable :: tmp_distances(:)
        integer(int32), allocatable :: tmp_perm(:)

        call set_ok(ierr)
        call validate_dimension_size(n_dimensions, ierr)
        call validate_dimension_size(n_vectors, ierr)
        call validate_in_range_real(mean_to_other_vecs_dist_quant, ierr, min=0.0_real64, max=1.0_real64)
        if (is_err(ierr)) return

        M_ALLOCATE(tmp_mean_vec(n_dimensions))
        M_ALLOCATE(tmp_distances(n_vectors))
        M_ALLOCATE(tmp_perm(n_vectors))

        call calculate_density_radius(vectors, n_dimensions, n_vectors, &
                                      tmp_mean_vec, tmp_distances, tmp_perm, &
                                      radius, mean_to_other_vecs_dist_quant, ierr)

    end subroutine calculate_density_radius_alloc

    !> Validated Entry Point for SKG `calculate_density_radius`.
    subroutine calculate_density_radius(vectors, n_dimensions, n_vectors, &
                                        mean_vec, distances, perm, radius, &
                                        mean_to_other_vecs_dist_quant, ierr)

        integer(int32), intent(in) :: n_dimensions
            !! Ambient dimension D
        integer(int32), intent(in) :: n_vectors
            !! Number of input vectors N
        real(real64),   intent(in) :: vectors(n_dimensions, n_vectors)
            !! Input data matrix
        real(real64),   intent(out) :: mean_vec(n_dimensions)
            !! Preallocated workspace for the global mean vector
        real(real64),   intent(out) :: distances(n_vectors)
            !! Preallocated workspace for mean-to-vector distances
        integer(int32), intent(out) :: perm(n_vectors)
            !! Preallocated workspace for the sort permutation
        real(real64),   intent(out) :: radius
            !! Resulting density search radius
        real(real64),   intent(in), optional :: mean_to_other_vecs_dist_quant
            !! Percentile (0.0 to 1.0) of mean-to-other-vector distances used
            !! as the radius; defaults to 0.15 (the 15th percentile)
        integer(int32), intent(out) :: ierr
            !! Error code: 0 = success

        call set_ok(ierr)
        call validate_dimension_size(n_dimensions, ierr)
        call validate_dimension_size(n_vectors, ierr)
        call validate_in_range_real(mean_to_other_vecs_dist_quant, ierr, min=0.0_real64, max=1.0_real64)
        if (is_err(ierr)) return

        call calculate_density_radius_helper(vectors, n_dimensions, n_vectors, &
                                             mean_vec, distances, perm, &
                                             mean_to_other_vecs_dist_quant, radius, ierr)

    end subroutine calculate_density_radius

    !> Core Implementation for SKG `calculate_density_radius`. Not `pure`:
    !| the per-vector distance loop uses `!$omp parallel do` rather than
    !| `do concurrent`, since it calls the external procedure
    !| `euclidean_distance` -- gfortran's `do concurrent` auto-parallelizer
    !| refuses to parallelize a loop whose body calls an external procedure,
    !| even a `pure` one (confirmed for this exact class of loop in
    !| `src/lomanle.F90`; see that module's docstring and
    !| `misc/STC_for_LoManLe.md` section 5), and OpenMP worksharing
    !| directives are not allowed inside a `pure` procedure.
    subroutine calculate_density_radius_helper(vectors, n_dimensions, n_vectors, &
                                               mean_vec, distances, perm, &
                                               mean_to_other_vecs_dist_quant, radius, ierr)

        integer(int32), intent(in) :: n_dimensions
            !! Ambient dimension D
        integer(int32), intent(in) :: n_vectors
            !! Number of input vectors N
        real(real64),   intent(in) :: vectors(n_dimensions, n_vectors)
            !! Input data matrix
        real(real64),   intent(out) :: mean_vec(n_dimensions)
            !! Preallocated workspace for the global mean vector
        real(real64),   intent(out) :: distances(n_vectors)
            !! Preallocated workspace for mean-to-vector distances
        integer(int32), intent(out) :: perm(n_vectors)
            !! Preallocated workspace for the sort permutation
        real(real64),   intent(in), optional :: mean_to_other_vecs_dist_quant
            !! Percentile (0.0 to 1.0) of mean-to-other-vector distances used
            !! as the radius; defaults to 0.15 (the 15th percentile)
        real(real64),   intent(out) :: radius
            !! Resulting density search radius
        integer(int32), intent(out) :: ierr
            !! Error code: 0 = success

        real(real64)   :: actual_quant, percentile
        integer(int32) :: i_vec

        call set_ok(ierr)
        M_DEFAULT_VAL(mean_to_other_vecs_dist_quant, actual_quant, 0.15_real64)

        mean_vec = sum(vectors, dim=2) / real(n_vectors, real64)

        !$omp parallel do default(shared)
        do i_vec = 1, n_vectors
            call euclidean_distance(mean_vec, vectors(:, i_vec), n_dimensions, distances(i_vec))
        end do
        !$omp end parallel do

        do i_vec = 1, n_vectors
            perm(i_vec) = i_vec
        end do

        call sort_real_heapsort(distances, perm)

        percentile = actual_quant * 100.0_real64
        call calc_percentile(distances, perm, percentile, radius, ierr)

    end subroutine calculate_density_radius_helper

    !> Allocating Wrapper for SKG `calc_ensemble_growth_radius` (`misc/mod_STC.md`).
    !| Builds no k-d tree itself: `kd_indices`/`dimension_order` are a tree
    !| already built once over the whole dataset (see `build_kd_index`,
    !| `kd_tree` module) and reused across every seed, rather than rebuilt
    !| per call.
    subroutine calc_ensemble_growth_radius_alloc(vectors, n_dimensions, n_vectors, &
                                                  kd_indices, dimension_order, seed_index, &
                                                  growth_radius, k_min, ierr)

        integer(int32), intent(in) :: n_dimensions
            !! Ambient dimension D
        integer(int32), intent(in) :: n_vectors
            !! Number of input vectors N
        real(real64),   intent(in) :: vectors(n_dimensions, n_vectors)
            !! Input data matrix
        integer(int32), intent(in) :: kd_indices(n_vectors)
            !! Pre-built k-d tree indices, see `kd_tree` module
        integer(int32), intent(in) :: dimension_order(n_dimensions)
            !! k-d tree dimension-splitting order, see `kd_tree` module
        integer(int32), intent(in) :: seed_index
            !! Index into `vectors`/`kd_indices` of the seed to compute the
            !! growth radius for
        real(real64),   intent(out) :: growth_radius
            !! Median distance among the seed's own k_min nearest neighbors
        integer(int32), intent(in), optional :: k_min
            !! Neighborhood size the median is taken over; defaults to 30
        integer(int32), intent(out) :: ierr
            !! Error code: 0 = success

        integer(int32) :: actual_k_min
        integer(int32), allocatable :: tmp_neighbors(:)
        real(real64),   allocatable :: tmp_distances(:)
        integer(int32), allocatable :: tmp_sort_perm(:)

        call set_ok(ierr)
        call validate_dimension_size(n_dimensions, ierr)
        call validate_dimension_size(n_vectors, ierr)
        call validate_in_range_int(seed_index, ierr, min=1_int32, max=n_vectors)
        if (is_err(ierr)) return

        M_DEFAULT_VAL(k_min, actual_k_min, 30_int32)
        call validate_in_range_int(actual_k_min, ierr, min=1_int32, max=n_vectors-1)
        if (is_err(ierr)) return

        M_ALLOCATE(tmp_neighbors(actual_k_min + 1))
        M_ALLOCATE(tmp_distances(actual_k_min + 1))
        M_ALLOCATE(tmp_sort_perm(actual_k_min))

        call calc_ensemble_growth_radius(vectors, n_dimensions, n_vectors, kd_indices, dimension_order, &
                                         seed_index, tmp_neighbors, tmp_distances, tmp_sort_perm, &
                                         growth_radius, k_min, ierr)

    end subroutine calc_ensemble_growth_radius_alloc

    !> Validated Entry Point for SKG `calc_ensemble_growth_radius`.
    subroutine calc_ensemble_growth_radius(vectors, n_dimensions, n_vectors, &
                                           kd_indices, dimension_order, seed_index, &
                                           tmp_neighbors, tmp_distances, tmp_sort_perm, &
                                           growth_radius, k_min, ierr)

        integer(int32), intent(in) :: n_dimensions
            !! Ambient dimension D
        integer(int32), intent(in) :: n_vectors
            !! Number of input vectors N
        real(real64),   intent(in) :: vectors(n_dimensions, n_vectors)
            !! Input data matrix
        integer(int32), intent(in) :: kd_indices(n_vectors)
            !! Pre-built k-d tree indices, see `kd_tree` module
        integer(int32), intent(in) :: dimension_order(n_dimensions)
            !! k-d tree dimension-splitting order, see `kd_tree` module
        integer(int32), intent(in) :: seed_index
            !! Index into `vectors`/`kd_indices` of the seed to compute the
            !! growth radius for
        integer(int32), intent(out) :: tmp_neighbors(:)
            !! Preallocated workspace, size >= (resolved k_min)+1
        real(real64),   intent(out) :: tmp_distances(:)
            !! Preallocated workspace, size >= (resolved k_min)+1
        integer(int32), intent(out) :: tmp_sort_perm(:)
            !! Preallocated workspace, size >= resolved k_min
        real(real64),   intent(out) :: growth_radius
            !! Median distance among the seed's own k_min nearest neighbors
        integer(int32), intent(in), optional :: k_min
            !! Neighborhood size the median is taken over; defaults to 30
        integer(int32), intent(out) :: ierr
            !! Error code: 0 = success

        integer(int32) :: actual_k_min

        call set_ok(ierr)
        call validate_dimension_size(n_dimensions, ierr)
        call validate_dimension_size(n_vectors, ierr)
        call validate_in_range_int(seed_index, ierr, min=1_int32, max=n_vectors)
        if (is_err(ierr)) return

        M_DEFAULT_VAL(k_min, actual_k_min, 30_int32)
        call validate_in_range_int(actual_k_min, ierr, min=1_int32, max=n_vectors-1)
        if (is_err(ierr)) return

        call calc_ensemble_growth_radius_helper(vectors, n_dimensions, n_vectors, kd_indices, dimension_order, &
                                                seed_index, actual_k_min, &
                                                tmp_neighbors(1:actual_k_min+1), tmp_distances(1:actual_k_min+1), &
                                                tmp_sort_perm(1:actual_k_min), growth_radius, ierr)

    end subroutine calc_ensemble_growth_radius

    !> Core Implementation for SKG `calc_ensemble_growth_radius`: the median
    !| distance among the seed's own k_min nearest neighbors (excluding the
    !| seed itself), matching LoManLe's `local_scale_i` exactly -- see
    !| `misc/STC_for_LoManLe.md` section 2.2 for why this must stay a
    !| per-seed, locally adaptive radius rather than a single dataset-wide
    !| one.
    pure subroutine calc_ensemble_growth_radius_helper(vectors, n_dimensions, n_vectors, &
                                                        kd_indices, dimension_order, seed_index, k_min, &
                                                        tmp_neighbors, tmp_distances, tmp_sort_perm, &
                                                        growth_radius, ierr)

        integer(int32), intent(in) :: n_dimensions
            !! Ambient dimension D
        integer(int32), intent(in) :: n_vectors
            !! Number of input vectors N
        real(real64),   intent(in) :: vectors(n_dimensions, n_vectors)
            !! Input data matrix
        integer(int32), intent(in) :: kd_indices(n_vectors)
            !! Pre-built k-d tree indices, see `kd_tree` module
        integer(int32), intent(in) :: dimension_order(n_dimensions)
            !! k-d tree dimension-splitting order, see `kd_tree` module
        integer(int32), intent(in) :: seed_index
            !! Index into `vectors`/`kd_indices` of the seed to compute the
            !! growth radius for
        integer(int32), intent(in) :: k_min
            !! Neighborhood size the median is taken over
        integer(int32), intent(out) :: tmp_neighbors(k_min+1)
            !! Workspace: k-NN query result, indices
        real(real64),   intent(out) :: tmp_distances(k_min+1)
            !! Workspace: k-NN query result, distances
        integer(int32), intent(out) :: tmp_sort_perm(k_min)
            !! Workspace: sort-permutation scratch
        real(real64),   intent(out) :: growth_radius
            !! Median distance among the seed's own k_min nearest neighbors
        integer(int32), intent(out) :: ierr
            !! Error code: 0 = success

        integer(int32) :: self_pos, j, k_query

        k_query = k_min + 1
        call kd_knn_query(vectors, kd_indices, n_dimensions, n_vectors, dimension_order, &
                          vectors(:, seed_index), k_query, tmp_neighbors, tmp_distances, ierr)
        if (is_err(ierr)) return

        ! Exclude the seed itself, at distance 0, from its own k_min+1
        ! nearest neighbors: swap it into the last slot, then only the
        ! first k_min entries are used from here on. kd_knn_query does not
        ! guarantee its output is sorted by distance (max-heap based
        ! internally), so this swap-with-last is safe irrespective of
        ! result order.
        self_pos = 0
        do j = 1, k_query
            if (tmp_neighbors(j) == seed_index) then
                self_pos = j
                exit
            end if
        end do
        if (self_pos > 0 .and. self_pos < k_query) then
            tmp_neighbors(self_pos) = tmp_neighbors(k_query)
            tmp_distances(self_pos) = tmp_distances(k_query)
        end if

        ! Sort the remaining k_min distances ascending so the median below
        ! is meaningful.
        do j = 1, k_min
            tmp_sort_perm(j) = j
        end do
        call sort_real_heapsort(tmp_distances(1:k_min), tmp_sort_perm(1:k_min))

        if (mod(k_min, 2) == 1) then
            growth_radius = tmp_distances(tmp_sort_perm((k_min + 1) / 2))
        else
            growth_radius = 0.5_real64 * ( &
                tmp_distances(tmp_sort_perm(k_min / 2)) + &
                tmp_distances(tmp_sort_perm(k_min / 2 + 1)))
        end if

    end subroutine calc_ensemble_growth_radius_helper

    !> Validated Entry Point for `density_labels` (`misc/mod_STC.md`,
    !| referenced under SKG `seeds`): for every vector, the count of
    !| vectors -- including itself -- within `radius`,
    !| $\rho_i=\sum_j \mathbf{1}(d(v_i,v_j)\le r)$.
    subroutine density_labels(vectors, n_dimensions, n_vectors, kd_indices, dimension_order, radius, labels, ierr)

        integer(int32), intent(in) :: n_dimensions
            !! Ambient dimension D
        integer(int32), intent(in) :: n_vectors
            !! Number of input vectors N
        real(real64),   intent(in) :: vectors(n_dimensions, n_vectors)
            !! Input data matrix
        integer(int32), intent(in) :: kd_indices(n_vectors)
            !! Pre-built k-d tree indices, see `kd_tree` module
        integer(int32), intent(in) :: dimension_order(n_dimensions)
            !! k-d tree dimension-splitting order, see `kd_tree` module
        real(real64),   intent(in) :: radius
            !! Density search radius, see `calculate_density_radius`
        real(real64),   intent(out) :: labels(n_vectors)
            !! Per-vector density label
        integer(int32), intent(out) :: ierr
            !! Error code: 0 = success

        call set_ok(ierr)
        call validate_dimension_size(n_dimensions, ierr)
        call validate_dimension_size(n_vectors, ierr)
        call validate_in_range_real(radius, ierr, min=0.0_real64)
        if (is_err(ierr)) return

        call density_labels_helper(vectors, n_dimensions, n_vectors, kd_indices, dimension_order, radius, labels)

    end subroutine density_labels

    !> Core Implementation for `density_labels`. Not `pure`: uses
    !| `!$omp parallel do`, not `do concurrent`, since the loop body calls
    !| the external procedure `kd_range_query_list` -- see
    !| `calculate_density_radius_helper`'s docstring for why. Each
    !| iteration's workspace is a `block`-local automatic array, private
    !| per iteration/thread by construction -- deliberately not a dense
    !| O(N^2) structure (see `misc/STC_for_LoManLe.md` section 5 for why
    !| that matters at real dataset sizes).
    subroutine density_labels_helper(vectors, n_dimensions, n_vectors, kd_indices, dimension_order, radius, labels)

        integer(int32), intent(in) :: n_dimensions
            !! Ambient dimension D
        integer(int32), intent(in) :: n_vectors
            !! Number of input vectors N
        real(real64),   intent(in) :: vectors(n_dimensions, n_vectors)
            !! Input data matrix
        integer(int32), intent(in) :: kd_indices(n_vectors)
            !! Pre-built k-d tree indices, see `kd_tree` module
        integer(int32), intent(in) :: dimension_order(n_dimensions)
            !! k-d tree dimension-splitting order, see `kd_tree` module
        real(real64),   intent(in) :: radius
            !! Density search radius, see `calculate_density_radius`
        real(real64),   intent(out) :: labels(n_vectors)
            !! Per-vector density label

        integer(int32) :: i_vec

        !$omp parallel do default(shared)
        do i_vec = 1, n_vectors
            block
                integer(int32) :: range_buf(n_vectors), n_found
                call kd_range_query_list(vectors, kd_indices, n_dimensions, n_vectors, dimension_order, &
                                         vectors(:, i_vec), radius, range_buf, n_found)
                labels(i_vec) = real(n_found, real64)
            end block
        end do
        !$omp end parallel do

    end subroutine density_labels_helper

    !> Allocating Wrapper for SKG `seeds` (`misc/mod_STC.md`). Builds its
    !| own k-d tree internally, matching the SKG's stated interface (input:
    !| just `data_vectors`) -- `seeds` below takes an already-built one
    !| instead, so a caller that already has one (e.g. the eventual
    !| ensemble-growth orchestration, which also needs a tree) is not
    !| forced to rebuild it.
    subroutine seeds_alloc(vectors, n_dimensions, n_vectors, is_seed_mask, &
                           mean_to_other_vecs_dist_quant, ierr)

        integer(int32), intent(in) :: n_dimensions
            !! Ambient dimension D
        integer(int32), intent(in) :: n_vectors
            !! Number of input vectors N
        real(real64),   intent(in) :: vectors(n_dimensions, n_vectors)
            !! Input data matrix
        logical,        intent(out) :: is_seed_mask(n_vectors)
            !! .true. for points selected as seeds
        real(real64),   intent(in), optional :: mean_to_other_vecs_dist_quant
            !! Density-radius percentile (0.0 to 1.0); defaults to 0.15
        integer(int32), intent(out) :: ierr
            !! Error code: 0 = success

        integer(int32), allocatable :: tmp_kd_indices(:), tmp_dim_order(:)
        integer(int32), allocatable :: tmp_kd_workspace(:), tmp_kd_perm(:)
        integer(int32), allocatable :: tmp_kd_lstack(:), tmp_kd_rstack(:), tmp_kd_recstack(:,:)
        real(real64),   allocatable :: tmp_kd_valbuf(:)
        real(real64),   allocatable :: tmp_density_mean_vec(:), tmp_density_distances(:)
        integer(int32), allocatable :: tmp_density_radius_perm(:)
        real(real64),   allocatable :: tmp_labels(:)
        integer(int32), allocatable :: tmp_rank_perm(:)
        logical,        allocatable :: tmp_visited_mask(:), tmp_newly_covered_mask(:)
        integer(int32) :: i_dim

        call set_ok(ierr)
        call validate_dimension_size(n_dimensions, ierr)
        call validate_dimension_size(n_vectors, ierr)
        call validate_in_range_real(mean_to_other_vecs_dist_quant, ierr, min=0.0_real64, max=1.0_real64)
        if (is_err(ierr)) return

        M_ALLOCATE(tmp_kd_indices(n_vectors))
        M_ALLOCATE(tmp_dim_order(n_dimensions))
        M_ALLOCATE(tmp_kd_workspace(n_vectors))
        M_ALLOCATE(tmp_kd_perm(n_vectors))
        M_ALLOCATE(tmp_kd_lstack(n_vectors))
        M_ALLOCATE(tmp_kd_rstack(n_vectors))
        M_ALLOCATE(tmp_kd_recstack(3, n_vectors))
        M_ALLOCATE(tmp_kd_valbuf(n_vectors))
        M_ALLOCATE(tmp_density_mean_vec(n_dimensions))
        M_ALLOCATE(tmp_density_distances(n_vectors))
        M_ALLOCATE(tmp_density_radius_perm(n_vectors))
        M_ALLOCATE(tmp_labels(n_vectors))
        M_ALLOCATE(tmp_rank_perm(n_vectors))
        M_ALLOCATE(tmp_visited_mask(n_vectors))
        M_ALLOCATE(tmp_newly_covered_mask(n_vectors))

        do i_dim = 1, n_dimensions
            tmp_dim_order(i_dim) = i_dim
        end do
        call build_kd_index(vectors, n_dimensions, n_vectors, tmp_kd_indices, tmp_dim_order, &
                            tmp_kd_workspace, tmp_kd_valbuf, tmp_kd_perm, tmp_kd_lstack, tmp_kd_rstack, &
                            tmp_kd_recstack, ierr)
        if (is_err(ierr)) return

        call seeds(vectors, n_dimensions, n_vectors, tmp_kd_indices, tmp_dim_order, &
                  tmp_density_mean_vec, tmp_density_distances, tmp_density_radius_perm, &
                  tmp_labels, tmp_rank_perm, tmp_visited_mask, tmp_newly_covered_mask, &
                  is_seed_mask, mean_to_other_vecs_dist_quant, ierr)

    end subroutine seeds_alloc

    !> Validated Entry Point for SKG `seeds`. Takes an already-built k-d
    !| tree (`kd_indices`, `dimension_order`) rather than building one, for
    !| the reason given on `seeds_alloc`.
    subroutine seeds(vectors, n_dimensions, n_vectors, kd_indices, dimension_order, &
                     density_mean_vec, density_distances, density_radius_perm, &
                     labels, rank_perm, visited_mask, newly_covered_mask, &
                     is_seed_mask, mean_to_other_vecs_dist_quant, ierr)

        integer(int32), intent(in) :: n_dimensions
            !! Ambient dimension D
        integer(int32), intent(in) :: n_vectors
            !! Number of input vectors N
        real(real64),   intent(in) :: vectors(n_dimensions, n_vectors)
            !! Input data matrix
        integer(int32), intent(in) :: kd_indices(n_vectors)
            !! Pre-built k-d tree indices, see `kd_tree` module
        integer(int32), intent(in) :: dimension_order(n_dimensions)
            !! k-d tree dimension-splitting order, see `kd_tree` module
        real(real64),   intent(inout) :: density_mean_vec(n_dimensions)
            !! Preallocated workspace for `calculate_density_radius_helper`
        real(real64),   intent(inout) :: density_distances(n_vectors)
            !! Preallocated workspace for `calculate_density_radius_helper`
        integer(int32), intent(inout) :: density_radius_perm(n_vectors)
            !! Preallocated workspace for `calculate_density_radius_helper`
        real(real64),   intent(inout) :: labels(n_vectors)
            !! Preallocated workspace: per-vector density labels
        integer(int32), intent(inout) :: rank_perm(n_vectors)
            !! Preallocated workspace: density-descending sort permutation
        logical,        intent(inout) :: visited_mask(n_vectors)
            !! Preallocated workspace: coverage tracker across the greedy loop
        logical,        intent(inout) :: newly_covered_mask(n_vectors)
            !! Preallocated workspace: per-candidate range-query result
        logical,        intent(out) :: is_seed_mask(n_vectors)
            !! .true. for points selected as seeds
        real(real64),   intent(in), optional :: mean_to_other_vecs_dist_quant
            !! Density-radius percentile (0.0 to 1.0); defaults to 0.15
        integer(int32), intent(out) :: ierr
            !! Error code: 0 = success

        call set_ok(ierr)
        call validate_dimension_size(n_dimensions, ierr)
        call validate_dimension_size(n_vectors, ierr)
        call validate_in_range_real(mean_to_other_vecs_dist_quant, ierr, min=0.0_real64, max=1.0_real64)
        if (is_err(ierr)) return

        call seeds_helper(vectors, n_dimensions, n_vectors, kd_indices, dimension_order, &
                          density_mean_vec, density_distances, density_radius_perm, &
                          labels, rank_perm, visited_mask, newly_covered_mask, &
                          is_seed_mask, mean_to_other_vecs_dist_quant, ierr)

    end subroutine seeds

    !> Core Implementation for SKG `seeds`: rank vectors by density,
    !| descending, and greedily select seeds -- each pick marks every
    !| vector within the density radius as covered, so only genuinely
    !| uncovered regions can seed a subsequent ensemble
    !| (`misc/STC_for_LoManLe.md` section 2.1: sequential, density-ranked,
    !| excluding already-covered points -- what avoids seeding many
    !| near-duplicate ensembles in the same dense region). Not `pure`:
    !| transitively calls `calculate_density_radius_helper` and
    !| `density_labels_helper`, neither of which is `pure` (both use
    !| `!$omp parallel do`); the greedy loop itself is inherently
    !| sequential -- each pick depends on the previous ones' coverage --
    !| so nothing here is independently parallelizable.
    subroutine seeds_helper(vectors, n_dimensions, n_vectors, kd_indices, dimension_order, &
                            density_mean_vec, density_distances, density_radius_perm, &
                            labels, rank_perm, visited_mask, newly_covered_mask, &
                            is_seed_mask, mean_to_other_vecs_dist_quant, ierr)

        integer(int32), intent(in) :: n_dimensions
            !! Ambient dimension D
        integer(int32), intent(in) :: n_vectors
            !! Number of input vectors N
        real(real64),   intent(in) :: vectors(n_dimensions, n_vectors)
            !! Input data matrix
        integer(int32), intent(in) :: kd_indices(n_vectors)
            !! Pre-built k-d tree indices, see `kd_tree` module
        integer(int32), intent(in) :: dimension_order(n_dimensions)
            !! k-d tree dimension-splitting order, see `kd_tree` module
        real(real64),   intent(inout) :: density_mean_vec(n_dimensions)
            !! Workspace for `calculate_density_radius_helper`
        real(real64),   intent(inout) :: density_distances(n_vectors)
            !! Workspace for `calculate_density_radius_helper`
        integer(int32), intent(inout) :: density_radius_perm(n_vectors)
            !! Workspace for `calculate_density_radius_helper`
        real(real64),   intent(inout) :: labels(n_vectors)
            !! Workspace: per-vector density labels
        integer(int32), intent(inout) :: rank_perm(n_vectors)
            !! Workspace: density-descending sort permutation
        logical,        intent(inout) :: visited_mask(n_vectors)
            !! Workspace: coverage tracker across the greedy loop
        logical,        intent(inout) :: newly_covered_mask(n_vectors)
            !! Workspace: per-candidate range-query result
        logical,        intent(out) :: is_seed_mask(n_vectors)
            !! .true. for points selected as seeds
        real(real64),   intent(in), optional :: mean_to_other_vecs_dist_quant
            !! Density-radius percentile (0.0 to 1.0); defaults to 0.15
        integer(int32), intent(out) :: ierr
            !! Error code: 0 = success

        real(real64)   :: radius
        integer(int32) :: i, rank, candidate, swap_tmp

        call calculate_density_radius_helper(vectors, n_dimensions, n_vectors, &
                                             density_mean_vec, density_distances, density_radius_perm, &
                                             mean_to_other_vecs_dist_quant, radius, ierr)
        if (is_err(ierr)) return

        call density_labels_helper(vectors, n_dimensions, n_vectors, kd_indices, dimension_order, radius, labels)

        ! Rank vectors by density, descending: heapsort gives ascending
        ! order, so reverse the resulting permutation.
        do i = 1, n_vectors
            rank_perm(i) = i
        end do
        call sort_real_heapsort(labels, rank_perm)
        do i = 1, n_vectors / 2
            swap_tmp = rank_perm(i)
            rank_perm(i) = rank_perm(n_vectors - i + 1)
            rank_perm(n_vectors - i + 1) = swap_tmp
        end do

        is_seed_mask = .false.
        visited_mask = .false.

        do rank = 1, n_vectors
            candidate = rank_perm(rank)
            if (visited_mask(candidate)) cycle
            is_seed_mask(candidate) = .true.
            call kd_range_query_mask(vectors, kd_indices, n_dimensions, n_vectors, dimension_order, &
                                     vectors(:, candidate), radius, newly_covered_mask)
            visited_mask = visited_mask .or. newly_covered_mask
        end do

    end subroutine seeds_helper

    !> Allocating Wrapper for SKG `grow_ensemble` (`misc/mod_STC.md`). Takes
    !| an already-built k-d tree (`kd_indices`, `dimension_order`), matching
    !| `calc_ensemble_growth_radius`'s and `seeds`' own pattern: the tree is
    !| built once over the whole dataset and reused across every seed and
    !| every growth step, not rebuilt per call.
    subroutine grow_ensemble_alloc(vectors, n_dimensions, n_vectors, kd_indices, dimension_order, &
                                   is_member_mask, growth_radius, is_member_mask_next, ierr)

        integer(int32), intent(in) :: n_dimensions
            !! Ambient dimension D
        integer(int32), intent(in) :: n_vectors
            !! Number of input vectors N
        real(real64),   intent(in) :: vectors(n_dimensions, n_vectors)
            !! Input data matrix
        integer(int32), intent(in) :: kd_indices(n_vectors)
            !! Pre-built k-d tree indices, see `kd_tree` module
        integer(int32), intent(in) :: dimension_order(n_dimensions)
            !! k-d tree dimension-splitting order, see `kd_tree` module
        logical,        intent(in) :: is_member_mask(n_vectors)
            !! Current ensemble membership; must have at least one member
        real(real64),   intent(in) :: growth_radius
            !! This ensemble's growth radius, see `calc_ensemble_growth_radius`
        logical,        intent(out) :: is_member_mask_next(n_vectors)
            !! Grown ensemble membership (superset of `is_member_mask`)
        integer(int32), intent(out) :: ierr
            !! Error code: 0 = success

        logical, allocatable :: tmp_member_mask_buf(:)

        call set_ok(ierr)
        call validate_dimension_size(n_dimensions, ierr)
        call validate_dimension_size(n_vectors, ierr)
        call validate_in_range_real(growth_radius, ierr, min=0.0_real64)
        if (is_err(ierr)) return
        if (count(is_member_mask) == 0) then
            call set_err(ierr, ERR_INVALID_INPUT)
            return
        end if

        M_ALLOCATE(tmp_member_mask_buf(n_vectors))

        call grow_ensemble(vectors, n_dimensions, n_vectors, kd_indices, dimension_order, &
                           is_member_mask, growth_radius, tmp_member_mask_buf, is_member_mask_next, ierr)

    end subroutine grow_ensemble_alloc

    !> Validated Entry Point for SKG `grow_ensemble`.
    pure subroutine grow_ensemble(vectors, n_dimensions, n_vectors, kd_indices, dimension_order, &
                                  is_member_mask, growth_radius, tmp_member_mask_buf, is_member_mask_next, ierr)

        integer(int32), intent(in) :: n_dimensions
            !! Ambient dimension D
        integer(int32), intent(in) :: n_vectors
            !! Number of input vectors N
        real(real64),   intent(in) :: vectors(n_dimensions, n_vectors)
            !! Input data matrix
        integer(int32), intent(in) :: kd_indices(n_vectors)
            !! Pre-built k-d tree indices, see `kd_tree` module
        integer(int32), intent(in) :: dimension_order(n_dimensions)
            !! k-d tree dimension-splitting order, see `kd_tree` module
        logical,        intent(in) :: is_member_mask(n_vectors)
            !! Current ensemble membership; must have at least one member
        real(real64),   intent(in) :: growth_radius
            !! This ensemble's growth radius, see `calc_ensemble_growth_radius`
        logical,        intent(inout) :: tmp_member_mask_buf(n_vectors)
            !! Preallocated workspace: per-member range-query result
        logical,        intent(out) :: is_member_mask_next(n_vectors)
            !! Grown ensemble membership (superset of `is_member_mask`)
        integer(int32), intent(out) :: ierr
            !! Error code: 0 = success

        call set_ok(ierr)
        call validate_dimension_size(n_dimensions, ierr)
        call validate_dimension_size(n_vectors, ierr)
        call validate_in_range_real(growth_radius, ierr, min=0.0_real64)
        if (is_err(ierr)) return
        if (count(is_member_mask) == 0) then
            call set_err(ierr, ERR_INVALID_INPUT)
            return
        end if

        call grow_ensemble_helper(vectors, n_dimensions, n_vectors, kd_indices, dimension_order, &
                                  is_member_mask, growth_radius, tmp_member_mask_buf, is_member_mask_next)

    end subroutine grow_ensemble

    !> Core Implementation for SKG `grow_ensemble`: the union, over every
    !| current member, of the points within `growth_radius` of it --
    !| $\mathcal{C}_{\mathcal{E}_{t+1}} = \{x_k \mid \|x_k-x_i\|\le
    !| r_{\mathcal{E}} \;\exists\, x_i\in\mathcal{E}\}$. Deliberately a
    !| plain sequential loop, not `!$omp parallel do`: this SKG runs once
    !| per growth iteration per ensemble, and "Cluster identification"
    !| (`misc/mod_STC.md`) already parallelizes across ensembles/seeds at
    !| the outer level -- adding a second, nested layer of parallelism
    !| here would risk nested-parallel-region complications for little
    !| benefit, since a single ensemble's member count is typically small,
    !| especially early in growth. Kept `pure` as a result: every call it
    !| makes (`kd_range_query_mask`) is `pure`, and range queries have no
    !| failure mode of their own (any-size result, 0 to N, is always
    !| well-defined) -- unlike `calculate_density_radius_helper` and
    !| `density_labels_helper`, `grow_ensemble_helper` never needed OpenMP
    !| in the first place.
    pure subroutine grow_ensemble_helper(vectors, n_dimensions, n_vectors, kd_indices, dimension_order, &
                                         is_member_mask, growth_radius, tmp_member_mask_buf, is_member_mask_next)

        integer(int32), intent(in) :: n_dimensions
            !! Ambient dimension D
        integer(int32), intent(in) :: n_vectors
            !! Number of input vectors N
        real(real64),   intent(in) :: vectors(n_dimensions, n_vectors)
            !! Input data matrix
        integer(int32), intent(in) :: kd_indices(n_vectors)
            !! Pre-built k-d tree indices, see `kd_tree` module
        integer(int32), intent(in) :: dimension_order(n_dimensions)
            !! k-d tree dimension-splitting order, see `kd_tree` module
        logical,        intent(in) :: is_member_mask(n_vectors)
            !! Current ensemble membership
        real(real64),   intent(in) :: growth_radius
            !! This ensemble's growth radius
        logical,        intent(inout) :: tmp_member_mask_buf(n_vectors)
            !! Workspace: per-member range-query result
        logical,        intent(out) :: is_member_mask_next(n_vectors)
            !! Grown ensemble membership (superset of `is_member_mask`)

        integer(int32) :: i

        is_member_mask_next = is_member_mask

        do i = 1, n_vectors
            if (.not. is_member_mask(i)) cycle
            call kd_range_query_mask(vectors, kd_indices, n_dimensions, n_vectors, dimension_order, &
                                     vectors(:, i), growth_radius, tmp_member_mask_buf)
            is_member_mask_next = is_member_mask_next .or. tmp_member_mask_buf
        end do

    end subroutine grow_ensemble_helper

    !> Allocating Wrapper for SKG `observable` (`misc/mod_STC.md`): computes
    !| the tuple (U, d, G, mu, normal_error, tangent_scales) for the given
    !| ensemble via the economy-mode singular value decomposition (LAPACK
    !| `dgesdd`) of its centered member vectors -- never an eigendecomposition
    !| of an explicitly formed covariance matrix, see "Numerical Linear
    !| Algebra" in `misc/mod_STC.md`. `U` and `eigenvalues` are returned
    !| zero-padded to the full ambient dimension `n_dimensions`: the economy
    !| SVD only yields `rank = min(n_dimensions, n_members)` genuine
    !| columns/values, which is less than `n_dimensions` whenever an ensemble
    !| is smaller than the ambient space (typical early in growth) -- padding
    !| keeps every call's output shape fixed regardless of ensemble size,
    !| satisfies `normal_error`/`tangent_scales`'s existing
    !| `eigenvalues(n_dimensions)` interface directly with no changes to
    !| either, and matches the fixed-D-width convention already decided for
    !| `cluster_identification`'s own per-iteration history storage
    !| (`misc/mod_STC.md`, "Cluster identification" > "Output").
    subroutine observable_alloc(vectors, n_dimensions, n_vectors, is_member_mask, &
                                U, eigenvalues, mu, d, G, normal_error_value, tangent_scales_value, ierr)

        integer(int32), intent(in)  :: n_dimensions
            !! Ambient dimension D
        integer(int32), intent(in)  :: n_vectors
            !! Number of input vectors N
        real(real64),   intent(in)  :: vectors(n_dimensions, n_vectors)
            !! Input data matrix
        logical,        intent(in)  :: is_member_mask(n_vectors)
            !! Ensemble membership; must have at least 2 members
        real(real64),   intent(out) :: U(n_dimensions, n_dimensions)
            !! Tangent+normal basis, zero-padded beyond rank
        real(real64),   intent(out) :: eigenvalues(n_dimensions)
            !! Covariance eigenvalues, descending, zero-padded beyond rank
        real(real64),   intent(out) :: mu(n_dimensions)
            !! Ensemble center
        integer(int32), intent(out) :: d
            !! Estimated intrinsic (tangent) dimension
        real(real64),   intent(out) :: G
            !! Spectral gap at d
        real(real64),   intent(out) :: normal_error_value
            !! Mean squared residual off the tangent subspace, see `normal_error`
        real(real64),   intent(out) :: tangent_scales_value(n_dimensions)
            !! Extent along each tangent direction, zero-padded beyond d, see `tangent_scales`
        integer(int32), intent(out) :: ierr
            !! Error code: 0 = success

        integer(int32) :: n_members, rank, lwork
        real(real64),   allocatable :: tmp_y(:,:), tmp_s(:), tmp_u_econ(:,:), tmp_vt_econ(:,:), tmp_work(:)
        integer(int32), allocatable :: tmp_iwork(:)

        call set_ok(ierr)
        call validate_dimension_size(n_dimensions, ierr)
        call validate_dimension_size(n_vectors, ierr)
        call validate_in_range_int(n_dimensions, ierr, min=2_int32)
        if (is_err(ierr)) return

        n_members = count(is_member_mask)
        call validate_in_range_int(n_members, ierr, min=2_int32, max=n_vectors)
        if (is_err(ierr)) return

        rank  = min(n_dimensions, n_members)
        lwork = 4 * rank * rank + 7 * rank
            ! LAPACK dgesdd's documented minimum workspace for JOBZ='S':
            ! LWORK >= 4*min(M,N)**2 + 7*min(M,N) (see `man dgesdd`).

        M_ALLOCATE(tmp_y(n_dimensions, n_members))
        M_ALLOCATE(tmp_s(rank))
        M_ALLOCATE(tmp_u_econ(n_dimensions, rank))
        M_ALLOCATE(tmp_vt_econ(rank, n_members))
        M_ALLOCATE(tmp_work(lwork))
        M_ALLOCATE(tmp_iwork(8 * rank))

        call observable(vectors, n_dimensions, n_vectors, is_member_mask, lwork, &
                       tmp_y, tmp_s, tmp_u_econ, tmp_vt_econ, tmp_work, tmp_iwork, &
                       U, eigenvalues, mu, d, G, normal_error_value, tangent_scales_value, ierr)

    end subroutine observable_alloc

    !> Validated Entry Point for SKG `observable`. Takes pre-allocated
    !| workspace (see `observable_alloc`). `n_members` and `rank` are not
    !| separate dummies here -- they're re-derived from `is_member_mask` via
    !| `count`/`min` directly in the workspace arrays' bound expressions
    !| below, so there is exactly one source of truth for them and nothing
    !| to fall out of sync. `lwork` alone is passed explicitly, mirroring
    !| `dsyev`'s own `lwork` dummy in `src/lomanle.F90`, since it comes from
    !| LAPACK's minimum-workspace formula rather than being trivially
    !| re-derivable inline.
    subroutine observable(vectors, n_dimensions, n_vectors, is_member_mask, lwork, &
                          tmp_y, tmp_s, tmp_u_econ, tmp_vt_econ, tmp_work, tmp_iwork, &
                          U, eigenvalues, mu, d, G, normal_error_value, tangent_scales_value, ierr)

        integer(int32), intent(in)    :: n_dimensions
            !! Ambient dimension D
        integer(int32), intent(in)    :: n_vectors
            !! Number of input vectors N
        real(real64),   intent(in)    :: vectors(n_dimensions, n_vectors)
            !! Input data matrix
        logical,        intent(in)    :: is_member_mask(n_vectors)
            !! Ensemble membership; must have at least 2 members
        integer(int32), intent(in)    :: lwork
            !! Size of `tmp_work`, see LAPACK `dgesdd`'s minimum-workspace formula
        real(real64),   intent(inout) :: tmp_y(n_dimensions, count(is_member_mask))
            !! Workspace: centered member matrix
        real(real64),   intent(inout) :: tmp_s(min(n_dimensions, count(is_member_mask)))
            !! Workspace: singular values
        real(real64),   intent(inout) :: tmp_u_econ(n_dimensions, min(n_dimensions, count(is_member_mask)))
            !! Workspace: economy-mode left singular vectors
        real(real64),   intent(inout) :: tmp_vt_econ(min(n_dimensions, count(is_member_mask)), count(is_member_mask))
            !! Workspace: economy-mode right singular vectors, transposed (unused beyond the SVD call)
        real(real64),   intent(inout) :: tmp_work(lwork)
            !! Workspace: LAPACK `dgesdd` scratch
        integer(int32), intent(inout) :: tmp_iwork(8 * min(n_dimensions, count(is_member_mask)))
            !! Workspace: LAPACK `dgesdd` integer scratch
        real(real64),   intent(out)   :: U(n_dimensions, n_dimensions)
            !! Tangent+normal basis, zero-padded beyond rank
        real(real64),   intent(out)   :: eigenvalues(n_dimensions)
            !! Covariance eigenvalues, descending, zero-padded beyond rank
        real(real64),   intent(out)   :: mu(n_dimensions)
            !! Ensemble center
        integer(int32), intent(out)   :: d
            !! Estimated intrinsic (tangent) dimension
        real(real64),   intent(out)   :: G
            !! Spectral gap at d
        real(real64),   intent(out)   :: normal_error_value
            !! Mean squared residual off the tangent subspace
        real(real64),   intent(out)   :: tangent_scales_value(n_dimensions)
            !! Extent along each tangent direction, zero-padded beyond d
        integer(int32), intent(out)   :: ierr
            !! Error code: 0 = success

        call set_ok(ierr)
        call validate_dimension_size(n_dimensions, ierr)
        call validate_dimension_size(n_vectors, ierr)
        call validate_in_range_int(n_dimensions, ierr, min=2_int32)
        if (is_err(ierr)) return

        call validate_in_range_int(count(is_member_mask), ierr, min=2_int32, max=n_vectors)
        if (is_err(ierr)) return

        call observable_helper(vectors, n_dimensions, n_vectors, is_member_mask, lwork, &
                               tmp_y, tmp_s, tmp_u_econ, tmp_vt_econ, tmp_work, tmp_iwork, &
                               U, eigenvalues, mu, d, G, normal_error_value, tangent_scales_value, ierr)

    end subroutine observable

    !> Core Implementation for SKG `observable`: economy-mode SVD (LAPACK
    !| `dgesdd`) of the ensemble's centered member vectors, spectral-gap
    !| intrinsic-dimension estimation, and the resulting `normal_error`/
    !| `tangent_scales` derived quantities -- see `misc/mod_STC.md`'s SKG
    !| `observable` for the full derivation. `pure` because `dgesdd` is
    !| declared `pure` above (documentation of its own determinism/
    !| thread-safety, matching `dsyev`'s precedent in `src/lomanle.F90`), and
    !| every other call here (`normal_error_helper`, `tangent_scales_helper`)
    !| already is.
    pure subroutine observable_helper(vectors, n_dimensions, n_vectors, is_member_mask, lwork, &
                                      tmp_y, tmp_s, tmp_u_econ, tmp_vt_econ, tmp_work, tmp_iwork, &
                                      U, eigenvalues, mu, d, G, normal_error_value, tangent_scales_value, ierr)

        integer(int32), intent(in)    :: n_dimensions
            !! Ambient dimension D
        integer(int32), intent(in)    :: n_vectors
            !! Number of input vectors N
        real(real64),   intent(in)    :: vectors(n_dimensions, n_vectors)
            !! Input data matrix
        logical,        intent(in)    :: is_member_mask(n_vectors)
            !! Ensemble membership
        integer(int32), intent(in)    :: lwork
            !! Size of `tmp_work`
        real(real64),   intent(inout) :: tmp_y(n_dimensions, count(is_member_mask))
            !! Workspace: centered member matrix
        real(real64),   intent(inout) :: tmp_s(min(n_dimensions, count(is_member_mask)))
            !! Workspace: singular values
        real(real64),   intent(inout) :: tmp_u_econ(n_dimensions, min(n_dimensions, count(is_member_mask)))
            !! Workspace: economy-mode left singular vectors
        real(real64),   intent(inout) :: tmp_vt_econ(min(n_dimensions, count(is_member_mask)), count(is_member_mask))
            !! Workspace: economy-mode right singular vectors, transposed
        real(real64),   intent(inout) :: tmp_work(lwork)
            !! Workspace: LAPACK `dgesdd` scratch
        integer(int32), intent(inout) :: tmp_iwork(8 * min(n_dimensions, count(is_member_mask)))
            !! Workspace: LAPACK `dgesdd` integer scratch
        real(real64),   intent(out)   :: U(n_dimensions, n_dimensions)
            !! Tangent+normal basis, zero-padded beyond rank
        real(real64),   intent(out)   :: eigenvalues(n_dimensions)
            !! Covariance eigenvalues, descending, zero-padded beyond rank
        real(real64),   intent(out)   :: mu(n_dimensions)
            !! Ensemble center
        integer(int32), intent(out)   :: d
            !! Estimated intrinsic (tangent) dimension
        real(real64),   intent(out)   :: G
            !! Spectral gap at d
        real(real64),   intent(out)   :: normal_error_value
            !! Mean squared residual off the tangent subspace
        real(real64),   intent(out)   :: tangent_scales_value(n_dimensions)
            !! Extent along each tangent direction, zero-padded beyond d
        integer(int32), intent(out)   :: ierr
            !! Error code: 0 = success

        integer(int32) :: n_members, rank, i_vec, i_member, j, r, info, d_local
        real(real64)   :: gap, best_gap

        n_members = count(is_member_mask)
        rank      = min(n_dimensions, n_members)

        i_member = 0
        do i_vec = 1, n_vectors
            if (.not. is_member_mask(i_vec)) cycle
            i_member = i_member + 1
            tmp_y(:, i_member) = vectors(:, i_vec)
        end do

        mu    = sum(tmp_y, dim=2) / real(n_members, real64)
        tmp_y = tmp_y - spread(mu, dim=2, ncopies=n_members)

        call dgesdd('S', n_dimensions, n_members, tmp_y, n_dimensions, tmp_s, &
                   tmp_u_econ, n_dimensions, tmp_vt_econ, rank, tmp_work, lwork, tmp_iwork, info)
        if (info /= 0) then
            call set_err(ierr, ERR_INTERNAL)
            return
        end if
        call set_ok(ierr)

        U = 0.0_real64
        U(:, 1:rank) = tmp_u_econ(:, 1:rank)

        eigenvalues = 0.0_real64
        do j = 1, rank
            eigenvalues(j) = tmp_s(j)**2 / real(n_members - 1, real64)
        end do

        ! Spectral gap G(r) = lambda_r / (lambda_{r+1} + eps), r = 1..D-1;
        ! d = argmax_r G(r). eigenvalues is zero-padded beyond `rank`, so any
        ! r >= rank+1 contributes G(r) = 0/eps = 0 -- never a spurious
        ! maximum -- while r = rank itself, whenever rank < n_dimensions,
        ! divides a genuine nonzero eigenvalue by (0 + eps) and correctly
        ! wins whenever the ensemble's data does not span the full ambient
        ! space (e.g. still few members, early in growth).
        best_gap = -1.0_real64
        d_local  = 1
        do r = 1, n_dimensions - 1
            gap = eigenvalues(r) / (eigenvalues(r + 1) + epsilon(1.0_real64))
            if (gap > best_gap) then
                best_gap = gap
                d_local  = r
            end if
        end do
        d = d_local
        G = best_gap

        call normal_error_helper(d_local, eigenvalues, n_dimensions, normal_error_value)

        tangent_scales_value = 0.0_real64
        block
            real(real64) :: tangent_scales_local(d_local)
            call tangent_scales_helper(d_local, eigenvalues, n_dimensions, tangent_scales_local)
            tangent_scales_value(1:d_local) = tangent_scales_local
        end block

    end subroutine observable_helper

    !> Allocating Wrapper for SKG `accept_ensemble` (`misc/mod_STC.md`):
    !| whether a grown ensemble at t+1 is still compatible with its own
    !| state at t, judged by three independent criteria (principal angle
    !| between the two tangent bases, change in intrinsic dimension, and
    !| relative change in spectral gap). This compares the SAME ensemble
    !| across one growth step -- not two different ensembles/anchors at a
    !| possible junction -- so, unlike `misc/STC_for_LoManLe.md` section 4's
    !| explicit "angle never gates a junction" rule, a principal-angle
    !| mismatch here legitimately contributes to rejection.
    subroutine accept_ensemble_alloc(n_dimensions, U_t, d_t, G_t, U_tp1, d_tp1, G_tp1, &
                                     alpha_max, d_max, G_max, is_accepted, ierr)

        integer(int32), intent(in)  :: n_dimensions
            !! Ambient dimension D
        real(real64),   intent(in)  :: U_t(n_dimensions, n_dimensions)
            !! Ensemble's tangent+normal basis at t, see `observable`
        integer(int32), intent(in)  :: d_t
            !! Ensemble's intrinsic dimension at t
        real(real64),   intent(in)  :: G_t
            !! Ensemble's spectral gap at t; must be > 0
        real(real64),   intent(in)  :: U_tp1(n_dimensions, n_dimensions)
            !! Ensemble's tangent+normal basis at t+1
        integer(int32), intent(in)  :: d_tp1
            !! Ensemble's intrinsic dimension at t+1
        real(real64),   intent(in)  :: G_tp1
            !! Ensemble's spectral gap at t+1; must be > 0
        real(real64),   intent(in)  :: alpha_max
            !! Maximum tolerated principal angle (radians), 0 <= alpha_max <= pi/2
        integer(int32), intent(in)  :: d_max
            !! Maximum tolerated change in intrinsic dimension, >= 0
        real(real64),   intent(in)  :: G_max
            !! Maximum tolerated |log(G_tp1/G_t)|, >= 0
        logical,        intent(out) :: is_accepted
            !! .true. if all three acceptance criteria are satisfied
        integer(int32), intent(out) :: ierr
            !! Error code: 0 = success

        integer(int32) :: lwork
        real(real64),   allocatable :: tmp_m(:,:), tmp_s(:), tmp_work(:)

        call set_ok(ierr)
        call validate_dimension_size(n_dimensions, ierr)
        call validate_in_range_int(d_t, ierr, min=0_int32, max=n_dimensions)
        call validate_in_range_int(d_tp1, ierr, min=0_int32, max=n_dimensions)
        call validate_in_range_real(alpha_max, ierr, min=0.0_real64, max=2.0_real64 * atan(1.0_real64))
        call validate_in_range_int(d_max, ierr, min=0_int32)
        call validate_in_range_real(G_max, ierr, min=0.0_real64)
        if (is_err(ierr)) return
        if (G_t <= 0.0_real64 .or. G_tp1 <= 0.0_real64) then
            call set_err(ierr, ERR_INVALID_INPUT)
            return
        end if

        lwork = max(1, 5 * min(d_t, d_tp1))
            ! LAPACK dgesvd's documented minimum workspace for a square
            ! input (M=N=min(d_t,d_tp1)) with JOBU='N', JOBVT='N':
            ! LWORK >= max(1, 3*min(M,N)+max(M,N), 5*min(M,N)), which for a
            ! square matrix reduces to max(1, 5*min(M,N)) (see `man dgesvd`).

        M_ALLOCATE(tmp_m(min(d_t, d_tp1), min(d_t, d_tp1)))
        M_ALLOCATE(tmp_s(min(d_t, d_tp1)))
        M_ALLOCATE(tmp_work(lwork))

        call accept_ensemble(n_dimensions, U_t, d_t, G_t, U_tp1, d_tp1, G_tp1, &
                             alpha_max, d_max, G_max, lwork, tmp_m, tmp_s, tmp_work, is_accepted, ierr)

    end subroutine accept_ensemble_alloc

    !> Validated Entry Point for SKG `accept_ensemble`. Takes pre-allocated
    !| workspace (see `accept_ensemble_alloc`); its shapes are direct
    !| specification expressions of `d_t`/`d_tp1` (both required,
    !| non-optional dummies, so no default-resolution timing issue -- unlike
    !| e.g. `calc_ensemble_growth_radius`'s optional `k_min`), matching
    !| `grow_ensemble`'s and `observable`'s own explicit-shape-at-every-tier
    !| convention.
    subroutine accept_ensemble(n_dimensions, U_t, d_t, G_t, U_tp1, d_tp1, G_tp1, &
                               alpha_max, d_max, G_max, lwork, tmp_m, tmp_s, tmp_work, is_accepted, ierr)

        integer(int32), intent(in)    :: n_dimensions
            !! Ambient dimension D
        real(real64),   intent(in)    :: U_t(n_dimensions, n_dimensions)
            !! Ensemble's tangent+normal basis at t
        integer(int32), intent(in)    :: d_t
            !! Ensemble's intrinsic dimension at t
        real(real64),   intent(in)    :: G_t
            !! Ensemble's spectral gap at t; must be > 0
        real(real64),   intent(in)    :: U_tp1(n_dimensions, n_dimensions)
            !! Ensemble's tangent+normal basis at t+1
        integer(int32), intent(in)    :: d_tp1
            !! Ensemble's intrinsic dimension at t+1
        real(real64),   intent(in)    :: G_tp1
            !! Ensemble's spectral gap at t+1; must be > 0
        real(real64),   intent(in)    :: alpha_max
            !! Maximum tolerated principal angle (radians)
        integer(int32), intent(in)    :: d_max
            !! Maximum tolerated change in intrinsic dimension
        real(real64),   intent(in)    :: G_max
            !! Maximum tolerated |log(G_tp1/G_t)|
        integer(int32), intent(in)    :: lwork
            !! Size of `tmp_work`, see LAPACK `dgesvd`'s minimum-workspace formula
        real(real64),   intent(inout) :: tmp_m(min(d_t, d_tp1), min(d_t, d_tp1))
            !! Workspace: M = U_t(:,1:d)^T U_tp1(:,1:d)
        real(real64),   intent(inout) :: tmp_s(min(d_t, d_tp1))
            !! Workspace: singular values of M (= cos of the principal angles)
        real(real64),   intent(inout) :: tmp_work(lwork)
            !! Workspace: LAPACK `dgesvd` scratch
        logical,        intent(out)   :: is_accepted
            !! .true. if all three acceptance criteria are satisfied
        integer(int32), intent(out)   :: ierr
            !! Error code: 0 = success

        call set_ok(ierr)
        call validate_dimension_size(n_dimensions, ierr)
        call validate_in_range_int(d_t, ierr, min=0_int32, max=n_dimensions)
        call validate_in_range_int(d_tp1, ierr, min=0_int32, max=n_dimensions)
        call validate_in_range_real(alpha_max, ierr, min=0.0_real64, max=2.0_real64 * atan(1.0_real64))
        call validate_in_range_int(d_max, ierr, min=0_int32)
        call validate_in_range_real(G_max, ierr, min=0.0_real64)
        if (is_err(ierr)) return
        if (G_t <= 0.0_real64 .or. G_tp1 <= 0.0_real64) then
            call set_err(ierr, ERR_INVALID_INPUT)
            return
        end if

        call accept_ensemble_helper(n_dimensions, U_t, d_t, G_t, U_tp1, d_tp1, G_tp1, &
                                    alpha_max, d_max, G_max, lwork, tmp_m, tmp_s, tmp_work, is_accepted, ierr)

    end subroutine accept_ensemble

    !> Core Implementation for SKG `accept_ensemble`. `pure` because
    !| `dgesvd` is declared `pure` above (same documentation-of-thread-
    !| safety rationale as `dgesdd` in `observable_helper`).
    !|
    !| The three criteria (`misc/mod_STC.md`):
    !| (1) principal angles between the d-dimensional tangent bases, via
    !|     `dgesvd` on M = U_t(:,1:d)^T U_tp1(:,1:d), whose singular values
    !|     are cos(alpha_i) directly -- but only when d_t == d_tp1: when the
    !|     estimated intrinsic dimension itself changed, the two tangent
    !|     bases don't share a common dimension to compare angles over at
    !|     all, so this criterion is skipped (no SVD is computed) and
    !|     treated as vacuously satisfied; criterion (2) below is what
    !|     actually judges whether that change in d is acceptable.
    !| (2) |d_tp1 - d_t| <= d_max.
    !| (3) |log(G_tp1/G_t)| <= G_max.
    pure subroutine accept_ensemble_helper(n_dimensions, U_t, d_t, G_t, U_tp1, d_tp1, G_tp1, &
                                           alpha_max, d_max, G_max, lwork, tmp_m, tmp_s, tmp_work, &
                                           is_accepted, ierr)

        integer(int32), intent(in)    :: n_dimensions
            !! Ambient dimension D
        real(real64),   intent(in)    :: U_t(n_dimensions, n_dimensions)
            !! Ensemble's tangent+normal basis at t
        integer(int32), intent(in)    :: d_t
            !! Ensemble's intrinsic dimension at t
        real(real64),   intent(in)    :: G_t
            !! Ensemble's spectral gap at t
        real(real64),   intent(in)    :: U_tp1(n_dimensions, n_dimensions)
            !! Ensemble's tangent+normal basis at t+1
        integer(int32), intent(in)    :: d_tp1
            !! Ensemble's intrinsic dimension at t+1
        real(real64),   intent(in)    :: G_tp1
            !! Ensemble's spectral gap at t+1
        real(real64),   intent(in)    :: alpha_max
            !! Maximum tolerated principal angle (radians)
        integer(int32), intent(in)    :: d_max
            !! Maximum tolerated change in intrinsic dimension
        real(real64),   intent(in)    :: G_max
            !! Maximum tolerated |log(G_tp1/G_t)|
        integer(int32), intent(in)    :: lwork
            !! Size of `tmp_work`
        real(real64),   intent(inout) :: tmp_m(min(d_t, d_tp1), min(d_t, d_tp1))
            !! Workspace: M = U_t(:,1:d)^T U_tp1(:,1:d)
        real(real64),   intent(inout) :: tmp_s(min(d_t, d_tp1))
            !! Workspace: singular values of M
        real(real64),   intent(inout) :: tmp_work(lwork)
            !! Workspace: LAPACK `dgesvd` scratch
        logical,        intent(out)   :: is_accepted
            !! .true. if all three acceptance criteria are satisfied
        integer(int32), intent(out)   :: ierr
            !! Error code: 0 = success

        real(real64)   :: u_dummy(1,1), vt_dummy(1,1), cos_alpha, alpha_i
        logical        :: angle_ok, d_diff_ok, g_ratio_ok
        integer(int32) :: d_common, i, info

        call set_ok(ierr)

        d_diff_ok  = abs(d_tp1 - d_t) <= d_max
        g_ratio_ok = abs(log(G_tp1 / G_t)) <= G_max

        angle_ok = .true.
        if (d_t == d_tp1) then
            d_common = d_t
            if (d_common > 0) then
                tmp_m = matmul(transpose(U_t(:, 1:d_common)), U_tp1(:, 1:d_common))
                call dgesvd('N', 'N', d_common, d_common, tmp_m, d_common, tmp_s, &
                           u_dummy, 1, vt_dummy, 1, tmp_work, lwork, info)
                if (info /= 0) then
                    call set_err(ierr, ERR_INTERNAL)
                    return
                end if
                do i = 1, d_common
                    cos_alpha = max(-1.0_real64, min(1.0_real64, tmp_s(i)))
                    alpha_i   = acos(cos_alpha)
                    if (alpha_i > alpha_max) angle_ok = .false.
                end do
            end if
        end if

        is_accepted = angle_ok .and. d_diff_ok .and. g_ratio_ok

    end subroutine accept_ensemble_helper

end module shape_truthful_clustering
