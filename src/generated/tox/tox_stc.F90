#include <src/macros.h>

!> Shape-truthful clustering with a global or per-seed growth radius.
!| Compute density labels, select seeds, grow their ensembles, and merge overlapping ensembles.
!| Both radius variants use the same density labels, seed ordering, growth rule and merge rule.
!|
!| Generated from [[tox_stc_impl(module)]]; do not edit -- regenerate instead.
module tox_stc
    use f42_safeguard
    use tox_stc_impl, only: MODE_GLOBAL_RADIUS, MODE_SEED_RADIUS, accept_ensemble_impl, calculate_density_radius_impl
    use tox_stc_impl, only: calculate_labels_as_density_impl, compute_ambient_density_stats_impl, compute_ensemble_observable_impl, grow_ensemble_impl
    use tox_stc_impl, only: identify_ensemble_seeds_impl, merge_ensembles_impl, obtain_ensembles_impl
    use, intrinsic :: iso_c_binding, only: c_bool
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use tox_errors, only: set_ok, is_err, ERR_ALLOC_FAIL, clear_err_arg_pos
    use tox_errors, only: set_err, validate_all_in_range_int, validate_all_in_range_real, validate_dimension_size
    use tox_errors, only: validate_in_range_int, validate_in_range_real
    M_IMPLICIT_NONE
    private

    public :: calculate_density_radius
    public :: calculate_density_radius_expert
    public :: calculate_labels_as_density
    public :: calculate_labels_as_density_expert
    public :: identify_ensemble_seeds
    public :: identify_ensemble_seeds_expert
    public :: compute_ambient_density_stats
    public :: compute_ambient_density_stats_expert
    public :: obtain_ensembles_global_radius
    public :: obtain_ensembles_global_radius_expert
    public :: obtain_ensembles_seed_radius
    public :: obtain_ensembles_seed_radius_expert
    public :: merge_ensembles
    public :: merge_ensembles_expert
    public :: grow_ensemble
    public :: grow_ensemble_expert
    public :: compute_ensemble_observable
    public :: accept_ensemble

contains

    !> summary: Validates its inputs, prepares what [[tox_stc_impl(module):calculate_density_radius_impl]] needs, then calls it. The entry point to reach for first; see [[tox_stc(module):calculate_density_radius_expert]] to prepare it yourself.
    pure subroutine calculate_density_radius(&
            vectors,&
            n_dimensions,&
            n_vectors,&
            mean_to_other_vecs_dist_quant,&
            radius,&
            ierr&
        )
        integer(int32), intent(in) :: n_dimensions
            !! Number of dimensions
        integer(int32), intent(in) :: n_vectors
            !! Number of vectors
        real(real64), dimension(n_dimensions, n_vectors), intent(in) :: vectors
            !! Input data matrix (n_dimensions x n_vectors)
        real(real64), intent(in), optional :: mean_to_other_vecs_dist_quant
            !! Quantile of the distances to the data centroid that becomes the density radius
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.15_real64`.
        real(real64), intent(out) :: radius
            !! The resulting density search radius
        integer(int32), intent(out) :: ierr
            !! Error code
        real(real64), dimension(:), allocatable :: tmp_mean_vec
        real(real64), dimension(:), allocatable :: tmp_distances
        integer(int32), dimension(:), allocatable :: tmp_perm

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_dimensions, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_vectors, ierr, arg_pos=3_int32)
        call validate_in_range_real(mean_to_other_vecs_dist_quant, ierr, arg_pos=4_int32, min=0.0_real64, max=1.0_real64)
        call validate_all_in_range_real(vectors, n_dimensions * n_vectors, ierr, arg_pos=1_int32)
        if (is_err(ierr)) return
#endif

        M_ALLOCATE(tmp_mean_vec(n_dimensions))
        M_ALLOCATE(tmp_distances(n_vectors))
        M_ALLOCATE(tmp_perm(n_vectors))

        call calculate_density_radius_impl(&
            vectors = vectors,&
            n_dimensions = n_dimensions,&
            n_vectors = n_vectors,&
            tmp_mean_vec = tmp_mean_vec,&
            tmp_distances = tmp_distances,&
            tmp_perm = tmp_perm,&
            mean_to_other_vecs_dist_quant = mean_to_other_vecs_dist_quant,&
            radius = radius,&
            ierr = ierr&
        )
        call clear_err_arg_pos(ierr)
    end subroutine calculate_density_radius

    !> summary: Validates its inputs, then calls [[tox_stc_impl(module):calculate_density_radius_impl]] with what you supply. The expert entry point: it allocates nothing and prepares nothing; [[tox_stc(module):calculate_density_radius]] does both.
    pure subroutine calculate_density_radius_expert(&
            vectors,&
            n_dimensions,&
            n_vectors,&
            tmp_mean_vec,&
            tmp_distances,&
            tmp_perm,&
            mean_to_other_vecs_dist_quant,&
            radius,&
            ierr&
        )
        integer(int32), intent(in) :: n_dimensions
            !! Number of dimensions
        integer(int32), intent(in) :: n_vectors
            !! Number of vectors
        real(real64), dimension(n_dimensions, n_vectors), intent(in) :: vectors
            !! Input data matrix (n_dimensions x n_vectors)
        real(real64), dimension(n_dimensions), intent(out) :: tmp_mean_vec
            !! Preallocated workspace for the mean vector
        real(real64), dimension(n_vectors), intent(out) :: tmp_distances
            !! Preallocated workspace for calculated distances
        integer(int32), dimension(n_vectors), intent(out) :: tmp_perm
            !! Preallocated workspace for permutation indices
        real(real64), intent(in), optional :: mean_to_other_vecs_dist_quant
            !! Quantile of the distances to the data centroid that becomes the density radius
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.15_real64`.
        real(real64), intent(out) :: radius
            !! The resulting density search radius
        integer(int32), intent(out) :: ierr
            !! Error code

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_dimensions, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_vectors, ierr, arg_pos=3_int32)
        call validate_in_range_real(mean_to_other_vecs_dist_quant, ierr, arg_pos=7_int32, min=0.0_real64, max=1.0_real64)
        call validate_all_in_range_real(vectors, n_dimensions * n_vectors, ierr, arg_pos=1_int32)
        if (is_err(ierr)) return
#endif

        call calculate_density_radius_impl(&
            vectors = vectors,&
            n_dimensions = n_dimensions,&
            n_vectors = n_vectors,&
            tmp_mean_vec = tmp_mean_vec,&
            tmp_distances = tmp_distances,&
            tmp_perm = tmp_perm,&
            mean_to_other_vecs_dist_quant = mean_to_other_vecs_dist_quant,&
            radius = radius,&
            ierr = ierr&
        )
        call clear_err_arg_pos(ierr)
    end subroutine calculate_density_radius_expert

    !> summary: Validates its inputs, prepares what [[tox_stc_impl(module):calculate_labels_as_density_impl]] needs, then calls it. The entry point to reach for first; see [[tox_stc(module):calculate_labels_as_density_expert]] to prepare it yourself.
    pure subroutine calculate_labels_as_density(&
            vectors,&
            n_dimensions,&
            n_vectors,&
            r,&
            dimension_order,&
            kd_indices,&
            n_tiles,&
            label_densities,&
            ierr&
        )
        integer(int32), intent(in) :: n_dimensions
            !! Number of dimensions
        integer(int32), intent(in) :: n_vectors
            !! Number of vectors
        real(real64), dimension(n_dimensions, n_vectors), intent(in) :: vectors
            !! Input data matrix (n_dimensions x n_vectors)
        real(real64), intent(in) :: r
            !! Label-sphere search radius
            !! The minimum valid value is `0.0_real64`.
        integer(int32), dimension(n_dimensions), intent(in) :: dimension_order
            !! Sequence array tracking tree split axes by variance
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_dimensions`.
        integer(int32), dimension(n_vectors), intent(in) :: kd_indices
            !! KD-tree index sequence array computed via [[f42_kd_tree(module):build_kd_index(subroutine)]].
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_vectors`.
        integer(int32), intent(in) :: n_tiles
            !! Number of concurrent vector tiles, at least 1. Capped at the number of points.
            !! The minimum valid value is `1_int32`.
        real(real64), dimension(n_vectors), intent(out) :: label_densities
            !! Output array storing generated density scalars
        integer(int32), intent(out) :: ierr
            !! Error code
        integer(int32), dimension(:, :, :), allocatable :: tmp_stack

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_dimensions, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_vectors, ierr, arg_pos=3_int32)
        call validate_in_range_real(r, ierr, arg_pos=4_int32, min=0.0_real64)
        call validate_in_range_int(n_tiles, ierr, arg_pos=7_int32, min=1_int32)
        call validate_all_in_range_real(vectors, n_dimensions * n_vectors, ierr, arg_pos=1_int32)
        call validate_all_in_range_int(dimension_order, n_dimensions, ierr, arg_pos=5_int32, min=1_int32, max=n_dimensions)
        call validate_all_in_range_int(kd_indices, n_vectors, ierr, arg_pos=6_int32, min=1_int32, max=n_vectors)
        if (is_err(ierr)) return
#endif

        M_ALLOCATE(tmp_stack(3, 64, min(n_tiles,n_vectors)))

        call calculate_labels_as_density_impl(&
            vectors = vectors,&
            n_dimensions = n_dimensions,&
            n_vectors = n_vectors,&
            r = r,&
            dimension_order = dimension_order,&
            kd_indices = kd_indices,&
            n_tiles = n_tiles,&
            tmp_stack = tmp_stack,&
            label_densities = label_densities,&
            ierr = ierr&
        )
        call clear_err_arg_pos(ierr)
    end subroutine calculate_labels_as_density

    !> summary: Validates its inputs, then calls [[tox_stc_impl(module):calculate_labels_as_density_impl]] with what you supply. The expert entry point: it allocates nothing and prepares nothing; [[tox_stc(module):calculate_labels_as_density]] does both.
    pure subroutine calculate_labels_as_density_expert(&
            vectors,&
            n_dimensions,&
            n_vectors,&
            r,&
            dimension_order,&
            kd_indices,&
            n_tiles,&
            tmp_stack,&
            label_densities,&
            ierr&
        )
        integer(int32), intent(in) :: n_dimensions
            !! Number of dimensions
        integer(int32), intent(in) :: n_vectors
            !! Number of vectors
        integer(int32), intent(in) :: n_tiles
            !! Number of concurrent vector tiles, at least 1. Capped at the number of points.
            !! The minimum valid value is `1_int32`.
        real(real64), dimension(n_dimensions, n_vectors), intent(in) :: vectors
            !! Input data matrix (n_dimensions x n_vectors)
        real(real64), intent(in) :: r
            !! Label-sphere search radius
            !! The minimum valid value is `0.0_real64`.
        integer(int32), dimension(n_dimensions), intent(in) :: dimension_order
            !! Sequence array tracking tree split axes by variance
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_dimensions`.
        integer(int32), dimension(n_vectors), intent(in) :: kd_indices
            !! KD-tree index sequence array computed via [[f42_kd_tree(module):build_kd_index(subroutine)]].
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_vectors`.
        integer(int32), dimension(3, 64, min(n_tiles,n_vectors)), intent(inout) :: tmp_stack
            !! Preallocated traversal stack per tile, reused serially across its vectors
        real(real64), dimension(n_vectors), intent(out) :: label_densities
            !! Output array storing generated density scalars
        integer(int32), intent(out) :: ierr
            !! Error code

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_dimensions, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_vectors, ierr, arg_pos=3_int32)
        call validate_in_range_real(r, ierr, arg_pos=4_int32, min=0.0_real64)
        call validate_in_range_int(n_tiles, ierr, arg_pos=7_int32, min=1_int32)
        call validate_all_in_range_real(vectors, n_dimensions * n_vectors, ierr, arg_pos=1_int32)
        call validate_all_in_range_int(dimension_order, n_dimensions, ierr, arg_pos=5_int32, min=1_int32, max=n_dimensions)
        call validate_all_in_range_int(kd_indices, n_vectors, ierr, arg_pos=6_int32, min=1_int32, max=n_vectors)
        if (is_err(ierr)) return
#endif

        call calculate_labels_as_density_impl(&
            vectors = vectors,&
            n_dimensions = n_dimensions,&
            n_vectors = n_vectors,&
            r = r,&
            dimension_order = dimension_order,&
            kd_indices = kd_indices,&
            n_tiles = n_tiles,&
            tmp_stack = tmp_stack,&
            label_densities = label_densities,&
            ierr = ierr&
        )
        call clear_err_arg_pos(ierr)
    end subroutine calculate_labels_as_density_expert

    !> summary: Validates its inputs, prepares what [[tox_stc_impl(module):identify_ensemble_seeds_impl]] needs, then calls it. The entry point to reach for first; see [[tox_stc(module):identify_ensemble_seeds_expert]] to prepare it yourself.
    pure subroutine identify_ensemble_seeds(&
            vectors,&
            n_dimensions,&
            n_vectors,&
            density_labels,&
            dimension_order,&
            kd_indices,&
            k_seeding,&
            sorted_perm,&
            n_seeds,&
            seed_mask,&
            seed_radii,&
            coverage_quant,&
            ierr&
        )
        integer(int32), intent(in) :: n_dimensions
            !! Number of dimensions
        integer(int32), intent(in) :: n_vectors
            !! Number of vectors
            !! The minimum valid value is `2_int32`.
        real(real64), dimension(n_dimensions, n_vectors), intent(in) :: vectors
            !! Input data matrix (n_dimensions x n_vectors)
        real(real64), dimension(n_vectors), intent(in) :: density_labels
            !! Precalculated density labels for all ambient vectors
            !! The minimum valid value is `1.0_real64`.
        integer(int32), dimension(n_dimensions), intent(in) :: dimension_order
            !! Dimension split order array tracking the K-D tree structure
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_dimensions`.
        integer(int32), dimension(n_vectors), intent(in) :: kd_indices
            !! K-D tree index sequence array
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_vectors`.
        integer(int32), intent(in) :: k_seeding
            !! Number of nearest neighbors used to estimate each seed's local coverage radius
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_vectors - 1_int32`.
        integer(int32), dimension(n_vectors), intent(out) :: sorted_perm
            !! Vector indices sorted by density descending
        integer(int32), intent(out) :: n_seeds
            !! Number of selected seed vectors
        logical(c_bool), dimension(n_vectors), intent(out) :: seed_mask
            !! Logical mask identifying selected seed vectors
        real(real64), dimension(n_vectors), intent(out) :: seed_radii
            !! Per-seed growth radius, stored at each seed's own vector index; zero at non-seeds
        real(real64), intent(in), optional :: coverage_quant
            !! Quantile of each seed's k-nearest-neighbor distances used as its coverage and growth radius
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.5_real64`.
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.
        integer(int32), dimension(:), allocatable :: tmp_perm
        real(real64), dimension(:), allocatable :: tmp_distances
        integer(int32), dimension(:, :), allocatable :: tmp_stack
        logical(c_bool), dimension(:), allocatable :: tmp_visited_mask
        logical(c_bool), dimension(:), allocatable :: tmp_newly_covered_mask

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_dimensions, ierr, arg_pos=2_int32)
        call validate_in_range_int(n_vectors, ierr, arg_pos=3_int32, min=2_int32)
        call validate_in_range_int(k_seeding, ierr, arg_pos=7_int32, min=1_int32, max=n_vectors - 1_int32)
        call validate_in_range_real(coverage_quant, ierr, arg_pos=12_int32, min=0.0_real64, max=1.0_real64)
        call validate_all_in_range_real(vectors, n_dimensions * n_vectors, ierr, arg_pos=1_int32)
        call validate_all_in_range_real(density_labels, n_vectors, ierr, arg_pos=4_int32, min=1.0_real64)
        call validate_all_in_range_int(dimension_order, n_dimensions, ierr, arg_pos=5_int32, min=1_int32, max=n_dimensions)
        call validate_all_in_range_int(kd_indices, n_vectors, ierr, arg_pos=6_int32, min=1_int32, max=n_vectors)
        if (is_err(ierr)) return
#endif

        M_ALLOCATE(tmp_perm(n_vectors))
        M_ALLOCATE(tmp_distances(n_vectors))
        M_ALLOCATE(tmp_stack(3, 64))
        M_ALLOCATE(tmp_visited_mask(n_vectors))
        M_ALLOCATE(tmp_newly_covered_mask(n_vectors))

        call identify_ensemble_seeds_impl(&
            vectors = vectors,&
            n_dimensions = n_dimensions,&
            n_vectors = n_vectors,&
            density_labels = density_labels,&
            dimension_order = dimension_order,&
            kd_indices = kd_indices,&
            k_seeding = k_seeding,&
            tmp_perm = tmp_perm,&
            tmp_distances = tmp_distances,&
            tmp_stack = tmp_stack,&
            tmp_visited_mask = tmp_visited_mask,&
            tmp_newly_covered_mask = tmp_newly_covered_mask,&
            sorted_perm = sorted_perm,&
            n_seeds = n_seeds,&
            seed_mask = seed_mask,&
            seed_radii = seed_radii,&
            coverage_quant = coverage_quant&
        )
    end subroutine identify_ensemble_seeds

    !> summary: Validates its inputs, then calls [[tox_stc_impl(module):identify_ensemble_seeds_impl]] with what you supply. The expert entry point: it allocates nothing and prepares nothing; [[tox_stc(module):identify_ensemble_seeds]] does both.
    pure subroutine identify_ensemble_seeds_expert(&
            vectors,&
            n_dimensions,&
            n_vectors,&
            density_labels,&
            dimension_order,&
            kd_indices,&
            k_seeding,&
            tmp_perm,&
            tmp_distances,&
            tmp_stack,&
            tmp_visited_mask,&
            tmp_newly_covered_mask,&
            sorted_perm,&
            n_seeds,&
            seed_mask,&
            seed_radii,&
            coverage_quant,&
            ierr&
        )
        integer(int32), intent(in) :: n_dimensions
            !! Number of dimensions
        integer(int32), intent(in) :: n_vectors
            !! Number of vectors
            !! The minimum valid value is `2_int32`.
        real(real64), dimension(n_dimensions, n_vectors), intent(in) :: vectors
            !! Input data matrix (n_dimensions x n_vectors)
        real(real64), dimension(n_vectors), intent(in) :: density_labels
            !! Precalculated density labels for all ambient vectors
            !! The minimum valid value is `1.0_real64`.
        integer(int32), dimension(n_dimensions), intent(in) :: dimension_order
            !! Dimension split order array tracking the K-D tree structure
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_dimensions`.
        integer(int32), dimension(n_vectors), intent(in) :: kd_indices
            !! K-D tree index sequence array
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_vectors`.
        integer(int32), intent(in) :: k_seeding
            !! Number of nearest neighbors used to estimate each seed's local coverage radius
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_vectors - 1_int32`.
        integer(int32), dimension(n_vectors), intent(out) :: tmp_perm
            !! Preallocated permutation workspace
        real(real64), dimension(n_vectors), intent(out) :: tmp_distances
            !! Preallocated candidate-to-vector distance workspace
        integer(int32), dimension(3, 64), intent(inout) :: tmp_stack
            !! Preallocated K-D tree traversal stack
        logical(c_bool), dimension(n_vectors), intent(out) :: tmp_visited_mask
            !! Preallocated workspace tracking vectors already covered during seeding
        logical(c_bool), dimension(n_vectors), intent(out) :: tmp_newly_covered_mask
            !! Preallocated workspace for the current seed's coverage query
        integer(int32), dimension(n_vectors), intent(out) :: sorted_perm
            !! Vector indices sorted by density descending
        integer(int32), intent(out) :: n_seeds
            !! Number of selected seed vectors
        logical(c_bool), dimension(n_vectors), intent(out) :: seed_mask
            !! Logical mask identifying selected seed vectors
        real(real64), dimension(n_vectors), intent(out) :: seed_radii
            !! Per-seed growth radius, stored at each seed's own vector index; zero at non-seeds
        real(real64), intent(in), optional :: coverage_quant
            !! Quantile of each seed's k-nearest-neighbor distances used as its coverage and growth radius
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.5_real64`.
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_dimensions, ierr, arg_pos=2_int32)
        call validate_in_range_int(n_vectors, ierr, arg_pos=3_int32, min=2_int32)
        call validate_in_range_int(k_seeding, ierr, arg_pos=7_int32, min=1_int32, max=n_vectors - 1_int32)
        call validate_in_range_real(coverage_quant, ierr, arg_pos=17_int32, min=0.0_real64, max=1.0_real64)
        call validate_all_in_range_real(vectors, n_dimensions * n_vectors, ierr, arg_pos=1_int32)
        call validate_all_in_range_real(density_labels, n_vectors, ierr, arg_pos=4_int32, min=1.0_real64)
        call validate_all_in_range_int(dimension_order, n_dimensions, ierr, arg_pos=5_int32, min=1_int32, max=n_dimensions)
        call validate_all_in_range_int(kd_indices, n_vectors, ierr, arg_pos=6_int32, min=1_int32, max=n_vectors)
        if (is_err(ierr)) return
#endif

        call identify_ensemble_seeds_impl(&
            vectors = vectors,&
            n_dimensions = n_dimensions,&
            n_vectors = n_vectors,&
            density_labels = density_labels,&
            dimension_order = dimension_order,&
            kd_indices = kd_indices,&
            k_seeding = k_seeding,&
            tmp_perm = tmp_perm,&
            tmp_distances = tmp_distances,&
            tmp_stack = tmp_stack,&
            tmp_visited_mask = tmp_visited_mask,&
            tmp_newly_covered_mask = tmp_newly_covered_mask,&
            sorted_perm = sorted_perm,&
            n_seeds = n_seeds,&
            seed_mask = seed_mask,&
            seed_radii = seed_radii,&
            coverage_quant = coverage_quant&
        )
    end subroutine identify_ensemble_seeds_expert

    !> summary: Validates its inputs, prepares what [[tox_stc_impl(module):compute_ambient_density_stats_impl]] needs, then calls it. The entry point to reach for first; see [[tox_stc(module):compute_ambient_density_stats_expert]] to prepare it yourself.
    pure subroutine compute_ambient_density_stats(&
            density_labels,&
            n_vectors,&
            median_ambient,&
            mad_ambient,&
            ierr&
        )
        integer(int32), intent(in) :: n_vectors
            !! Number of vectors
        real(real64), dimension(n_vectors), intent(in) :: density_labels
            !! Density labels for all vectors, from `calculate_labels_as_density`.
            !! The minimum valid value is `1.0_real64`.
        real(real64), intent(out) :: median_ambient
            !! Ambient median density median_{i=1..N}(rho_i)
        real(real64), intent(out) :: mad_ambient
            !! Ambient median absolute deviation median_{i=1..N}(|rho_i - median_ambient|)
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.
        integer(int32), dimension(:), allocatable :: tmp_perm
        real(real64), dimension(:), allocatable :: tmp_abs_diff

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_vectors, ierr, arg_pos=2_int32)
        call validate_all_in_range_real(density_labels, n_vectors, ierr, arg_pos=1_int32, min=1.0_real64)
        if (is_err(ierr)) return
#endif

        M_ALLOCATE(tmp_perm(n_vectors))
        M_ALLOCATE(tmp_abs_diff(n_vectors))

        call compute_ambient_density_stats_impl(&
            density_labels = density_labels,&
            n_vectors = n_vectors,&
            tmp_perm = tmp_perm,&
            tmp_abs_diff = tmp_abs_diff,&
            median_ambient = median_ambient,&
            mad_ambient = mad_ambient&
        )
    end subroutine compute_ambient_density_stats

    !> summary: Validates its inputs, then calls [[tox_stc_impl(module):compute_ambient_density_stats_impl]] with what you supply. The expert entry point: it allocates nothing and prepares nothing; [[tox_stc(module):compute_ambient_density_stats]] does both.
    pure subroutine compute_ambient_density_stats_expert(&
            density_labels,&
            n_vectors,&
            tmp_perm,&
            tmp_abs_diff,&
            median_ambient,&
            mad_ambient,&
            ierr&
        )
        integer(int32), intent(in) :: n_vectors
            !! Number of vectors
        real(real64), dimension(n_vectors), intent(in) :: density_labels
            !! Density labels for all vectors, from `calculate_labels_as_density`.
            !! The minimum valid value is `1.0_real64`.
        integer(int32), dimension(n_vectors), intent(inout) :: tmp_perm
            !! Preallocated workspace array for sorting and percentiles
        real(real64), dimension(n_vectors), intent(out) :: tmp_abs_diff
            !! Preallocated workspace array for the absolute deviations
        real(real64), intent(out) :: median_ambient
            !! Ambient median density median_{i=1..N}(rho_i)
        real(real64), intent(out) :: mad_ambient
            !! Ambient median absolute deviation median_{i=1..N}(|rho_i - median_ambient|)
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_vectors, ierr, arg_pos=2_int32)
        call validate_all_in_range_real(density_labels, n_vectors, ierr, arg_pos=1_int32, min=1.0_real64)
        if (is_err(ierr)) return
#endif

        call compute_ambient_density_stats_impl(&
            density_labels = density_labels,&
            n_vectors = n_vectors,&
            tmp_perm = tmp_perm,&
            tmp_abs_diff = tmp_abs_diff,&
            median_ambient = median_ambient,&
            mad_ambient = mad_ambient&
        )
    end subroutine compute_ambient_density_stats_expert

    !> summary: Validates its inputs, prepares what [[tox_stc_impl(module):obtain_ensembles_impl]] needs, then calls it. The entry point to reach for first; see [[tox_stc(module):obtain_ensembles_global_radius_expert]] to prepare it yourself.
    pure subroutine obtain_ensembles_global_radius(&
            vectors,&
            n_dimensions,&
            n_vectors,&
            dimension_order,&
            kd_indices,&
            density_labels,&
            seed_indices,&
            n_seeds,&
            r,&
            alpha_mad,&
            alpha_accept,&
            t_observables,&
            n_tiles,&
            ensemble_matrix,&
            stop_reasons,&
            mad_ambient,&
            ierr&
        )
        integer(int32), intent(in) :: n_dimensions
            !! Number of dimensions
        integer(int32), intent(in) :: n_vectors
            !! Number of vectors
        integer(int32), intent(in) :: n_seeds
            !! Number of seed vectors, from `identify_ensemble_seeds`.
            !! The minimum valid value is `0_int32`.
            !! The maximum valid value is `n_vectors`.
        real(real64), dimension(n_dimensions, n_vectors), intent(in) :: vectors
            !! Input data matrix (n_dimensions x n_vectors)
        integer(int32), dimension(n_dimensions), intent(in) :: dimension_order
            !! Sequence array tracking split dimensions for KD-tree
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_dimensions`.
        integer(int32), dimension(n_vectors), intent(in) :: kd_indices
            !! KD-tree index sequence array computed via [[f42_kd_tree(module):build_kd_index(subroutine)]].
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_vectors`.
        real(real64), dimension(n_vectors), intent(in) :: density_labels
            !! Density labels for all vectors, from `calculate_labels_as_density`.
            !! The minimum valid value is `1.0_real64`.
        integer(int32), dimension(n_seeds), intent(in) :: seed_indices
            !! Indices of the seed vectors to grow from, from `identify_ensemble_seeds`.
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_vectors`.
        real(real64), intent(in) :: r
            !! Uniform growth radius for every seed.
            !! The minimum valid value is `0.0_real64`.
        real(real64), intent(in), optional :: alpha_mad
            !! MAD multiplier for density compatibility threshold
            !! The minimum valid value is `0.0_real64`.
            !! The default value is `0.5_real64`.
        real(real64), intent(in), optional :: alpha_accept
            !! Acceptance log2 fold-change threshold
            !! The minimum valid value is `0.0_real64`.
            !! The default value is `0.5_real64`.
        integer(int32), intent(in) :: t_observables
            !! History window in iterations, at least two because acceptance compares two consecutive
            !! steps. Pass n_vectors to retain the full history.
            !! The minimum valid value is `2_int32`.
        integer(int32), intent(in) :: n_tiles
            !! Number of concurrent seed tiles, at least 1. Capped at the number of seeds when nonzero.
            !! The minimum valid value is `1_int32`.
        logical(c_bool), dimension(n_vectors, n_seeds), intent(out) :: ensemble_matrix
            !! Output matrix storing unmerged grown ensemble masks
        integer(int32), dimension(n_seeds), intent(out) :: stop_reasons
            !! Per-seed stop code. 1: fixed point, 2: rejection after acceptance,
            !! 3: first step rejected, 4: no candidates.
        real(real64), intent(out) :: mad_ambient
            !! Median absolute deviation of `density_labels`, which scales the `alpha_mad` density tolerance.
            !! Zero when there are no seeds; otherwise positive on success, because a zero value is an error.
        integer(int32), intent(out) :: ierr
            !! Invalid input when the median absolute deviation of `density_labels` is zero (more than
            !! half the labels are equal), which leaves no density tolerance for growth.
        integer(int32), dimension(:, :, :), allocatable :: tmp_stack
        logical(c_bool), dimension(:, :), allocatable :: tmp_vicinity_mask
        logical(c_bool), dimension(:, :), allocatable :: tmp_surface_mask
        integer(int32), dimension(:, :), allocatable :: tmp_perm
        real(real64), dimension(:, :), allocatable :: tmp_abs_diff
        real(real64), dimension(:, :, :), allocatable :: tmp_observables
        logical(c_bool), dimension(:, :), allocatable :: tmp_current_mask

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_dimensions, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_vectors, ierr, arg_pos=3_int32)
        call validate_in_range_int(n_seeds, ierr, arg_pos=8_int32, min=0_int32, max=n_vectors)
        call validate_in_range_real(r, ierr, arg_pos=9_int32, min=0.0_real64)
        call validate_in_range_real(alpha_mad, ierr, arg_pos=10_int32, min=0.0_real64)
        call validate_in_range_real(alpha_accept, ierr, arg_pos=11_int32, min=0.0_real64)
        call validate_in_range_int(t_observables, ierr, arg_pos=12_int32, min=2_int32)
        call validate_in_range_int(n_tiles, ierr, arg_pos=13_int32, min=1_int32)
        call validate_all_in_range_real(vectors, n_dimensions * n_vectors, ierr, arg_pos=1_int32)
        call validate_all_in_range_int(dimension_order, n_dimensions, ierr, arg_pos=4_int32, min=1_int32, max=n_dimensions)
        call validate_all_in_range_int(kd_indices, n_vectors, ierr, arg_pos=5_int32, min=1_int32, max=n_vectors)
        call validate_all_in_range_real(density_labels, n_vectors, ierr, arg_pos=6_int32, min=1.0_real64)
        call validate_all_in_range_int(seed_indices, n_seeds, ierr, arg_pos=7_int32, min=1_int32, max=n_vectors)
        if (is_err(ierr)) return
#endif

        M_ALLOCATE(tmp_stack(3, 64, min(n_tiles,max(1_int32,n_seeds))))
        M_ALLOCATE(tmp_vicinity_mask(n_vectors, min(n_tiles,max(1_int32,n_seeds))))
        M_ALLOCATE(tmp_surface_mask(n_vectors, min(n_tiles,max(1_int32,n_seeds))))
        M_ALLOCATE(tmp_perm(n_vectors, min(n_tiles,max(1_int32,n_seeds))))
        M_ALLOCATE(tmp_abs_diff(n_vectors, min(n_tiles,max(1_int32,n_seeds))))
        M_ALLOCATE(tmp_observables(5, max(2_int32,min(t_observables,n_vectors)), min(n_tiles,max(1_int32,n_seeds))))
        M_ALLOCATE(tmp_current_mask(n_vectors, min(n_tiles,max(1_int32,n_seeds))))

        call obtain_ensembles_impl(&
            vectors = vectors,&
            n_dimensions = n_dimensions,&
            n_vectors = n_vectors,&
            dimension_order = dimension_order,&
            kd_indices = kd_indices,&
            density_labels = density_labels,&
            seed_indices = seed_indices,&
            n_seeds = n_seeds,&
            growth_mode = MODE_GLOBAL_RADIUS,&
            r = r,&
            alpha_mad = alpha_mad,&
            alpha_accept = alpha_accept,&
            t_observables = t_observables,&
            n_tiles = n_tiles,&
            tmp_stack = tmp_stack,&
            tmp_vicinity_mask = tmp_vicinity_mask,&
            tmp_surface_mask = tmp_surface_mask,&
            tmp_perm = tmp_perm,&
            tmp_abs_diff = tmp_abs_diff,&
            tmp_observables = tmp_observables,&
            tmp_current_mask = tmp_current_mask,&
            ensemble_matrix = ensemble_matrix,&
            stop_reasons = stop_reasons,&
            mad_ambient = mad_ambient,&
            ierr = ierr&
        )
        call clear_err_arg_pos(ierr)
    end subroutine obtain_ensembles_global_radius

    !> summary: Validates its inputs, then calls [[tox_stc_impl(module):obtain_ensembles_impl]] with what you supply. The expert entry point: it allocates nothing and prepares nothing; [[tox_stc(module):obtain_ensembles_global_radius]] does both.
    pure subroutine obtain_ensembles_global_radius_expert(&
            vectors,&
            n_dimensions,&
            n_vectors,&
            dimension_order,&
            kd_indices,&
            density_labels,&
            seed_indices,&
            n_seeds,&
            r,&
            alpha_mad,&
            alpha_accept,&
            t_observables,&
            n_tiles,&
            tmp_stack,&
            tmp_vicinity_mask,&
            tmp_surface_mask,&
            tmp_perm,&
            tmp_abs_diff,&
            tmp_observables,&
            tmp_current_mask,&
            ensemble_matrix,&
            stop_reasons,&
            mad_ambient,&
            ierr&
        )
        integer(int32), intent(in) :: n_dimensions
            !! Number of dimensions
        integer(int32), intent(in) :: n_vectors
            !! Number of vectors
        integer(int32), intent(in) :: n_seeds
            !! Number of seed vectors, from `identify_ensemble_seeds`.
            !! The minimum valid value is `0_int32`.
            !! The maximum valid value is `n_vectors`.
        integer(int32), intent(in) :: t_observables
            !! History window in iterations, at least two because acceptance compares two consecutive
            !! steps. Pass n_vectors to retain the full history.
            !! The minimum valid value is `2_int32`.
        integer(int32), intent(in) :: n_tiles
            !! Number of concurrent seed tiles, at least 1. Capped at the number of seeds when nonzero.
            !! The minimum valid value is `1_int32`.
        real(real64), dimension(n_dimensions, n_vectors), intent(in) :: vectors
            !! Input data matrix (n_dimensions x n_vectors)
        integer(int32), dimension(n_dimensions), intent(in) :: dimension_order
            !! Sequence array tracking split dimensions for KD-tree
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_dimensions`.
        integer(int32), dimension(n_vectors), intent(in) :: kd_indices
            !! KD-tree index sequence array computed via [[f42_kd_tree(module):build_kd_index(subroutine)]].
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_vectors`.
        real(real64), dimension(n_vectors), intent(in) :: density_labels
            !! Density labels for all vectors, from `calculate_labels_as_density`.
            !! The minimum valid value is `1.0_real64`.
        integer(int32), dimension(n_seeds), intent(in) :: seed_indices
            !! Indices of the seed vectors to grow from, from `identify_ensemble_seeds`.
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_vectors`.
        real(real64), intent(in) :: r
            !! Uniform growth radius for every seed.
            !! The minimum valid value is `0.0_real64`.
        real(real64), intent(in), optional :: alpha_mad
            !! MAD multiplier for density compatibility threshold
            !! The minimum valid value is `0.0_real64`.
            !! The default value is `0.5_real64`.
        real(real64), intent(in), optional :: alpha_accept
            !! Acceptance log2 fold-change threshold
            !! The minimum valid value is `0.0_real64`.
            !! The default value is `0.5_real64`.
        integer(int32), dimension(3, 64, min(n_tiles,max(1_int32,n_seeds))), intent(inout) :: tmp_stack
            !! Workspace stack for tree traversal per tile
        logical(c_bool), dimension(n_vectors, min(n_tiles,max(1_int32,n_seeds))), intent(out) :: tmp_vicinity_mask
            !! Workspace mask holding one member's neighborhood query result per tile
        logical(c_bool), dimension(n_vectors, min(n_tiles,max(1_int32,n_seeds))), intent(out) :: tmp_surface_mask
            !! Cached spatial neighborhood union per tile, reset for each seed
        integer(int32), dimension(n_vectors, min(n_tiles,max(1_int32,n_seeds))), intent(inout) :: tmp_perm
            !! Sorting workspace, reused for newly accepted member indices between growth rounds per tile
        real(real64), dimension(n_vectors, min(n_tiles,max(1_int32,n_seeds))), intent(out) :: tmp_abs_diff
            !! Workspace array for deviation calculations per tile
        real(real64), dimension(5, max(2_int32,min(t_observables,n_vectors)), min(n_tiles,max(1_int32,n_seeds))), intent(out) :: tmp_observables
            !! Observable history per tile. Growth adds at least one member per iteration, so more than
            !! n_vectors columns are never needed; at least two are kept so acceptance can compare steps.
        logical(c_bool), dimension(n_vectors, min(n_tiles,max(1_int32,n_seeds))), intent(out) :: tmp_current_mask
            !! Workspace holding accepted seed ensemble membership per tile
        logical(c_bool), dimension(n_vectors, n_seeds), intent(out) :: ensemble_matrix
            !! Output matrix storing unmerged grown ensemble masks
        integer(int32), dimension(n_seeds), intent(out) :: stop_reasons
            !! Per-seed stop code. 1: fixed point, 2: rejection after acceptance,
            !! 3: first step rejected, 4: no candidates.
        real(real64), intent(out) :: mad_ambient
            !! Median absolute deviation of `density_labels`, which scales the `alpha_mad` density tolerance.
            !! Zero when there are no seeds; otherwise positive on success, because a zero value is an error.
        integer(int32), intent(out) :: ierr
            !! Invalid input when the median absolute deviation of `density_labels` is zero (more than
            !! half the labels are equal), which leaves no density tolerance for growth.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_dimensions, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_vectors, ierr, arg_pos=3_int32)
        call validate_in_range_int(n_seeds, ierr, arg_pos=8_int32, min=0_int32, max=n_vectors)
        call validate_in_range_real(r, ierr, arg_pos=9_int32, min=0.0_real64)
        call validate_in_range_real(alpha_mad, ierr, arg_pos=10_int32, min=0.0_real64)
        call validate_in_range_real(alpha_accept, ierr, arg_pos=11_int32, min=0.0_real64)
        call validate_in_range_int(t_observables, ierr, arg_pos=12_int32, min=2_int32)
        call validate_in_range_int(n_tiles, ierr, arg_pos=13_int32, min=1_int32)
        call validate_all_in_range_real(vectors, n_dimensions * n_vectors, ierr, arg_pos=1_int32)
        call validate_all_in_range_int(dimension_order, n_dimensions, ierr, arg_pos=4_int32, min=1_int32, max=n_dimensions)
        call validate_all_in_range_int(kd_indices, n_vectors, ierr, arg_pos=5_int32, min=1_int32, max=n_vectors)
        call validate_all_in_range_real(density_labels, n_vectors, ierr, arg_pos=6_int32, min=1.0_real64)
        call validate_all_in_range_int(seed_indices, n_seeds, ierr, arg_pos=7_int32, min=1_int32, max=n_vectors)
        if (is_err(ierr)) return
#endif

        call obtain_ensembles_impl(&
            vectors = vectors,&
            n_dimensions = n_dimensions,&
            n_vectors = n_vectors,&
            dimension_order = dimension_order,&
            kd_indices = kd_indices,&
            density_labels = density_labels,&
            seed_indices = seed_indices,&
            n_seeds = n_seeds,&
            growth_mode = MODE_GLOBAL_RADIUS,&
            r = r,&
            alpha_mad = alpha_mad,&
            alpha_accept = alpha_accept,&
            t_observables = t_observables,&
            n_tiles = n_tiles,&
            tmp_stack = tmp_stack,&
            tmp_vicinity_mask = tmp_vicinity_mask,&
            tmp_surface_mask = tmp_surface_mask,&
            tmp_perm = tmp_perm,&
            tmp_abs_diff = tmp_abs_diff,&
            tmp_observables = tmp_observables,&
            tmp_current_mask = tmp_current_mask,&
            ensemble_matrix = ensemble_matrix,&
            stop_reasons = stop_reasons,&
            mad_ambient = mad_ambient,&
            ierr = ierr&
        )
        call clear_err_arg_pos(ierr)
    end subroutine obtain_ensembles_global_radius_expert

    !> summary: Validates its inputs, prepares what [[tox_stc_impl(module):obtain_ensembles_impl]] needs, then calls it. The entry point to reach for first; see [[tox_stc(module):obtain_ensembles_seed_radius_expert]] to prepare it yourself.
    pure subroutine obtain_ensembles_seed_radius(&
            vectors,&
            n_dimensions,&
            n_vectors,&
            dimension_order,&
            kd_indices,&
            density_labels,&
            seed_indices,&
            n_seeds,&
            seed_radii,&
            alpha_mad,&
            alpha_accept,&
            t_observables,&
            n_tiles,&
            ensemble_matrix,&
            stop_reasons,&
            mad_ambient,&
            ierr&
        )
        integer(int32), intent(in) :: n_dimensions
            !! Number of dimensions
        integer(int32), intent(in) :: n_vectors
            !! Number of vectors
        integer(int32), intent(in) :: n_seeds
            !! Number of seed vectors, from `identify_ensemble_seeds`.
            !! The minimum valid value is `0_int32`.
            !! The maximum valid value is `n_vectors`.
        real(real64), dimension(n_dimensions, n_vectors), intent(in) :: vectors
            !! Input data matrix (n_dimensions x n_vectors)
        integer(int32), dimension(n_dimensions), intent(in) :: dimension_order
            !! Sequence array tracking split dimensions for KD-tree
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_dimensions`.
        integer(int32), dimension(n_vectors), intent(in) :: kd_indices
            !! KD-tree index sequence array computed via [[f42_kd_tree(module):build_kd_index(subroutine)]].
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_vectors`.
        real(real64), dimension(n_vectors), intent(in) :: density_labels
            !! Density labels for all vectors, from `calculate_labels_as_density`.
            !! The minimum valid value is `1.0_real64`.
        integer(int32), dimension(n_seeds), intent(in) :: seed_indices
            !! Indices of the seed vectors to grow from, from `identify_ensemble_seeds`.
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_vectors`.
        real(real64), dimension(n_seeds), intent(in) :: seed_radii
            !! Per-seed growth radius, allowing tighter growth in dense regions and wider growth in sparse regions.
            !! The minimum valid value is `0.0_real64`.
        real(real64), intent(in), optional :: alpha_mad
            !! MAD multiplier for density compatibility threshold
            !! The minimum valid value is `0.0_real64`.
            !! The default value is `0.5_real64`.
        real(real64), intent(in), optional :: alpha_accept
            !! Acceptance log2 fold-change threshold
            !! The minimum valid value is `0.0_real64`.
            !! The default value is `0.5_real64`.
        integer(int32), intent(in) :: t_observables
            !! History window in iterations, at least two because acceptance compares two consecutive
            !! steps. Pass n_vectors to retain the full history.
            !! The minimum valid value is `2_int32`.
        integer(int32), intent(in) :: n_tiles
            !! Number of concurrent seed tiles, at least 1. Capped at the number of seeds when nonzero.
            !! The minimum valid value is `1_int32`.
        logical(c_bool), dimension(n_vectors, n_seeds), intent(out) :: ensemble_matrix
            !! Output matrix storing unmerged grown ensemble masks
        integer(int32), dimension(n_seeds), intent(out) :: stop_reasons
            !! Per-seed stop code. 1: fixed point, 2: rejection after acceptance,
            !! 3: first step rejected, 4: no candidates.
        real(real64), intent(out) :: mad_ambient
            !! Median absolute deviation of `density_labels`, which scales the `alpha_mad` density tolerance.
            !! Zero when there are no seeds; otherwise positive on success, because a zero value is an error.
        integer(int32), intent(out) :: ierr
            !! Invalid input when the median absolute deviation of `density_labels` is zero (more than
            !! half the labels are equal), which leaves no density tolerance for growth.
        integer(int32), dimension(:, :, :), allocatable :: tmp_stack
        logical(c_bool), dimension(:, :), allocatable :: tmp_vicinity_mask
        logical(c_bool), dimension(:, :), allocatable :: tmp_surface_mask
        integer(int32), dimension(:, :), allocatable :: tmp_perm
        real(real64), dimension(:, :), allocatable :: tmp_abs_diff
        real(real64), dimension(:, :, :), allocatable :: tmp_observables
        logical(c_bool), dimension(:, :), allocatable :: tmp_current_mask

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_dimensions, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_vectors, ierr, arg_pos=3_int32)
        call validate_in_range_int(n_seeds, ierr, arg_pos=8_int32, min=0_int32, max=n_vectors)
        call validate_in_range_real(alpha_mad, ierr, arg_pos=10_int32, min=0.0_real64)
        call validate_in_range_real(alpha_accept, ierr, arg_pos=11_int32, min=0.0_real64)
        call validate_in_range_int(t_observables, ierr, arg_pos=12_int32, min=2_int32)
        call validate_in_range_int(n_tiles, ierr, arg_pos=13_int32, min=1_int32)
        call validate_all_in_range_real(vectors, n_dimensions * n_vectors, ierr, arg_pos=1_int32)
        call validate_all_in_range_int(dimension_order, n_dimensions, ierr, arg_pos=4_int32, min=1_int32, max=n_dimensions)
        call validate_all_in_range_int(kd_indices, n_vectors, ierr, arg_pos=5_int32, min=1_int32, max=n_vectors)
        call validate_all_in_range_real(density_labels, n_vectors, ierr, arg_pos=6_int32, min=1.0_real64)
        call validate_all_in_range_int(seed_indices, n_seeds, ierr, arg_pos=7_int32, min=1_int32, max=n_vectors)
        call validate_all_in_range_real(seed_radii, n_seeds, ierr, arg_pos=9_int32, min=0.0_real64)
        if (is_err(ierr)) return
#endif

        M_ALLOCATE(tmp_stack(3, 64, min(n_tiles,max(1_int32,n_seeds))))
        M_ALLOCATE(tmp_vicinity_mask(n_vectors, min(n_tiles,max(1_int32,n_seeds))))
        M_ALLOCATE(tmp_surface_mask(n_vectors, min(n_tiles,max(1_int32,n_seeds))))
        M_ALLOCATE(tmp_perm(n_vectors, min(n_tiles,max(1_int32,n_seeds))))
        M_ALLOCATE(tmp_abs_diff(n_vectors, min(n_tiles,max(1_int32,n_seeds))))
        M_ALLOCATE(tmp_observables(5, max(2_int32,min(t_observables,n_vectors)), min(n_tiles,max(1_int32,n_seeds))))
        M_ALLOCATE(tmp_current_mask(n_vectors, min(n_tiles,max(1_int32,n_seeds))))

        call obtain_ensembles_impl(&
            vectors = vectors,&
            n_dimensions = n_dimensions,&
            n_vectors = n_vectors,&
            dimension_order = dimension_order,&
            kd_indices = kd_indices,&
            density_labels = density_labels,&
            seed_indices = seed_indices,&
            n_seeds = n_seeds,&
            growth_mode = MODE_SEED_RADIUS,&
            seed_radii = seed_radii,&
            alpha_mad = alpha_mad,&
            alpha_accept = alpha_accept,&
            t_observables = t_observables,&
            n_tiles = n_tiles,&
            tmp_stack = tmp_stack,&
            tmp_vicinity_mask = tmp_vicinity_mask,&
            tmp_surface_mask = tmp_surface_mask,&
            tmp_perm = tmp_perm,&
            tmp_abs_diff = tmp_abs_diff,&
            tmp_observables = tmp_observables,&
            tmp_current_mask = tmp_current_mask,&
            ensemble_matrix = ensemble_matrix,&
            stop_reasons = stop_reasons,&
            mad_ambient = mad_ambient,&
            ierr = ierr&
        )
        call clear_err_arg_pos(ierr)
    end subroutine obtain_ensembles_seed_radius

    !> summary: Validates its inputs, then calls [[tox_stc_impl(module):obtain_ensembles_impl]] with what you supply. The expert entry point: it allocates nothing and prepares nothing; [[tox_stc(module):obtain_ensembles_seed_radius]] does both.
    pure subroutine obtain_ensembles_seed_radius_expert(&
            vectors,&
            n_dimensions,&
            n_vectors,&
            dimension_order,&
            kd_indices,&
            density_labels,&
            seed_indices,&
            n_seeds,&
            seed_radii,&
            alpha_mad,&
            alpha_accept,&
            t_observables,&
            n_tiles,&
            tmp_stack,&
            tmp_vicinity_mask,&
            tmp_surface_mask,&
            tmp_perm,&
            tmp_abs_diff,&
            tmp_observables,&
            tmp_current_mask,&
            ensemble_matrix,&
            stop_reasons,&
            mad_ambient,&
            ierr&
        )
        integer(int32), intent(in) :: n_dimensions
            !! Number of dimensions
        integer(int32), intent(in) :: n_vectors
            !! Number of vectors
        integer(int32), intent(in) :: n_seeds
            !! Number of seed vectors, from `identify_ensemble_seeds`.
            !! The minimum valid value is `0_int32`.
            !! The maximum valid value is `n_vectors`.
        integer(int32), intent(in) :: t_observables
            !! History window in iterations, at least two because acceptance compares two consecutive
            !! steps. Pass n_vectors to retain the full history.
            !! The minimum valid value is `2_int32`.
        integer(int32), intent(in) :: n_tiles
            !! Number of concurrent seed tiles, at least 1. Capped at the number of seeds when nonzero.
            !! The minimum valid value is `1_int32`.
        real(real64), dimension(n_dimensions, n_vectors), intent(in) :: vectors
            !! Input data matrix (n_dimensions x n_vectors)
        integer(int32), dimension(n_dimensions), intent(in) :: dimension_order
            !! Sequence array tracking split dimensions for KD-tree
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_dimensions`.
        integer(int32), dimension(n_vectors), intent(in) :: kd_indices
            !! KD-tree index sequence array computed via [[f42_kd_tree(module):build_kd_index(subroutine)]].
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_vectors`.
        real(real64), dimension(n_vectors), intent(in) :: density_labels
            !! Density labels for all vectors, from `calculate_labels_as_density`.
            !! The minimum valid value is `1.0_real64`.
        integer(int32), dimension(n_seeds), intent(in) :: seed_indices
            !! Indices of the seed vectors to grow from, from `identify_ensemble_seeds`.
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_vectors`.
        real(real64), dimension(n_seeds), intent(in) :: seed_radii
            !! Per-seed growth radius, allowing tighter growth in dense regions and wider growth in sparse regions.
            !! The minimum valid value is `0.0_real64`.
        real(real64), intent(in), optional :: alpha_mad
            !! MAD multiplier for density compatibility threshold
            !! The minimum valid value is `0.0_real64`.
            !! The default value is `0.5_real64`.
        real(real64), intent(in), optional :: alpha_accept
            !! Acceptance log2 fold-change threshold
            !! The minimum valid value is `0.0_real64`.
            !! The default value is `0.5_real64`.
        integer(int32), dimension(3, 64, min(n_tiles,max(1_int32,n_seeds))), intent(inout) :: tmp_stack
            !! Workspace stack for tree traversal per tile
        logical(c_bool), dimension(n_vectors, min(n_tiles,max(1_int32,n_seeds))), intent(out) :: tmp_vicinity_mask
            !! Workspace mask holding one member's neighborhood query result per tile
        logical(c_bool), dimension(n_vectors, min(n_tiles,max(1_int32,n_seeds))), intent(out) :: tmp_surface_mask
            !! Cached spatial neighborhood union per tile, reset for each seed
        integer(int32), dimension(n_vectors, min(n_tiles,max(1_int32,n_seeds))), intent(inout) :: tmp_perm
            !! Sorting workspace, reused for newly accepted member indices between growth rounds per tile
        real(real64), dimension(n_vectors, min(n_tiles,max(1_int32,n_seeds))), intent(out) :: tmp_abs_diff
            !! Workspace array for deviation calculations per tile
        real(real64), dimension(5, max(2_int32,min(t_observables,n_vectors)), min(n_tiles,max(1_int32,n_seeds))), intent(out) :: tmp_observables
            !! Observable history per tile. Growth adds at least one member per iteration, so more than
            !! n_vectors columns are never needed; at least two are kept so acceptance can compare steps.
        logical(c_bool), dimension(n_vectors, min(n_tiles,max(1_int32,n_seeds))), intent(out) :: tmp_current_mask
            !! Workspace holding accepted seed ensemble membership per tile
        logical(c_bool), dimension(n_vectors, n_seeds), intent(out) :: ensemble_matrix
            !! Output matrix storing unmerged grown ensemble masks
        integer(int32), dimension(n_seeds), intent(out) :: stop_reasons
            !! Per-seed stop code. 1: fixed point, 2: rejection after acceptance,
            !! 3: first step rejected, 4: no candidates.
        real(real64), intent(out) :: mad_ambient
            !! Median absolute deviation of `density_labels`, which scales the `alpha_mad` density tolerance.
            !! Zero when there are no seeds; otherwise positive on success, because a zero value is an error.
        integer(int32), intent(out) :: ierr
            !! Invalid input when the median absolute deviation of `density_labels` is zero (more than
            !! half the labels are equal), which leaves no density tolerance for growth.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_dimensions, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_vectors, ierr, arg_pos=3_int32)
        call validate_in_range_int(n_seeds, ierr, arg_pos=8_int32, min=0_int32, max=n_vectors)
        call validate_in_range_real(alpha_mad, ierr, arg_pos=10_int32, min=0.0_real64)
        call validate_in_range_real(alpha_accept, ierr, arg_pos=11_int32, min=0.0_real64)
        call validate_in_range_int(t_observables, ierr, arg_pos=12_int32, min=2_int32)
        call validate_in_range_int(n_tiles, ierr, arg_pos=13_int32, min=1_int32)
        call validate_all_in_range_real(vectors, n_dimensions * n_vectors, ierr, arg_pos=1_int32)
        call validate_all_in_range_int(dimension_order, n_dimensions, ierr, arg_pos=4_int32, min=1_int32, max=n_dimensions)
        call validate_all_in_range_int(kd_indices, n_vectors, ierr, arg_pos=5_int32, min=1_int32, max=n_vectors)
        call validate_all_in_range_real(density_labels, n_vectors, ierr, arg_pos=6_int32, min=1.0_real64)
        call validate_all_in_range_int(seed_indices, n_seeds, ierr, arg_pos=7_int32, min=1_int32, max=n_vectors)
        call validate_all_in_range_real(seed_radii, n_seeds, ierr, arg_pos=9_int32, min=0.0_real64)
        if (is_err(ierr)) return
#endif

        call obtain_ensembles_impl(&
            vectors = vectors,&
            n_dimensions = n_dimensions,&
            n_vectors = n_vectors,&
            dimension_order = dimension_order,&
            kd_indices = kd_indices,&
            density_labels = density_labels,&
            seed_indices = seed_indices,&
            n_seeds = n_seeds,&
            growth_mode = MODE_SEED_RADIUS,&
            seed_radii = seed_radii,&
            alpha_mad = alpha_mad,&
            alpha_accept = alpha_accept,&
            t_observables = t_observables,&
            n_tiles = n_tiles,&
            tmp_stack = tmp_stack,&
            tmp_vicinity_mask = tmp_vicinity_mask,&
            tmp_surface_mask = tmp_surface_mask,&
            tmp_perm = tmp_perm,&
            tmp_abs_diff = tmp_abs_diff,&
            tmp_observables = tmp_observables,&
            tmp_current_mask = tmp_current_mask,&
            ensemble_matrix = ensemble_matrix,&
            stop_reasons = stop_reasons,&
            mad_ambient = mad_ambient,&
            ierr = ierr&
        )
        call clear_err_arg_pos(ierr)
    end subroutine obtain_ensembles_seed_radius_expert

    !> summary: Validates its inputs, prepares what [[tox_stc_impl(module):merge_ensembles_impl]] needs, then calls it. The entry point to reach for first; see [[tox_stc(module):merge_ensembles_expert]] to prepare it yourself.
    pure subroutine merge_ensembles(&
            raw_masks,&
            n_vectors,&
            n_seeds,&
            min_overlap_coefficient,&
            merged_masks,&
            n_ensembles,&
            ierr&
        )
        integer(int32), intent(in) :: n_vectors
            !! Total number of ambient vectors
        integer(int32), intent(in) :: n_seeds
            !! Total number of grown raw ensemble columns
            !! The minimum valid value is `0_int32`.
            !! The maximum valid value is `n_vectors`.
        logical(c_bool), dimension(n_vectors, n_seeds), intent(in) :: raw_masks
            !! Matrix of raw unmerged boolean ensemble masks
        real(real64), intent(in), optional :: min_overlap_coefficient
            !! Minimum Overlap Coefficient required to merge two ensembles
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.9_real64`.
        logical(c_bool), dimension(n_vectors, n_seeds), intent(out) :: merged_masks
            !! Output matrix storing merged boolean ensemble masks
            !! The first `n_ensembles` elements will hold the results.
        integer(int32), intent(out) :: n_ensembles
            !! Final count of merged non-singleton ensembles
        integer(int32), intent(out) :: ierr
            !! Error status flag
        logical(c_bool), dimension(:), allocatable :: tmp_active_flag
        integer(int32), dimension(:), allocatable :: tmp_parent

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_vectors, ierr, arg_pos=2_int32)
        call validate_in_range_int(n_seeds, ierr, arg_pos=3_int32, min=0_int32, max=n_vectors)
        call validate_in_range_real(min_overlap_coefficient, ierr, arg_pos=4_int32, min=0.0_real64, max=1.0_real64)
        if (is_err(ierr)) return
#endif

        M_ALLOCATE(tmp_active_flag(n_seeds))
        M_ALLOCATE(tmp_parent(n_seeds))

        call merge_ensembles_impl(&
            raw_masks = raw_masks,&
            n_vectors = n_vectors,&
            n_seeds = n_seeds,&
            min_overlap_coefficient = min_overlap_coefficient,&
            merged_masks = merged_masks,&
            tmp_active_flag = tmp_active_flag,&
            tmp_parent = tmp_parent,&
            n_ensembles = n_ensembles,&
            ierr = ierr&
        )
        call clear_err_arg_pos(ierr)
    end subroutine merge_ensembles

    !> summary: Validates its inputs, then calls [[tox_stc_impl(module):merge_ensembles_impl]] with what you supply. The expert entry point: it allocates nothing and prepares nothing; [[tox_stc(module):merge_ensembles]] does both.
    pure subroutine merge_ensembles_expert(&
            raw_masks,&
            n_vectors,&
            n_seeds,&
            min_overlap_coefficient,&
            merged_masks,&
            tmp_active_flag,&
            tmp_parent,&
            n_ensembles,&
            ierr&
        )
        integer(int32), intent(in) :: n_vectors
            !! Total number of ambient vectors
        integer(int32), intent(in) :: n_seeds
            !! Total number of grown raw ensemble columns
            !! The minimum valid value is `0_int32`.
            !! The maximum valid value is `n_vectors`.
        logical(c_bool), dimension(n_vectors, n_seeds), intent(in) :: raw_masks
            !! Matrix of raw unmerged boolean ensemble masks
        real(real64), intent(in), optional :: min_overlap_coefficient
            !! Minimum Overlap Coefficient required to merge two ensembles
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.9_real64`.
        logical(c_bool), dimension(n_vectors, n_seeds), intent(out) :: merged_masks
            !! Output matrix storing merged boolean ensemble masks
            !! The first `n_ensembles` elements will hold the results.
        logical(c_bool), dimension(n_seeds), intent(out) :: tmp_active_flag
            !! Preallocated workspace marking the component roots that survive as ensembles
        integer(int32), dimension(n_seeds), intent(inout) :: tmp_parent
            !! Preallocated union-find parent pointers
        integer(int32), intent(out) :: n_ensembles
            !! Final count of merged non-singleton ensembles
        integer(int32), intent(out) :: ierr
            !! Error status flag

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_vectors, ierr, arg_pos=2_int32)
        call validate_in_range_int(n_seeds, ierr, arg_pos=3_int32, min=0_int32, max=n_vectors)
        call validate_in_range_real(min_overlap_coefficient, ierr, arg_pos=4_int32, min=0.0_real64, max=1.0_real64)
        if (is_err(ierr)) return
#endif

        call merge_ensembles_impl(&
            raw_masks = raw_masks,&
            n_vectors = n_vectors,&
            n_seeds = n_seeds,&
            min_overlap_coefficient = min_overlap_coefficient,&
            merged_masks = merged_masks,&
            tmp_active_flag = tmp_active_flag,&
            tmp_parent = tmp_parent,&
            n_ensembles = n_ensembles,&
            ierr = ierr&
        )
        call clear_err_arg_pos(ierr)
    end subroutine merge_ensembles_expert

    !> summary: Validates its inputs, prepares what [[tox_stc_impl(module):grow_ensemble_impl]] needs, then calls it. The entry point to reach for first; see [[tox_stc(module):grow_ensemble_expert]] to prepare it yourself.
    pure subroutine grow_ensemble(&
            vectors,&
            n_dimensions,&
            n_vectors,&
            ensemble_mask,&
            r,&
            dimension_order,&
            kd_indices,&
            density_labels,&
            alpha_mad,&
            grown_mask,&
            ierr&
        )
        integer(int32), intent(in) :: n_dimensions
            !! Number of dimensions
        integer(int32), intent(in) :: n_vectors
            !! Number of vectors
        real(real64), dimension(n_dimensions, n_vectors), intent(in) :: vectors
            !! Input data matrix (n_dimensions x n_vectors)
        logical(c_bool), dimension(n_vectors), intent(in) :: ensemble_mask
            !! Membership before the growth step.
        real(real64), intent(in) :: r
            !! Search radius threshold for surface growth
            !! The minimum valid value is `0.0_real64`.
        integer(int32), dimension(n_dimensions), intent(in) :: dimension_order
            !! Sequence array tracking split dimensions
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_dimensions`.
        integer(int32), dimension(n_vectors), intent(in) :: kd_indices
            !! KD-tree index sequence array computed via [[f42_kd_tree(module):build_kd_index(subroutine)]].
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_vectors`.
        real(real64), dimension(n_vectors), intent(in) :: density_labels
            !! Density labels for all vectors, from `calculate_labels_as_density`.
            !! The minimum valid value is `1.0_real64`.
        real(real64), intent(in), optional :: alpha_mad
            !! Multiplier factor for MAD density compatibility threshold
            !! The minimum valid value is `0.0_real64`.
            !! The default value is `0.5_real64`.
        logical(c_bool), dimension(n_vectors), intent(out) :: grown_mask
            !! Membership after the growth step: every input member plus each density-compatible neighbor.
        integer(int32), intent(out) :: ierr
            !! Invalid input when the median absolute deviation of `density_labels` is zero (more than
            !! half the labels are equal), which leaves no density tolerance for growth.
        integer(int32), dimension(:, :), allocatable :: tmp_stack
        logical(c_bool), dimension(:), allocatable :: tmp_vicinity_mask
        logical(c_bool), dimension(:), allocatable :: tmp_surface_mask
        integer(int32), dimension(:), allocatable :: tmp_perm
        real(real64), dimension(:), allocatable :: tmp_abs_diff

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_dimensions, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_vectors, ierr, arg_pos=3_int32)
        call validate_in_range_real(r, ierr, arg_pos=5_int32, min=0.0_real64)
        call validate_in_range_real(alpha_mad, ierr, arg_pos=9_int32, min=0.0_real64)
        call validate_all_in_range_real(vectors, n_dimensions * n_vectors, ierr, arg_pos=1_int32)
        call validate_all_in_range_int(dimension_order, n_dimensions, ierr, arg_pos=6_int32, min=1_int32, max=n_dimensions)
        call validate_all_in_range_int(kd_indices, n_vectors, ierr, arg_pos=7_int32, min=1_int32, max=n_vectors)
        call validate_all_in_range_real(density_labels, n_vectors, ierr, arg_pos=8_int32, min=1.0_real64)
        if (is_err(ierr)) return
#endif

        M_ALLOCATE(tmp_stack(3, 64))
        M_ALLOCATE(tmp_vicinity_mask(n_vectors))
        M_ALLOCATE(tmp_surface_mask(n_vectors))
        M_ALLOCATE(tmp_perm(n_vectors))
        M_ALLOCATE(tmp_abs_diff(n_vectors))

        call grow_ensemble_impl(&
            vectors = vectors,&
            n_dimensions = n_dimensions,&
            n_vectors = n_vectors,&
            ensemble_mask = ensemble_mask,&
            r = r,&
            dimension_order = dimension_order,&
            kd_indices = kd_indices,&
            density_labels = density_labels,&
            alpha_mad = alpha_mad,&
            tmp_stack = tmp_stack,&
            tmp_vicinity_mask = tmp_vicinity_mask,&
            tmp_surface_mask = tmp_surface_mask,&
            tmp_perm = tmp_perm,&
            tmp_abs_diff = tmp_abs_diff,&
            grown_mask = grown_mask,&
            ierr = ierr&
        )
        call clear_err_arg_pos(ierr)
    end subroutine grow_ensemble

    !> summary: Validates its inputs, then calls [[tox_stc_impl(module):grow_ensemble_impl]] with what you supply. The expert entry point: it allocates nothing and prepares nothing; [[tox_stc(module):grow_ensemble]] does both.
    pure subroutine grow_ensemble_expert(&
            vectors,&
            n_dimensions,&
            n_vectors,&
            ensemble_mask,&
            r,&
            dimension_order,&
            kd_indices,&
            density_labels,&
            alpha_mad,&
            tmp_stack,&
            tmp_vicinity_mask,&
            tmp_surface_mask,&
            tmp_perm,&
            tmp_abs_diff,&
            grown_mask,&
            ierr&
        )
        integer(int32), intent(in) :: n_dimensions
            !! Number of dimensions
        integer(int32), intent(in) :: n_vectors
            !! Number of vectors
        real(real64), dimension(n_dimensions, n_vectors), intent(in) :: vectors
            !! Input data matrix (n_dimensions x n_vectors)
        logical(c_bool), dimension(n_vectors), intent(in) :: ensemble_mask
            !! Membership before the growth step.
        real(real64), intent(in) :: r
            !! Search radius threshold for surface growth
            !! The minimum valid value is `0.0_real64`.
        integer(int32), dimension(n_dimensions), intent(in) :: dimension_order
            !! Sequence array tracking split dimensions
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_dimensions`.
        integer(int32), dimension(n_vectors), intent(in) :: kd_indices
            !! KD-tree index sequence array computed via [[f42_kd_tree(module):build_kd_index(subroutine)]].
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_vectors`.
        real(real64), dimension(n_vectors), intent(in) :: density_labels
            !! Density labels for all vectors, from `calculate_labels_as_density`.
            !! The minimum valid value is `1.0_real64`.
        real(real64), intent(in), optional :: alpha_mad
            !! Multiplier factor for MAD density compatibility threshold
            !! The minimum valid value is `0.0_real64`.
            !! The default value is `0.5_real64`.
        integer(int32), dimension(3, 64), intent(inout) :: tmp_stack
            !! Preallocated workspace stack for a single tree traversal
        logical(c_bool), dimension(n_vectors), intent(out) :: tmp_vicinity_mask
            !! Workspace for one member's neighborhood query
        logical(c_bool), dimension(n_vectors), intent(out) :: tmp_surface_mask
            !! Geometric union of accepted member neighborhoods before this growth step
        integer(int32), dimension(n_vectors), intent(inout) :: tmp_perm
            !! Preallocated workspace array for sorting and percentiles
        real(real64), dimension(n_vectors), intent(out) :: tmp_abs_diff
            !! Preallocated workspace array holding the active member densities
        logical(c_bool), dimension(n_vectors), intent(out) :: grown_mask
            !! Membership after the growth step: every input member plus each density-compatible neighbor.
        integer(int32), intent(out) :: ierr
            !! Invalid input when the median absolute deviation of `density_labels` is zero (more than
            !! half the labels are equal), which leaves no density tolerance for growth.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_dimensions, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_vectors, ierr, arg_pos=3_int32)
        call validate_in_range_real(r, ierr, arg_pos=5_int32, min=0.0_real64)
        call validate_in_range_real(alpha_mad, ierr, arg_pos=9_int32, min=0.0_real64)
        call validate_all_in_range_real(vectors, n_dimensions * n_vectors, ierr, arg_pos=1_int32)
        call validate_all_in_range_int(dimension_order, n_dimensions, ierr, arg_pos=6_int32, min=1_int32, max=n_dimensions)
        call validate_all_in_range_int(kd_indices, n_vectors, ierr, arg_pos=7_int32, min=1_int32, max=n_vectors)
        call validate_all_in_range_real(density_labels, n_vectors, ierr, arg_pos=8_int32, min=1.0_real64)
        if (is_err(ierr)) return
#endif

        call grow_ensemble_impl(&
            vectors = vectors,&
            n_dimensions = n_dimensions,&
            n_vectors = n_vectors,&
            ensemble_mask = ensemble_mask,&
            r = r,&
            dimension_order = dimension_order,&
            kd_indices = kd_indices,&
            density_labels = density_labels,&
            alpha_mad = alpha_mad,&
            tmp_stack = tmp_stack,&
            tmp_vicinity_mask = tmp_vicinity_mask,&
            tmp_surface_mask = tmp_surface_mask,&
            tmp_perm = tmp_perm,&
            tmp_abs_diff = tmp_abs_diff,&
            grown_mask = grown_mask,&
            ierr = ierr&
        )
        call clear_err_arg_pos(ierr)
    end subroutine grow_ensemble_expert

    !> summary: Validates its inputs, then calls [[tox_stc_impl(module):compute_ensemble_observable_impl]].
    pure subroutine compute_ensemble_observable(&
            ensemble_mask,&
            density_labels,&
            n_vectors,&
            candidate_count,&
            current_iter,&
            observables,&
            n_history,&
            n_observables,&
            ierr&
        )
        integer(int32), intent(in) :: n_vectors
            !! Number of vectors
        integer(int32), intent(in) :: n_history
            !! Number of columns in the observable history. This is the history window: once
            !! `current_iter` exceeds it, the columns shift left and the newest step is stored last.
            !! Acceptance compares two consecutive steps, so at least two columns are required.
            !! The minimum valid value is `2_int32`.
        integer(int32), intent(in) :: n_observables
            !! Number of observable rows, which must be exactly five.
            !! The minimum valid value is `5`.
            !! The maximum valid value is `5`.
        logical(c_bool), dimension(n_vectors), intent(in) :: ensemble_mask
            !! Selection mask for active ensemble members
        real(real64), dimension(n_vectors), intent(in) :: density_labels
            !! Density labels for all vectors, from `calculate_labels_as_density`.
            !! The minimum valid value is `1.0_real64`.
        integer(int32), intent(in) :: candidate_count
            !! Number of candidates added in current growth step
            !! The minimum valid value is `0_int32`.
        integer(int32), intent(in) :: current_iter
            !! Current growth iteration index (1-based)
            !! The minimum valid value is `1_int32`.
        real(real64), dimension(n_observables, n_history), intent(inout) :: observables
            !! Observable history, one column per growth iteration within the window
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_vectors, ierr, arg_pos=3_int32)
        call validate_in_range_int(candidate_count, ierr, arg_pos=4_int32, min=0_int32)
        call validate_in_range_int(current_iter, ierr, arg_pos=5_int32, min=1_int32)
        call validate_in_range_int(n_history, ierr, arg_pos=7_int32, min=2_int32)
        call validate_in_range_int(n_observables, ierr, arg_pos=8_int32, min=5, max=5)
        call validate_all_in_range_real(density_labels, n_vectors, ierr, arg_pos=2_int32, min=1.0_real64)
        call validate_all_in_range_real(observables, n_observables * n_history, ierr, arg_pos=6_int32)
        if (is_err(ierr)) return
#endif

        call compute_ensemble_observable_impl(&
            ensemble_mask = ensemble_mask,&
            density_labels = density_labels,&
            n_vectors = n_vectors,&
            candidate_count = candidate_count,&
            current_iter = current_iter,&
            observables = observables,&
            n_history = n_history,&
            n_observables = n_observables&
        )
    end subroutine compute_ensemble_observable

    !> summary: Validates its inputs, then calls [[tox_stc_impl(module):accept_ensemble_impl]].
    pure subroutine accept_ensemble(&
            observables,&
            n_history,&
            n_observables,&
            current_iter,&
            alpha_accept,&
            is_accepted,&
            ierr&
        )
        integer(int32), intent(in) :: n_history
            !! Number of columns in the observable history, the same window the history was written with.
            !! Acceptance compares two consecutive steps, so at least two columns are required.
            !! The minimum valid value is `2_int32`.
        integer(int32), intent(in) :: n_observables
            !! Number of observable rows, which must be exactly five.
            !! The minimum valid value is `5`.
            !! The maximum valid value is `5`.
        real(real64), dimension(n_observables, n_history), intent(in) :: observables
            !! Observable history written by [[tox_stc(module):compute_ensemble_observable]]
        integer(int32), intent(in) :: current_iter
            !! Current growth iteration index (1-based)
            !! The minimum valid value is `1_int32`.
        real(real64), intent(in), optional :: alpha_accept
            !! Log2 fold-change threshold; a step whose density heterogeneity changes by at least this much is rejected
            !! The minimum valid value is `0.0_real64`.
            !! The default value is `0.5_real64`.
        logical(c_bool), intent(out) :: is_accepted
            !! Output flag indicating step acceptance decision
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_in_range_int(n_history, ierr, arg_pos=2_int32, min=2_int32)
        call validate_in_range_int(n_observables, ierr, arg_pos=3_int32, min=5, max=5)
        call validate_in_range_int(current_iter, ierr, arg_pos=4_int32, min=1_int32)
        call validate_in_range_real(alpha_accept, ierr, arg_pos=5_int32, min=0.0_real64)
        call validate_all_in_range_real(observables, n_observables * n_history, ierr, arg_pos=1_int32)
        if (is_err(ierr)) return
#endif

        call accept_ensemble_impl(&
            observables = observables,&
            n_history = n_history,&
            n_observables = n_observables,&
            current_iter = current_iter,&
            alpha_accept = alpha_accept,&
            is_accepted = is_accepted&
        )
    end subroutine accept_ensemble

end module tox_stc
