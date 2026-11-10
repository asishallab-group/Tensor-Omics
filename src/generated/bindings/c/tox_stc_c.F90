#ifndef NO_C_BINDING
#include <src/macros.h>

!> summary: C-wrappers for [[tox_stc(module)]]
!| Shape-truthful clustering with a global or per-seed growth radius.
!| Compute density labels, select seeds, grow their ensembles, and merge overlapping ensembles.
!| Both radius variants use the same density labels, seed ordering, growth rule and merge rule.
module tox_stc_c
    use f42_safeguard
    use, intrinsic :: iso_c_binding, only: c_associated, c_bool, c_double, c_int, c_loc
    use tox_errors, only: set_ok, set_err, ERR_POINTER_NULL
    M_IMPLICIT_NONE
    private

    public :: calculate_density_radius_c
    public :: calculate_density_radius_expert_c
    public :: calculate_labels_as_density_c
    public :: calculate_labels_as_density_expert_c
    public :: identify_ensemble_seeds_c
    public :: identify_ensemble_seeds_expert_c
    public :: compute_ambient_density_stats_c
    public :: compute_ambient_density_stats_expert_c
    public :: merge_ensembles_c
    public :: merge_ensembles_expert_c
    public :: grow_ensemble_c
    public :: grow_ensemble_expert_c
    public :: compute_ensemble_observable_c
    public :: accept_ensemble_c
    public :: obtain_ensembles_global_radius_c
    public :: obtain_ensembles_global_radius_expert_c
    public :: obtain_ensembles_seed_radius_c
    public :: obtain_ensembles_seed_radius_expert_c

contains

    !> summary: C-wrapper for [[tox_stc(module):calculate_density_radius(subroutine)]]
    subroutine calculate_density_radius_c(&
            vectors,&
            n_dimensions,&
            n_vectors,&
            mean_to_other_vecs_dist_quant,&
            radius,&
            ierr&
        ) bind(C, name="calculate_density_radius_c")
        use tox_stc, only: calculate_density_radius

        integer(c_int), intent(in), target :: n_dimensions
            !! Number of dimensions
        integer(c_int), intent(in), target :: n_vectors
            !! Number of vectors
        real(c_double), dimension(n_dimensions, n_vectors), intent(in), target :: vectors
            !! Input data matrix (n_dimensions x n_vectors)
        real(c_double), intent(in), target :: mean_to_other_vecs_dist_quant
            !! Optional quantile fraction (0.0 to 1.0), defaults to 0.15
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.15_real64`.
        real(c_double), intent(out), target :: radius
            !! The resulting density search radius
        integer(c_int), intent(out), target :: ierr
            !! Error code

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_dimensions)
        M_CHECK_NON_NULL(n_vectors)
        M_CHECK_NON_NULL(mean_to_other_vecs_dist_quant)
        M_CHECK_NON_NULL(radius)
        M_CHECK_ARRAY_NON_NULL(vectors, n_dimensions * n_vectors)

        call calculate_density_radius(&
            vectors = vectors,&
            n_dimensions = n_dimensions,&
            n_vectors = n_vectors,&
            mean_to_other_vecs_dist_quant = mean_to_other_vecs_dist_quant,&
            radius = radius,&
            ierr = ierr&
        )
    end subroutine calculate_density_radius_c

    !> summary: C-wrapper for [[tox_stc(module):calculate_density_radius_expert(subroutine)]]
    subroutine calculate_density_radius_expert_c(&
            vectors,&
            n_dimensions,&
            n_vectors,&
            tmp_mean_vec,&
            tmp_distances,&
            tmp_perm,&
            mean_to_other_vecs_dist_quant,&
            radius,&
            ierr&
        ) bind(C, name="calculate_density_radius_expert_c")
        use tox_stc, only: calculate_density_radius_expert

        integer(c_int), intent(in), target :: n_dimensions
            !! Number of dimensions
        integer(c_int), intent(in), target :: n_vectors
            !! Number of vectors
        real(c_double), dimension(n_dimensions, n_vectors), intent(in), target :: vectors
            !! Input data matrix (n_dimensions x n_vectors)
        real(c_double), dimension(n_dimensions), intent(out), target :: tmp_mean_vec
            !! Preallocated workspace for the mean vector
        real(c_double), dimension(n_vectors), intent(out), target :: tmp_distances
            !! Preallocated workspace for calculated distances
        integer(c_int), dimension(n_vectors), intent(out), target :: tmp_perm
            !! Preallocated workspace for permutation indices
        real(c_double), intent(in), target :: mean_to_other_vecs_dist_quant
            !! Optional quantile fraction (0.0 to 1.0), defaults to 0.15
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.15_real64`.
        real(c_double), intent(out), target :: radius
            !! The resulting density search radius
        integer(c_int), intent(out), target :: ierr
            !! Error code

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_dimensions)
        M_CHECK_NON_NULL(n_vectors)
        M_CHECK_NON_NULL(mean_to_other_vecs_dist_quant)
        M_CHECK_NON_NULL(radius)
        M_CHECK_ARRAY_NON_NULL(vectors, n_dimensions * n_vectors)
        M_CHECK_ARRAY_NON_NULL(tmp_mean_vec, n_dimensions)
        M_CHECK_ARRAY_NON_NULL(tmp_distances, n_vectors)
        M_CHECK_ARRAY_NON_NULL(tmp_perm, n_vectors)

        call calculate_density_radius_expert(&
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
    end subroutine calculate_density_radius_expert_c

    !> summary: C-wrapper for [[tox_stc(module):calculate_labels_as_density(subroutine)]]
    subroutine calculate_labels_as_density_c(&
            vectors,&
            n_dimensions,&
            n_vectors,&
            r,&
            dimension_order,&
            kd_indices,&
            n_tiles,&
            label_densities,&
            ierr&
        ) bind(C, name="calculate_labels_as_density_c")
        use tox_stc, only: calculate_labels_as_density

        integer(c_int), intent(in), target :: n_dimensions
            !! Number of dimensions
        integer(c_int), intent(in), target :: n_vectors
            !! Number of vectors
        real(c_double), dimension(n_dimensions, n_vectors), intent(in), target :: vectors
            !! Input data matrix (n_dimensions x n_vectors)
        real(c_double), intent(in), target :: r
            !! Label-sphere search radius
            !! The minimum valid value is `0.0_real64`.
        integer(c_int), dimension(n_dimensions), intent(in), target :: dimension_order
            !! Sequence array tracking tree split axes by variance
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_dimensions`.
        integer(c_int), dimension(n_vectors), intent(in), target :: kd_indices
            !! KD-tree index sequence array computed via [[f42_kd_tree(module):build_kd_index(subroutine)]].
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_vectors`.
        integer(c_int), intent(in), target :: n_tiles
            !! Number of concurrent vector tiles, at least 1. Capped at the number of points.
            !! The minimum valid value is `1_int32`.
        real(c_double), dimension(n_vectors), intent(out), target :: label_densities
            !! Output array storing generated density scalars
        integer(c_int), intent(out), target :: ierr
            !! Error code

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_dimensions)
        M_CHECK_NON_NULL(n_vectors)
        M_CHECK_NON_NULL(r)
        M_CHECK_NON_NULL(n_tiles)
        M_CHECK_ARRAY_NON_NULL(vectors, n_dimensions * n_vectors)
        M_CHECK_ARRAY_NON_NULL(dimension_order, n_dimensions)
        M_CHECK_ARRAY_NON_NULL(kd_indices, n_vectors)
        M_CHECK_ARRAY_NON_NULL(label_densities, n_vectors)

        call calculate_labels_as_density(&
            vectors = vectors,&
            n_dimensions = n_dimensions,&
            n_vectors = n_vectors,&
            r = r,&
            dimension_order = dimension_order,&
            kd_indices = kd_indices,&
            n_tiles = n_tiles,&
            label_densities = label_densities,&
            ierr = ierr&
        )
    end subroutine calculate_labels_as_density_c

    !> summary: C-wrapper for [[tox_stc(module):calculate_labels_as_density_expert(subroutine)]]
    subroutine calculate_labels_as_density_expert_c(&
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
        ) bind(C, name="calculate_labels_as_density_expert_c")
        use tox_stc, only: calculate_labels_as_density_expert

        integer(c_int), intent(in), target :: n_dimensions
            !! Number of dimensions
        integer(c_int), intent(in), target :: n_vectors
            !! Number of vectors
        integer(c_int), intent(in), target :: n_tiles
            !! Number of concurrent vector tiles, at least 1. Capped at the number of points.
            !! The minimum valid value is `1_int32`.
        real(c_double), dimension(n_dimensions, n_vectors), intent(in), target :: vectors
            !! Input data matrix (n_dimensions x n_vectors)
        real(c_double), intent(in), target :: r
            !! Label-sphere search radius
            !! The minimum valid value is `0.0_real64`.
        integer(c_int), dimension(n_dimensions), intent(in), target :: dimension_order
            !! Sequence array tracking tree split axes by variance
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_dimensions`.
        integer(c_int), dimension(n_vectors), intent(in), target :: kd_indices
            !! KD-tree index sequence array computed via [[f42_kd_tree(module):build_kd_index(subroutine)]].
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_vectors`.
        integer(c_int), dimension(3, 64, min(n_tiles,n_vectors)), intent(inout), target :: tmp_stack
            !! Preallocated traversal stack per tile, reused serially across its vectors
        real(c_double), dimension(n_vectors), intent(out), target :: label_densities
            !! Output array storing generated density scalars
        integer(c_int), intent(out), target :: ierr
            !! Error code

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_dimensions)
        M_CHECK_NON_NULL(n_vectors)
        M_CHECK_NON_NULL(r)
        M_CHECK_NON_NULL(n_tiles)
        M_CHECK_ARRAY_NON_NULL(vectors, n_dimensions * n_vectors)
        M_CHECK_ARRAY_NON_NULL(dimension_order, n_dimensions)
        M_CHECK_ARRAY_NON_NULL(kd_indices, n_vectors)
        M_CHECK_ARRAY_NON_NULL(tmp_stack, 3 * 64 * (min(n_tiles,n_vectors)))
        M_CHECK_ARRAY_NON_NULL(label_densities, n_vectors)

        call calculate_labels_as_density_expert(&
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
    end subroutine calculate_labels_as_density_expert_c

    !> summary: C-wrapper for [[tox_stc(module):identify_ensemble_seeds(subroutine)]]
    subroutine identify_ensemble_seeds_c(&
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
        ) bind(C, name="identify_ensemble_seeds_c")
        use tox_stc, only: identify_ensemble_seeds

        integer(c_int), intent(in), target :: n_dimensions
            !! Number of dimensions
        integer(c_int), intent(in), target :: n_vectors
            !! Number of vectors
            !! The minimum valid value is `2_int32`.
        real(c_double), dimension(n_dimensions, n_vectors), intent(in), target :: vectors
            !! Input data matrix (n_dimensions x n_vectors)
        real(c_double), dimension(n_vectors), intent(in), target :: density_labels
            !! Precalculated density labels for all ambient vectors
            !! The minimum valid value is `1.0_real64`.
        integer(c_int), dimension(n_dimensions), intent(in), target :: dimension_order
            !! Dimension split order array tracking the K-D tree structure
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_dimensions`.
        integer(c_int), dimension(n_vectors), intent(in), target :: kd_indices
            !! K-D tree index sequence array
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_vectors`.
        integer(c_int), intent(in), target :: k_seeding
            !! Number of nearest neighbors used to estimate each seed's local coverage radius
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_vectors - 1_int32`.
        integer(c_int), dimension(n_vectors), intent(out), target :: sorted_perm
            !! Vector indices sorted by density descending
        integer(c_int), intent(out), target :: n_seeds
            !! Number of selected seed vectors
        logical(c_bool), dimension(n_vectors), intent(out), target :: seed_mask
            !! Logical mask identifying selected seed vectors
        real(c_double), dimension(n_vectors), intent(out), target :: seed_radii
            !! Per-seed growth radius, stored at each seed's own vector index; zero at non-seeds
        real(c_double), intent(in), target :: coverage_quant
            !! Optional percentile of k-nearest-neighbor distances used as the seed coverage and growth radius; defaults to 0.5.
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.5_real64`.
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_dimensions)
        M_CHECK_NON_NULL(n_vectors)
        M_CHECK_NON_NULL(k_seeding)
        M_CHECK_NON_NULL(n_seeds)
        M_CHECK_NON_NULL(coverage_quant)
        M_CHECK_ARRAY_NON_NULL(vectors, n_dimensions * n_vectors)
        M_CHECK_ARRAY_NON_NULL(density_labels, n_vectors)
        M_CHECK_ARRAY_NON_NULL(dimension_order, n_dimensions)
        M_CHECK_ARRAY_NON_NULL(kd_indices, n_vectors)
        M_CHECK_ARRAY_NON_NULL(sorted_perm, n_vectors)
        M_CHECK_ARRAY_NON_NULL(seed_mask, n_vectors)
        M_CHECK_ARRAY_NON_NULL(seed_radii, n_vectors)

        call identify_ensemble_seeds(&
            vectors = vectors,&
            n_dimensions = n_dimensions,&
            n_vectors = n_vectors,&
            density_labels = density_labels,&
            dimension_order = dimension_order,&
            kd_indices = kd_indices,&
            k_seeding = k_seeding,&
            sorted_perm = sorted_perm,&
            n_seeds = n_seeds,&
            seed_mask = seed_mask,&
            seed_radii = seed_radii,&
            coverage_quant = coverage_quant,&
            ierr = ierr&
        )
    end subroutine identify_ensemble_seeds_c

    !> summary: C-wrapper for [[tox_stc(module):identify_ensemble_seeds_expert(subroutine)]]
    subroutine identify_ensemble_seeds_expert_c(&
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
        ) bind(C, name="identify_ensemble_seeds_expert_c")
        use tox_stc, only: identify_ensemble_seeds_expert

        integer(c_int), intent(in), target :: n_dimensions
            !! Number of dimensions
        integer(c_int), intent(in), target :: n_vectors
            !! Number of vectors
            !! The minimum valid value is `2_int32`.
        real(c_double), dimension(n_dimensions, n_vectors), intent(in), target :: vectors
            !! Input data matrix (n_dimensions x n_vectors)
        real(c_double), dimension(n_vectors), intent(in), target :: density_labels
            !! Precalculated density labels for all ambient vectors
            !! The minimum valid value is `1.0_real64`.
        integer(c_int), dimension(n_dimensions), intent(in), target :: dimension_order
            !! Dimension split order array tracking the K-D tree structure
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_dimensions`.
        integer(c_int), dimension(n_vectors), intent(in), target :: kd_indices
            !! K-D tree index sequence array
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_vectors`.
        integer(c_int), intent(in), target :: k_seeding
            !! Number of nearest neighbors used to estimate each seed's local coverage radius
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_vectors - 1_int32`.
        integer(c_int), dimension(n_vectors), intent(out), target :: tmp_perm
            !! Preallocated permutation workspace
        real(c_double), dimension(n_vectors), intent(out), target :: tmp_distances
            !! Preallocated candidate-to-vector distance workspace
        integer(c_int), dimension(3, 64), intent(inout), target :: tmp_stack
            !! Preallocated K-D tree traversal stack
        logical(c_bool), dimension(n_vectors), intent(out), target :: tmp_visited_mask
            !! Preallocated workspace tracking vectors already covered during seeding
        logical(c_bool), dimension(n_vectors), intent(out), target :: tmp_newly_covered_mask
            !! Preallocated workspace for the current seed's coverage query
        integer(c_int), dimension(n_vectors), intent(out), target :: sorted_perm
            !! Vector indices sorted by density descending
        integer(c_int), intent(out), target :: n_seeds
            !! Number of selected seed vectors
        logical(c_bool), dimension(n_vectors), intent(out), target :: seed_mask
            !! Logical mask identifying selected seed vectors
        real(c_double), dimension(n_vectors), intent(out), target :: seed_radii
            !! Per-seed growth radius, stored at each seed's own vector index; zero at non-seeds
        real(c_double), intent(in), target :: coverage_quant
            !! Optional percentile of k-nearest-neighbor distances used as the seed coverage and growth radius; defaults to 0.5.
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.5_real64`.
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_dimensions)
        M_CHECK_NON_NULL(n_vectors)
        M_CHECK_NON_NULL(k_seeding)
        M_CHECK_NON_NULL(n_seeds)
        M_CHECK_NON_NULL(coverage_quant)
        M_CHECK_ARRAY_NON_NULL(vectors, n_dimensions * n_vectors)
        M_CHECK_ARRAY_NON_NULL(density_labels, n_vectors)
        M_CHECK_ARRAY_NON_NULL(dimension_order, n_dimensions)
        M_CHECK_ARRAY_NON_NULL(kd_indices, n_vectors)
        M_CHECK_ARRAY_NON_NULL(tmp_perm, n_vectors)
        M_CHECK_ARRAY_NON_NULL(tmp_distances, n_vectors)
        M_CHECK_ARRAY_NON_NULL(tmp_stack, 3 * 64)
        M_CHECK_ARRAY_NON_NULL(tmp_visited_mask, n_vectors)
        M_CHECK_ARRAY_NON_NULL(tmp_newly_covered_mask, n_vectors)
        M_CHECK_ARRAY_NON_NULL(sorted_perm, n_vectors)
        M_CHECK_ARRAY_NON_NULL(seed_mask, n_vectors)
        M_CHECK_ARRAY_NON_NULL(seed_radii, n_vectors)

        call identify_ensemble_seeds_expert(&
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
            coverage_quant = coverage_quant,&
            ierr = ierr&
        )
    end subroutine identify_ensemble_seeds_expert_c

    !> summary: C-wrapper for [[tox_stc(module):compute_ambient_density_stats(subroutine)]]
    subroutine compute_ambient_density_stats_c(&
            density_labels,&
            n_vectors,&
            median_ambient,&
            mad_ambient,&
            ierr&
        ) bind(C, name="compute_ambient_density_stats_c")
        use tox_stc, only: compute_ambient_density_stats

        integer(c_int), intent(in), target :: n_vectors
            !! Number of vectors
        real(c_double), dimension(n_vectors), intent(in), target :: density_labels
            !! Density labels for all vectors the density-label routine.
            !! The minimum valid value is `1.0_real64`.
        real(c_double), intent(out), target :: median_ambient
            !! Ambient median density median_{i=1..N}(rho_i)
        real(c_double), intent(out), target :: mad_ambient
            !! Ambient median absolute deviation median_{i=1..N}(|rho_i - median_ambient|)
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_vectors)
        M_CHECK_NON_NULL(median_ambient)
        M_CHECK_NON_NULL(mad_ambient)
        M_CHECK_ARRAY_NON_NULL(density_labels, n_vectors)

        call compute_ambient_density_stats(&
            density_labels = density_labels,&
            n_vectors = n_vectors,&
            median_ambient = median_ambient,&
            mad_ambient = mad_ambient,&
            ierr = ierr&
        )
    end subroutine compute_ambient_density_stats_c

    !> summary: C-wrapper for [[tox_stc(module):compute_ambient_density_stats_expert(subroutine)]]
    subroutine compute_ambient_density_stats_expert_c(&
            density_labels,&
            n_vectors,&
            tmp_perm,&
            tmp_abs_diff,&
            median_ambient,&
            mad_ambient,&
            ierr&
        ) bind(C, name="compute_ambient_density_stats_expert_c")
        use tox_stc, only: compute_ambient_density_stats_expert

        integer(c_int), intent(in), target :: n_vectors
            !! Number of vectors
        real(c_double), dimension(n_vectors), intent(in), target :: density_labels
            !! Density labels for all vectors the density-label routine.
            !! The minimum valid value is `1.0_real64`.
        integer(c_int), dimension(n_vectors), intent(inout), target :: tmp_perm
            !! Preallocated workspace array for sorting and percentiles
        real(c_double), dimension(n_vectors), intent(out), target :: tmp_abs_diff
            !! Preallocated workspace array for the absolute deviations
        real(c_double), intent(out), target :: median_ambient
            !! Ambient median density median_{i=1..N}(rho_i)
        real(c_double), intent(out), target :: mad_ambient
            !! Ambient median absolute deviation median_{i=1..N}(|rho_i - median_ambient|)
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_vectors)
        M_CHECK_NON_NULL(median_ambient)
        M_CHECK_NON_NULL(mad_ambient)
        M_CHECK_ARRAY_NON_NULL(density_labels, n_vectors)
        M_CHECK_ARRAY_NON_NULL(tmp_perm, n_vectors)
        M_CHECK_ARRAY_NON_NULL(tmp_abs_diff, n_vectors)

        call compute_ambient_density_stats_expert(&
            density_labels = density_labels,&
            n_vectors = n_vectors,&
            tmp_perm = tmp_perm,&
            tmp_abs_diff = tmp_abs_diff,&
            median_ambient = median_ambient,&
            mad_ambient = mad_ambient,&
            ierr = ierr&
        )
    end subroutine compute_ambient_density_stats_expert_c

    !> summary: C-wrapper for [[tox_stc(module):merge_ensembles(subroutine)]]
    subroutine merge_ensembles_c(&
            raw_masks,&
            n_vectors,&
            n_seeds,&
            min_overlap_coefficient,&
            merged_masks,&
            n_ensembles,&
            ierr&
        ) bind(C, name="merge_ensembles_c")
        use tox_stc, only: merge_ensembles

        integer(c_int), intent(in), target :: n_vectors
            !! Total number of ambient vectors
        integer(c_int), intent(in), target :: n_seeds
            !! Total number of grown raw ensemble columns
            !! The minimum valid value is `0_int32`.
            !! The maximum valid value is `n_vectors`.
        logical(c_bool), dimension(n_vectors, n_seeds), intent(in), target :: raw_masks
            !! Matrix of raw unmerged boolean ensemble masks
        real(c_double), intent(in), target :: min_overlap_coefficient
            !! Minimum Overlap Coefficient required to merge two ensembles
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.9_real64`.
        logical(c_bool), dimension(n_vectors, n_seeds), intent(out), target :: merged_masks
            !! Output matrix storing merged boolean ensemble masks
            !! The first `n_ensembles` elements will hold the results.
        integer(c_int), intent(out), target :: n_ensembles
            !! Final count of merged non-singleton ensembles
        integer(c_int), intent(out), target :: ierr
            !! Error status flag

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_vectors)
        M_CHECK_NON_NULL(n_seeds)
        M_CHECK_NON_NULL(min_overlap_coefficient)
        M_CHECK_NON_NULL(n_ensembles)
        M_CHECK_ARRAY_NON_NULL(raw_masks, n_vectors * n_seeds)
        M_CHECK_ARRAY_NON_NULL(merged_masks, n_vectors * n_seeds)

        call merge_ensembles(&
            raw_masks = raw_masks,&
            n_vectors = n_vectors,&
            n_seeds = n_seeds,&
            min_overlap_coefficient = min_overlap_coefficient,&
            merged_masks = merged_masks,&
            n_ensembles = n_ensembles,&
            ierr = ierr&
        )
    end subroutine merge_ensembles_c

    !> summary: C-wrapper for [[tox_stc(module):merge_ensembles_expert(subroutine)]]
    subroutine merge_ensembles_expert_c(&
            raw_masks,&
            n_vectors,&
            n_seeds,&
            min_overlap_coefficient,&
            merged_masks,&
            tmp_active_flag,&
            tmp_parent,&
            n_ensembles,&
            ierr&
        ) bind(C, name="merge_ensembles_expert_c")
        use tox_stc, only: merge_ensembles_expert

        integer(c_int), intent(in), target :: n_vectors
            !! Total number of ambient vectors
        integer(c_int), intent(in), target :: n_seeds
            !! Total number of grown raw ensemble columns
            !! The minimum valid value is `0_int32`.
            !! The maximum valid value is `n_vectors`.
        logical(c_bool), dimension(n_vectors, n_seeds), intent(in), target :: raw_masks
            !! Matrix of raw unmerged boolean ensemble masks
        real(c_double), intent(in), target :: min_overlap_coefficient
            !! Minimum Overlap Coefficient required to merge two ensembles
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.9_real64`.
        logical(c_bool), dimension(n_vectors, n_seeds), intent(out), target :: merged_masks
            !! Output matrix storing merged boolean ensemble masks
            !! The first `n_ensembles` elements will hold the results.
        logical(c_bool), dimension(n_seeds), intent(out), target :: tmp_active_flag
            !! Preallocated workspace marking the component roots that survive as ensembles
        integer(c_int), dimension(n_seeds), intent(inout), target :: tmp_parent
            !! Preallocated union-find parent pointers
        integer(c_int), intent(out), target :: n_ensembles
            !! Final count of merged non-singleton ensembles
        integer(c_int), intent(out), target :: ierr
            !! Error status flag

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_vectors)
        M_CHECK_NON_NULL(n_seeds)
        M_CHECK_NON_NULL(min_overlap_coefficient)
        M_CHECK_NON_NULL(n_ensembles)
        M_CHECK_ARRAY_NON_NULL(raw_masks, n_vectors * n_seeds)
        M_CHECK_ARRAY_NON_NULL(merged_masks, n_vectors * n_seeds)
        M_CHECK_ARRAY_NON_NULL(tmp_active_flag, n_seeds)
        M_CHECK_ARRAY_NON_NULL(tmp_parent, n_seeds)

        call merge_ensembles_expert(&
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
    end subroutine merge_ensembles_expert_c

    !> summary: C-wrapper for [[tox_stc(module):grow_ensemble(subroutine)]]
    subroutine grow_ensemble_c(&
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
        ) bind(C, name="grow_ensemble_c")
        use tox_stc, only: grow_ensemble

        integer(c_int), intent(in), target :: n_dimensions
            !! Number of dimensions
        integer(c_int), intent(in), target :: n_vectors
            !! Number of vectors
        real(c_double), dimension(n_dimensions, n_vectors), intent(in), target :: vectors
            !! Input data matrix (n_dimensions x n_vectors)
        logical(c_bool), dimension(n_vectors), intent(in), target :: ensemble_mask
            !! Membership before the growth step.
        real(c_double), intent(in), target :: r
            !! Search radius threshold for surface growth
            !! The minimum valid value is `0.0_real64`.
        integer(c_int), dimension(n_dimensions), intent(in), target :: dimension_order
            !! Sequence array tracking split dimensions
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_dimensions`.
        integer(c_int), dimension(n_vectors), intent(in), target :: kd_indices
            !! KD-tree index sequence array computed via [[f42_kd_tree(module):build_kd_index(subroutine)]].
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_vectors`.
        real(c_double), dimension(n_vectors), intent(in), target :: density_labels
            !! Density labels for all vectors the density-label routine.
            !! The minimum valid value is `1.0_real64`.
        real(c_double), intent(in), target :: alpha_mad
            !! Multiplier factor for MAD density compatibility threshold
            !! The minimum valid value is `0.0_real64`.
            !! The default value is `0.5_real64`.
        logical(c_bool), dimension(n_vectors), intent(out), target :: grown_mask
            !! Membership after the growth step.
        integer(c_int), intent(out), target :: ierr
            !! Runtime status for the density contrast check.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_dimensions)
        M_CHECK_NON_NULL(n_vectors)
        M_CHECK_NON_NULL(r)
        M_CHECK_NON_NULL(alpha_mad)
        M_CHECK_ARRAY_NON_NULL(vectors, n_dimensions * n_vectors)
        M_CHECK_ARRAY_NON_NULL(ensemble_mask, n_vectors)
        M_CHECK_ARRAY_NON_NULL(dimension_order, n_dimensions)
        M_CHECK_ARRAY_NON_NULL(kd_indices, n_vectors)
        M_CHECK_ARRAY_NON_NULL(density_labels, n_vectors)
        M_CHECK_ARRAY_NON_NULL(grown_mask, n_vectors)

        call grow_ensemble(&
            vectors = vectors,&
            n_dimensions = n_dimensions,&
            n_vectors = n_vectors,&
            ensemble_mask = ensemble_mask,&
            r = r,&
            dimension_order = dimension_order,&
            kd_indices = kd_indices,&
            density_labels = density_labels,&
            alpha_mad = alpha_mad,&
            grown_mask = grown_mask,&
            ierr = ierr&
        )
    end subroutine grow_ensemble_c

    !> summary: C-wrapper for [[tox_stc(module):grow_ensemble_expert(subroutine)]]
    subroutine grow_ensemble_expert_c(&
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
        ) bind(C, name="grow_ensemble_expert_c")
        use tox_stc, only: grow_ensemble_expert

        integer(c_int), intent(in), target :: n_dimensions
            !! Number of dimensions
        integer(c_int), intent(in), target :: n_vectors
            !! Number of vectors
        real(c_double), dimension(n_dimensions, n_vectors), intent(in), target :: vectors
            !! Input data matrix (n_dimensions x n_vectors)
        logical(c_bool), dimension(n_vectors), intent(in), target :: ensemble_mask
            !! Membership before the growth step.
        real(c_double), intent(in), target :: r
            !! Search radius threshold for surface growth
            !! The minimum valid value is `0.0_real64`.
        integer(c_int), dimension(n_dimensions), intent(in), target :: dimension_order
            !! Sequence array tracking split dimensions
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_dimensions`.
        integer(c_int), dimension(n_vectors), intent(in), target :: kd_indices
            !! KD-tree index sequence array computed via [[f42_kd_tree(module):build_kd_index(subroutine)]].
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_vectors`.
        real(c_double), dimension(n_vectors), intent(in), target :: density_labels
            !! Density labels for all vectors the density-label routine.
            !! The minimum valid value is `1.0_real64`.
        real(c_double), intent(in), target :: alpha_mad
            !! Multiplier factor for MAD density compatibility threshold
            !! The minimum valid value is `0.0_real64`.
            !! The default value is `0.5_real64`.
        integer(c_int), dimension(3, 64), intent(inout), target :: tmp_stack
            !! Preallocated workspace stack for a single tree traversal
        logical(c_bool), dimension(n_vectors), intent(out), target :: tmp_vicinity_mask
            !! Workspace for one neighborhood query, then the complete trial membership
        logical(c_bool), dimension(n_vectors), intent(out), target :: tmp_surface_mask
            !! Geometric union of accepted member neighborhoods before this growth step
        integer(c_int), dimension(n_vectors), intent(inout), target :: tmp_perm
            !! Preallocated workspace array for sorting and percentiles
        real(c_double), dimension(n_vectors), intent(out), target :: tmp_abs_diff
            !! Preallocated workspace array holding the active member densities
        logical(c_bool), dimension(n_vectors), intent(out), target :: grown_mask
            !! Membership after the growth step.
        integer(c_int), intent(out), target :: ierr
            !! Runtime status for the density contrast check.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_dimensions)
        M_CHECK_NON_NULL(n_vectors)
        M_CHECK_NON_NULL(r)
        M_CHECK_NON_NULL(alpha_mad)
        M_CHECK_ARRAY_NON_NULL(vectors, n_dimensions * n_vectors)
        M_CHECK_ARRAY_NON_NULL(ensemble_mask, n_vectors)
        M_CHECK_ARRAY_NON_NULL(dimension_order, n_dimensions)
        M_CHECK_ARRAY_NON_NULL(kd_indices, n_vectors)
        M_CHECK_ARRAY_NON_NULL(density_labels, n_vectors)
        M_CHECK_ARRAY_NON_NULL(tmp_stack, 3 * 64)
        M_CHECK_ARRAY_NON_NULL(tmp_vicinity_mask, n_vectors)
        M_CHECK_ARRAY_NON_NULL(tmp_surface_mask, n_vectors)
        M_CHECK_ARRAY_NON_NULL(tmp_perm, n_vectors)
        M_CHECK_ARRAY_NON_NULL(tmp_abs_diff, n_vectors)
        M_CHECK_ARRAY_NON_NULL(grown_mask, n_vectors)

        call grow_ensemble_expert(&
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
    end subroutine grow_ensemble_expert_c

    !> summary: C-wrapper for [[tox_stc(module):compute_ensemble_observable(subroutine)]]
    subroutine compute_ensemble_observable_c(&
            ensemble_mask,&
            density_labels,&
            n_vectors,&
            candidate_count,&
            current_iter,&
            observables,&
            n_history,&
            n_observables,&
            t_observables,&
            ierr&
        ) bind(C, name="compute_ensemble_observable_c")
        use tox_stc, only: compute_ensemble_observable

        integer(c_int), intent(in), target :: n_vectors
            !! Number of vectors
        integer(c_int), intent(in), target :: n_history
            !! Number of columns in the observable history.
        integer(c_int), intent(in), target :: n_observables
            !! Number of observable rows, which must be exactly five.
            !! The minimum valid value is `5`.
            !! The maximum valid value is `5`.
        logical(c_bool), dimension(n_vectors), intent(in), target :: ensemble_mask
            !! Selection mask for active ensemble members
        real(c_double), dimension(n_vectors), intent(in), target :: density_labels
            !! Density labels for all vectors the density-label routine.
            !! The minimum valid value is `1.0_real64`.
        integer(c_int), intent(in), target :: candidate_count
            !! Number of candidates added in current growth step
            !! The minimum valid value is `0_int32`.
        integer(c_int), intent(in), target :: current_iter
            !! Current growth iteration index (1-based)
            !! The minimum valid value is `1_int32`.
        real(c_double), dimension(n_observables, n_history), intent(inout), target :: observables
            !! 2D observable trajectory matrix where columns match iteration indices
        integer(c_int), intent(in), target :: t_observables
            !! History tracking depth (defaults to 10; <= 0 represents infinity / all)
            !! The default value is `10_int32`.
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_vectors)
        M_CHECK_NON_NULL(candidate_count)
        M_CHECK_NON_NULL(current_iter)
        M_CHECK_NON_NULL(n_history)
        M_CHECK_NON_NULL(n_observables)
        M_CHECK_NON_NULL(t_observables)
        M_CHECK_ARRAY_NON_NULL(ensemble_mask, n_vectors)
        M_CHECK_ARRAY_NON_NULL(density_labels, n_vectors)
        M_CHECK_ARRAY_NON_NULL(observables, n_observables * n_history)

        call compute_ensemble_observable(&
            ensemble_mask = ensemble_mask,&
            density_labels = density_labels,&
            n_vectors = n_vectors,&
            candidate_count = candidate_count,&
            current_iter = current_iter,&
            observables = observables,&
            n_history = n_history,&
            n_observables = n_observables,&
            t_observables = t_observables,&
            ierr = ierr&
        )
    end subroutine compute_ensemble_observable_c

    !> summary: C-wrapper for [[tox_stc(module):accept_ensemble(subroutine)]]
    subroutine accept_ensemble_c(&
            observables,&
            n_history,&
            n_observables,&
            current_iter,&
            alpha_accept,&
            is_accepted,&
            ierr&
        ) bind(C, name="accept_ensemble_c")
        use tox_stc, only: accept_ensemble

        integer(c_int), intent(in), target :: n_history
            !! Number of columns in the observable history.
            !! The minimum valid value is `0_int32`.
        integer(c_int), intent(in), target :: n_observables
            !! Number of observable rows, which must be exactly five.
            !! The minimum valid value is `5`.
            !! The maximum valid value is `5`.
        real(c_double), dimension(n_observables, n_history), intent(in), target :: observables
            !! 2D observables trajectory matrix [5 x window_capacity]
        integer(c_int), intent(in), target :: current_iter
            !! Current growth iteration index (1-based)
            !! The minimum valid value is `1_int32`.
        real(c_double), intent(in), target :: alpha_accept
            !! Optional log2 fold change acceptance threshold (defaults to 0.5)
            !! The minimum valid value is `0.0_real64`.
            !! The default value is `0.5_real64`.
        logical(c_bool), intent(out), target :: is_accepted
            !! Output flag indicating step acceptance decision
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_history)
        M_CHECK_NON_NULL(n_observables)
        M_CHECK_NON_NULL(current_iter)
        M_CHECK_NON_NULL(alpha_accept)
        M_CHECK_NON_NULL(is_accepted)
        M_CHECK_ARRAY_NON_NULL(observables, n_observables * n_history)

        call accept_ensemble(&
            observables = observables,&
            n_history = n_history,&
            n_observables = n_observables,&
            current_iter = current_iter,&
            alpha_accept = alpha_accept,&
            is_accepted = is_accepted,&
            ierr = ierr&
        )
    end subroutine accept_ensemble_c

    !> summary: C-wrapper for [[tox_stc(module):obtain_ensembles_global_radius(subroutine)]]
    subroutine obtain_ensembles_global_radius_c(&
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
            n_ensembles,&
            ierr&
        ) bind(C, name="obtain_ensembles_global_radius_c")
        use tox_stc, only: obtain_ensembles_global_radius

        integer(c_int), intent(in), target :: n_dimensions
            !! Number of dimensions
        integer(c_int), intent(in), target :: n_vectors
            !! Number of vectors
        integer(c_int), intent(in), target :: n_seeds
            !! Total number of precalculated seed vectors the seed-selection routine.
            !! The minimum valid value is `0_int32`.
            !! The maximum valid value is `n_vectors`.
        real(c_double), dimension(n_dimensions, n_vectors), intent(in), target :: vectors
            !! Input data matrix (n_dimensions x n_vectors)
        integer(c_int), dimension(n_dimensions), intent(in), target :: dimension_order
            !! Sequence array tracking split dimensions for KD-tree
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_dimensions`.
        integer(c_int), dimension(n_vectors), intent(in), target :: kd_indices
            !! KD-tree index sequence array computed via [[f42_kd_tree(module):build_kd_index(subroutine)]].
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_vectors`.
        real(c_double), dimension(n_vectors), intent(in), target :: density_labels
            !! Precalculated density labels for ambient vectors the density-label routine.
            !! The minimum valid value is `1.0_real64`.
        integer(c_int), dimension(n_seeds), intent(in), target :: seed_indices
            !! Array of starting seed vector indices the seed-selection routine.
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_vectors`.
        real(c_double), intent(in), target :: r
            !! Uniform growth radius for every seed.
            !! The minimum valid value is `0.0_real64`.
        real(c_double), intent(in), target :: alpha_mad
            !! MAD multiplier for density compatibility threshold
            !! The minimum valid value is `0.0_real64`.
            !! The default value is `0.5_real64`.
        real(c_double), intent(in), target :: alpha_accept
            !! Acceptance log2 fold-change threshold
            !! The minimum valid value is `0.0_real64`.
            !! The default value is `0.5_real64`.
        integer(c_int), intent(in), target :: t_observables
            !! History window in iterations. Pass n_vectors to retain the full history.
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: n_tiles
            !! Number of concurrent seed tiles, at least 1. Capped at the number of seeds when nonzero.
            !! The minimum valid value is `1_int32`.
        logical(c_bool), dimension(n_vectors, n_seeds), intent(out), target :: ensemble_matrix
            !! Output matrix storing unmerged grown ensemble masks
        integer(c_int), dimension(n_seeds), intent(out), target :: stop_reasons
            !! Per-seed stop code. 1: fixed point, 2: rejection after acceptance,
            !! 3: first step rejected, 4: no candidates.
        real(c_double), intent(out), target :: mad_ambient
            !! Ambient MAD of density labels; zero makes alpha_mad ineffective by allowing only exactly equal densities.
        integer(c_int), intent(out), target :: n_ensembles
            !! Total count of extracted raw ensembles (equals n_seeds)
        integer(c_int), intent(out), target :: ierr
            !! Runtime status for the radius selection and density contrast checks.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_dimensions)
        M_CHECK_NON_NULL(n_vectors)
        M_CHECK_NON_NULL(n_seeds)
        M_CHECK_NON_NULL(r)
        M_CHECK_NON_NULL(alpha_mad)
        M_CHECK_NON_NULL(alpha_accept)
        M_CHECK_NON_NULL(t_observables)
        M_CHECK_NON_NULL(n_tiles)
        M_CHECK_NON_NULL(mad_ambient)
        M_CHECK_NON_NULL(n_ensembles)
        M_CHECK_ARRAY_NON_NULL(vectors, n_dimensions * n_vectors)
        M_CHECK_ARRAY_NON_NULL(dimension_order, n_dimensions)
        M_CHECK_ARRAY_NON_NULL(kd_indices, n_vectors)
        M_CHECK_ARRAY_NON_NULL(density_labels, n_vectors)
        M_CHECK_ARRAY_NON_NULL(seed_indices, n_seeds)
        M_CHECK_ARRAY_NON_NULL(ensemble_matrix, n_vectors * n_seeds)
        M_CHECK_ARRAY_NON_NULL(stop_reasons, n_seeds)

        call obtain_ensembles_global_radius(&
            vectors = vectors,&
            n_dimensions = n_dimensions,&
            n_vectors = n_vectors,&
            dimension_order = dimension_order,&
            kd_indices = kd_indices,&
            density_labels = density_labels,&
            seed_indices = seed_indices,&
            n_seeds = n_seeds,&
            r = r,&
            alpha_mad = alpha_mad,&
            alpha_accept = alpha_accept,&
            t_observables = t_observables,&
            n_tiles = n_tiles,&
            ensemble_matrix = ensemble_matrix,&
            stop_reasons = stop_reasons,&
            mad_ambient = mad_ambient,&
            n_ensembles = n_ensembles,&
            ierr = ierr&
        )
    end subroutine obtain_ensembles_global_radius_c

    !> summary: C-wrapper for [[tox_stc(module):obtain_ensembles_global_radius_expert(subroutine)]]
    subroutine obtain_ensembles_global_radius_expert_c(&
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
            tmp_growth_radii,&
            ensemble_matrix,&
            stop_reasons,&
            mad_ambient,&
            n_ensembles,&
            ierr&
        ) bind(C, name="obtain_ensembles_global_radius_expert_c")
        use tox_stc, only: obtain_ensembles_global_radius_expert

        integer(c_int), intent(in), target :: n_dimensions
            !! Number of dimensions
        integer(c_int), intent(in), target :: n_vectors
            !! Number of vectors
        integer(c_int), intent(in), target :: n_seeds
            !! Total number of precalculated seed vectors the seed-selection routine.
            !! The minimum valid value is `0_int32`.
            !! The maximum valid value is `n_vectors`.
        integer(c_int), intent(in), target :: t_observables
            !! History window in iterations. Pass n_vectors to retain the full history.
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: n_tiles
            !! Number of concurrent seed tiles, at least 1. Capped at the number of seeds when nonzero.
            !! The minimum valid value is `1_int32`.
        real(c_double), dimension(n_dimensions, n_vectors), intent(in), target :: vectors
            !! Input data matrix (n_dimensions x n_vectors)
        integer(c_int), dimension(n_dimensions), intent(in), target :: dimension_order
            !! Sequence array tracking split dimensions for KD-tree
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_dimensions`.
        integer(c_int), dimension(n_vectors), intent(in), target :: kd_indices
            !! KD-tree index sequence array computed via [[f42_kd_tree(module):build_kd_index(subroutine)]].
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_vectors`.
        real(c_double), dimension(n_vectors), intent(in), target :: density_labels
            !! Precalculated density labels for ambient vectors the density-label routine.
            !! The minimum valid value is `1.0_real64`.
        integer(c_int), dimension(n_seeds), intent(in), target :: seed_indices
            !! Array of starting seed vector indices the seed-selection routine.
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_vectors`.
        real(c_double), intent(in), target :: r
            !! Uniform growth radius for every seed.
            !! The minimum valid value is `0.0_real64`.
        real(c_double), intent(in), target :: alpha_mad
            !! MAD multiplier for density compatibility threshold
            !! The minimum valid value is `0.0_real64`.
            !! The default value is `0.5_real64`.
        real(c_double), intent(in), target :: alpha_accept
            !! Acceptance log2 fold-change threshold
            !! The minimum valid value is `0.0_real64`.
            !! The default value is `0.5_real64`.
        integer(c_int), dimension(3, 64, min(n_tiles,max(1_int32,n_seeds))), intent(inout), target :: tmp_stack
            !! Workspace stack for tree traversal per tile
        logical(c_bool), dimension(n_vectors, min(n_tiles,max(1_int32,n_seeds))), intent(out), target :: tmp_vicinity_mask
            !! Workspace mask holding one member's neighborhood query result per tile
        logical(c_bool), dimension(n_vectors, min(n_tiles,max(1_int32,n_seeds))), intent(out), target :: tmp_surface_mask
            !! Cached spatial neighborhood union per tile, reset for each seed
        integer(c_int), dimension(n_vectors, min(n_tiles,max(1_int32,n_seeds))), intent(inout), target :: tmp_perm
            !! Sorting workspace, reused for newly accepted member indices between growth rounds per tile
        real(c_double), dimension(n_vectors, min(n_tiles,max(1_int32,n_seeds))), intent(out), target :: tmp_abs_diff
            !! Workspace array for deviation calculations per tile
        real(c_double), dimension(5, min(t_observables,n_vectors), min(n_tiles,max(1_int32,n_seeds))), intent(out), target :: tmp_observables
            !! Workspace matrix for storing observable history per tile
        logical(c_bool), dimension(n_vectors, min(n_tiles,max(1_int32,n_seeds))), intent(out), target :: tmp_current_mask
            !! Workspace holding accepted seed ensemble membership per tile
        real(c_double), dimension(n_seeds), intent(out), target :: tmp_growth_radii
            !! Growth radius in seed order, shared by both modes.
        logical(c_bool), dimension(n_vectors, n_seeds), intent(out), target :: ensemble_matrix
            !! Output matrix storing unmerged grown ensemble masks
        integer(c_int), dimension(n_seeds), intent(out), target :: stop_reasons
            !! Per-seed stop code. 1: fixed point, 2: rejection after acceptance,
            !! 3: first step rejected, 4: no candidates.
        real(c_double), intent(out), target :: mad_ambient
            !! Ambient MAD of density labels; zero makes alpha_mad ineffective by allowing only exactly equal densities.
        integer(c_int), intent(out), target :: n_ensembles
            !! Total count of extracted raw ensembles (equals n_seeds)
        integer(c_int), intent(out), target :: ierr
            !! Runtime status for the radius selection and density contrast checks.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_dimensions)
        M_CHECK_NON_NULL(n_vectors)
        M_CHECK_NON_NULL(n_seeds)
        M_CHECK_NON_NULL(r)
        M_CHECK_NON_NULL(alpha_mad)
        M_CHECK_NON_NULL(alpha_accept)
        M_CHECK_NON_NULL(t_observables)
        M_CHECK_NON_NULL(n_tiles)
        M_CHECK_NON_NULL(mad_ambient)
        M_CHECK_NON_NULL(n_ensembles)
        M_CHECK_ARRAY_NON_NULL(vectors, n_dimensions * n_vectors)
        M_CHECK_ARRAY_NON_NULL(dimension_order, n_dimensions)
        M_CHECK_ARRAY_NON_NULL(kd_indices, n_vectors)
        M_CHECK_ARRAY_NON_NULL(density_labels, n_vectors)
        M_CHECK_ARRAY_NON_NULL(seed_indices, n_seeds)
        M_CHECK_ARRAY_NON_NULL(tmp_stack, 3 * 64 * (min(n_tiles,max(1_int32,n_seeds))))
        M_CHECK_ARRAY_NON_NULL(tmp_vicinity_mask, n_vectors * (min(n_tiles,max(1_int32,n_seeds))))
        M_CHECK_ARRAY_NON_NULL(tmp_surface_mask, n_vectors * (min(n_tiles,max(1_int32,n_seeds))))
        M_CHECK_ARRAY_NON_NULL(tmp_perm, n_vectors * (min(n_tiles,max(1_int32,n_seeds))))
        M_CHECK_ARRAY_NON_NULL(tmp_abs_diff, n_vectors * (min(n_tiles,max(1_int32,n_seeds))))
        M_CHECK_ARRAY_NON_NULL(tmp_observables, 5 * (min(t_observables,n_vectors)) * (min(n_tiles,max(1_int32,n_seeds))))
        M_CHECK_ARRAY_NON_NULL(tmp_current_mask, n_vectors * (min(n_tiles,max(1_int32,n_seeds))))
        M_CHECK_ARRAY_NON_NULL(tmp_growth_radii, n_seeds)
        M_CHECK_ARRAY_NON_NULL(ensemble_matrix, n_vectors * n_seeds)
        M_CHECK_ARRAY_NON_NULL(stop_reasons, n_seeds)

        call obtain_ensembles_global_radius_expert(&
            vectors = vectors,&
            n_dimensions = n_dimensions,&
            n_vectors = n_vectors,&
            dimension_order = dimension_order,&
            kd_indices = kd_indices,&
            density_labels = density_labels,&
            seed_indices = seed_indices,&
            n_seeds = n_seeds,&
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
            tmp_growth_radii = tmp_growth_radii,&
            ensemble_matrix = ensemble_matrix,&
            stop_reasons = stop_reasons,&
            mad_ambient = mad_ambient,&
            n_ensembles = n_ensembles,&
            ierr = ierr&
        )
    end subroutine obtain_ensembles_global_radius_expert_c

    !> summary: C-wrapper for [[tox_stc(module):obtain_ensembles_seed_radius(subroutine)]]
    subroutine obtain_ensembles_seed_radius_c(&
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
            n_ensembles,&
            ierr&
        ) bind(C, name="obtain_ensembles_seed_radius_c")
        use tox_stc, only: obtain_ensembles_seed_radius

        integer(c_int), intent(in), target :: n_dimensions
            !! Number of dimensions
        integer(c_int), intent(in), target :: n_vectors
            !! Number of vectors
        integer(c_int), intent(in), target :: n_seeds
            !! Total number of precalculated seed vectors the seed-selection routine.
            !! The minimum valid value is `0_int32`.
            !! The maximum valid value is `n_vectors`.
        real(c_double), dimension(n_dimensions, n_vectors), intent(in), target :: vectors
            !! Input data matrix (n_dimensions x n_vectors)
        integer(c_int), dimension(n_dimensions), intent(in), target :: dimension_order
            !! Sequence array tracking split dimensions for KD-tree
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_dimensions`.
        integer(c_int), dimension(n_vectors), intent(in), target :: kd_indices
            !! KD-tree index sequence array computed via [[f42_kd_tree(module):build_kd_index(subroutine)]].
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_vectors`.
        real(c_double), dimension(n_vectors), intent(in), target :: density_labels
            !! Precalculated density labels for ambient vectors the density-label routine.
            !! The minimum valid value is `1.0_real64`.
        integer(c_int), dimension(n_seeds), intent(in), target :: seed_indices
            !! Array of starting seed vector indices the seed-selection routine.
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_vectors`.
        real(c_double), dimension(n_seeds), intent(in), target :: seed_radii
            !! Per-seed growth radius, allowing tighter growth in dense regions and wider growth in sparse regions.
            !! The minimum valid value is `0.0_real64`.
        real(c_double), intent(in), target :: alpha_mad
            !! MAD multiplier for density compatibility threshold
            !! The minimum valid value is `0.0_real64`.
            !! The default value is `0.5_real64`.
        real(c_double), intent(in), target :: alpha_accept
            !! Acceptance log2 fold-change threshold
            !! The minimum valid value is `0.0_real64`.
            !! The default value is `0.5_real64`.
        integer(c_int), intent(in), target :: t_observables
            !! History window in iterations. Pass n_vectors to retain the full history.
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: n_tiles
            !! Number of concurrent seed tiles, at least 1. Capped at the number of seeds when nonzero.
            !! The minimum valid value is `1_int32`.
        logical(c_bool), dimension(n_vectors, n_seeds), intent(out), target :: ensemble_matrix
            !! Output matrix storing unmerged grown ensemble masks
        integer(c_int), dimension(n_seeds), intent(out), target :: stop_reasons
            !! Per-seed stop code. 1: fixed point, 2: rejection after acceptance,
            !! 3: first step rejected, 4: no candidates.
        real(c_double), intent(out), target :: mad_ambient
            !! Ambient MAD of density labels; zero makes alpha_mad ineffective by allowing only exactly equal densities.
        integer(c_int), intent(out), target :: n_ensembles
            !! Total count of extracted raw ensembles (equals n_seeds)
        integer(c_int), intent(out), target :: ierr
            !! Runtime status for the radius selection and density contrast checks.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_dimensions)
        M_CHECK_NON_NULL(n_vectors)
        M_CHECK_NON_NULL(n_seeds)
        M_CHECK_NON_NULL(alpha_mad)
        M_CHECK_NON_NULL(alpha_accept)
        M_CHECK_NON_NULL(t_observables)
        M_CHECK_NON_NULL(n_tiles)
        M_CHECK_NON_NULL(mad_ambient)
        M_CHECK_NON_NULL(n_ensembles)
        M_CHECK_ARRAY_NON_NULL(vectors, n_dimensions * n_vectors)
        M_CHECK_ARRAY_NON_NULL(dimension_order, n_dimensions)
        M_CHECK_ARRAY_NON_NULL(kd_indices, n_vectors)
        M_CHECK_ARRAY_NON_NULL(density_labels, n_vectors)
        M_CHECK_ARRAY_NON_NULL(seed_indices, n_seeds)
        M_CHECK_ARRAY_NON_NULL(seed_radii, n_seeds)
        M_CHECK_ARRAY_NON_NULL(ensemble_matrix, n_vectors * n_seeds)
        M_CHECK_ARRAY_NON_NULL(stop_reasons, n_seeds)

        call obtain_ensembles_seed_radius(&
            vectors = vectors,&
            n_dimensions = n_dimensions,&
            n_vectors = n_vectors,&
            dimension_order = dimension_order,&
            kd_indices = kd_indices,&
            density_labels = density_labels,&
            seed_indices = seed_indices,&
            n_seeds = n_seeds,&
            seed_radii = seed_radii,&
            alpha_mad = alpha_mad,&
            alpha_accept = alpha_accept,&
            t_observables = t_observables,&
            n_tiles = n_tiles,&
            ensemble_matrix = ensemble_matrix,&
            stop_reasons = stop_reasons,&
            mad_ambient = mad_ambient,&
            n_ensembles = n_ensembles,&
            ierr = ierr&
        )
    end subroutine obtain_ensembles_seed_radius_c

    !> summary: C-wrapper for [[tox_stc(module):obtain_ensembles_seed_radius_expert(subroutine)]]
    subroutine obtain_ensembles_seed_radius_expert_c(&
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
            tmp_growth_radii,&
            ensemble_matrix,&
            stop_reasons,&
            mad_ambient,&
            n_ensembles,&
            ierr&
        ) bind(C, name="obtain_ensembles_seed_radius_expert_c")
        use tox_stc, only: obtain_ensembles_seed_radius_expert

        integer(c_int), intent(in), target :: n_dimensions
            !! Number of dimensions
        integer(c_int), intent(in), target :: n_vectors
            !! Number of vectors
        integer(c_int), intent(in), target :: n_seeds
            !! Total number of precalculated seed vectors the seed-selection routine.
            !! The minimum valid value is `0_int32`.
            !! The maximum valid value is `n_vectors`.
        integer(c_int), intent(in), target :: t_observables
            !! History window in iterations. Pass n_vectors to retain the full history.
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: n_tiles
            !! Number of concurrent seed tiles, at least 1. Capped at the number of seeds when nonzero.
            !! The minimum valid value is `1_int32`.
        real(c_double), dimension(n_dimensions, n_vectors), intent(in), target :: vectors
            !! Input data matrix (n_dimensions x n_vectors)
        integer(c_int), dimension(n_dimensions), intent(in), target :: dimension_order
            !! Sequence array tracking split dimensions for KD-tree
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_dimensions`.
        integer(c_int), dimension(n_vectors), intent(in), target :: kd_indices
            !! KD-tree index sequence array computed via [[f42_kd_tree(module):build_kd_index(subroutine)]].
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_vectors`.
        real(c_double), dimension(n_vectors), intent(in), target :: density_labels
            !! Precalculated density labels for ambient vectors the density-label routine.
            !! The minimum valid value is `1.0_real64`.
        integer(c_int), dimension(n_seeds), intent(in), target :: seed_indices
            !! Array of starting seed vector indices the seed-selection routine.
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_vectors`.
        real(c_double), dimension(n_seeds), intent(in), target :: seed_radii
            !! Per-seed growth radius, allowing tighter growth in dense regions and wider growth in sparse regions.
            !! The minimum valid value is `0.0_real64`.
        real(c_double), intent(in), target :: alpha_mad
            !! MAD multiplier for density compatibility threshold
            !! The minimum valid value is `0.0_real64`.
            !! The default value is `0.5_real64`.
        real(c_double), intent(in), target :: alpha_accept
            !! Acceptance log2 fold-change threshold
            !! The minimum valid value is `0.0_real64`.
            !! The default value is `0.5_real64`.
        integer(c_int), dimension(3, 64, min(n_tiles,max(1_int32,n_seeds))), intent(inout), target :: tmp_stack
            !! Workspace stack for tree traversal per tile
        logical(c_bool), dimension(n_vectors, min(n_tiles,max(1_int32,n_seeds))), intent(out), target :: tmp_vicinity_mask
            !! Workspace mask holding one member's neighborhood query result per tile
        logical(c_bool), dimension(n_vectors, min(n_tiles,max(1_int32,n_seeds))), intent(out), target :: tmp_surface_mask
            !! Cached spatial neighborhood union per tile, reset for each seed
        integer(c_int), dimension(n_vectors, min(n_tiles,max(1_int32,n_seeds))), intent(inout), target :: tmp_perm
            !! Sorting workspace, reused for newly accepted member indices between growth rounds per tile
        real(c_double), dimension(n_vectors, min(n_tiles,max(1_int32,n_seeds))), intent(out), target :: tmp_abs_diff
            !! Workspace array for deviation calculations per tile
        real(c_double), dimension(5, min(t_observables,n_vectors), min(n_tiles,max(1_int32,n_seeds))), intent(out), target :: tmp_observables
            !! Workspace matrix for storing observable history per tile
        logical(c_bool), dimension(n_vectors, min(n_tiles,max(1_int32,n_seeds))), intent(out), target :: tmp_current_mask
            !! Workspace holding accepted seed ensemble membership per tile
        real(c_double), dimension(n_seeds), intent(out), target :: tmp_growth_radii
            !! Growth radius in seed order, shared by both modes.
        logical(c_bool), dimension(n_vectors, n_seeds), intent(out), target :: ensemble_matrix
            !! Output matrix storing unmerged grown ensemble masks
        integer(c_int), dimension(n_seeds), intent(out), target :: stop_reasons
            !! Per-seed stop code. 1: fixed point, 2: rejection after acceptance,
            !! 3: first step rejected, 4: no candidates.
        real(c_double), intent(out), target :: mad_ambient
            !! Ambient MAD of density labels; zero makes alpha_mad ineffective by allowing only exactly equal densities.
        integer(c_int), intent(out), target :: n_ensembles
            !! Total count of extracted raw ensembles (equals n_seeds)
        integer(c_int), intent(out), target :: ierr
            !! Runtime status for the radius selection and density contrast checks.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_dimensions)
        M_CHECK_NON_NULL(n_vectors)
        M_CHECK_NON_NULL(n_seeds)
        M_CHECK_NON_NULL(alpha_mad)
        M_CHECK_NON_NULL(alpha_accept)
        M_CHECK_NON_NULL(t_observables)
        M_CHECK_NON_NULL(n_tiles)
        M_CHECK_NON_NULL(mad_ambient)
        M_CHECK_NON_NULL(n_ensembles)
        M_CHECK_ARRAY_NON_NULL(vectors, n_dimensions * n_vectors)
        M_CHECK_ARRAY_NON_NULL(dimension_order, n_dimensions)
        M_CHECK_ARRAY_NON_NULL(kd_indices, n_vectors)
        M_CHECK_ARRAY_NON_NULL(density_labels, n_vectors)
        M_CHECK_ARRAY_NON_NULL(seed_indices, n_seeds)
        M_CHECK_ARRAY_NON_NULL(seed_radii, n_seeds)
        M_CHECK_ARRAY_NON_NULL(tmp_stack, 3 * 64 * (min(n_tiles,max(1_int32,n_seeds))))
        M_CHECK_ARRAY_NON_NULL(tmp_vicinity_mask, n_vectors * (min(n_tiles,max(1_int32,n_seeds))))
        M_CHECK_ARRAY_NON_NULL(tmp_surface_mask, n_vectors * (min(n_tiles,max(1_int32,n_seeds))))
        M_CHECK_ARRAY_NON_NULL(tmp_perm, n_vectors * (min(n_tiles,max(1_int32,n_seeds))))
        M_CHECK_ARRAY_NON_NULL(tmp_abs_diff, n_vectors * (min(n_tiles,max(1_int32,n_seeds))))
        M_CHECK_ARRAY_NON_NULL(tmp_observables, 5 * (min(t_observables,n_vectors)) * (min(n_tiles,max(1_int32,n_seeds))))
        M_CHECK_ARRAY_NON_NULL(tmp_current_mask, n_vectors * (min(n_tiles,max(1_int32,n_seeds))))
        M_CHECK_ARRAY_NON_NULL(tmp_growth_radii, n_seeds)
        M_CHECK_ARRAY_NON_NULL(ensemble_matrix, n_vectors * n_seeds)
        M_CHECK_ARRAY_NON_NULL(stop_reasons, n_seeds)

        call obtain_ensembles_seed_radius_expert(&
            vectors = vectors,&
            n_dimensions = n_dimensions,&
            n_vectors = n_vectors,&
            dimension_order = dimension_order,&
            kd_indices = kd_indices,&
            density_labels = density_labels,&
            seed_indices = seed_indices,&
            n_seeds = n_seeds,&
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
            tmp_growth_radii = tmp_growth_radii,&
            ensemble_matrix = ensemble_matrix,&
            stop_reasons = stop_reasons,&
            mad_ambient = mad_ambient,&
            n_ensembles = n_ensembles,&
            ierr = ierr&
        )
    end subroutine obtain_ensembles_seed_radius_expert_c

end module tox_stc_c
#endif
