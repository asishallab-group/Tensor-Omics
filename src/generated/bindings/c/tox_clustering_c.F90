#ifndef NO_C_BINDING
#include <src/macros.h>

!> summary: C-wrappers for [[tox_clustering(module)]]
!| Generated from the implementation; do not edit -- regenerate instead.
module tox_clustering_c
    use f42_safeguard
    use, intrinsic :: iso_c_binding, only: c_associated, c_char, c_double, c_int, c_loc
    use tox_conversions, only: c_char_1d_as_string
    use tox_errors, only: set_ok, set_err, is_err, ERR_POINTER_NULL, ERR_INVALID_INPUT
    M_IMPLICIT_NONE
    private

    public :: cluster_factor_trajectories_k_means_c
    public :: k_means_clustering_c
    public :: linkage_clustering_c

contains

    !> summary: C-wrapper for [[tox_clustering(module):cluster_factor_trajectories_k_means(subroutine)]]
    subroutine cluster_factor_trajectories_k_means_c(&
            n_clusters,&
            trajectories,&
            n_factors,&
            n_samples,&
            n_timepoints,&
            centroids,&
            labels,&
            label_counts,&
            max_iterations,&
            ierr&
        ) bind(C, name="cluster_factor_trajectories_k_means_c")
        use tox_clustering, only: cluster_factor_trajectories_k_means

        integer(c_int), intent(in), target :: n_clusters
            !! number (`k`) of clusters
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_samples*n_timepoints`.
        integer(c_int), intent(in), target :: n_factors
            !! number of factors
        integer(c_int), intent(in), target :: n_samples
            !! number of samples
        integer(c_int), intent(in), target :: n_timepoints
            !! number of timepoints
        real(c_double), dimension(n_factors, n_samples, n_timepoints), intent(in), target :: trajectories
            !! matrix with data points to cluster
        real(c_double), dimension(n_factors, n_clusters), intent(inout), target :: centroids
            !! matrix with initial centroids of the clusters, could be random data or actual points or unassigned garbage.
            !! The centroids should be unique. This is not checked in this routine.
            !!
            !! The final values will be the final centroids of the clusters
        integer(c_int), dimension(n_samples*n_timepoints), intent(out), target :: labels
            !! array of labels, each index corresponds to the respective point's index, so first label is first point's label.
            !!
            !! each label is the index of its related cluster -> `1<=label<=n_clusters=k`
        integer(c_int), dimension(n_clusters), intent(out), target :: label_counts
            !! holds the number of points having the respective label assigned
        integer(c_int), intent(in), target :: max_iterations
            !! number of maximum iterations of the clustering
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_clusters)
        M_CHECK_NON_NULL(n_factors)
        M_CHECK_NON_NULL(n_samples)
        M_CHECK_NON_NULL(n_timepoints)
        M_CHECK_NON_NULL(max_iterations)
        M_CHECK_ARRAY_NON_NULL(trajectories, n_factors * n_samples * n_timepoints)
        M_CHECK_ARRAY_NON_NULL(centroids, n_factors * n_clusters)
        M_CHECK_ARRAY_NON_NULL(labels, (n_samples*n_timepoints))
        M_CHECK_ARRAY_NON_NULL(label_counts, n_clusters)

        call cluster_factor_trajectories_k_means(&
            n_clusters = n_clusters,&
            trajectories = trajectories,&
            n_factors = n_factors,&
            n_samples = n_samples,&
            n_timepoints = n_timepoints,&
            centroids = centroids,&
            labels = labels,&
            label_counts = label_counts,&
            max_iterations = max_iterations,&
            ierr = ierr&
        )
    end subroutine cluster_factor_trajectories_k_means_c

    !> summary: C-wrapper for [[tox_clustering(module):k_means_clustering(subroutine)]]
    !| 1. Assigns each data point to one of `k` clusters whose centroid is clostest
    !| 2. Recalculates the centroids using the mean of its assigned points
    !| 3. repeat 1-2 until assignment remains unchanged
    subroutine k_means_clustering_c(&
            n_clusters,&
            data_points,&
            n_points,&
            n_dims,&
            centroids,&
            labels,&
            label_counts,&
            max_iterations,&
            ierr&
        ) bind(C, name="k_means_clustering_c")
        use tox_clustering, only: k_means_clustering

        integer(c_int), intent(in), target :: n_clusters
            !! number (`k`) of clusters
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_points`.
        integer(c_int), intent(in), target :: n_points
            !! number of points to cluster
        integer(c_int), intent(in), target :: n_dims
            !! number of elements a point has
        real(c_double), dimension(n_dims, n_points), intent(in), target :: data_points
            !! matrix with data points to cluster
        real(c_double), dimension(n_dims, n_clusters), intent(inout), target :: centroids
            !! matrix with initial centroids of the clusters, could be random data or actual points or unassigned garbage.
            !! The centroids should be unique. This is not checked in this routine.
            !!
            !! The final values will be the final centroids of the clusters
        integer(c_int), dimension(n_points), intent(out), target :: labels
            !! array of labels, each index corresponds to the respective point's index, so first label is first point's label.
            !!
            !! each label is the index of its related cluster -> `1<=label<=n_clusters=k`
        integer(c_int), dimension(n_clusters), intent(out), target :: label_counts
            !! holds the number of points having the respective label assigned
        integer(c_int), intent(in), target :: max_iterations
            !! number of maximum iterations of the clustering.
            !! The default value is `300_int32`.
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_clusters)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(n_dims)
        M_CHECK_NON_NULL(max_iterations)
        M_CHECK_ARRAY_NON_NULL(data_points, n_dims * n_points)
        M_CHECK_ARRAY_NON_NULL(centroids, n_dims * n_clusters)
        M_CHECK_ARRAY_NON_NULL(labels, n_points)
        M_CHECK_ARRAY_NON_NULL(label_counts, n_clusters)

        call k_means_clustering(&
            n_clusters = n_clusters,&
            data_points = data_points,&
            n_points = n_points,&
            n_dims = n_dims,&
            centroids = centroids,&
            labels = labels,&
            label_counts = label_counts,&
            max_iterations = max_iterations,&
            ierr = ierr&
        )
    end subroutine k_means_clustering_c

    !> summary: C-wrapper for [[tox_clustering(module):linkage_clustering(subroutine)]]
    !| @note
    !| The bottom triangle is used as scratch and restored from the top triangle before
    !| returning, on success or on error, so the matrix comes back unchanged. There is no
    !| need to copy it before calling.
    !| @endnote
    subroutine linkage_clustering_c(&
            distances,&
            n_points,&
            merge_i,&
            merge_j,&
            heights,&
            cluster_sizes,&
            method,&
            ierr&
        ) bind(C, name="linkage_clustering_c")
        use tox_clustering, only: linkage_clustering
        use tox_clustering_impl, only: METHOD_AVERAGE, METHOD_WARD, METHOD_WEIGHTED

        integer(c_int), intent(in), target :: n_points
            !! number of points to cluster
        real(c_double), dimension(n_points, n_points), intent(inout), target :: distances
            !! symmetric distance matrix, holding the positive distances between points. Distance of X->X is always zero.
            !!
            !! @note
            !! The bottom triangle is used as scratch and restored from the top triangle before
            !! returning, on success or on error, so the matrix comes back unchanged. There is no
            !! need to copy it before calling.
            !! @endnote
            !!
            !! Its structure (symmetry, non-negativity, zero diagonal) is validated by the
            !! distance-matrix naming convention in the generated wrapper.
        integer(c_int), dimension(n_points - 1), intent(out), target :: merge_i
            !! holds cluster labels of the merged node pair at iteration k -> positives relate to leafs/data point indices, negatives to inner nodes
        integer(c_int), dimension(n_points - 1), intent(out), target :: merge_j
            !! holds cluster labels of the merged node pair at iteration k -> positives relate to leafs/data point indices, negatives to inner nodes
        real(c_double), dimension(n_points - 1), intent(out), target :: heights
            !! height of the shorter branch of the merge, e.g. if (A,B)+(C) merges to ((A,B),C), the branch to (A,B) is shorter
        integer(c_int), dimension(n_points - 1), intent(out), target :: cluster_sizes
            !! size of cluster at iteration k
        character(len=1, kind=c_char), dimension(8), intent(in), target :: method
            !! used algorithm
            !! The minimum valid value is `0_int32`.
            !! The maximum valid value is `2_int32`.
            !!
            !! | Method           | Value                                                     |
            !! |------------------|-----------------------------------------------------------|
            !! | Average / UPGMA  | [[tox_clustering_impl(module):METHOD_AVERAGE(variable)]]  |
            !! | Weighted / WPGMA | [[tox_clustering_impl(module):METHOD_WEIGHTED(variable)]] |
            !! | Ward             | [[tox_clustering_impl(module):METHOD_WARD(variable)]]     |
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.
        integer(int32) :: method_mode_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_ARRAY_NON_NULL(distances, n_points * n_points)
        M_CHECK_ARRAY_NON_NULL(merge_i, (n_points - 1))
        M_CHECK_ARRAY_NON_NULL(merge_j, (n_points - 1))
        M_CHECK_ARRAY_NON_NULL(heights, (n_points - 1))
        M_CHECK_ARRAY_NON_NULL(cluster_sizes, (n_points - 1))
        M_CHECK_ARRAY_NON_NULL(method, 8)

        block
            character(len=:), allocatable :: method_f
            call c_char_1d_as_string(method, method_f, ierr)
            if (is_err(ierr)) return

            select case (method_f)
                case ("average")
                    method_mode_f = METHOD_AVERAGE
                case ("weighted")
                    method_mode_f = METHOD_WEIGHTED
                case ("ward")
                    method_mode_f = METHOD_WARD
                case default
                    call set_err(ierr, ERR_INVALID_INPUT)
                    return
            end select
        end block

        call linkage_clustering(&
            distances = distances,&
            n_points = n_points,&
            merge_i = merge_i,&
            merge_j = merge_j,&
            heights = heights,&
            cluster_sizes = cluster_sizes,&
            method = method_mode_f,&
            ierr = ierr&
        )
    end subroutine linkage_clustering_c

end module tox_clustering_c
#endif
