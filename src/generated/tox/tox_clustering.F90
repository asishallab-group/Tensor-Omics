#include <src/macros.h>

!> Clustering for tensor omics: k-means over factors and trajectories, and hierarchical
!| (agglomerative) linkage clustering over a precomputed distance matrix.
!|
!| `k_means_clustering` partitions points in n dimensions; `cluster_factor_trajectories_k_means`
!| applies the same to whole time series, treating each factor's trajectory as one point.
!| `linkage_clustering` takes the distances already computed and merges under the linkage
!| criterion asked for, so the same matrix can be re-clustered without recomputing it.
!|
!| Generated from [[tox_clustering_impl(module)]]; do not edit -- regenerate instead.
module tox_clustering
    use tox_clustering_impl, only: METHOD_AVERAGE, METHOD_WARD, METHOD_WEIGHTED, cluster_factor_trajectories_k_means_impl
    use tox_clustering_impl, only: k_means_clustering_impl, linkage_clustering_impl
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use tox_errors, only: set_ok, is_err, ERR_INVALID_INPUT, set_err_once
    use tox_errors, only: validate_all_in_range_real, validate_dimension_size, validate_distance_matrix, validate_in_range_int
    M_IMPLICIT_NONE
    private

    public :: cluster_factor_trajectories_k_means
    public :: k_means_clustering
    public :: linkage_clustering

contains

    !> summary: Validates its inputs, then calls [[tox_clustering_impl(module):cluster_factor_trajectories_k_means_impl]].
    pure subroutine cluster_factor_trajectories_k_means(&
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
        )
        integer(int32), intent(in) :: n_clusters
            !! number (`k`) of clusters
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_samples*n_timepoints`.
        integer(int32), intent(in) :: n_factors
            !! number of factors
        integer(int32), intent(in) :: n_samples
            !! number of samples
        integer(int32), intent(in) :: n_timepoints
            !! number of timepoints
        real(real64), dimension(n_factors, n_samples, n_timepoints), intent(in) :: trajectories
            !! matrix with data points to cluster
        real(real64), dimension(n_factors, n_clusters), intent(inout) :: centroids
            !! matrix with initial centroids of the clusters, could be random data or actual points or unassigned garbage.
            !! The centroids should be unique. This is not checked in this routine.
            !!
            !! The final values will be the final centroids of the clusters
        integer(int32), dimension(n_samples*n_timepoints), intent(out) :: labels
            !! array of labels, each index corresponds to the respective point's index, so first label is first point's label.
            !!
            !! each label is the index of its related cluster -> `1<=label<=n_clusters=k`
        integer(int32), dimension(n_clusters), intent(out) :: label_counts
            !! holds the number of points having the respective label assigned
        integer(int32), intent(in) :: max_iterations
            !! number of maximum iterations of the clustering
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_in_range_int(n_clusters, ierr, arg_pos=1_int32, min=1_int32, max=n_samples*n_timepoints)
        call validate_dimension_size(n_factors, ierr, arg_pos=3_int32)
        call validate_dimension_size(n_samples, ierr, arg_pos=4_int32)
        call validate_dimension_size(n_timepoints, ierr, arg_pos=5_int32)
        call validate_all_in_range_real(trajectories, n_factors * n_samples * n_timepoints, ierr, arg_pos=2_int32)
        call validate_all_in_range_real(centroids, n_factors * n_clusters, ierr, arg_pos=6_int32)
        if (is_err(ierr)) return
#endif

        call cluster_factor_trajectories_k_means_impl(&
            n_clusters = n_clusters,&
            trajectories = trajectories,&
            n_factors = n_factors,&
            n_samples = n_samples,&
            n_timepoints = n_timepoints,&
            centroids = centroids,&
            labels = labels,&
            label_counts = label_counts,&
            max_iterations = max_iterations&
        )
    end subroutine cluster_factor_trajectories_k_means

    !> summary: Validates its inputs, then calls [[tox_clustering_impl(module):k_means_clustering_impl]].
    !| 1. Assigns each data point to one of `k` clusters whose centroid is clostest
    !| 2. Recalculates the centroids using the mean of its assigned points
    !| 3. repeat 1-2 until assignment remains unchanged
    pure subroutine k_means_clustering(&
            n_clusters,&
            data_points,&
            n_points,&
            n_dims,&
            centroids,&
            labels,&
            label_counts,&
            max_iterations,&
            ierr&
        )
        integer(int32), intent(in) :: n_clusters
            !! number (`k`) of clusters
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_points`.
        integer(int32), intent(in) :: n_points
            !! number of points to cluster
        integer(int32), intent(in) :: n_dims
            !! number of elements a point has
        real(real64), dimension(n_dims, n_points), intent(in) :: data_points
            !! matrix with data points to cluster
        real(real64), dimension(n_dims, n_clusters), intent(inout) :: centroids
            !! matrix with initial centroids of the clusters, could be random data or actual points or unassigned garbage.
            !! The centroids should be unique. This is not checked in this routine.
            !!
            !! The final values will be the final centroids of the clusters
        integer(int32), dimension(n_points), intent(out) :: labels
            !! array of labels, each index corresponds to the respective point's index, so first label is first point's label.
            !!
            !! each label is the index of its related cluster -> `1<=label<=n_clusters=k`
        integer(int32), dimension(n_clusters), intent(out) :: label_counts
            !! holds the number of points having the respective label assigned
        integer(int32), intent(in), optional :: max_iterations
            !! number of maximum iterations of the clustering.
            !! The default value is `300_int32`.
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_in_range_int(n_clusters, ierr, arg_pos=1_int32, min=1_int32, max=n_points)
        call validate_dimension_size(n_points, ierr, arg_pos=3_int32)
        call validate_dimension_size(n_dims, ierr, arg_pos=4_int32)
        call validate_all_in_range_real(data_points, n_dims * n_points, ierr, arg_pos=2_int32)
        call validate_all_in_range_real(centroids, n_dims * n_clusters, ierr, arg_pos=5_int32)
        if (is_err(ierr)) return
#endif

        call k_means_clustering_impl(&
            n_clusters = n_clusters,&
            data_points = data_points,&
            n_points = n_points,&
            n_dims = n_dims,&
            centroids = centroids,&
            labels = labels,&
            label_counts = label_counts,&
            max_iterations = max_iterations&
        )
    end subroutine k_means_clustering

    !> summary: Validates its inputs, then calls [[tox_clustering_impl(module):linkage_clustering_impl]].
    !| @note
    !| The bottom triangle is used as scratch and restored from the top triangle before
    !| returning, on success or on error, so the matrix comes back unchanged. There is no
    !| need to copy it before calling.
    !| @endnote
    pure subroutine linkage_clustering(&
            distances,&
            n_points,&
            merge_i,&
            merge_j,&
            heights,&
            cluster_sizes,&
            method,&
            ierr&
        )
        integer(int32), intent(in) :: n_points
            !! number of points to cluster
        real(real64), dimension(n_points, n_points), intent(inout) :: distances
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
        integer(int32), dimension(n_points - 1), intent(out) :: merge_i
            !! holds cluster labels of the merged node pair at iteration k -> positives relate to leafs/data point indices, negatives to inner nodes
        integer(int32), dimension(n_points - 1), intent(out) :: merge_j
            !! holds cluster labels of the merged node pair at iteration k -> positives relate to leafs/data point indices, negatives to inner nodes
        real(real64), dimension(n_points - 1), intent(out) :: heights
            !! height of the shorter branch of the merge, e.g. if (A,B)+(C) merges to ((A,B),C), the branch to (A,B) is shorter
        integer(int32), dimension(n_points - 1), intent(out) :: cluster_sizes
            !! size of cluster at iteration k
        integer(int32), intent(in) :: method
            !! used algorithm
            !! The minimum valid value is `0_int32`.
            !! The maximum valid value is `2_int32`.
            !!
            !! | Method           | Value                                                     |
            !! |------------------|-----------------------------------------------------------|
            !! | Average / UPGMA  | [[tox_clustering_impl(module):METHOD_AVERAGE(variable)]]  |
            !! | Weighted / WPGMA | [[tox_clustering_impl(module):METHOD_WEIGHTED(variable)]] |
            !! | Ward             | [[tox_clustering_impl(module):METHOD_WARD(variable)]]     |
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_points, ierr, arg_pos=2_int32)
        call validate_in_range_int(method, ierr, arg_pos=7_int32, min=0_int32, max=2_int32)
        call validate_distance_matrix(distances, n_points, ierr, arg_pos=1_int32)
        if (method /= METHOD_AVERAGE .and. method /= METHOD_WEIGHTED .and. method /= METHOD_WARD) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=7_int32)
        if (is_err(ierr)) return
#endif

        call linkage_clustering_impl(&
            distances = distances,&
            n_points = n_points,&
            merge_i = merge_i,&
            merge_j = merge_j,&
            heights = heights,&
            cluster_sizes = cluster_sizes,&
            method = method&
        )
    end subroutine linkage_clustering

end module tox_clustering
