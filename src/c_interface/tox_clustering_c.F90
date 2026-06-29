#ifndef NO_C_INTERFACE
#include <src/macros.h>

!> summary: Module for C-wrappers for [[tox_clustering(module)]]
module tox_clustering_c
    use safeguard
    use, intrinsic :: iso_c_binding, only: c_int, c_double, c_char, c_double_complex
    use, intrinsic :: iso_c_binding, only: c_loc, c_associated

    use tox_conversions, only: logical_as_c_int, c_int_as_logical
    use tox_conversions, only: c_char_as_char, char_as_c_char
    use tox_conversions, only: string_as_c_char_1d, c_char_1d_as_string
    use tox_conversions, only: string_as_c_char_2d, c_char_2d_as_string

    use tox_errors, only: ERR_POINTER_NULL, is_err, set_err, ERR_ALLOC_FAIL, ERR_INVALID_INPUT
    implicit none
contains

    !> summary: C-wrapper for [[tox_clustering(module):cluster_factor_trajectories_k_means(subroutine)]]
    !| Performs k-means clustering on factor trajectories, so factor evolution over time
    subroutine cluster_factor_trajectories_k_means_c(&
            n_clusters,&
            trajectories,&
            n_factors,&
            n_samples,&
            n_timepoints,&
            centroids,&
            labels,&
            label_counts,&
            ierr,&
            max_iterations&
            ) bind(C, name="cluster_factor_trajectories_k_means_c")
        use tox_clustering, only: cluster_factor_trajectories_k_means
        use tox_clustering
        integer(c_int), intent(in), target :: n_clusters
            !! number (`k`) of clusters
        integer(c_int), intent(in), target :: n_factors
            !! number of factors
        integer(c_int), intent(in), target :: n_samples
            !! number of samples
        integer(c_int), intent(in), target :: n_timepoints
            !! number of timepoints
        real(c_double), intent(in), dimension(n_factors, n_samples, n_timepoints), target :: trajectories
            !! matrix with data points to cluster
        real(c_double), intent(inout), dimension(n_factors, n_clusters), target :: centroids
            !! matrix with initial centroids of the clusters, could be random data or actual points or unassigned garbage.
            !! The centroids should be unique. This is not checked in this routine.
            !! 
            !! The final values will be the final centroids of the clusters
        integer(c_int), intent(out), dimension(n_samples * n_timepoints), target :: labels
            !! array of labels, each index corresponds to the respective point's index, so first label is first point's label.
            !! 
            !! each label is the index of its related cluster -> `1<=label<=n_clusters=k`
        integer(c_int), intent(out), dimension(n_clusters), target :: label_counts
            !! holds the number of points having the respective label assigned
        integer(c_int), intent(out), target :: ierr
            !! Error code
        integer(c_int), intent(in), target :: max_iterations
            !! number of maximum iterations of the clustering
            !! The default value is `300_int32`.
        M_CHECK_IERR_NON_NULL
        M_CHECK_NON_NULL(n_clusters)
        M_CHECK_NON_NULL(trajectories)
        M_CHECK_NON_NULL(n_factors)
        M_CHECK_NON_NULL(n_samples)
        M_CHECK_NON_NULL(n_timepoints)
        M_CHECK_NON_NULL(centroids)
        M_CHECK_NON_NULL(labels)
        M_CHECK_NON_NULL(label_counts)
        call cluster_factor_trajectories_k_means(&
            n_clusters = n_clusters,&
            trajectories = trajectories,&
            n_factors = n_factors,&
            n_samples = n_samples,&
            n_timepoints = n_timepoints,&
            centroids = centroids,&
            labels = labels,&
            label_counts = label_counts,&
            ierr = ierr&
        )
    end subroutine cluster_factor_trajectories_k_means_c

    !> summary: C-wrapper for [[tox_clustering(module):k_means_clustering(subroutine)]]
    !| k-means clustering algorithm:
    !| 
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
            ierr,&
            max_iterations&
            ) bind(C, name="k_means_clustering_c")
        use tox_clustering, only: k_means_clustering
        use tox_clustering
        integer(c_int), intent(in), target :: n_clusters
            !! number (`k`) of clusters
        integer(c_int), intent(in), target :: n_points
            !! number of points to cluster
        integer(c_int), intent(in), target :: n_dims
            !! number of elements a point has
        real(c_double), intent(in), dimension(n_dims, n_points), target :: data_points
            !! matrix with data points to cluster
        real(c_double), intent(inout), dimension(n_dims, n_clusters), target :: centroids
            !! matrix with initial centroids of the clusters, could be random data or actual points or unassigned garbage.
            !! The centroids should be unique. This is not checked in this routine.
            !! 
            !! The final values will be the final centroids of the clusters
        integer(c_int), intent(out), dimension(n_points), target :: labels
            !! array of labels, each index corresponds to the respective point's index, so first label is first point's label.
            !! 
            !! each label is the index of its related cluster -> `1<=label<=n_clusters=k`
        integer(c_int), intent(out), dimension(n_clusters), target :: label_counts
            !! holds the number of points having the respective label assigned
        integer(c_int), intent(out), target :: ierr
            !! Error code
        integer(c_int), intent(in), target :: max_iterations
            !! number of maximum iterations of the clustering
            !! The default value is `300_int32`.
        M_CHECK_IERR_NON_NULL
        M_CHECK_NON_NULL(n_clusters)
        M_CHECK_NON_NULL(data_points)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(n_dims)
        M_CHECK_NON_NULL(centroids)
        M_CHECK_NON_NULL(labels)
        M_CHECK_NON_NULL(label_counts)
        call k_means_clustering(&
            n_clusters = n_clusters,&
            data_points = data_points,&
            n_points = n_points,&
            n_dims = n_dims,&
            centroids = centroids,&
            labels = labels,&
            label_counts = label_counts,&
            ierr = ierr&
        )
    end subroutine k_means_clustering_c

    !> summary: C-wrapper for [[tox_clustering(module):linkage_clustering(subroutine)]]
    !| @note
    !| This subroutine operates in-place in the bottom triangle of the distance matrix and recovers it using the top triangle once done or on error.
    !| So there is no need to copy an existing distance matrix, just pass the original.
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
        use tox_clustering
        integer(c_int), intent(in), target :: n_points
            !! number of points to cluster
        real(c_double), intent(inout), dimension(n_points, n_points), target :: distances
            !! symmetric distance matrix, holding the positive distances between points. Distance of X->X is always zero.
            !! 
            !! @note
            !! This subroutine operates in-place in the bottom triangle of the distance matrix and recovers it using the top triangle once done or on error.
            !! So there is no need to copy an existing distance matrix, just pass the original.
            !! @endnote
        integer(c_int), intent(out), dimension(n_points - 1), target :: merge_i
            !! holds cluster labels of the merged node pair at iteration k -> positives relate to leafs/data point indices, negatives to inner nodes
        integer(c_int), intent(out), dimension(n_points - 1), target :: merge_j
            !! holds cluster labels of the merged node pair at iteration k -> positives relate to leafs/data point indices, negatives to inner nodes
        real(c_double), intent(out), dimension(n_points - 1), target :: heights
            !! height of the shorter branch of the merge, e.g. if (A,B)+(C) merges to ((A,B),C), the branch to (A,B) is shorter
        integer(c_int), intent(out), dimension(n_points - 1), target :: cluster_sizes
            !! size of cluster at iteration k
        character(len=1, kind=c_char), intent(in), dimension(8), target :: method
            !! used algorithm
            !! 
            !! |      Method      |   Value    |
            !! |------------------|------------|
            !! | Average / UPGMA  | "average"  |
            !! | Weighted / WPGMA | "weighted" |
            !! |       Ward       |   "ward"   |
        integer(c_int), intent(out), target :: ierr
            !! Error code
        character(len=:), allocatable :: method_f
        integer(int32) :: method_int_f
        M_CHECK_IERR_NON_NULL
        M_CHECK_NON_NULL(distances)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(merge_i)
        M_CHECK_NON_NULL(merge_j)
        M_CHECK_NON_NULL(heights)
        M_CHECK_NON_NULL(cluster_sizes)
        M_CHECK_NON_NULL(method)
        call c_char_1d_as_string(method, method_f, ierr)
        if (is_err(ierr)) return
        select case (method_f)
            case ("average")
                    method_int_f = METHOD_AVERAGE
            case ("weighted")
                    method_int_f = METHOD_WEIGHTED
            case ("ward")
                    method_int_f = METHOD_WARD
            case default
                call set_err(ierr, ERR_INVALID_INPUT)
                return
        end select
        call linkage_clustering(&
            distances = distances,&
            n_points = n_points,&
            merge_i = merge_i,&
            merge_j = merge_j,&
            heights = heights,&
            cluster_sizes = cluster_sizes,&
            method = method_int_f,&
            ierr = ierr&
        )
    end subroutine linkage_clustering_c

end module tox_clustering_c
#endif