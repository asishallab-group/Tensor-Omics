#include <src/macros.h>

!> Clustering for tensor omics: k-means over factors and trajectories, and hierarchical
!| (agglomerative) linkage clustering over a precomputed distance matrix.
!|
!| `k_means_clustering` partitions points in n dimensions; `cluster_factor_trajectories_k_means`
!| applies the same to whole time series, treating each sample's trajectory, over all its factors
!| and timepoints, as one point.
!| `linkage_clustering` takes the distances already computed and merges under the linkage
!| criterion asked for, so the same matrix can be re-clustered without recomputing it.
!|
!| Generated from [[tox_clustering_impl(module)]]; do not edit -- regenerate instead.
module tox_clustering
    use tox_clustering_impl, only: METHOD_AVERAGE, METHOD_WARD, METHOD_WEIGHTED, cluster_factor_trajectories_k_means_impl
    use tox_clustering_impl, only: k_means_clustering_impl, linkage_clustering_impl
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use tox_errors, only: set_ok, is_err, ERR_ALLOC_FAIL, ERR_INVALID_INPUT
    use tox_errors, only: clear_err_arg_pos, get_err_arg_pos, set_err, set_err_once
    use tox_errors, only: validate_all_in_range_real, validate_dimension_size, validate_distance_matrix, validate_in_range_int
    M_IMPLICIT_NONE
    private

    public :: cluster_factor_trajectories_k_means
    public :: cluster_factor_trajectories_k_means_expert
    public :: k_means_clustering
    public :: k_means_clustering_expert
    public :: linkage_clustering

contains

    !> summary: Validates its inputs, prepares what [[tox_clustering_impl(module):cluster_factor_trajectories_k_means_impl]] needs, then calls it. The entry point to reach for first; see [[tox_clustering(module):cluster_factor_trajectories_k_means_expert]] to prepare it yourself.
    !| One sample is one point: its trajectory is flattened over factors and time into a single
    !| vector, `factor1_t1, factor2_t1, ..., factor1_t2, ...`, so the clustering answers which
    !| samples follow a similar course. Clustering the individual `(sample, timepoint)` states
    !| instead is [[tox_clustering_impl(module):k_means_clustering_impl(subroutine)]] on the same
    !| array read as `(n_factors, n_samples*n_timepoints)`, which is a different question.
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
            !! The maximum valid value is `n_samples`.
        integer(int32), intent(in) :: n_factors
            !! number of factors
        integer(int32), intent(in) :: n_samples
            !! number of samples
        integer(int32), intent(in) :: n_timepoints
            !! number of timepoints
        real(real64), dimension(n_factors, n_samples, n_timepoints), intent(in) :: trajectories
            !! matrix with the trajectories to cluster
        real(real64), dimension(n_factors*n_timepoints, n_clusters), intent(inout) :: centroids
            !! matrix with initial centroids of the clusters, could be random data or actual points or unassigned garbage.
            !! The centroids should be unique. This is not checked in this routine.
            !!
            !! One centroid is a whole flattened trajectory, so it has `n_factors*n_timepoints` entries.
            !!
            !! The final values will be the final centroids of the clusters. A cluster that ends up
            !! without samples (its `label_counts` is 0) keeps the centroid it had when it lost its
            !! last sample, or its starting centroid if it never had one.
        integer(int32), dimension(n_samples), intent(out) :: labels
            !! array of labels, one per sample, so the first label is the first sample's.
            !!
            !! each label is the index of its related cluster -> `1<=label<=n_clusters=k`
        integer(int32), dimension(n_clusters), intent(out) :: label_counts
            !! holds the number of samples having the respective label assigned
        integer(int32), intent(in) :: max_iterations
            !! number of maximum iterations of the clustering. Unlike
            !! [[tox_clustering_impl(module):k_means_clustering_impl(subroutine)]], which defaults to
            !! 300, this one has to be given.
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.
        real(real64), dimension(:, :), allocatable :: tmp_flattened_trajectories
        real(real64), dimension(:, :), allocatable :: tmp_centroid_factors

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_in_range_int(n_clusters, ierr, arg_pos=1_int32, min=1_int32, max=n_samples)
        call validate_dimension_size(n_factors, ierr, arg_pos=3_int32)
        call validate_dimension_size(n_samples, ierr, arg_pos=4_int32)
        call validate_dimension_size(n_timepoints, ierr, arg_pos=5_int32)
        call validate_in_range_int(max_iterations, ierr, arg_pos=9_int32, min=1_int32)
        call validate_all_in_range_real(trajectories, n_factors * n_samples * n_timepoints, ierr, arg_pos=2_int32)
        call validate_all_in_range_real(centroids, (n_factors*n_timepoints) * n_clusters, ierr, arg_pos=6_int32)
        if (is_err(ierr)) return
#endif

        M_ALLOCATE(tmp_flattened_trajectories(n_factors*n_timepoints, n_samples))
        M_ALLOCATE(tmp_centroid_factors(n_factors*n_timepoints, n_clusters))

        call cluster_factor_trajectories_k_means_impl(&
            n_clusters = n_clusters,&
            trajectories = trajectories,&
            n_factors = n_factors,&
            n_samples = n_samples,&
            n_timepoints = n_timepoints,&
            centroids = centroids,&
            labels = labels,&
            label_counts = label_counts,&
            tmp_flattened_trajectories = tmp_flattened_trajectories,&
            tmp_centroid_factors = tmp_centroid_factors,&
            max_iterations = max_iterations&
        )
    end subroutine cluster_factor_trajectories_k_means

    !> summary: Validates its inputs, then calls [[tox_clustering_impl(module):cluster_factor_trajectories_k_means_impl]] with what you supply. The expert entry point: it allocates nothing and prepares nothing; [[tox_clustering(module):cluster_factor_trajectories_k_means]] does both.
    !| One sample is one point: its trajectory is flattened over factors and time into a single
    !| vector, `factor1_t1, factor2_t1, ..., factor1_t2, ...`, so the clustering answers which
    !| samples follow a similar course. Clustering the individual `(sample, timepoint)` states
    !| instead is [[tox_clustering_impl(module):k_means_clustering_impl(subroutine)]] on the same
    !| array read as `(n_factors, n_samples*n_timepoints)`, which is a different question.
    pure subroutine cluster_factor_trajectories_k_means_expert(&
            n_clusters,&
            trajectories,&
            n_factors,&
            n_samples,&
            n_timepoints,&
            centroids,&
            labels,&
            label_counts,&
            tmp_flattened_trajectories,&
            tmp_centroid_factors,&
            max_iterations,&
            ierr&
        )
        integer(int32), intent(in) :: n_clusters
            !! number (`k`) of clusters
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_samples`.
        integer(int32), intent(in) :: n_factors
            !! number of factors
        integer(int32), intent(in) :: n_samples
            !! number of samples
        integer(int32), intent(in) :: n_timepoints
            !! number of timepoints
        real(real64), dimension(n_factors, n_samples, n_timepoints), intent(in) :: trajectories
            !! matrix with the trajectories to cluster
        real(real64), dimension(n_factors*n_timepoints, n_clusters), intent(inout) :: centroids
            !! matrix with initial centroids of the clusters, could be random data or actual points or unassigned garbage.
            !! The centroids should be unique. This is not checked in this routine.
            !!
            !! One centroid is a whole flattened trajectory, so it has `n_factors*n_timepoints` entries.
            !!
            !! The final values will be the final centroids of the clusters. A cluster that ends up
            !! without samples (its `label_counts` is 0) keeps the centroid it had when it lost its
            !! last sample, or its starting centroid if it never had one.
        integer(int32), dimension(n_samples), intent(out) :: labels
            !! array of labels, one per sample, so the first label is the first sample's.
            !!
            !! each label is the index of its related cluster -> `1<=label<=n_clusters=k`
        integer(int32), dimension(n_clusters), intent(out) :: label_counts
            !! holds the number of samples having the respective label assigned
        real(real64), dimension(n_factors*n_timepoints, n_samples), intent(out) :: tmp_flattened_trajectories
            !! Work matrix holding each sample's trajectory as one column
        real(real64), dimension(n_factors*n_timepoints, n_clusters), intent(out) :: tmp_centroid_factors
            !! Work matrix: for each cluster and each entry of a flattened trajectory, the power of
            !! two its values are summed with when the centroids are recomputed
        integer(int32), intent(in) :: max_iterations
            !! number of maximum iterations of the clustering. Unlike
            !! [[tox_clustering_impl(module):k_means_clustering_impl(subroutine)]], which defaults to
            !! 300, this one has to be given.
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_in_range_int(n_clusters, ierr, arg_pos=1_int32, min=1_int32, max=n_samples)
        call validate_dimension_size(n_factors, ierr, arg_pos=3_int32)
        call validate_dimension_size(n_samples, ierr, arg_pos=4_int32)
        call validate_dimension_size(n_timepoints, ierr, arg_pos=5_int32)
        call validate_in_range_int(max_iterations, ierr, arg_pos=11_int32, min=1_int32)
        call validate_all_in_range_real(trajectories, n_factors * n_samples * n_timepoints, ierr, arg_pos=2_int32)
        call validate_all_in_range_real(centroids, (n_factors*n_timepoints) * n_clusters, ierr, arg_pos=6_int32)
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
            tmp_flattened_trajectories = tmp_flattened_trajectories,&
            tmp_centroid_factors = tmp_centroid_factors,&
            max_iterations = max_iterations&
        )
    end subroutine cluster_factor_trajectories_k_means_expert

    !> summary: Validates its inputs, prepares what [[tox_clustering_impl(module):k_means_clustering_impl]] needs, then calls it. The entry point to reach for first; see [[tox_clustering(module):k_means_clustering_expert]] to prepare it yourself.
    !| 1. Assigns each data point to one of `k` clusters whose centroid is closest
    !| 2. Recalculates the centroids using the mean of its assigned points
    !| 3. repeat 1-2 until assignment remains unchanged
    !|
    !| A point exactly as close to several centroids goes to the last of them, the one with the
    !| highest cluster index.
    !|
    !| Each point's squared distances are compared scaled by a power of two picked for that point,
    !| by f42's rule for sums of squares, from its largest coordinate difference to the nearest
    !| centroid in that measure; ordinary data is not scaled at all.
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
            !! The final values will be the final centroids of the clusters. A cluster that ends up
            !! without points (its `label_counts` is 0) keeps the centroid it had when it lost its
            !! last point, or its starting centroid if it never had one.
        integer(int32), dimension(n_points), intent(out) :: labels
            !! array of labels, each index corresponds to the respective point's index, so first label is first point's label.
            !!
            !! each label is the index of its related cluster -> `1<=label<=n_clusters=k`
        integer(int32), dimension(n_clusters), intent(out) :: label_counts
            !! holds the number of points having the respective label assigned
        integer(int32), intent(in), optional :: max_iterations
            !! number of maximum iterations of the clustering. Defaults to 300, while
            !! [[tox_clustering_impl(module):cluster_factor_trajectories_k_means_impl(subroutine)]]
            !! requires it.
            !! The minimum valid value is `1_int32`.
            !! The default value is `300_int32`.
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.
        real(real64), dimension(:, :), allocatable :: tmp_centroid_factors

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_in_range_int(n_clusters, ierr, arg_pos=1_int32, min=1_int32, max=n_points)
        call validate_dimension_size(n_points, ierr, arg_pos=3_int32)
        call validate_dimension_size(n_dims, ierr, arg_pos=4_int32)
        call validate_in_range_int(max_iterations, ierr, arg_pos=8_int32, min=1_int32)
        call validate_all_in_range_real(data_points, n_dims * n_points, ierr, arg_pos=2_int32)
        call validate_all_in_range_real(centroids, n_dims * n_clusters, ierr, arg_pos=5_int32)
        if (is_err(ierr)) return
#endif

        M_ALLOCATE(tmp_centroid_factors(n_dims, n_clusters))

        call k_means_clustering_impl(&
            n_clusters = n_clusters,&
            data_points = data_points,&
            n_points = n_points,&
            n_dims = n_dims,&
            centroids = centroids,&
            labels = labels,&
            label_counts = label_counts,&
            tmp_centroid_factors = tmp_centroid_factors,&
            max_iterations = max_iterations&
        )
    end subroutine k_means_clustering

    !> summary: Validates its inputs, then calls [[tox_clustering_impl(module):k_means_clustering_impl]] with what you supply. The expert entry point: it allocates nothing and prepares nothing; [[tox_clustering(module):k_means_clustering]] does both.
    !| 1. Assigns each data point to one of `k` clusters whose centroid is closest
    !| 2. Recalculates the centroids using the mean of its assigned points
    !| 3. repeat 1-2 until assignment remains unchanged
    !|
    !| A point exactly as close to several centroids goes to the last of them, the one with the
    !| highest cluster index.
    !|
    !| Each point's squared distances are compared scaled by a power of two picked for that point,
    !| by f42's rule for sums of squares, from its largest coordinate difference to the nearest
    !| centroid in that measure; ordinary data is not scaled at all.
    pure subroutine k_means_clustering_expert(&
            n_clusters,&
            data_points,&
            n_points,&
            n_dims,&
            centroids,&
            labels,&
            label_counts,&
            tmp_centroid_factors,&
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
            !! The final values will be the final centroids of the clusters. A cluster that ends up
            !! without points (its `label_counts` is 0) keeps the centroid it had when it lost its
            !! last point, or its starting centroid if it never had one.
        integer(int32), dimension(n_points), intent(out) :: labels
            !! array of labels, each index corresponds to the respective point's index, so first label is first point's label.
            !!
            !! each label is the index of its related cluster -> `1<=label<=n_clusters=k`
        integer(int32), dimension(n_clusters), intent(out) :: label_counts
            !! holds the number of points having the respective label assigned
        real(real64), dimension(n_dims, n_clusters), intent(out) :: tmp_centroid_factors
            !! Work matrix: for each cluster and dimension, the power of two the values are summed
            !! with when the centroids are recomputed
        integer(int32), intent(in), optional :: max_iterations
            !! number of maximum iterations of the clustering. Defaults to 300, while
            !! [[tox_clustering_impl(module):cluster_factor_trajectories_k_means_impl(subroutine)]]
            !! requires it.
            !! The minimum valid value is `1_int32`.
            !! The default value is `300_int32`.
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_in_range_int(n_clusters, ierr, arg_pos=1_int32, min=1_int32, max=n_points)
        call validate_dimension_size(n_points, ierr, arg_pos=3_int32)
        call validate_dimension_size(n_dims, ierr, arg_pos=4_int32)
        call validate_in_range_int(max_iterations, ierr, arg_pos=9_int32, min=1_int32)
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
            tmp_centroid_factors = tmp_centroid_factors,&
            max_iterations = max_iterations&
        )
    end subroutine k_means_clustering_expert

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
            !! symmetric distance matrix, holding the non-negative distances between points. Distance of X->X is always zero.
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
            !! the first of the two clusters merged at merge k: a positive label is a point, its
            !! 1-based index; a negative label `-j` is the cluster made by merge j. Of the two, this
            !! is the one whose lowest point index is the higher.
            !!
            !! The signs are the reverse of R's `hclust`, whose `merge` holds points as negative
            !! and merged clusters as positive numbers.
        integer(int32), dimension(n_points - 1), intent(out) :: merge_j
            !! the second of the two clusters merged at merge k, labelled as in `merge_i`: the one
            !! whose lowest point index is the lower.
            !!
            !! Two pairs at the same distance merge in the order of their lowest point indices: the
            !! pair whose lower cluster has the lowest point index first, then the one whose higher
            !! cluster has.
        real(real64), dimension(n_points - 1), intent(out) :: heights
            !! the linkage distance at which merge k joins its two clusters, by the method's
            !! criterion. The heights never decrease from one merge to the next.
        integer(int32), dimension(n_points - 1), intent(out) :: cluster_sizes
            !! number of points in the cluster made by merge k
        integer(int32), intent(in) :: method
            !! used algorithm
            !!
            !! | Method           | Value                                                     |
            !! |------------------|-----------------------------------------------------------|
            !! | Average / UPGMA  | [[tox_clustering_impl(module):METHOD_AVERAGE(variable)]]  |
            !! | Weighted / WPGMA | [[tox_clustering_impl(module):METHOD_WEIGHTED(variable)]] |
            !! | Ward             | [[tox_clustering_impl(module):METHOD_WARD(variable)]]     |
        integer(int32), intent(out) :: ierr
            !! Error code: `ERR_NAN_INF` when a merge distance passes `huge`, which only Ward's can
            !! for finite input, as it grows with the clusters' sizes. The merges before it are
            !! written, the later ones are not, and the matrix still comes back unchanged.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_points, ierr, arg_pos=2_int32)
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
            method = method,&
            ierr = ierr&
        )
        if (is_err(ierr)) then
            select case (get_err_arg_pos(ierr))
            case (1_int32:7_int32)
                ! the same position here
            case default
                call clear_err_arg_pos(ierr)
            end select
            return
        end if
    end subroutine linkage_clustering

end module tox_clustering
