#include <src/macros.h>

!> Clustering for tensor omics: k-means over factors and trajectories, and hierarchical
!| (agglomerative) linkage clustering over a precomputed distance matrix.
!|
!| `k_means_clustering` partitions points in n dimensions; `cluster_factor_trajectories_k_means`
!| applies the same to whole time series, treating each sample's trajectory, over all its factors
!| and timepoints, as one point.
!| `linkage_clustering` takes the distances already computed and merges under the linkage
!| criterion asked for, so the same matrix can be re-clustered without recomputing it.
module tox_clustering_impl
    use f42_safeguard
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use, intrinsic :: iso_c_binding, only: c_bool
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_positive_inf, ieee_is_nan
    use f42_math_impl, only: scaling_exponent_for_largest, scaling_factors
    use tox_errors, only: set_ok, set_err, ERR_NAN_INF
    M_IMPLICIT_NONE

    integer(int32), parameter :: METHOD_AVERAGE = 0_int32
        !! Method code for average linkage clustering in [[tox_clustering(module):linkage_clustering(subroutine)]]
    integer(int32), parameter :: METHOD_WEIGHTED = 1_int32
        !! Method code for weighted linkage clustering in [[tox_clustering(module):linkage_clustering(subroutine)]]
    integer(int32), parameter :: METHOD_WARD = 2_int32
        !! Method code for ward linkage clustering in [[tox_clustering(module):linkage_clustering(subroutine)]]

#define DEFAULT_MAX_ITER_K_MEANS 300_int32
contains

    !> summary: Performs k-means clustering on whole factor trajectories, so which samples evolve alike
    !| AUTHOR_FRANZ_ERIC_SILL
    !| One sample is one point: its trajectory is flattened over factors and time into a single
    !| vector, `factor1_t1, factor2_t1, ..., factor1_t2, ...`, so the clustering answers which
    !| samples follow a similar course. Clustering the individual `(sample, timepoint)` states
    !| instead is [[tox_clustering_impl(module):k_means_clustering_impl(subroutine)]] on the same
    !| array read as `(n_factors, n_samples*n_timepoints)`, which is a different question.
    pure subroutine cluster_factor_trajectories_k_means_impl(n_clusters, trajectories, n_factors, n_samples, n_timepoints, centroids, labels, label_counts, tmp_flattened_trajectories, &
                                                         tmp_centroid_factors, max_iterations)
        integer(int32), intent(in) :: n_clusters
            !! number (`k`) of clusters
            !! DM_MIN(1_int32)
            !! DM_MAX(n_samples)
        integer(int32), intent(in) :: n_factors
            !! number of factors
        integer(int32), intent(in) :: n_timepoints
            !! number of timepoints
        integer(int32), intent(in) :: n_samples
            !! number of samples
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
            !! DM_MIN(1_int32)

        integer(int32) :: i_sample, i_time, i_factor

        ! A sample's values are not contiguous in `trajectories` -- factor varies fastest, then
        ! sample, then time -- so the flattening is a copy, not a reshape. Factor-major within each
        ! timepoint keeps the two axes of a trajectory adjacent in the order they are read in.
        do concurrent (i_sample = 1:n_samples) shared(trajectories, tmp_flattened_trajectories, n_timepoints, n_factors)
            do i_time = 1, n_timepoints
                do i_factor = 1, n_factors
                    tmp_flattened_trajectories((i_time - 1)*n_factors + i_factor, i_sample) = &
                        trajectories(i_factor, i_sample, i_time)
                end do
            end do
        end do

        call k_means_clustering_impl(n_clusters, tmp_flattened_trajectories, n_samples, n_factors*n_timepoints, &
                                     centroids, labels, label_counts, tmp_centroid_factors, max_iterations)
    end subroutine cluster_factor_trajectories_k_means_impl

    !> summary: k-means clustering algorithm
    !| AUTHOR_FRANZ_ERIC_SILL
    !|
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
    pure subroutine k_means_clustering_impl(n_clusters, data_points, n_points, n_dims, centroids, labels, label_counts, &
                                            tmp_centroid_factors, max_iterations)
        integer(int32), intent(in) :: n_clusters
            !! number (`k`) of clusters
            !! DM_MIN(1_int32)
            !! DM_MAX(n_points)
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
            !! DM_MIN(1_int32)
            !! DM_DEFAULT(DEFAULT_MAX_ITER_K_MEANS)

        integer(int32) :: label, iteration, i_point, max_iter
        logical(c_bool) :: labels_changed

        M_DEFAULT_VAL(max_iterations, max_iter, DEFAULT_MAX_ITER_K_MEANS)

        ! initialize labels to a value outside the valid [1, n_clusters] range so the first
        ! comparison below never reads undefined intent(out) memory and always triggers
        ! labels_changed on the first iteration
        labels = 0_int32

        iteration = 0_int32
        do while (iteration < max_iter)
            iteration = iteration + 1

            ! assign points to their closest cluster centroid
            labels_changed = .false.
            label_counts = 0
            do i_point = 1, n_points
                call k_means_assign_cluster_helper(i_point, n_clusters, data_points, n_points, n_dims, centroids, label)

                label_counts(label) = label_counts(label) + 1
                if (labels(i_point) /= label) then
                    labels(i_point) = label
                    labels_changed = .true.
                end if
            end do

            ! if the assignments did not change, the clustering is done
            if (.not. labels_changed) exit

            call k_means_recompute_cluster_centroids_helper(data_points, n_points, n_dims, centroids, n_clusters, labels, &
                                                            label_counts, tmp_centroid_factors)
        end do
    end subroutine k_means_clustering_impl
    !> AUTHOR_FRANZ_ERIC_SILL
    !| Helper to recompute the centroids in k-means by taking the mean of assigned points. Each
    !| coordinate of a centroid is the mean of one set, its points' values in that dimension, and
    !| is summed scaled by the power of two f42's rule picks for that set's largest magnitude, as
    !| f42's mean sums a set: the mean of two values at `2**1023` is `2**1023`, and a coordinate of
    !| small values beside a large one keeps its own. The factors are multiplied in, not applied
    !| with scale(), which is a library call per element; values whose set needs no scaling are
    !| summed exactly as before, in point order.
    pure subroutine k_means_recompute_cluster_centroids_helper(data_points, n_points, n_dims, centroids, n_clusters, labels, &
                                                               label_counts, factors)
        integer(int32), intent(in) :: n_clusters
            !! number (`k`) of clusters
        integer(int32), intent(in) :: n_points
            !! number of points to cluster
        integer(int32), intent(in) :: n_dims
            !! number of elements a point has
        real(real64), dimension(n_dims, n_points), intent(in) :: data_points
            !! matrix with data points to cluster
        real(real64), dimension(n_dims, n_clusters), intent(inout) :: centroids
            !! the centroids, recomputed in place; a cluster without points keeps its own
        integer(int32), dimension(n_points), intent(in) :: labels
            !! array of labels, each index corresponds to the respective point's index, so first label is first point's label.
            !!
            !! each label is the index of its related cluster -> `1<=label<=n_clusters=k`
        integer(int32), dimension(n_clusters), intent(in) :: label_counts
            !! holds the number of points having the respective label assigned
        real(real64), dimension(n_dims, n_clusters), intent(out) :: factors
            !! work: first each set's largest magnitude, then the factor its values are summed with

        integer(int32) :: i_point, label, i_dim, i_cluster
        real(real64) :: down, up

        ! the largest magnitude of each (cluster, dimension) set
        factors = 0.0_real64
        do i_point = 1, n_points
            label = labels(i_point)
            do i_dim = 1, n_dims
                factors(i_dim, label) = max(factors(i_dim, label), abs(data_points(i_dim, i_point)))
            end do
        end do

        ! A cluster without points keeps the centroid it has: it has no mean, and any value made
        ! up for it (the origin, say) would make the result depend on the coordinate system.
        do i_cluster = 1, n_clusters
            if (label_counts(i_cluster) == 0) cycle
            do i_dim = 1, n_dims
                ! one of the two factors is 1, so their product is the power of two itself
                call scaling_factors(scaling_exponent_for_largest(factors(i_dim, i_cluster)), down, up)
                factors(i_dim, i_cluster) = down*up
                centroids(i_dim, i_cluster) = 0.0_real64
            end do
        end do

        do i_point = 1, n_points
            label = labels(i_point)
            do i_dim = 1, n_dims
                centroids(i_dim, label) = centroids(i_dim, label) + data_points(i_dim, i_point)*factors(i_dim, label)
            end do
        end do

        do i_cluster = 1, n_clusters
            if (label_counts(i_cluster) == 0) cycle
            do i_dim = 1, n_dims
                ! dividing by a power of two is exact, and brings the mean back to its values' scale
                centroids(i_dim, i_cluster) = (centroids(i_dim, i_cluster)/real(label_counts(i_cluster), real64)) &
                                              /factors(i_dim, i_cluster)
            end do
        end do
    end subroutine k_means_recompute_cluster_centroids_helper

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Helper routine to assign a cluster's label to a point in k-means
    pure subroutine k_means_assign_cluster_helper(i_point, n_clusters, data_points, n_points, n_dims, centroids, label)
        integer(int32), intent(in) :: i_point
            !! index of the point that should be assigned a cluster
        integer(int32), intent(in) :: n_clusters
            !! number (`k`) of clusters
        integer(int32), intent(in) :: n_points
            !! number of points to cluster
        integer(int32), intent(in) :: n_dims
            !! number of elements a point has
        real(real64), dimension(n_dims, n_points), intent(in) :: data_points
            !! matrix with data points to cluster
        real(real64), dimension(n_dims, n_clusters), intent(in) :: centroids
            !! matrix with initial centroids of the clusters, could be random data or actual points or unassigned garbage.
            !! The centroids should be unique. This is not checked in this routine.
            !!
            !! The final values will be the final centroids of the clusters
        integer(int32), intent(out) :: label
            !! assigned label

        integer(int32) :: i_cluster, i_dim
        real(real64) :: squared_dist_to_centroid, smallest_squared_dist, nearest_largest_difference, largest_difference
        real(real64) :: down, up

        ! The squared distances are compared scaled by a power of two picked for this point: from the
        ! largest coordinate difference to the centroid whose largest difference is smallest (an
        ! overflowing difference counts as +Inf). Scaled so, the nearest centroid's squares neither
        ! overflow nor all vanish; a farther one's may overflow, but it loses anyway. A common factor
        ! changes no comparison, and an ordinary point gets the exponent 0, which leaves its
        ! arithmetic as it was. The factors are multiplied in, not applied with scale(), which is a
        ! library call per element.
        nearest_largest_difference = ieee_value(1.0_real64, ieee_positive_inf)
        do i_cluster = 1, n_clusters
            largest_difference = 0.0_real64
            do i_dim = 1, n_dims
                largest_difference = max(largest_difference, abs(data_points(i_dim, i_point) - centroids(i_dim, i_cluster)))
            end do
            nearest_largest_difference = min(nearest_largest_difference, largest_difference)
        end do
        call scaling_factors(scaling_exponent_for_largest(nearest_largest_difference), down, up)

        label = 0_int32
        smallest_squared_dist = ieee_value(1.0_real64, ieee_positive_inf)

        ! calculate distance to each centroid and assign to closest
        do i_cluster = 1, n_clusters
            squared_dist_to_centroid = 0.0_real64
            ! GFORTRAN BUG: do concurrent (i_dim = 1:n_dims) shared(data_points, i_point, i_cluster) reduce(+:squared_dist_to_centroid)
            do i_dim = 1, n_dims
                squared_dist_to_centroid = squared_dist_to_centroid &
                                           + ((data_points(i_dim, i_point)*down - centroids(i_dim, i_cluster)*down)*up)**2
            end do

            if (squared_dist_to_centroid <= smallest_squared_dist) then
                smallest_squared_dist = squared_dist_to_centroid
                label = i_cluster
            end if
        end do
    end subroutine k_means_assign_cluster_helper

    !> summary: Perform linkage clustering on a distance matrix.
    !| AUTHOR_FRANZ_ERIC_SILL
    !|
    !| @note
    !| The bottom triangle is used as scratch and restored from the top triangle before
    !| returning, on success or on error, so the matrix comes back unchanged. There is no
    !| need to copy it before calling.
    !| @endnote
    pure subroutine linkage_clustering_impl(distances, n_points, merge_i, merge_j, heights, cluster_sizes, method, ierr)
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
            !! | Method | Value |
            !! |--------|-------|
            !! | Average / UPGMA | [[tox_clustering_impl(module):METHOD_AVERAGE(variable)]] |
            !! | Weighted / WPGMA | [[tox_clustering_impl(module):METHOD_WEIGHTED(variable)]] |
            !! | Ward | [[tox_clustering_impl(module):METHOD_WARD(variable)]] |
        integer(int32), intent(out) :: ierr
            !! Error code: `ERR_NAN_INF` when a merge distance passes `huge`, which only Ward's can
            !! for finite input, as it grows with the clusters' sizes. The merges before it are
            !! written, the later ones are not, and the matrix still comes back unchanged.

        integer(int32) :: i, idx_A, idx_B, size_B, size_A, cluster_label, storage_exponent, i_row, i_col
        real(real64) :: dist_AB, height
        logical(c_bool) :: found
        ! Ward's distances grow with the clusters' sizes, but stay below sqrt(n/2) times the largest
        ! input for n points, so below 2**15.5 times it for any int32 count: stored scaled by 2**-16,
        ! no Ward distance of an input up to huge can pass huge in the scratch triangle
        real(real64), parameter :: LARGEST_STORED_UNSCALED = 2.0_real64**1007
        integer(int32), parameter :: WARD_STORAGE_EXPONENT = -16_int32

        !TODO optimize: get_min_distance_indices_helper does a full O(n^2) scan of the lower triangle on every
        !               one of the n-1 merge iterations, making this O(n^3) overall. Classic agglomerative
        !               clustering can be done in O(n^2) total with the nearest-neighbor-chain algorithm, which
        !               would matter a lot for larger n_points.
        ! The scratch triangle holds Ward's distances scaled by 2**WARD_STORAGE_EXPONENT when an
        ! input exceeds 2**1007, so that a distance whose true value passes huge on the way, and
        ! may come back below it once more points join, stays finite and comparable; +Inf remains
        ! the marker of a retired entry alone. Every height is scaled back as it is written.
        storage_exponent = 0_int32
        if (method == METHOD_WARD) then
            if (maxval(distances) > LARGEST_STORED_UNSCALED) storage_exponent = WARD_STORAGE_EXPONENT
        end if
        if (storage_exponent /= 0_int32) then
            do i_col = 1, n_points - 1
                do i_row = i_col + 1, n_points
                    distances(i_row, i_col) = scale(distances(i_row, i_col), storage_exponent)
                end do
            end do
        end if

        call set_ok(ierr)
        do i = 1, n_points - 1
            call get_min_distance_indices_helper(distances, n_points, idx_A, idx_B, dist_AB, found)
            ! +Inf marks a retired entry, so a merge distance that passed huge cannot be told from
            ! one; with no entry below it the merge is not found
            if (.not. found) then
                call set_err(ierr, ERR_NAN_INF)
                exit
            end if
            ! a merge whose true height passes huge cannot be written; the merges before it are right
            height = scale(dist_AB, -storage_exponent)
            if (height > huge(height)) then
                call set_err(ierr, ERR_NAN_INF)
                exit
            end if
            heights(i) = height

            ! Get Weight and cluster index/label
            call get_cluster_data_linkage_helper(distances, n_points, cluster_sizes, idx_B, size_B, cluster_label)
            merge_j(i) = cluster_label
            call get_cluster_data_linkage_helper(distances, n_points, cluster_sizes, idx_A, size_A, cluster_label)
            merge_i(i) = cluster_label

            cluster_sizes(i) = size_B + size_A

            select case (method)
            case (METHOD_AVERAGE)
                call merge_distances_xPGMA_linkage_helper(distances, n_points, idx_A, idx_B, size_A, size_B, i, dist_AB)
            case (METHOD_WEIGHTED)
                call merge_distances_xPGMA_linkage_helper(distances, n_points, idx_A, idx_B, 1_int32, 1_int32, i, dist_AB)
            case (METHOD_WARD)
                call merge_distances_ward_linkage_helper(distances, n_points, idx_A, idx_B, size_A, size_B, i, cluster_sizes, &
                                                         dist_AB)
            end select
        end do

        call recover_distance_matrix_helper(distances, n_points)
    end subroutine linkage_clustering_impl

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Helper to recover a distance matrix from its top triangle.
    pure subroutine recover_distance_matrix_helper(distances, n_points)
        integer(int32), intent(in) :: n_points
            !! number of points to cluster
        real(real64), dimension(n_points, n_points), intent(inout) :: distances
            !! top triangle holds distances, bottom triangle will be recovered

        integer(int32) :: i_row, i_col

        do concurrent (i_col = 1:n_points) shared(distances)
            ! recover self distance
            distances(i_col, i_col) = 0.0_real64
            ! recover triangle
            do concurrent (i_row = i_col + 1:n_points) shared(distances, i_col)
                distances(i_row, i_col) = distances(i_col, i_row)
            end do
        end do
    end subroutine recover_distance_matrix_helper

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Helper for linkage clustering to get the size/weight and label of a cluster
    !| by using the self-distance as indicator in which iteration the cluster was merged
    pure subroutine get_cluster_data_linkage_helper(distances, n_points, cluster_sizes, idx, weight, cluster_label)
        integer(int32), intent(in) :: n_points
            !! number of points to cluster
        real(real64), dimension(n_points, n_points), intent(in) :: distances
            !! distance matrix, holding the distances between points
        integer(int32), dimension(n_points - 1), intent(in) :: cluster_sizes
            !! size of cluster at iteration k
        integer(int32), intent(in) :: idx
            !! index of cluster in `distances`
        integer(int32), intent(out) :: weight
            !! number of leafs of cluster
        integer(int32), intent(out) :: cluster_label
            !! index/label of cluster

        integer(int32) :: iteration_k

        iteration_k = int(distances(idx, idx), int32)
        if (iteration_k == 0) then
            weight = 1
            cluster_label = idx
        else
            weight = cluster_sizes(iteration_k)
            cluster_label = -iteration_k
        end if
    end subroutine get_cluster_data_linkage_helper

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Helper for *PGMA algorithms for updating the distance matrix on merge of two clusters
    pure subroutine merge_distances_xPGMA_linkage_helper(distances, n_points, idx_A, idx_B, size_A, size_B, iteration_k, dist_AB)
        integer(int32), intent(in) :: n_points
            !! number of points to cluster
        real(real64), dimension(n_points, n_points), intent(inout) :: distances
            !! distance matrix, holding the distances between points
        integer(int32), intent(in) :: idx_A
            !! row index of min distance in `distances`
        integer(int32), intent(in) :: idx_B
            !! column index of min distance in `distances`
        integer(int32), intent(in) :: size_A
            !! cluster size of cluster at `idx_A`
        integer(int32), intent(in) :: size_B
            !! cluster size of cluster at `idx_B`
        integer(int32), intent(in) :: iteration_k
            !! iteration of current merge
        real(real64), intent(in) :: dist_AB
            !! distance between cluster A and B, the merge's height; no updated distance is below it

        integer(int32) :: i_node
        real(real64) :: new_dist, weight_A, weight_B, new_dist_part_AC, new_dist_part_BC

        ! Merge distances
        ! Update merged nodes with arithmetic mean distances
        ! Use column index for merged distances -> fill idx_A with infinity values
        ! update bottom triangle
        ! Each distance is weighted before the two are added, so the sum is an average of finite
        ! distances and cannot overflow; forming size_A*d(A,C) + size_B*d(B,C) first could.
        weight_A = real(size_A, real64)/real(size_A + size_B, real64)
        weight_B = real(size_B, real64)/real(size_A + size_B, real64)
        ! Lance-Williams update formula for UPGMA/WPGMA: new d(AB,C) = (size_A*d(A,C) + size_B*d(B,C)) / (size_A+size_B).
        ! For WPGMA (METHOD_WEIGHTED) the caller passes size_A=size_B=1 instead of the true cluster
        ! sizes, which reduces this to the unweighted average (d(A,C)+d(B,C))/2, giving every
        ! previous merge equal weight regardless of how many leaves it represents.
        do i_node = 1, n_points
            if (distances(i_node, i_node) >= 0.0_real64) then
                ! i_node < idx_B < idx_A -> fill rows in both
                if (i_node < idx_B) then
                    new_dist_part_AC = distances(idx_A, i_node)*weight_A
                    new_dist_part_BC = distances(idx_B, i_node)*weight_B
                    ! an average of distances that are all at least dist_AB is too, but rounding can
                    ! land it an ulp below; that would make a later height lower than this one
                    new_dist = max(new_dist_part_AC + new_dist_part_BC, dist_AB)
                    distances(idx_B, i_node) = new_dist
                    distances(idx_A, i_node) = ieee_value(1.0_real64, ieee_positive_inf)
                    ! idx_B <= i_node < idx_A -> fill col of idx_B, row of idx_A
                else if (i_node < idx_A) then
                    new_dist_part_AC = distances(idx_A, i_node)*weight_A
                    new_dist_part_BC = distances(i_node, idx_B)*weight_B
                    ! an average of distances that are all at least dist_AB is too, but rounding can
                    ! land it an ulp below; that would make a later height lower than this one
                    new_dist = max(new_dist_part_AC + new_dist_part_BC, dist_AB)
                    distances(i_node, idx_B) = new_dist
                    distances(idx_A, i_node) = ieee_value(1.0_real64, ieee_positive_inf)
                    ! idx_B < idx_A < i_node -> fill columns in both
                else if (i_node > idx_A) then
                    new_dist_part_AC = distances(i_node, idx_A)*weight_A
                    new_dist_part_BC = distances(i_node, idx_B)*weight_B
                    ! an average of distances that are all at least dist_AB is too, but rounding can
                    ! land it an ulp below; that would make a later height lower than this one
                    new_dist = max(new_dist_part_AC + new_dist_part_BC, dist_AB)
                    distances(i_node, idx_B) = new_dist
                    ! don't fill idx_A with infinity, as whole column will be marked inactive
                end if
            end if
        end do

        call post_merge_distances_helper(distances, n_points, idx_A, idx_B, iteration_k)
    end subroutine merge_distances_xPGMA_linkage_helper

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Helper for ward algorithm for updating the distance matrix on merge of two clusters
    pure subroutine merge_distances_ward_linkage_helper(distances, n_points, idx_A, idx_B, size_A, size_B, iteration_k, cluster_sizes, &
                                                        dist_AB)
        integer(int32), intent(in) :: n_points
            !! number of points to cluster
        real(real64), dimension(n_points, n_points), intent(inout) :: distances
            !! distance matrix, holding the distances between points
        integer(int32), intent(in) :: idx_A
            !! row index of min distance in `distances`
        integer(int32), intent(in) :: idx_B
            !! column index of min distance in `distances`
        integer(int32), intent(in) :: size_A
            !! cluster size of cluster at `idx_A`
        integer(int32), intent(in) :: size_B
            !! cluster size of cluster at `idx_B`
        integer(int32), intent(in) :: iteration_k
            !! iteration of current merge
        integer(int32), dimension(n_points - 1), intent(in) :: cluster_sizes
            !! size of cluster at iteration k, needed for ward algorithm
        real(real64), intent(in) :: dist_AB
            !! distance between cluster A and B

        integer(int32) :: i_node, size_C, cluster_label_C, exponent
        real(real64) :: new_dist, new_dist_part_AC, new_dist_part_BC, new_dist_part_AB, size_AC_real, size_BC_real, size_ABC_real

        ! Merge distances
        ! Lance-Williams update formula for Ward's minimum-variance method:
        !   d(AB,C) = sqrt( ((size_A+size_C)*d(A,C)^2 + (size_B+size_C)*d(B,C)^2 - size_C*d(A,B)^2) / (size_A+size_B+size_C) )
        ! This is the increase in within-cluster variance caused by merging A and B, expressed in
        ! distance units via the square root (unlike xPGMA, the components are combined on squared
        ! distances, not the distances themselves).
        ! Use column index for merged distances -> fill idx_A with infinity values
        ! update bottom triangle
        do i_node = 1, n_points
            if (distances(i_node, i_node) >= 0.0_real64) then
                ! get number of leafs for current node/cluster
                call get_cluster_data_linkage_helper(distances, n_points, cluster_sizes, i_node, size_C, cluster_label_C)
                size_AC_real = real(size_A + size_C, real64)
                size_BC_real = real(size_B + size_C, real64)
                size_ABC_real = real(size_A + size_B + size_C, real64)

                ! i_node < idx_B < idx_A -> fill rows in both
                if (i_node < idx_B) then
                    ! scaled for the three distances this update combines
                    exponent = scaling_exponent_for_largest(max(distances(idx_A, i_node), distances(idx_B, i_node), dist_AB))
                    new_dist_part_AB = scale(dist_AB, exponent)**2*(real(size_C, real64)/size_ABC_real)
                    new_dist_part_AC = scale(distances(idx_A, i_node), exponent)**2*(size_AC_real/size_ABC_real)
                    new_dist_part_BC = scale(distances(idx_B, i_node), exponent)**2*(size_BC_real/size_ABC_real)
                    ! at least dist_AB, as d(A,C) and d(B,C) are; rounding can land it below
                    new_dist = max(scale(sqrt(new_dist_part_AC + new_dist_part_BC - new_dist_part_AB), -exponent), dist_AB)
                    distances(idx_B, i_node) = new_dist
                    distances(idx_A, i_node) = ieee_value(1.0_real64, ieee_positive_inf)
                    ! idx_B <= i_node < idx_A -> fill col of idx_B, row of idx_A
                else if (i_node < idx_A) then
                    ! scaled for the three distances this update combines
                    exponent = scaling_exponent_for_largest(max(distances(idx_A, i_node), distances(i_node, idx_B), dist_AB))
                    new_dist_part_AB = scale(dist_AB, exponent)**2*(real(size_C, real64)/size_ABC_real)
                    new_dist_part_AC = scale(distances(idx_A, i_node), exponent)**2*(size_AC_real/size_ABC_real)
                    new_dist_part_BC = scale(distances(i_node, idx_B), exponent)**2*(size_BC_real/size_ABC_real)
                    ! at least dist_AB, as d(A,C) and d(B,C) are; rounding can land it below
                    new_dist = max(scale(sqrt(new_dist_part_AC + new_dist_part_BC - new_dist_part_AB), -exponent), dist_AB)
                    distances(i_node, idx_B) = new_dist
                    distances(idx_A, i_node) = ieee_value(1.0_real64, ieee_positive_inf)
                    ! idx_B < idx_A < i_node -> fill columns in both
                else if (i_node > idx_A) then
                    ! scaled for the three distances this update combines
                    exponent = scaling_exponent_for_largest(max(distances(i_node, idx_A), distances(i_node, idx_B), dist_AB))
                    new_dist_part_AB = scale(dist_AB, exponent)**2*(real(size_C, real64)/size_ABC_real)
                    new_dist_part_AC = scale(distances(i_node, idx_A), exponent)**2*(size_AC_real/size_ABC_real)
                    new_dist_part_BC = scale(distances(i_node, idx_B), exponent)**2*(size_BC_real/size_ABC_real)
                    ! at least dist_AB, as d(A,C) and d(B,C) are; rounding can land it below
                    new_dist = max(scale(sqrt(new_dist_part_AC + new_dist_part_BC - new_dist_part_AB), -exponent), dist_AB)
                    distances(i_node, idx_B) = new_dist
                    ! don't fill idx_A with infinity, as whole column will be marked inactive
                end if
            end if
        end do

        call post_merge_distances_helper(distances, n_points, idx_A, idx_B, iteration_k)
    end subroutine merge_distances_ward_linkage_helper

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Helper for doing last updates after merge: disable merged cluster, update self-distance to current iteration
    pure subroutine post_merge_distances_helper(distances, n_points, idx_A, idx_B, iteration_k)
        integer(int32), intent(in) :: n_points
            !! number of points to cluster
        real(real64), dimension(n_points, n_points), intent(inout) :: distances
            !! distance matrix, holding the distances between points
        integer(int32), intent(in) :: idx_A
            !! row index of min distance in `distances`
        integer(int32), intent(in) :: idx_B
            !! column index of min distance in `distances`
        integer(int32), intent(in) :: iteration_k
            !! iteration of current merge

        ! Update index and cluster size
        distances(idx_B, idx_B) = real(iteration_k, real64)
        distances(idx_A, idx_B) = ieee_value(1.0_real64, ieee_positive_inf) ! column of B becomes new cluster -> A must not be future minimum
        distances(idx_A, idx_A) = -1.0_real64 ! mark column of A inactive
    end subroutine post_merge_distances_helper

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Helper routine to find the indices of the minimum value a distance matrix.
    !|
    !| - Expects `n>1`
    !| - Ignores columns with self-distance (i,i)<0
    !| - Searches in bottom triangle (excluding self-distance (i, i))
    !| - bottom triangle -> `row_idx > col_idx`
    pure subroutine get_min_distance_indices_helper(distances, n, row_idx, col_idx, min_dist, found)
        integer(int32), intent(in) :: n
            !! number of rows and columns of square matrix `distances`
        real(real64), dimension(n, n), intent(in) :: distances
            !! distance matrix, holding the distances between points
        integer(int32), intent(out) :: row_idx
            !! row index of min distance in `distances`
        integer(int32), intent(out) :: col_idx
            !! column index of min distance in `distances`
        real(real64), intent(out) :: min_dist
            !! min distance
        logical(c_bool), intent(out) :: found
            !! whether an entry below +Inf was found; if not, `row_idx` and `col_idx` are 0

        integer(int32) :: i_row, i_col

        min_dist = ieee_value(1.0_real64, ieee_positive_inf)
        row_idx = 0_int32
        col_idx = 0_int32

        do i_col = 1, n - 1
            if (distances(i_col, i_col) >= 0.0_real64) then
                do i_row = i_col + 1, n
                    if (distances(i_row, i_col) < min_dist) then
                        min_dist = distances(i_row, i_col)
                        row_idx = i_row
                        col_idx = i_col
                    end if
                end do
            end if
        end do
        found = row_idx /= 0_int32
    end subroutine get_min_distance_indices_helper
end module tox_clustering_impl
