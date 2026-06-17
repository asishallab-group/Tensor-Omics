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
    subroutine cluster_factor_trajectories_k_means_c(n_clusters, trajectories, n_factors, n_samples, n_timepoints, centroids, labels, label_counts, ierr, max_iterations) bind(C, name="cluster_factor_trajectories_k_means_c")
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
        M_CHECK_IERR_NON_NULL
        M_CHECK_NON_NULL(n_clusters)
        M_CHECK_NON_NULL(trajectories)
        M_CHECK_NON_NULL(n_factors)
        M_CHECK_NON_NULL(n_samples)
        M_CHECK_NON_NULL(n_timepoints)
        M_CHECK_NON_NULL(centroids)
        M_CHECK_NON_NULL(labels)
        M_CHECK_NON_NULL(label_counts)
        M_CHECK_NON_NULL(max_iterations)
        call cluster_factor_trajectories_k_means(n_clusters = n_clusters, trajectories = trajectories, n_factors = n_factors, n_samples = n_samples, n_timepoints = n_timepoints, centroids = centroids, labels = labels, label_counts = label_counts, ierr = ierr, max_iterations = max_iterations)
    end subroutine cluster_factor_trajectories_k_means_c

end module tox_clustering_c
#endif