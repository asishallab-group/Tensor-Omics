#include <Rcpp.h>

using namespace Rcpp;

static_assert(sizeof(Rcomplex) == sizeof(double _Complex) && alignof(Rcomplex) == alignof(double _Complex), "Rcomplex layout is incompatible with C's 'double _Complex' and thus Fortran's 'c_double_complex'");

extern "C" {
    void cluster_factor_trajectories_k_means_c(
        const int* n_clusters,
        const double* trajectories,
        const int* n_factors,
        const int* n_samples,
        const int* n_timepoints,
        double* centroids,
        int* labels,
        int* label_counts,
        int* ierr,
        const int* max_iterations
    );

    void k_means_clustering_c(
        const int* n_clusters,
        const double* data_points,
        const int* n_points,
        const int* n_dims,
        double* centroids,
        int* labels,
        int* label_counts,
        int* ierr,
        const int* max_iterations
    );

    void linkage_clustering_c(
        double* distances,
        const int* n_points,
        int* merge_i,
        int* merge_j,
        double* heights,
        int* cluster_sizes,
        const char* method,
        int* ierr
    );
}

List cluster_factor_trajectories_k_means_rcpp(
    const NumericVector& trajectories,
    NumericMatrix& centroids,
    const int max_iterations = 300
) {


    IntegerVector trajectories_shape = trajectories.attr("dim");
    n_factors = trajectories_shape[0]
    n_samples = trajectories_shape[1]
    n_timepoints = trajectories_shape[2]
    IntegerVector centroids_shape = centroids.attr("dim");
    n_clusters = centroids_shape[1]
    IntegerVector labels_shape = labels.attr("dim");
    n_samples * n_timepoints = labels_shape[0]



    IntegerVector labels(n_samples * n_timepoints);
    IntegerVector label_counts(n_clusters);
    int ierr = 0;

    cluster_factor_trajectories_k_means_c(
        &n_clusters,
        trajectories.begin(),
        &n_factors,
        &n_samples,
        &n_timepoints,
        centroids.begin(),
        labels.begin(),
        label_counts.begin(),
        &ierr,
        &max_iterations
    );

//{format(self.arguments, "type_conversion_outputs")}

    return List::create(
        Named("labels") = labels,
        Named("label_counts") = label_counts,
        Named("ierr") = ierr
    );
}

List k_means_clustering_rcpp(
    const NumericMatrix& data_points,
    NumericMatrix& centroids,
    const int max_iterations = 300
) {


    IntegerVector data_points_shape = data_points.attr("dim");
    n_dims = data_points_shape[0]
    n_points = data_points_shape[1]
    IntegerVector centroids_shape = centroids.attr("dim");
    n_clusters = centroids_shape[1]



    IntegerVector labels(n_points);
    IntegerVector label_counts(n_clusters);
    int ierr = 0;

    k_means_clustering_c(
        &n_clusters,
        data_points.begin(),
        &n_points,
        &n_dims,
        centroids.begin(),
        labels.begin(),
        label_counts.begin(),
        &ierr,
        &max_iterations
    );

//{format(self.arguments, "type_conversion_outputs")}

    return List::create(
        Named("labels") = labels,
        Named("label_counts") = label_counts,
        Named("ierr") = ierr
    );
}

List linkage_clustering_rcpp(
    NumericMatrix& distances,
    const String& method
) {


    IntegerVector distances_shape = distances.attr("dim");
    n_points = distances_shape[0]
    IntegerVector merge_i_shape = merge_i.attr("dim");
    n_points - 1 = merge_i_shape[0]

    if (method == NA_STRING) {
        return List::create(Named("ierr") = 204);
    }
    const char* method_p = method.get_cstring();

    IntegerVector merge_i(n_points - 1);
    IntegerVector merge_j(n_points - 1);
    NumericVector heights(n_points - 1);
    IntegerVector cluster_sizes(n_points - 1);
    int ierr = 0;

    linkage_clustering_c(
        distances.begin(),
        &n_points,
        merge_i.begin(),
        merge_j.begin(),
        heights.begin(),
        cluster_sizes.begin(),
        method_p,
        &ierr
    );

//{format(self.arguments, "type_conversion_outputs")}

    return List::create(
        Named("merge_i") = merge_i,
        Named("merge_j") = merge_j,
        Named("heights") = heights,
        Named("cluster_sizes") = cluster_sizes,
        Named("ierr") = ierr
    );
}