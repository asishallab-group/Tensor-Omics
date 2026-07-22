// Generated. Do not edit.
#include <Rcpp.h>
#include <numeric>
#include <functional>
#include "tox_marshal.h"

using namespace Rcpp;

extern "C" {
    void cluster_factor_trajectories_k_means_c(const int*, const double*, const int*, const int*, const int*, double*, int*, int*, int*, const int*);
    void k_means_clustering_c(const int*, const double*, const int*, const int*, double*, int*, int*, int*, const int*);
    void linkage_clustering_c(double*, const int*, int*, int*, double*, int*, const char*, int*);
}

// [[Rcpp::export(.cluster_factor_trajectories_k_means_rcpp)]]
List cluster_factor_trajectories_k_means_rcpp(NumericVector trajectories, NumericVector centroids, int max_iterations) {
    // derived from the inputs, not asked of the caller
    int n_clusters = (int) IntegerVector(centroids.attr("dim"))[1];
    int n_factors = (int) IntegerVector(trajectories.attr("dim"))[0];
    int n_samples = (int) IntegerVector(trajectories.attr("dim"))[1];
    int n_timepoints = (int) IntegerVector(trajectories.attr("dim"))[2];

    // copy what is modified in place, so the caller's stays intact
    NumericVector centroids_out = clone(centroids);

    // outputs and work space
    IntegerVector labels((n_samples*n_timepoints));
    IntegerVector label_counts(n_clusters);
    int ierr = 0;

    cluster_factor_trajectories_k_means_c(
        &n_clusters,
        trajectories.begin(),
        &n_factors,
        &n_samples,
        &n_timepoints,
        centroids_out.begin(),
        labels.begin(),
        label_counts.begin(),
        &ierr,
        &max_iterations
    );

    return List::create(
        _["centroids"] = centroids_out,
        _["labels"] = labels,
        _["label_counts"] = label_counts,
        _["ierr"] = ierr
    );
}

// [[Rcpp::export(.k_means_clustering_rcpp)]]
List k_means_clustering_rcpp(NumericVector data_points, NumericVector centroids, int max_iterations) {
    // derived from the inputs, not asked of the caller
    int n_clusters = (int) IntegerVector(centroids.attr("dim"))[1];
    int n_points = (int) IntegerVector(data_points.attr("dim"))[1];
    int n_dims = (int) IntegerVector(data_points.attr("dim"))[0];

    // copy what is modified in place, so the caller's stays intact
    NumericVector centroids_out = clone(centroids);

    // outputs and work space
    IntegerVector labels(n_points);
    IntegerVector label_counts(n_clusters);
    int ierr = 0;

    k_means_clustering_c(
        &n_clusters,
        data_points.begin(),
        &n_points,
        &n_dims,
        centroids_out.begin(),
        labels.begin(),
        label_counts.begin(),
        &ierr,
        &max_iterations
    );

    return List::create(
        _["centroids"] = centroids_out,
        _["labels"] = labels,
        _["label_counts"] = label_counts,
        _["ierr"] = ierr
    );
}

// [[Rcpp::export(.linkage_clustering_rcpp)]]
List linkage_clustering_rcpp(NumericVector distances, CharacterVector method) {
    // derived from the inputs, not asked of the caller
    int n_points = (int) IntegerVector(distances.attr("dim"))[0];

    // copy what is modified in place, so the caller's stays intact
    NumericVector distances_out = clone(distances);

    // convert what C cannot take directly
    tox::CharBuffer method_c(method, 8);

    // outputs and work space
    IntegerVector merge_i((n_points - 1));
    IntegerVector merge_j((n_points - 1));
    NumericVector heights((n_points - 1));
    IntegerVector cluster_sizes((n_points - 1));
    int ierr = 0;

    linkage_clustering_c(
        distances_out.begin(),
        &n_points,
        merge_i.begin(),
        merge_j.begin(),
        heights.begin(),
        cluster_sizes.begin(),
        method_c.data(),
        &ierr
    );

    return List::create(
        _["distances"] = distances_out,
        _["merge_i"] = merge_i,
        _["merge_j"] = merge_j,
        _["heights"] = heights,
        _["cluster_sizes"] = cluster_sizes,
        _["ierr"] = ierr
    );
}
