// Generated. Do not edit.
#include <Rcpp.h>
#include <numeric>
#include <functional>
#include "tox_marshal.h"

using namespace Rcpp;

extern "C" {
    void euclidean_distance_c(const double*, const double*, const int*, double*, int*);
    void distance_to_centroid_c(const int*, const int*, const double*, const double*, const int*, double*, const int*, int*);
}

// [[Rcpp::export(.euclidean_distance_rcpp)]]
List euclidean_distance_rcpp(NumericVector vec1, NumericVector vec2) {
    // derived from the inputs, not asked of the caller
    int n_elements = (int) vec1.size();

    // outputs and work space
    double result = 0;
    int ierr = 0;

    euclidean_distance_c(
        vec1.begin(),
        vec2.begin(),
        &n_elements,
        &result,
        &ierr
    );

    return List::create(
        _["result"] = result,
        _["ierr"] = ierr
    );
}

// [[Rcpp::export(.distance_to_centroid_rcpp)]]
List distance_to_centroid_rcpp(NumericVector genes, NumericVector centroids, IntegerVector gene_to_fam) {
    // derived from the inputs, not asked of the caller
    int n_genes = (int) IntegerVector(genes.attr("dim"))[1];
    int n_families = (int) IntegerVector(centroids.attr("dim"))[1];
    int n_tissues = (int) IntegerVector(genes.attr("dim"))[0];

    // outputs and work space
    NumericVector distances(n_genes);
    int ierr = 0;

    distance_to_centroid_c(
        &n_genes,
        &n_families,
        genes.begin(),
        centroids.begin(),
        gene_to_fam.begin(),
        distances.begin(),
        &n_tissues,
        &ierr
    );

    return List::create(
        _["distances"] = distances,
        _["ierr"] = ierr
    );
}
