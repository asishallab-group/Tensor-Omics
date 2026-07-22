// Generated. Do not edit.
#include <Rcpp.h>
#include <numeric>
#include <functional>
#include "tox_marshal.h"

using namespace Rcpp;

extern "C" {
    void compute_shift_vector_field_c(const int*, const int*, const int*, const double*, const double*, const int*, double*, int*);
}

// [[Rcpp::export(.compute_shift_vector_field_rcpp)]]
List compute_shift_vector_field_rcpp(NumericVector expression_vectors, NumericVector family_centroids, IntegerVector gene_to_fam) {
    // derived from the inputs, not asked of the caller
    int n_tissues = (int) IntegerVector(expression_vectors.attr("dim"))[0];
    int n_genes = (int) IntegerVector(expression_vectors.attr("dim"))[1];
    int n_families = (int) IntegerVector(family_centroids.attr("dim"))[1];

    // outputs and work space
    NumericVector shift_vectors(n_tissues * 2 * n_genes);
    shift_vectors.attr("dim") = IntegerVector::create(n_tissues, 2, n_genes);
    int ierr = 0;

    compute_shift_vector_field_c(
        &n_tissues,
        &n_genes,
        &n_families,
        expression_vectors.begin(),
        family_centroids.begin(),
        gene_to_fam.begin(),
        shift_vectors.begin(),
        &ierr
    );

    return List::create(
        _["shift_vectors"] = shift_vectors,
        _["ierr"] = ierr
    );
}
