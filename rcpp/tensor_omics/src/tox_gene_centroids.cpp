// Generated. Do not edit.
#include <Rcpp.h>
#include <numeric>
#include <functional>
#include "tox_marshal.h"

using namespace Rcpp;

extern "C" {
    void mean_vector_c(const double*, const int*, const int*, const int*, const int*, double*, int*);
    void group_centroid_c(const double*, const int*, const int*, const int*, const int*, double*, const char*, int*, int*, const bool*);
}

// [[Rcpp::export(.mean_vector_rcpp)]]
List mean_vector_rcpp(NumericVector expression_vectors, IntegerVector gene_indices) {
    // derived from the inputs, not asked of the caller
    int n_axes = (int) IntegerVector(expression_vectors.attr("dim"))[0];
    int n_genes = (int) IntegerVector(expression_vectors.attr("dim"))[1];
    int n_selected_genes = (int) gene_indices.size();

    // outputs and work space
    NumericVector centroid(n_axes);
    int ierr = 0;

    mean_vector_c(
        expression_vectors.begin(),
        &n_axes,
        &n_genes,
        gene_indices.begin(),
        &n_selected_genes,
        centroid.begin(),
        &ierr
    );

    return List::create(
        _["centroid"] = centroid,
        _["ierr"] = ierr
    );
}

// [[Rcpp::export(.group_centroid_rcpp)]]
List group_centroid_rcpp(NumericVector expression_vectors, IntegerVector gene_to_family, int n_families, CharacterVector mode, Nullable<LogicalVector> ortholog_set = R_NilValue) {
    // optionals: a null pointer and size 0 when the caller omits them
    int ortholog_set_size = 0;
    LogicalVector ortholog_set_val;
    if (ortholog_set.isNotNull()) {
        ortholog_set_val = ortholog_set.get();
        ortholog_set_size = ortholog_set_val.size();
    }

    // derived from the inputs, not asked of the caller
    int n_axes = (int) IntegerVector(expression_vectors.attr("dim"))[0];
    int n_genes = (int) IntegerVector(expression_vectors.attr("dim"))[1];

    // convert what C cannot take directly
    tox::CharBuffer mode_c(mode, 15);
    tox::BoolBuffer ortholog_set_c(ortholog_set_val);

    // outputs and work space
    NumericVector centroid_matrix(n_axes * n_families);
    centroid_matrix.attr("dim") = IntegerVector::create(n_axes, n_families);
    std::vector<int> tmp_group_indices(n_genes);
    int ierr = 0;

    group_centroid_c(
        expression_vectors.begin(),
        &n_axes,
        &n_genes,
        gene_to_family.begin(),
        &n_families,
        centroid_matrix.begin(),
        mode_c.data(),
        tmp_group_indices.data(),
        &ierr,
        ortholog_set.isNotNull() ? ortholog_set_c.data() : nullptr
    );

    return List::create(
        _["centroid_matrix"] = centroid_matrix,
        _["ierr"] = ierr
    );
}
