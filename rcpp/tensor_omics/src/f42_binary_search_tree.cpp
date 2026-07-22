// Generated. Do not edit.
#include <Rcpp.h>
#include <numeric>
#include <functional>
#include "tox_marshal.h"

using namespace Rcpp;

extern "C" {
    void build_bst_index_c(const double*, const int*, int*, int*);
    void bst_range_query_c(const double*, const int*, const int*, const double*, const double*, int*, int*, int*);
}

// [[Rcpp::export(.build_bst_index_rcpp)]]
List build_bst_index_rcpp(NumericVector values) {
    // derived from the inputs, not asked of the caller
    int n_values = (int) values.size();

    // outputs and work space
    IntegerVector sorted_indices(n_values);
    int ierr = 0;

    build_bst_index_c(
        values.begin(),
        &n_values,
        sorted_indices.begin(),
        &ierr
    );

    return List::create(
        _["sorted_indices"] = sorted_indices,
        _["ierr"] = ierr
    );
}

// [[Rcpp::export(.bst_range_query_rcpp)]]
List bst_range_query_rcpp(NumericVector values, IntegerVector sorted_indices, double lower_bound, double upper_bound) {
    // derived from the inputs, not asked of the caller
    int n_values = (int) values.size();

    // outputs and work space
    IntegerVector output_indices(n_values);
    int n_matches = 0;
    int ierr = 0;

    bst_range_query_c(
        values.begin(),
        sorted_indices.begin(),
        &n_values,
        &lower_bound,
        &upper_bound,
        output_indices.begin(),
        &n_matches,
        &ierr
    );

    return List::create(
        _["output_indices"] = output_indices,
        _["n_matches"] = n_matches,
        _["ierr"] = ierr
    );
}
