// Generated. Do not edit.
#include <Rcpp.h>
#include <numeric>
#include <functional>
#include "tox_marshal.h"

using namespace Rcpp;

extern "C" {
    void serialize_logical_helper_c(const bool*, const int*, const int*, const int*, const char*, const int*, int*);
}

// [[Rcpp::export(.serialize_logical_helper_rcpp)]]
List serialize_logical_helper_rcpp(LogicalVector arr, CharacterVector filename) {
    // derived from the inputs, not asked of the caller
    IntegerVector arr_shape = Rf_isNull(arr.attr("dim")) ? IntegerVector::create(arr.size()) : IntegerVector(arr.attr("dim"));
    int n_elements = (int) arr.size();
    int n_arr_shape_elements = (int) arr_shape.size();
    int filename_strlen = tox::max_strlen(filename);

    // convert what C cannot take directly
    tox::BoolBuffer arr_c(arr);
    tox::CharBuffer filename_c(filename, filename_strlen);

    // outputs and work space
    int ierr = 0;

    serialize_logical_helper_c(
        arr_c.data(),
        &n_elements,
        arr_shape.begin(),
        &n_arr_shape_elements,
        filename_c.data(),
        &filename_strlen,
        &ierr
    );

    return List::create(
        _["ierr"] = ierr
    );
}
