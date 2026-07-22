// Generated. Do not edit.
#include <Rcpp.h>
#include <numeric>
#include <functional>
#include "tox_marshal.h"

using namespace Rcpp;

extern "C" {
    void get_array_metadata_c(const char*, const int*, int*, const int*, int*, int*, int*);
}

// [[Rcpp::export(.get_array_metadata_rcpp)]]
List get_array_metadata_rcpp(CharacterVector filename, int dims_out_capacity) {
    // derived from the inputs, not asked of the caller
    int filename_strlen = tox::max_strlen(filename);

    // convert what C cannot take directly
    tox::CharBuffer filename_c(filename, filename_strlen);

    // outputs and work space
    IntegerVector dims_out(dims_out_capacity);
    int ndims = 0;
    int type_code = 0;
    int ierr = 0;

    get_array_metadata_c(
        filename_c.data(),
        &filename_strlen,
        dims_out.begin(),
        &dims_out_capacity,
        &ndims,
        &type_code,
        &ierr
    );

    return List::create(
        _["dims_out"] = dims_out,
        _["ndims"] = ndims,
        _["type_code"] = type_code,
        _["ierr"] = ierr
    );
}
