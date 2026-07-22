// Generated. Do not edit.
#include <Rcpp.h>
#include <numeric>
#include <functional>
#include "tox_marshal.h"

using namespace Rcpp;

extern "C" {
    void deserialize_char_helper_c(char*, const int*, const int*, const int*, const int*, const char*, const int*, int*);
}

// [[Rcpp::export(.deserialize_char_helper_rcpp)]]
List deserialize_char_helper_rcpp(int strlen, IntegerVector arr_shape, CharacterVector filename) {
    // derived from the inputs, not asked of the caller
    int n_strings = (int) std::accumulate(arr_shape.begin(), arr_shape.end(), 1, std::multiplies<int>());
    int n_arr_shape_elements = (int) arr_shape.size();
    int filename_strlen = tox::max_strlen(filename);

    // convert what C cannot take directly
    tox::CharBuffer filename_c(filename, filename_strlen);

    // outputs and work space
    tox::CharBuffer arr_c(strlen, (int) std::accumulate(arr_shape.begin(), arr_shape.end(), 1, std::multiplies<int>()));
    int ierr = 0;

    deserialize_char_helper_c(
        arr_c.data(),
        &n_strings,
        &strlen,
        arr_shape.begin(),
        &n_arr_shape_elements,
        filename_c.data(),
        &filename_strlen,
        &ierr
    );

    // convert the outputs back
    CharacterVector arr = arr_c.to_r();
    arr.attr("dim") = arr_shape;

    return List::create(
        _["arr"] = arr,
        _["ierr"] = ierr
    );
}
