// Generated. Do not edit.
#include <Rcpp.h>
#include <numeric>
#include <functional>
#include "tox_marshal.h"

using namespace Rcpp;

extern "C" {
    void deserialize_complex_helper_c(double _Complex*, const int*, const int*, const int*, const char*, const int*, int*);
}

// [[Rcpp::export(.deserialize_complex_helper_rcpp)]]
List deserialize_complex_helper_rcpp(IntegerVector arr_shape, CharacterVector filename) {
    // derived from the inputs, not asked of the caller
    int n_elements = (int) std::accumulate(arr_shape.begin(), arr_shape.end(), 1, std::multiplies<int>());
    int n_arr_shape_elements = (int) arr_shape.size();
    int filename_strlen = tox::max_strlen(filename);

    // convert what C cannot take directly
    tox::CharBuffer filename_c(filename, filename_strlen);

    // outputs and work space
    ComplexVector arr((int) std::accumulate(arr_shape.begin(), arr_shape.end(), 1, std::multiplies<int>()));
    arr.attr("dim") = arr_shape;
    int ierr = 0;

    deserialize_complex_helper_c(
        reinterpret_cast<double _Complex*>(arr.begin()),
        &n_elements,
        arr_shape.begin(),
        &n_arr_shape_elements,
        filename_c.data(),
        &filename_strlen,
        &ierr
    );

    return List::create(
        _["arr"] = arr,
        _["ierr"] = ierr
    );
}
