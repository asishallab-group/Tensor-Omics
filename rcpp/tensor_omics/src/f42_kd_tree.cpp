// Generated. Do not edit.
#include <Rcpp.h>
#include <numeric>
#include <functional>
#include "tox_marshal.h"

using namespace Rcpp;

extern "C" {
    void build_kd_index_c(const double*, const int*, const int*, int*, const int*, int*);
    void build_spherical_kd_c(const double*, const int*, const int*, int*, const int*, int*);
}

// [[Rcpp::export(.build_kd_index_rcpp)]]
List build_kd_index_rcpp(NumericVector points, IntegerVector dimension_order) {
    // derived from the inputs, not asked of the caller
    int n_dimensions = (int) IntegerVector(points.attr("dim"))[0];
    int n_points = (int) IntegerVector(points.attr("dim"))[1];

    // outputs and work space
    IntegerVector kd_indices(n_points);
    int ierr = 0;

    build_kd_index_c(
        points.begin(),
        &n_dimensions,
        &n_points,
        kd_indices.begin(),
        dimension_order.begin(),
        &ierr
    );

    return List::create(
        _["kd_indices"] = kd_indices,
        _["ierr"] = ierr
    );
}

// [[Rcpp::export(.build_spherical_kd_rcpp)]]
List build_spherical_kd_rcpp(NumericVector points, IntegerVector dimension_order) {
    // derived from the inputs, not asked of the caller
    int n_dimensions = (int) IntegerVector(points.attr("dim"))[0];
    int n_points = (int) IntegerVector(points.attr("dim"))[1];

    // outputs and work space
    IntegerVector kd_indices(n_points);
    int ierr = 0;

    build_spherical_kd_c(
        points.begin(),
        &n_dimensions,
        &n_points,
        kd_indices.begin(),
        dimension_order.begin(),
        &ierr
    );

    return List::create(
        _["kd_indices"] = kd_indices,
        _["ierr"] = ierr
    );
}
