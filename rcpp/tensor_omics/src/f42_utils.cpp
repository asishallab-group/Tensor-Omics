// Generated. Do not edit.
#include <Rcpp.h>
#include <numeric>
#include <functional>
#include "tox_marshal.h"

using namespace Rcpp;

extern "C" {
    void loess_smooth_2d_c(const int*, const int*, const double*, const double*, const int*, const int*, const double*, const double*, const double*, double*, int*);
    void compute_edf_expert_c(const double*, const int*, const int*, double*, double*, int*, int*);
    void compute_edf_c(const double*, const int*, double*, double*, int*, int*);
    void compute_scaled_distance_quantile_c(const int*, const double*, const double*, const int*, double*, const double*, int*);
}

// [[Rcpp::export(.loess_smooth_2d_rcpp)]]
List loess_smooth_2d_rcpp(NumericVector x_ref, NumericVector y_ref, IntegerVector indices_used, NumericVector x_query, double kernel_sigma, double kernel_cutoff) {
    // derived from the inputs, not asked of the caller
    int n_total = (int) x_ref.size();
    int n_target = (int) x_query.size();
    int n_used = (int) indices_used.size();

    // outputs and work space
    NumericVector y_out(n_target);
    int ierr = 0;

    loess_smooth_2d_c(
        &n_total,
        &n_target,
        x_ref.begin(),
        y_ref.begin(),
        indices_used.begin(),
        &n_used,
        x_query.begin(),
        &kernel_sigma,
        &kernel_cutoff,
        y_out.begin(),
        &ierr
    );

    return List::create(
        _["y_out"] = y_out,
        _["ierr"] = ierr
    );
}

// [[Rcpp::export(.compute_edf_expert_rcpp)]]
List compute_edf_expert_rcpp(NumericVector values, IntegerVector perm) {
    // derived from the inputs, not asked of the caller
    int n_values = (int) values.size();

    // outputs and work space
    NumericVector unique_values(n_values);
    NumericVector cdf_values(n_values);
    int n_unique = 0;
    int ierr = 0;

    compute_edf_expert_c(
        values.begin(),
        &n_values,
        perm.begin(),
        unique_values.begin(),
        cdf_values.begin(),
        &n_unique,
        &ierr
    );

    return List::create(
        _["unique_values"] = unique_values,
        _["cdf_values"] = cdf_values,
        _["n_unique"] = n_unique,
        _["ierr"] = ierr
    );
}

// [[Rcpp::export(.compute_edf_rcpp)]]
List compute_edf_rcpp(NumericVector values) {
    // derived from the inputs, not asked of the caller
    int n_values = (int) values.size();

    // outputs and work space
    NumericVector unique_values(n_values);
    NumericVector cdf_values(n_values);
    int n_unique = 0;
    int ierr = 0;

    compute_edf_c(
        values.begin(),
        &n_values,
        unique_values.begin(),
        cdf_values.begin(),
        &n_unique,
        &ierr
    );

    return List::create(
        _["unique_values"] = unique_values,
        _["cdf_values"] = cdf_values,
        _["n_unique"] = n_unique,
        _["ierr"] = ierr
    );
}

// [[Rcpp::export(.compute_scaled_distance_quantile_rcpp)]]
List compute_scaled_distance_quantile_rcpp(NumericVector rdi, NumericVector sorted_rdi, IntegerVector perm, double c_const) {
    // derived from the inputs, not asked of the caller
    int n_genes = (int) rdi.size();

    // outputs and work space
    NumericVector quantile(n_genes);
    int ierr = 0;

    compute_scaled_distance_quantile_c(
        &n_genes,
        rdi.begin(),
        sorted_rdi.begin(),
        perm.begin(),
        quantile.begin(),
        &c_const,
        &ierr
    );

    return List::create(
        _["quantile"] = quantile,
        _["ierr"] = ierr
    );
}
