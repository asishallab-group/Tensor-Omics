// Generated. Do not edit.
#include <Rcpp.h>
#include <numeric>
#include <functional>
#include "tox_marshal.h"

using namespace Rcpp;

extern "C" {
    void tox_loess_required_workspace_c(const int*, const int*, int*, int*, const bool*, int*);
    void loess_fit_plain_c(const int*, const double*, const double*, const double*, const double*, const double*, const int*, const int*, const bool*, const bool*, int*, const int*, double*, const int*, double*, double*, int*);
    void loess_fit_robust_c(const int*, const double*, const double*, const double*, const double*, const double*, const int*, const int*, const bool*, const bool*, const int*, int*, const int*, double*, const int*, double*, double*, double*, double*, int*, double*, int*);
    void loess_c(const double*, const int*, const double*, const int*, const double*, const int*, double*, const char*, const int*, int*);
}

// [[Rcpp::export(.tox_loess_required_workspace_rcpp)]]
List tox_loess_required_workspace_rcpp(int n_dim, int max_neighborhood_size, bool save_factorization) {
    // outputs and work space
    int int_workspace_size = 0;
    int real_workspace_size = 0;
    int ierr = 0;

    tox_loess_required_workspace_c(
        &n_dim,
        &max_neighborhood_size,
        &int_workspace_size,
        &real_workspace_size,
        &save_factorization,
        &ierr
    );

    return List::create(
        _["int_workspace_size"] = int_workspace_size,
        _["real_workspace_size"] = real_workspace_size,
        _["ierr"] = ierr
    );
}

// [[Rcpp::export(.loess_fit_plain_rcpp)]]
List loess_fit_plain_rcpp(NumericVector x, NumericVector y, NumericVector weights, NumericVector eval_points, double span, int degree, int max_neighborhood_size, bool compute_influence, bool save_factorization, IntegerVector int_workspace, NumericVector real_workspace, NumericVector hat_diag) {
    // derived from the inputs, not asked of the caller
    int n = (int) x.size();
    int int_workspace_size = (int) int_workspace.size();
    int real_workspace_size = (int) real_workspace.size();

    // copy what is modified in place, so the caller's stays intact
    IntegerVector int_workspace_out = clone(int_workspace);
    NumericVector real_workspace_out = clone(real_workspace);
    NumericVector hat_diag_out = clone(hat_diag);

    // outputs and work space
    NumericVector fitted_values(n);
    int ierr = 0;

    loess_fit_plain_c(
        &n,
        x.begin(),
        y.begin(),
        weights.begin(),
        eval_points.begin(),
        &span,
        &degree,
        &max_neighborhood_size,
        &compute_influence,
        &save_factorization,
        int_workspace_out.begin(),
        &int_workspace_size,
        real_workspace_out.begin(),
        &real_workspace_size,
        hat_diag_out.begin(),
        fitted_values.begin(),
        &ierr
    );

    return List::create(
        _["int_workspace"] = int_workspace_out,
        _["real_workspace"] = real_workspace_out,
        _["hat_diag"] = hat_diag_out,
        _["fitted_values"] = fitted_values,
        _["ierr"] = ierr
    );
}

// [[Rcpp::export(.loess_fit_robust_rcpp)]]
List loess_fit_robust_rcpp(NumericVector x, NumericVector y, NumericVector weights, NumericVector eval_points, double span, int degree, int max_neighborhood_size, bool compute_influence, bool save_factorization, int n_iters, IntegerVector int_workspace, NumericVector real_workspace, NumericVector hat_diag, NumericVector robust_weights, NumericVector combined_weights, NumericVector residuals, IntegerVector permutation_indices) {
    // derived from the inputs, not asked of the caller
    int n = (int) x.size();
    int int_workspace_size = (int) int_workspace.size();
    int real_workspace_size = (int) real_workspace.size();

    // copy what is modified in place, so the caller's stays intact
    IntegerVector int_workspace_out = clone(int_workspace);
    NumericVector real_workspace_out = clone(real_workspace);
    NumericVector hat_diag_out = clone(hat_diag);
    NumericVector robust_weights_out = clone(robust_weights);
    NumericVector combined_weights_out = clone(combined_weights);
    NumericVector residuals_out = clone(residuals);
    IntegerVector permutation_indices_out = clone(permutation_indices);

    // outputs and work space
    NumericVector fitted_values(n);
    int ierr = 0;

    loess_fit_robust_c(
        &n,
        x.begin(),
        y.begin(),
        weights.begin(),
        eval_points.begin(),
        &span,
        &degree,
        &max_neighborhood_size,
        &compute_influence,
        &save_factorization,
        &n_iters,
        int_workspace_out.begin(),
        &int_workspace_size,
        real_workspace_out.begin(),
        &real_workspace_size,
        hat_diag_out.begin(),
        robust_weights_out.begin(),
        combined_weights_out.begin(),
        residuals_out.begin(),
        permutation_indices_out.begin(),
        fitted_values.begin(),
        &ierr
    );

    return List::create(
        _["int_workspace"] = int_workspace_out,
        _["real_workspace"] = real_workspace_out,
        _["hat_diag"] = hat_diag_out,
        _["robust_weights"] = robust_weights_out,
        _["combined_weights"] = combined_weights_out,
        _["residuals"] = residuals_out,
        _["permutation_indices"] = permutation_indices_out,
        _["fitted_values"] = fitted_values,
        _["ierr"] = ierr
    );
}

// [[Rcpp::export(.loess_rcpp)]]
List loess_rcpp(NumericVector x, NumericVector y, double span, int degree, CharacterVector mode, int n_iters) {
    // derived from the inputs, not asked of the caller
    int n_x_elements = (int) x.size();
    int n_y_elements = (int) y.size();

    // convert what C cannot take directly
    tox::CharBuffer mode_c(mode, 6);

    // outputs and work space
    NumericVector fitted_values(((int) y.size()));
    int ierr = 0;

    loess_c(
        x.begin(),
        &n_x_elements,
        y.begin(),
        &n_y_elements,
        &span,
        &degree,
        fitted_values.begin(),
        mode_c.data(),
        &n_iters,
        &ierr
    );

    return List::create(
        _["fitted_values"] = fitted_values,
        _["ierr"] = ierr
    );
}
