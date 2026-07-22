// Generated. Do not edit.
#include <Rcpp.h>
#include <numeric>
#include <functional>
#include "tox_marshal.h"

using namespace Rcpp;

extern "C" {
    void normalize_variable_timeseries_c(const double*, double*, const int*, int*, int*);
    void normalize_single_trajectory_c(const double*, double*, const int*, const int*, int*, int*);
    void normalize_all_trajectories_c(const double*, double*, const int*, const int*, const int*, int*, int*);
}

// [[Rcpp::export(.normalize_variable_timeseries_rcpp)]]
List normalize_variable_timeseries_rcpp(NumericVector v) {
    // derived from the inputs, not asked of the caller
    int n_points = (int) v.size();

    // outputs and work space
    NumericVector v_norm(n_points);
    int ierr = 0;
    int status = 0;

    normalize_variable_timeseries_c(
        v.begin(),
        v_norm.begin(),
        &n_points,
        &ierr,
        &status
    );

    return List::create(
        _["v_norm"] = v_norm,
        _["status"] = status,
        _["ierr"] = ierr
    );
}

// [[Rcpp::export(.normalize_single_trajectory_rcpp)]]
List normalize_single_trajectory_rcpp(NumericVector trajectory) {
    // derived from the inputs, not asked of the caller
    int n_factors = (int) IntegerVector(trajectory.attr("dim"))[1];
    int n_timepoints = (int) IntegerVector(trajectory.attr("dim"))[0];

    // outputs and work space
    NumericVector trajectory_norm(n_timepoints * n_factors);
    trajectory_norm.attr("dim") = IntegerVector::create(n_timepoints, n_factors);
    int ierr = 0;
    IntegerVector status(n_factors);

    normalize_single_trajectory_c(
        trajectory.begin(),
        trajectory_norm.begin(),
        &n_factors,
        &n_timepoints,
        &ierr,
        status.begin()
    );

    return List::create(
        _["trajectory_norm"] = trajectory_norm,
        _["status"] = status,
        _["ierr"] = ierr
    );
}

// [[Rcpp::export(.normalize_all_trajectories_rcpp)]]
List normalize_all_trajectories_rcpp(NumericVector trajectories) {
    // derived from the inputs, not asked of the caller
    int n_factors = (int) IntegerVector(trajectories.attr("dim"))[0];
    int n_samples = (int) IntegerVector(trajectories.attr("dim"))[1];
    int n_timepoints = (int) IntegerVector(trajectories.attr("dim"))[2];

    // outputs and work space
    NumericVector trajectories_norm(n_factors * n_samples * n_timepoints);
    trajectories_norm.attr("dim") = IntegerVector::create(n_factors, n_samples, n_timepoints);
    int ierr = 0;
    IntegerVector status(n_factors * n_samples);
    status.attr("dim") = IntegerVector::create(n_factors, n_samples);

    normalize_all_trajectories_c(
        trajectories.begin(),
        trajectories_norm.begin(),
        &n_factors,
        &n_samples,
        &n_timepoints,
        &ierr,
        status.begin()
    );

    return List::create(
        _["trajectories_norm"] = trajectories_norm,
        _["status"] = status,
        _["ierr"] = ierr
    );
}
