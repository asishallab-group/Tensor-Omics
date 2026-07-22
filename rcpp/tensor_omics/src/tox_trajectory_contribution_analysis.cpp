// Generated. Do not edit.
#include <Rcpp.h>
#include <numeric>
#include <functional>
#include "tox_marshal.h"

using namespace Rcpp;

extern "C" {
    void perform_permutation_test_c(const double*, const int*, const int*, const int*, const int*, const int*, const int*, const char*, const int*, double*, double*, double*, double*, int*, const int*);
    void compute_p_values_c(const double*, const double*, const double*, const double*, const int*, const int*, double*, double*, int*);
    void compute_contributions_c(const double*, const double*, const int*, const char*, double*, double*, int*);
    void compute_all_contributions_c(const double*, const int*, const int*, const int*, const int*, const int*, const int*, const int*, const char*, double*, double*, double*, double*, int*);
    void compute_baselines_factor_dependent_c(const int*, const double*, const double*, const char*, double*, double*, int*);
    void compute_velocity_trajectory_c(const double*, double*, const int*, int*);
    void compute_acceleration_from_velocity_trajectory_c(const double*, double*, const int*, int*);
    void compute_velocity_trajectories_c(const double*, double*, const int*, const int*, const int*, int*);
    void compute_acceleration_from_velocity_c(const double*, double*, const int*, const int*, const int*, int*);
    void compute_velocity_acceleration_contributions_expert_c(const double*, const int*, const int*, const int*, const char*, double*, double*, double*, double*, double*, double*, double*, int*);
    void compute_velocity_acceleration_contributions_c(const double*, const int*, const int*, const int*, const char*, double*, double*, double*, double*, int*);
}

// [[Rcpp::export(.perform_permutation_test_rcpp)]]
List perform_permutation_test_rcpp(NumericVector trajectories, int factor_idx, int dependent_idx, int sample_idx, CharacterVector baseline_mode, int n_permutations, Nullable<IntegerVector> random_seed = R_NilValue) {
    // optionals: a null pointer and size 0 when the caller omits them
    const int* random_seed_p = nullptr;
    int random_seed_size = 0;
    IntegerVector random_seed_val;
    if (random_seed.isNotNull()) {
        random_seed_val = random_seed.get();
        random_seed_p = random_seed_val.begin();
        random_seed_size = random_seed_val.size();
    }

    // derived from the inputs, not asked of the caller
    int n_factors = (int) IntegerVector(trajectories.attr("dim"))[0];
    int n_samples = (int) IntegerVector(trajectories.attr("dim"))[1];
    int n_timepoints = (int) IntegerVector(trajectories.attr("dim"))[2];

    // convert what C cannot take directly
    tox::CharBuffer baseline_mode_c(baseline_mode, 4);

    // outputs and work space
    NumericVector local_contributions(n_timepoints * n_permutations);
    local_contributions.attr("dim") = IntegerVector::create(n_timepoints, n_permutations);
    NumericVector total_contributions(n_permutations);
    std::vector<double> tmp_factor(n_timepoints);
    std::vector<double> tmp_dependent(n_timepoints);
    int ierr = 0;

    perform_permutation_test_c(
        trajectories.begin(),
        &n_factors,
        &n_samples,
        &n_timepoints,
        &factor_idx,
        &dependent_idx,
        &sample_idx,
        baseline_mode_c.data(),
        &n_permutations,
        local_contributions.begin(),
        total_contributions.begin(),
        tmp_factor.data(),
        tmp_dependent.data(),
        &ierr,
        random_seed_p
    );

    return List::create(
        _["local_contributions"] = local_contributions,
        _["total_contributions"] = total_contributions,
        _["ierr"] = ierr
    );
}

// [[Rcpp::export(.compute_p_values_rcpp)]]
List compute_p_values_rcpp(NumericVector local_contributions_observed, double total_contribution_observed, NumericVector local_contributions_perm_test, NumericVector total_contributions_perm_test) {
    // derived from the inputs, not asked of the caller
    int n_timepoints = (int) local_contributions_observed.size();
    int n_permutations = (int) IntegerVector(local_contributions_perm_test.attr("dim"))[1];

    // outputs and work space
    NumericVector local_p_values(n_timepoints);
    double total_p_value = 0;
    int ierr = 0;

    compute_p_values_c(
        local_contributions_observed.begin(),
        &total_contribution_observed,
        local_contributions_perm_test.begin(),
        total_contributions_perm_test.begin(),
        &n_timepoints,
        &n_permutations,
        local_p_values.begin(),
        &total_p_value,
        &ierr
    );

    return List::create(
        _["local_p_values"] = local_p_values,
        _["total_p_value"] = total_p_value,
        _["ierr"] = ierr
    );
}

// [[Rcpp::export(.compute_contributions_rcpp)]]
List compute_contributions_rcpp(NumericVector factor, NumericVector dependent, CharacterVector baseline_mode) {
    // derived from the inputs, not asked of the caller
    int n_dims = (int) factor.size();

    // convert what C cannot take directly
    tox::CharBuffer baseline_mode_c(baseline_mode, 4);

    // outputs and work space
    NumericVector local_contributions(n_dims);
    double total_contribution = 0;
    int ierr = 0;

    compute_contributions_c(
        factor.begin(),
        dependent.begin(),
        &n_dims,
        baseline_mode_c.data(),
        local_contributions.begin(),
        &total_contribution,
        &ierr
    );

    return List::create(
        _["local_contributions"] = local_contributions,
        _["total_contribution"] = total_contribution,
        _["ierr"] = ierr
    );
}

// [[Rcpp::export(.compute_all_contributions_rcpp)]]
List compute_all_contributions_rcpp(NumericVector trajectories, IntegerVector factor_indices, IntegerVector dependent_indices, CharacterVector baseline_mode) {
    // derived from the inputs, not asked of the caller
    int n_factors = (int) IntegerVector(trajectories.attr("dim"))[0];
    int n_samples = (int) IntegerVector(trajectories.attr("dim"))[1];
    int n_timepoints = (int) IntegerVector(trajectories.attr("dim"))[2];
    int n_selected_factors = (int) factor_indices.size();
    int n_selected_dependents = (int) dependent_indices.size();

    // convert what C cannot take directly
    tox::CharBuffer baseline_mode_c(baseline_mode, 4);

    // outputs and work space
    NumericVector local_contributions(n_timepoints * n_selected_factors * n_selected_dependents * n_samples);
    local_contributions.attr("dim") = IntegerVector::create(n_timepoints, n_selected_factors, n_selected_dependents, n_samples);
    NumericVector total_contributions(n_selected_factors * n_selected_dependents * n_samples);
    total_contributions.attr("dim") = IntegerVector::create(n_selected_factors, n_selected_dependents, n_samples);
    std::vector<double> tmp_factors(n_timepoints * n_selected_factors);
    std::vector<double> tmp_dependent(n_timepoints);
    int ierr = 0;

    compute_all_contributions_c(
        trajectories.begin(),
        &n_factors,
        &n_samples,
        &n_timepoints,
        factor_indices.begin(),
        &n_selected_factors,
        dependent_indices.begin(),
        &n_selected_dependents,
        baseline_mode_c.data(),
        local_contributions.begin(),
        total_contributions.begin(),
        tmp_factors.data(),
        tmp_dependent.data(),
        &ierr
    );

    return List::create(
        _["local_contributions"] = local_contributions,
        _["total_contributions"] = total_contributions,
        _["ierr"] = ierr
    );
}

// [[Rcpp::export(.compute_baselines_factor_dependent_rcpp)]]
List compute_baselines_factor_dependent_rcpp(NumericVector factor, NumericVector dependent, CharacterVector baseline_mode) {
    // derived from the inputs, not asked of the caller
    int n_timepoints = (int) factor.size();

    // convert what C cannot take directly
    tox::CharBuffer baseline_mode_c(baseline_mode, 4);

    // outputs and work space
    double factor_baseline = 0;
    double dependent_baseline = 0;
    int ierr = 0;

    compute_baselines_factor_dependent_c(
        &n_timepoints,
        factor.begin(),
        dependent.begin(),
        baseline_mode_c.data(),
        &factor_baseline,
        &dependent_baseline,
        &ierr
    );

    return List::create(
        _["factor_baseline"] = factor_baseline,
        _["dependent_baseline"] = dependent_baseline,
        _["ierr"] = ierr
    );
}

// [[Rcpp::export(.compute_velocity_trajectory_rcpp)]]
List compute_velocity_trajectory_rcpp(NumericVector trajectory) {
    // derived from the inputs, not asked of the caller
    int n_timepoints = (int) trajectory.size();

    // outputs and work space
    NumericVector velocity((std::max(0, n_timepoints - 1)));
    int ierr = 0;

    compute_velocity_trajectory_c(
        trajectory.begin(),
        velocity.begin(),
        &n_timepoints,
        &ierr
    );

    return List::create(
        _["velocity"] = velocity,
        _["ierr"] = ierr
    );
}

// [[Rcpp::export(.compute_acceleration_from_velocity_trajectory_rcpp)]]
List compute_acceleration_from_velocity_trajectory_rcpp(NumericVector velocity, int n_timepoints) {
    // outputs and work space
    NumericVector acceleration((std::max(0, n_timepoints - 2)));
    int ierr = 0;

    compute_acceleration_from_velocity_trajectory_c(
        velocity.begin(),
        acceleration.begin(),
        &n_timepoints,
        &ierr
    );

    return List::create(
        _["acceleration"] = acceleration,
        _["ierr"] = ierr
    );
}

// [[Rcpp::export(.compute_velocity_trajectories_rcpp)]]
List compute_velocity_trajectories_rcpp(NumericVector trajectories) {
    // derived from the inputs, not asked of the caller
    int n_factors = (int) IntegerVector(trajectories.attr("dim"))[0];
    int n_samples = (int) IntegerVector(trajectories.attr("dim"))[1];
    int n_timepoints = (int) IntegerVector(trajectories.attr("dim"))[2];

    // outputs and work space
    NumericVector velocity((std::max(0, n_timepoints - 1)) * n_factors * n_samples);
    velocity.attr("dim") = IntegerVector::create(std::max(0, n_timepoints - 1), n_factors, n_samples);
    int ierr = 0;

    compute_velocity_trajectories_c(
        trajectories.begin(),
        velocity.begin(),
        &n_factors,
        &n_samples,
        &n_timepoints,
        &ierr
    );

    return List::create(
        _["velocity"] = velocity,
        _["ierr"] = ierr
    );
}

// [[Rcpp::export(.compute_acceleration_from_velocity_rcpp)]]
List compute_acceleration_from_velocity_rcpp(NumericVector velocity, int n_timepoints) {
    // derived from the inputs, not asked of the caller
    int n_factors = (int) IntegerVector(velocity.attr("dim"))[1];
    int n_samples = (int) IntegerVector(velocity.attr("dim"))[2];

    // outputs and work space
    NumericVector acceleration((std::max(0, n_timepoints - 2)) * n_factors * n_samples);
    acceleration.attr("dim") = IntegerVector::create(std::max(0, n_timepoints - 2), n_factors, n_samples);
    int ierr = 0;

    compute_acceleration_from_velocity_c(
        velocity.begin(),
        acceleration.begin(),
        &n_factors,
        &n_samples,
        &n_timepoints,
        &ierr
    );

    return List::create(
        _["acceleration"] = acceleration,
        _["ierr"] = ierr
    );
}

// [[Rcpp::export(.compute_velocity_acceleration_contributions_expert_rcpp)]]
List compute_velocity_acceleration_contributions_expert_rcpp(NumericVector trajectories, CharacterVector baseline_mode) {
    // derived from the inputs, not asked of the caller
    int n_factors = (int) IntegerVector(trajectories.attr("dim"))[0];
    int n_samples = (int) IntegerVector(trajectories.attr("dim"))[1];
    int n_timepoints = (int) IntegerVector(trajectories.attr("dim"))[2];

    // convert what C cannot take directly
    tox::CharBuffer baseline_mode_c(baseline_mode, 4);

    // outputs and work space
    std::vector<double> tmp_factors((n_timepoints - 1) * n_factors);
    std::vector<double> tmp_dependent((n_timepoints - 1));
    std::vector<double> tmp_contributions((n_timepoints - 1));
    NumericVector contrib_velocity(n_factors * n_factors * n_samples);
    contrib_velocity.attr("dim") = IntegerVector::create(n_factors, n_factors, n_samples);
    NumericVector velocity_contribution_series(n_timepoints * n_factors * n_factors * n_samples);
    velocity_contribution_series.attr("dim") = IntegerVector::create(n_timepoints, n_factors, n_factors, n_samples);
    NumericVector contrib_acceleration(n_factors * n_factors * n_samples);
    contrib_acceleration.attr("dim") = IntegerVector::create(n_factors, n_factors, n_samples);
    NumericVector acceleration_contribution_series(n_timepoints * n_factors * n_factors * n_samples);
    acceleration_contribution_series.attr("dim") = IntegerVector::create(n_timepoints, n_factors, n_factors, n_samples);
    int ierr = 0;

    compute_velocity_acceleration_contributions_expert_c(
        trajectories.begin(),
        &n_factors,
        &n_samples,
        &n_timepoints,
        baseline_mode_c.data(),
        tmp_factors.data(),
        tmp_dependent.data(),
        tmp_contributions.data(),
        contrib_velocity.begin(),
        velocity_contribution_series.begin(),
        contrib_acceleration.begin(),
        acceleration_contribution_series.begin(),
        &ierr
    );

    return List::create(
        _["contrib_velocity"] = contrib_velocity,
        _["velocity_contribution_series"] = velocity_contribution_series,
        _["contrib_acceleration"] = contrib_acceleration,
        _["acceleration_contribution_series"] = acceleration_contribution_series,
        _["ierr"] = ierr
    );
}

// [[Rcpp::export(.compute_velocity_acceleration_contributions_rcpp)]]
List compute_velocity_acceleration_contributions_rcpp(NumericVector trajectories, CharacterVector baseline_mode) {
    // derived from the inputs, not asked of the caller
    int n_factors = (int) IntegerVector(trajectories.attr("dim"))[0];
    int n_samples = (int) IntegerVector(trajectories.attr("dim"))[1];
    int n_timepoints = (int) IntegerVector(trajectories.attr("dim"))[2];

    // convert what C cannot take directly
    tox::CharBuffer baseline_mode_c(baseline_mode, 4);

    // outputs and work space
    NumericVector contrib_velocity(n_factors * n_factors * n_samples);
    contrib_velocity.attr("dim") = IntegerVector::create(n_factors, n_factors, n_samples);
    NumericVector velocity_contribution_series(n_timepoints * n_factors * n_factors * n_samples);
    velocity_contribution_series.attr("dim") = IntegerVector::create(n_timepoints, n_factors, n_factors, n_samples);
    NumericVector contrib_acceleration(n_factors * n_factors * n_samples);
    contrib_acceleration.attr("dim") = IntegerVector::create(n_factors, n_factors, n_samples);
    NumericVector acceleration_contribution_series(n_timepoints * n_factors * n_factors * n_samples);
    acceleration_contribution_series.attr("dim") = IntegerVector::create(n_timepoints, n_factors, n_factors, n_samples);
    int ierr = 0;

    compute_velocity_acceleration_contributions_c(
        trajectories.begin(),
        &n_factors,
        &n_samples,
        &n_timepoints,
        baseline_mode_c.data(),
        contrib_velocity.begin(),
        velocity_contribution_series.begin(),
        contrib_acceleration.begin(),
        acceleration_contribution_series.begin(),
        &ierr
    );

    return List::create(
        _["contrib_velocity"] = contrib_velocity,
        _["velocity_contribution_series"] = velocity_contribution_series,
        _["contrib_acceleration"] = contrib_acceleration,
        _["acceleration_contribution_series"] = acceleration_contribution_series,
        _["ierr"] = ierr
    );
}
