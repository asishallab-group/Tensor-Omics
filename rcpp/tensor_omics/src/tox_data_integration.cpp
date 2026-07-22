// Generated. Do not edit.
#include <Rcpp.h>
#include <numeric>
#include <functional>
#include "tox_marshal.h"

using namespace Rcpp;

extern "C" {
    void compute_gene_means_c(const int*, const int*, const double*, double*, int*);
    void compute_residuals_c(const int*, const int*, const double*, const double*, double*, int*);
    void pool_means_c(const int*, const double*, const int*, const double*, const int*, int*, double*, int*);
    void pool_means_expert_c(const double*, const int*, const int*, const int*, int*, double*, int*);
    void calc_neighborhood_size_c(const int*, const int*, const int*, const double*, const int*, int*, int*);
    void construct_neighborhoods_c(const int*, const double*, const int*, const double*, const int*, const double*, double*, int*, const int*, int*);
    void gjct_permutation_test_c(const double*, const double*, const int*, const int*, const int*, const int*, const double*, const int*, const double*, const int*, double*, double*, int*, const int*, const bool*, const bool*);
    void gjct_permutation_test_expert_c(double*, double*, const int*, const int*, const int*, const int*, const double*, const int*, const double*, const int*, double*, double*, double*, double*, double*, int*, int*, int*, double*, double*, int*, const int*, const bool*, const bool*);
    void determine_shared_residual_range_expert_c(const double*, const int*, const int*, double*, int*, const double*);
    void determine_shared_residual_range_c(const double*, const double*, const int*, const int*, const int*, const int*, double*, int*, const double*);
    void build_residual_histograms_c(const double*, const int*, const int*, const int*, const double*, const int*, int*, double*, int*, int*, const bool*);
    void compute_divergence_per_reference_point_c(const double*, const double*, const int*, const int*, double*, int*);
    void compute_weighted_global_divergence_c(const double*, const int*, const int*, const int*, double*, double*, int*);
    void fjct_compute_jsd_c(const int*, const int*, const int*, const int*, const int*, const double*, const double*, const int*, const int*, const int*, const int*, const int*, const int*, const int*, const double*, double*, int*, int*, int*, double*, double*, int*);
    void fjct_compute_jsd_expert_c(const double*, const double*, const int*, const int*, const int*, const int*, const bool*, const bool*, const int*, const double*, double*, int*, int*, int*, double*, double*, double*, double*, int*, int*);
    void fjct_compute_contribution_scores_c(const double*, const int*, const int*, double*, double*, int*);
}

// [[Rcpp::export(.compute_gene_means_rcpp)]]
List compute_gene_means_rcpp(NumericVector expr) {
    // derived from the inputs, not asked of the caller
    int n_genes = (int) IntegerVector(expr.attr("dim"))[1];
    int n_reps = (int) IntegerVector(expr.attr("dim"))[0];

    // outputs and work space
    NumericVector means(n_genes);
    int ierr = 0;

    compute_gene_means_c(
        &n_genes,
        &n_reps,
        expr.begin(),
        means.begin(),
        &ierr
    );

    return List::create(
        _["means"] = means,
        _["ierr"] = ierr
    );
}

// [[Rcpp::export(.compute_residuals_rcpp)]]
List compute_residuals_rcpp(NumericVector expr, NumericVector means) {
    // derived from the inputs, not asked of the caller
    int n_genes = (int) IntegerVector(expr.attr("dim"))[1];
    int n_reps = (int) IntegerVector(expr.attr("dim"))[0];

    // outputs and work space
    NumericVector resid(n_reps * n_genes);
    resid.attr("dim") = IntegerVector::create(n_reps, n_genes);
    int ierr = 0;

    compute_residuals_c(
        &n_genes,
        &n_reps,
        expr.begin(),
        means.begin(),
        resid.begin(),
        &ierr
    );

    return List::create(
        _["resid"] = resid,
        _["ierr"] = ierr
    );
}

// [[Rcpp::export(.pool_means_rcpp)]]
List pool_means_rcpp(NumericVector mean_S1, NumericVector mean_S2, int n_points) {
    // derived from the inputs, not asked of the caller
    int n_genes_S1 = (int) mean_S1.size();
    int n_genes_S2 = (int) mean_S2.size();

    // outputs and work space
    int n_pool = 0;
    NumericVector x_star(n_points);
    int ierr = 0;

    pool_means_c(
        &n_genes_S1,
        mean_S1.begin(),
        &n_genes_S2,
        mean_S2.begin(),
        &n_points,
        &n_pool,
        x_star.begin(),
        &ierr
    );

    return List::create(
        _["n_pool"] = n_pool,
        _["x_star"] = x_star,
        _["ierr"] = ierr
    );
}

// [[Rcpp::export(.pool_means_expert_rcpp)]]
List pool_means_expert_rcpp(NumericVector pooled_means, IntegerVector pooled_means_perm, int n_points) {
    // derived from the inputs, not asked of the caller
    int pool_size = (int) pooled_means.size();

    // outputs and work space
    int n_pool = 0;
    NumericVector x_star(n_points);
    int ierr = 0;

    pool_means_expert_c(
        pooled_means.begin(),
        pooled_means_perm.begin(),
        &pool_size,
        &n_points,
        &n_pool,
        x_star.begin(),
        &ierr
    );

    return List::create(
        _["n_pool"] = n_pool,
        _["x_star"] = x_star,
        _["ierr"] = ierr
    );
}

// [[Rcpp::export(.calc_neighborhood_size_rcpp)]]
List calc_neighborhood_size_rcpp(int n_pool, int n_points, NumericVector mean_S, int desired_size) {
    // derived from the inputs, not asked of the caller
    int n_genes_S = (int) mean_S.size();

    // outputs and work space
    int n_neighbors = 0;
    int ierr = 0;

    calc_neighborhood_size_c(
        &n_pool,
        &n_points,
        &n_genes_S,
        mean_S.begin(),
        &desired_size,
        &n_neighbors,
        &ierr
    );

    return List::create(
        _["n_neighbors"] = n_neighbors,
        _["ierr"] = ierr
    );
}

// [[Rcpp::export(.construct_neighborhoods_rcpp)]]
List construct_neighborhoods_rcpp(NumericVector x_star, NumericVector mean_S, NumericVector resid_S, int n_neighbors) {
    // derived from the inputs, not asked of the caller
    int n_points = (int) x_star.size();
    int n_genes_S = (int) mean_S.size();
    int n_reps_S = (int) IntegerVector(resid_S.attr("dim"))[0];

    // outputs and work space
    NumericVector neighborhood_residuals(n_reps_S * n_neighbors * n_points);
    neighborhood_residuals.attr("dim") = IntegerVector::create(n_reps_S, n_neighbors, n_points);
    IntegerVector neighborhood_indices(n_neighbors * n_points);
    neighborhood_indices.attr("dim") = IntegerVector::create(n_neighbors, n_points);
    int ierr = 0;

    construct_neighborhoods_c(
        &n_points,
        x_star.begin(),
        &n_genes_S,
        mean_S.begin(),
        &n_reps_S,
        resid_S.begin(),
        neighborhood_residuals.begin(),
        neighborhood_indices.begin(),
        &n_neighbors,
        &ierr
    );

    return List::create(
        _["neighborhood_residuals"] = neighborhood_residuals,
        _["neighborhood_indices"] = neighborhood_indices,
        _["ierr"] = ierr
    );
}

// [[Rcpp::export(.gjct_permutation_test_rcpp)]]
List gjct_permutation_test_rcpp(NumericVector neighborhood_residuals_S1, NumericVector neighborhood_residuals_S2, double global_jsd_observed, int n_bins, double shared_residual_range, int n_permutations, Nullable<IntegerVector> random_seed = R_NilValue, Nullable<LogicalVector> neighbor_mask_S1 = R_NilValue, Nullable<LogicalVector> neighbor_mask_S2 = R_NilValue) {
    // optionals: a null pointer and size 0 when the caller omits them
    const int* random_seed_p = nullptr;
    int random_seed_size = 0;
    IntegerVector random_seed_val;
    if (random_seed.isNotNull()) {
        random_seed_val = random_seed.get();
        random_seed_p = random_seed_val.begin();
        random_seed_size = random_seed_val.size();
    }
    int neighbor_mask_S1_size = 0;
    std::vector<int> neighbor_mask_S1_dim(2, 0);
    LogicalVector neighbor_mask_S1_val;
    if (neighbor_mask_S1.isNotNull()) {
        neighbor_mask_S1_val = neighbor_mask_S1.get();
        neighbor_mask_S1_size = neighbor_mask_S1_val.size();
        if (!Rf_isNull(neighbor_mask_S1_val.attr("dim"))) {
            IntegerVector neighbor_mask_S1_d(neighbor_mask_S1_val.attr("dim"));
            for (int i = 0; i < 2 && i < neighbor_mask_S1_d.size(); ++i) neighbor_mask_S1_dim[i] = neighbor_mask_S1_d[i];
        }
    }
    int neighbor_mask_S2_size = 0;
    std::vector<int> neighbor_mask_S2_dim(2, 0);
    LogicalVector neighbor_mask_S2_val;
    if (neighbor_mask_S2.isNotNull()) {
        neighbor_mask_S2_val = neighbor_mask_S2.get();
        neighbor_mask_S2_size = neighbor_mask_S2_val.size();
        if (!Rf_isNull(neighbor_mask_S2_val.attr("dim"))) {
            IntegerVector neighbor_mask_S2_d(neighbor_mask_S2_val.attr("dim"));
            for (int i = 0; i < 2 && i < neighbor_mask_S2_d.size(); ++i) neighbor_mask_S2_dim[i] = neighbor_mask_S2_d[i];
        }
    }

    // derived from the inputs, not asked of the caller
    int n_reps_S1 = (int) IntegerVector(neighborhood_residuals_S1.attr("dim"))[0];
    int n_reps_S2 = (int) IntegerVector(neighborhood_residuals_S2.attr("dim"))[0];
    int n_neighbors = (int) IntegerVector(neighborhood_residuals_S1.attr("dim"))[1];
    int n_points = (int) IntegerVector(neighborhood_residuals_S1.attr("dim"))[2];

    // convert what C cannot take directly
    tox::BoolBuffer neighbor_mask_S1_c(neighbor_mask_S1_val);
    tox::BoolBuffer neighbor_mask_S2_c(neighbor_mask_S2_val);

    // outputs and work space
    NumericVector jsd_null(n_permutations);
    double p_value = 0;
    int ierr = 0;

    gjct_permutation_test_c(
        neighborhood_residuals_S1.begin(),
        neighborhood_residuals_S2.begin(),
        &n_reps_S1,
        &n_reps_S2,
        &n_neighbors,
        &n_points,
        &global_jsd_observed,
        &n_bins,
        &shared_residual_range,
        &n_permutations,
        jsd_null.begin(),
        &p_value,
        &ierr,
        random_seed_p,
        neighbor_mask_S1.isNotNull() ? neighbor_mask_S1_c.data() : nullptr,
        neighbor_mask_S2.isNotNull() ? neighbor_mask_S2_c.data() : nullptr
    );

    return List::create(
        _["jsd_null"] = jsd_null,
        _["p_value"] = p_value,
        _["ierr"] = ierr
    );
}

// [[Rcpp::export(.gjct_permutation_test_expert_rcpp)]]
List gjct_permutation_test_expert_rcpp(NumericVector neighborhood_residuals_S1_copy, NumericVector neighborhood_residuals_S2_copy, double global_jsd_observed, int n_bins, double shared_residual_range, int n_permutations, Nullable<IntegerVector> random_seed = R_NilValue, Nullable<LogicalVector> neighbor_mask_S1 = R_NilValue, Nullable<LogicalVector> neighbor_mask_S2 = R_NilValue) {
    // optionals: a null pointer and size 0 when the caller omits them
    const int* random_seed_p = nullptr;
    int random_seed_size = 0;
    IntegerVector random_seed_val;
    if (random_seed.isNotNull()) {
        random_seed_val = random_seed.get();
        random_seed_p = random_seed_val.begin();
        random_seed_size = random_seed_val.size();
    }
    int neighbor_mask_S1_size = 0;
    std::vector<int> neighbor_mask_S1_dim(2, 0);
    LogicalVector neighbor_mask_S1_val;
    if (neighbor_mask_S1.isNotNull()) {
        neighbor_mask_S1_val = neighbor_mask_S1.get();
        neighbor_mask_S1_size = neighbor_mask_S1_val.size();
        if (!Rf_isNull(neighbor_mask_S1_val.attr("dim"))) {
            IntegerVector neighbor_mask_S1_d(neighbor_mask_S1_val.attr("dim"));
            for (int i = 0; i < 2 && i < neighbor_mask_S1_d.size(); ++i) neighbor_mask_S1_dim[i] = neighbor_mask_S1_d[i];
        }
    }
    int neighbor_mask_S2_size = 0;
    std::vector<int> neighbor_mask_S2_dim(2, 0);
    LogicalVector neighbor_mask_S2_val;
    if (neighbor_mask_S2.isNotNull()) {
        neighbor_mask_S2_val = neighbor_mask_S2.get();
        neighbor_mask_S2_size = neighbor_mask_S2_val.size();
        if (!Rf_isNull(neighbor_mask_S2_val.attr("dim"))) {
            IntegerVector neighbor_mask_S2_d(neighbor_mask_S2_val.attr("dim"));
            for (int i = 0; i < 2 && i < neighbor_mask_S2_d.size(); ++i) neighbor_mask_S2_dim[i] = neighbor_mask_S2_d[i];
        }
    }

    // derived from the inputs, not asked of the caller
    int n_reps_S1 = (int) IntegerVector(neighborhood_residuals_S1_copy.attr("dim"))[0];
    int n_reps_S2 = (int) IntegerVector(neighborhood_residuals_S2_copy.attr("dim"))[0];
    int n_neighbors = (int) IntegerVector(neighborhood_residuals_S1_copy.attr("dim"))[1];
    int n_points = (int) IntegerVector(neighborhood_residuals_S1_copy.attr("dim"))[2];

    // copy what is modified in place, so the caller's stays intact
    NumericVector neighborhood_residuals_S1_copy_out = clone(neighborhood_residuals_S1_copy);
    NumericVector neighborhood_residuals_S2_copy_out = clone(neighborhood_residuals_S2_copy);

    // convert what C cannot take directly
    tox::BoolBuffer neighbor_mask_S1_c(neighbor_mask_S1_val);
    tox::BoolBuffer neighbor_mask_S2_c(neighbor_mask_S2_val);

    // outputs and work space
    NumericVector jsd_null(n_permutations);
    double p_value = 0;
    std::vector<double> tmp_pool((n_reps_S1 + n_reps_S2) * n_neighbors);
    std::vector<double> tmp_pmf_S1(n_points * n_bins);
    std::vector<double> tmp_pmf_S2(n_points * n_bins);
    std::vector<int> tmp_counts(n_points * n_bins);
    std::vector<int> tmp_included_n_reps_S1(n_points);
    std::vector<int> tmp_included_n_reps_S2(n_points);
    std::vector<double> tmp_js_divergences(n_points);
    std::vector<double> tmp_weights(n_points);
    int ierr = 0;

    gjct_permutation_test_expert_c(
        neighborhood_residuals_S1_copy_out.begin(),
        neighborhood_residuals_S2_copy_out.begin(),
        &n_reps_S1,
        &n_reps_S2,
        &n_neighbors,
        &n_points,
        &global_jsd_observed,
        &n_bins,
        &shared_residual_range,
        &n_permutations,
        jsd_null.begin(),
        &p_value,
        tmp_pool.data(),
        tmp_pmf_S1.data(),
        tmp_pmf_S2.data(),
        tmp_counts.data(),
        tmp_included_n_reps_S1.data(),
        tmp_included_n_reps_S2.data(),
        tmp_js_divergences.data(),
        tmp_weights.data(),
        &ierr,
        random_seed_p,
        neighbor_mask_S1.isNotNull() ? neighbor_mask_S1_c.data() : nullptr,
        neighbor_mask_S2.isNotNull() ? neighbor_mask_S2_c.data() : nullptr
    );

    return List::create(
        _["neighborhood_residuals_S1_copy"] = neighborhood_residuals_S1_copy_out,
        _["neighborhood_residuals_S2_copy"] = neighborhood_residuals_S2_copy_out,
        _["jsd_null"] = jsd_null,
        _["p_value"] = p_value,
        _["ierr"] = ierr
    );
}

// [[Rcpp::export(.determine_shared_residual_range_expert_rcpp)]]
List determine_shared_residual_range_expert_rcpp(NumericVector abs_residual_pool, IntegerVector abs_residual_pool_perm, double residual_range_quantile) {
    // derived from the inputs, not asked of the caller
    int pool_size = (int) abs_residual_pool.size();

    // outputs and work space
    double shared_residual_range = 0;
    int ierr = 0;

    determine_shared_residual_range_expert_c(
        abs_residual_pool.begin(),
        abs_residual_pool_perm.begin(),
        &pool_size,
        &shared_residual_range,
        &ierr,
        &residual_range_quantile
    );

    return List::create(
        _["shared_residual_range"] = shared_residual_range,
        _["ierr"] = ierr
    );
}

// [[Rcpp::export(.determine_shared_residual_range_rcpp)]]
List determine_shared_residual_range_rcpp(NumericVector neighborhood_residuals_S1, NumericVector neighborhood_residuals_S2, double residual_range_quantile) {
    // derived from the inputs, not asked of the caller
    int n_reps_S1 = (int) IntegerVector(neighborhood_residuals_S1.attr("dim"))[0];
    int n_reps_S2 = (int) IntegerVector(neighborhood_residuals_S2.attr("dim"))[0];
    int n_neighbors = (int) IntegerVector(neighborhood_residuals_S1.attr("dim"))[1];
    int n_points = (int) IntegerVector(neighborhood_residuals_S1.attr("dim"))[2];

    // outputs and work space
    double shared_residual_range = 0;
    int ierr = 0;

    determine_shared_residual_range_c(
        neighborhood_residuals_S1.begin(),
        neighborhood_residuals_S2.begin(),
        &n_reps_S1,
        &n_reps_S2,
        &n_neighbors,
        &n_points,
        &shared_residual_range,
        &ierr,
        &residual_range_quantile
    );

    return List::create(
        _["shared_residual_range"] = shared_residual_range,
        _["ierr"] = ierr
    );
}

// [[Rcpp::export(.build_residual_histograms_rcpp)]]
List build_residual_histograms_rcpp(NumericVector neighborhood_residuals, double shared_residual_range, int n_bins, Nullable<LogicalVector> neighbor_mask = R_NilValue) {
    // optionals: a null pointer and size 0 when the caller omits them
    int neighbor_mask_size = 0;
    std::vector<int> neighbor_mask_dim(2, 0);
    LogicalVector neighbor_mask_val;
    if (neighbor_mask.isNotNull()) {
        neighbor_mask_val = neighbor_mask.get();
        neighbor_mask_size = neighbor_mask_val.size();
        if (!Rf_isNull(neighbor_mask_val.attr("dim"))) {
            IntegerVector neighbor_mask_d(neighbor_mask_val.attr("dim"));
            for (int i = 0; i < 2 && i < neighbor_mask_d.size(); ++i) neighbor_mask_dim[i] = neighbor_mask_d[i];
        }
    }

    // derived from the inputs, not asked of the caller
    int n_reps = (int) IntegerVector(neighborhood_residuals.attr("dim"))[0];
    int n_neighbors = (int) IntegerVector(neighborhood_residuals.attr("dim"))[1];
    int n_points = (int) IntegerVector(neighborhood_residuals.attr("dim"))[2];

    // convert what C cannot take directly
    tox::BoolBuffer neighbor_mask_c(neighbor_mask_val);

    // outputs and work space
    IntegerVector counts(n_points * n_bins);
    counts.attr("dim") = IntegerVector::create(n_points, n_bins);
    NumericVector pmf(n_points * n_bins);
    pmf.attr("dim") = IntegerVector::create(n_points, n_bins);
    IntegerVector included_n_reps(n_points);
    int ierr = 0;

    build_residual_histograms_c(
        neighborhood_residuals.begin(),
        &n_reps,
        &n_neighbors,
        &n_points,
        &shared_residual_range,
        &n_bins,
        counts.begin(),
        pmf.begin(),
        included_n_reps.begin(),
        &ierr,
        neighbor_mask.isNotNull() ? neighbor_mask_c.data() : nullptr
    );

    return List::create(
        _["counts"] = counts,
        _["pmf"] = pmf,
        _["included_n_reps"] = included_n_reps,
        _["ierr"] = ierr
    );
}

// [[Rcpp::export(.compute_divergence_per_reference_point_rcpp)]]
List compute_divergence_per_reference_point_rcpp(NumericVector pmf_S1, NumericVector pmf_S2) {
    // derived from the inputs, not asked of the caller
    int n_points = (int) IntegerVector(pmf_S1.attr("dim"))[0];
    int n_bins = (int) IntegerVector(pmf_S1.attr("dim"))[1];

    // outputs and work space
    NumericVector js_divergences(n_points);
    int ierr = 0;

    compute_divergence_per_reference_point_c(
        pmf_S1.begin(),
        pmf_S2.begin(),
        &n_points,
        &n_bins,
        js_divergences.begin(),
        &ierr
    );

    return List::create(
        _["js_divergences"] = js_divergences,
        _["ierr"] = ierr
    );
}

// [[Rcpp::export(.compute_weighted_global_divergence_rcpp)]]
List compute_weighted_global_divergence_rcpp(NumericVector js_divergences, IntegerVector included_n_reps_S1, IntegerVector included_n_reps_S2) {
    // derived from the inputs, not asked of the caller
    int n_points = (int) js_divergences.size();

    // outputs and work space
    double global_js_divergence = 0;
    NumericVector weights(n_points);
    int ierr = 0;

    compute_weighted_global_divergence_c(
        js_divergences.begin(),
        &n_points,
        included_n_reps_S1.begin(),
        included_n_reps_S2.begin(),
        &global_js_divergence,
        weights.begin(),
        &ierr
    );

    return List::create(
        _["global_js_divergence"] = global_js_divergence,
        _["weights"] = weights,
        _["ierr"] = ierr
    );
}

// [[Rcpp::export(.fjct_compute_jsd_rcpp)]]
List fjct_compute_jsd_rcpp(int family_idx, IntegerVector gene_to_family_S1, IntegerVector gene_to_family_S2, NumericVector neighborhood_residuals_S1, NumericVector neighborhood_residuals_S2, IntegerVector neighborhood_genes_S1, IntegerVector neighborhood_genes_S2, int n_bins, double shared_residual_range) {
    // derived from the inputs, not asked of the caller
    int n_genes_S1 = (int) gene_to_family_S1.size();
    int n_genes_S2 = (int) gene_to_family_S2.size();
    int n_reps_S1 = (int) IntegerVector(neighborhood_residuals_S1.attr("dim"))[0];
    int n_reps_S2 = (int) IntegerVector(neighborhood_residuals_S2.attr("dim"))[0];
    int n_neighbors = (int) IntegerVector(neighborhood_residuals_S1.attr("dim"))[1];
    int n_points = (int) IntegerVector(neighborhood_residuals_S1.attr("dim"))[2];

    // outputs and work space
    NumericVector js_divergences(n_points);
    IntegerVector included_n_reps_S1(n_points);
    IntegerVector included_n_reps_S2(n_points);
    int total_included_n_reps = 0;
    double global_js_divergence = 0;
    NumericVector weights(n_points);
    int ierr = 0;

    fjct_compute_jsd_c(
        &family_idx,
        gene_to_family_S1.begin(),
        gene_to_family_S2.begin(),
        &n_genes_S1,
        &n_genes_S2,
        neighborhood_residuals_S1.begin(),
        neighborhood_residuals_S2.begin(),
        neighborhood_genes_S1.begin(),
        neighborhood_genes_S2.begin(),
        &n_reps_S1,
        &n_reps_S2,
        &n_neighbors,
        &n_points,
        &n_bins,
        &shared_residual_range,
        js_divergences.begin(),
        included_n_reps_S1.begin(),
        included_n_reps_S2.begin(),
        &total_included_n_reps,
        &global_js_divergence,
        weights.begin(),
        &ierr
    );

    return List::create(
        _["js_divergences"] = js_divergences,
        _["included_n_reps_S1"] = included_n_reps_S1,
        _["included_n_reps_S2"] = included_n_reps_S2,
        _["total_included_n_reps"] = total_included_n_reps,
        _["global_js_divergence"] = global_js_divergence,
        _["weights"] = weights,
        _["ierr"] = ierr
    );
}

// [[Rcpp::export(.fjct_compute_jsd_expert_rcpp)]]
List fjct_compute_jsd_expert_rcpp(NumericVector neighborhood_residuals_S1, NumericVector neighborhood_residuals_S2, LogicalVector neighbor_mask_S1, LogicalVector neighbor_mask_S2, int n_bins, double shared_residual_range) {
    // derived from the inputs, not asked of the caller
    int n_reps_S1 = (int) IntegerVector(neighborhood_residuals_S1.attr("dim"))[0];
    int n_reps_S2 = (int) IntegerVector(neighborhood_residuals_S2.attr("dim"))[0];
    int n_neighbors = (int) IntegerVector(neighborhood_residuals_S1.attr("dim"))[1];
    int n_points = (int) IntegerVector(neighborhood_residuals_S1.attr("dim"))[2];

    // convert what C cannot take directly
    tox::BoolBuffer neighbor_mask_S1_c(neighbor_mask_S1);
    tox::BoolBuffer neighbor_mask_S2_c(neighbor_mask_S2);

    // outputs and work space
    NumericVector js_divergences(n_points);
    IntegerVector included_n_reps_S1(n_points);
    IntegerVector included_n_reps_S2(n_points);
    int total_included_n_reps = 0;
    double global_js_divergence = 0;
    NumericVector weights(n_points);
    NumericVector pmf_S1(n_points * n_bins);
    pmf_S1.attr("dim") = IntegerVector::create(n_points, n_bins);
    NumericVector pmf_S2(n_points * n_bins);
    pmf_S2.attr("dim") = IntegerVector::create(n_points, n_bins);
    std::vector<int> tmp_counts(n_points * n_bins);
    int ierr = 0;

    fjct_compute_jsd_expert_c(
        neighborhood_residuals_S1.begin(),
        neighborhood_residuals_S2.begin(),
        &n_reps_S1,
        &n_reps_S2,
        &n_neighbors,
        &n_points,
        neighbor_mask_S1_c.data(),
        neighbor_mask_S2_c.data(),
        &n_bins,
        &shared_residual_range,
        js_divergences.begin(),
        included_n_reps_S1.begin(),
        included_n_reps_S2.begin(),
        &total_included_n_reps,
        &global_js_divergence,
        weights.begin(),
        pmf_S1.begin(),
        pmf_S2.begin(),
        tmp_counts.data(),
        &ierr
    );

    return List::create(
        _["js_divergences"] = js_divergences,
        _["included_n_reps_S1"] = included_n_reps_S1,
        _["included_n_reps_S2"] = included_n_reps_S2,
        _["total_included_n_reps"] = total_included_n_reps,
        _["global_js_divergence"] = global_js_divergence,
        _["weights"] = weights,
        _["pmf_S1"] = pmf_S1,
        _["pmf_S2"] = pmf_S2,
        _["ierr"] = ierr
    );
}

// [[Rcpp::export(.fjct_compute_contribution_scores_rcpp)]]
List fjct_compute_contribution_scores_rcpp(NumericVector global_js_divergences, IntegerVector total_included_n_reps_per_f) {
    // derived from the inputs, not asked of the caller
    int k_families = (int) global_js_divergences.size();

    // outputs and work space
    NumericVector support_weights(k_families);
    NumericVector contribution_scores(k_families);
    int ierr = 0;

    fjct_compute_contribution_scores_c(
        global_js_divergences.begin(),
        total_included_n_reps_per_f.begin(),
        &k_families,
        support_weights.begin(),
        contribution_scores.begin(),
        &ierr
    );

    return List::create(
        _["support_weights"] = support_weights,
        _["contribution_scores"] = contribution_scores,
        _["ierr"] = ierr
    );
}
