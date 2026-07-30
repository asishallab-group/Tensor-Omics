// Generated. Do not edit.
#include <Rcpp.h>
#include <numeric>
#include <functional>
#include "tox_marshal.h"

using namespace Rcpp;

extern "C" {
    void compute_family_scaling_expert_c(const int*, const int*, const double*, const int*, double*, double*, double*, int*, int*, int*, int*, int*, const int*, double*, const int*, double*, double*, double*, double*, double*, double*, int*, double*, const double*, const int*, const char*, const int*, double*, int*, double*, int*);
    void compute_family_scaling_c(const int*, const int*, const double*, const int*, double*, double*, double*, int*, int*);
    void compute_rdi_c(const int*, const double*, const int*, const double*, const int*, double*, double*, int*, int*, int*, int*);
    void identify_outliers_c(const int*, const double*, const double*, const int*, bool*, double*, double*, const double*, int*);
    void detect_outliers_c(const int*, const int*, const double*, const int*, double*, int*, int*, int*, bool*, double*, double*, int*, double*, int*, const double*);
}

// [[Rcpp::export(.compute_family_scaling_expert_rcpp)]]
List compute_family_scaling_expert_rcpp(int n_families, NumericVector distances, IntegerVector gene_to_fam, int int_workspace_size, int real_workspace_size, double span, int degree, CharacterVector mode, int n_iters) {
    // derived from the inputs, not asked of the caller
    int n_genes = (int) distances.size();

    // convert what C cannot take directly
    tox::CharBuffer mode_c(mode, 6);

    // outputs and work space
    NumericVector dscale(n_families);
    NumericVector loess_x(n_families);
    NumericVector loess_y(n_families);
    IntegerVector indices_used(n_families);
    std::vector<int> tmp_perm(n_genes);
    std::vector<int> tmp_stack_left(n_genes);
    std::vector<int> tmp_stack_right(n_genes);
    std::vector<int> tmp_int_workspace(int_workspace_size);
    std::vector<double> tmp_real_workspace(real_workspace_size);
    std::vector<double> tmp_diagl(n_families);
    std::vector<double> tmp_weights(n_families);
    std::vector<double> tmp_eval_points(n_families * 1);
    std::vector<double> tmp_robust_weights(n_families);
    std::vector<double> tmp_combined_weights(n_families);
    std::vector<double> tmp_residuals(n_families);
    std::vector<int> tmp_permutation_indices(n_families);
    std::vector<double> tmp_fitted_values(n_families);
    double low_sd_cutoff = 0;
    IntegerVector excluded_low_sd(n_families);
    std::vector<double> tmp_means_aux(n_families);
    int ierr = 0;

    compute_family_scaling_expert_c(
        &n_genes,
        &n_families,
        distances.begin(),
        gene_to_fam.begin(),
        dscale.begin(),
        loess_x.begin(),
        loess_y.begin(),
        indices_used.begin(),
        tmp_perm.data(),
        tmp_stack_left.data(),
        tmp_stack_right.data(),
        tmp_int_workspace.data(),
        &int_workspace_size,
        tmp_real_workspace.data(),
        &real_workspace_size,
        tmp_diagl.data(),
        tmp_weights.data(),
        tmp_eval_points.data(),
        tmp_robust_weights.data(),
        tmp_combined_weights.data(),
        tmp_residuals.data(),
        tmp_permutation_indices.data(),
        tmp_fitted_values.data(),
        &span,
        &degree,
        mode_c.data(),
        &n_iters,
        &low_sd_cutoff,
        excluded_low_sd.begin(),
        tmp_means_aux.data(),
        &ierr
    );

    return List::create(
        _["dscale"] = dscale,
        _["loess_x"] = loess_x,
        _["loess_y"] = loess_y,
        _["indices_used"] = indices_used,
        _["low_sd_cutoff"] = low_sd_cutoff,
        _["excluded_low_sd"] = excluded_low_sd,
        _["ierr"] = ierr
    );
}

// [[Rcpp::export(.compute_family_scaling_rcpp)]]
List compute_family_scaling_rcpp(int n_families, NumericVector distances, IntegerVector gene_to_fam) {
    // derived from the inputs, not asked of the caller
    int n_genes = (int) distances.size();

    // outputs and work space
    NumericVector dscale(n_families);
    NumericVector loess_x(n_families);
    NumericVector loess_y(n_families);
    IntegerVector indices_used(n_families);
    int ierr = 0;

    compute_family_scaling_c(
        &n_genes,
        &n_families,
        distances.begin(),
        gene_to_fam.begin(),
        dscale.begin(),
        loess_x.begin(),
        loess_y.begin(),
        indices_used.begin(),
        &ierr
    );

    return List::create(
        _["dscale"] = dscale,
        _["loess_x"] = loess_x,
        _["loess_y"] = loess_y,
        _["indices_used"] = indices_used,
        _["ierr"] = ierr
    );
}

// [[Rcpp::export(.compute_rdi_rcpp)]]
List compute_rdi_rcpp(NumericVector distances, IntegerVector gene_to_fam, NumericVector dscale) {
    // derived from the inputs, not asked of the caller
    int n_genes = (int) distances.size();
    int n_dscale_elements = (int) dscale.size();

    // outputs and work space
    NumericVector rdi(n_genes);
    NumericVector sorted_rdi(n_genes);
    IntegerVector perm(n_genes);
    std::vector<int> tmp_stack_left(n_genes);
    std::vector<int> tmp_stack_right(n_genes);
    int ierr = 0;

    compute_rdi_c(
        &n_genes,
        distances.begin(),
        gene_to_fam.begin(),
        dscale.begin(),
        &n_dscale_elements,
        rdi.begin(),
        sorted_rdi.begin(),
        perm.begin(),
        tmp_stack_left.data(),
        tmp_stack_right.data(),
        &ierr
    );

    return List::create(
        _["rdi"] = rdi,
        _["sorted_rdi"] = sorted_rdi,
        _["perm"] = perm,
        _["ierr"] = ierr
    );
}

// [[Rcpp::export(.identify_outliers_rcpp)]]
List identify_outliers_rcpp(NumericVector rdi, NumericVector sorted_rdi, IntegerVector perm, double percentile) {
    // derived from the inputs, not asked of the caller
    int n_genes = (int) rdi.size();

    // outputs and work space
    tox::BoolBuffer is_outlier_c(n_genes);
    double threshold = 0;
    NumericVector quantile(n_genes);
    int ierr = 0;

    identify_outliers_c(
        &n_genes,
        rdi.begin(),
        sorted_rdi.begin(),
        perm.begin(),
        is_outlier_c.data(),
        &threshold,
        quantile.begin(),
        &percentile,
        &ierr
    );

    // convert the outputs back
    LogicalVector is_outlier = is_outlier_c.to_r();

    return List::create(
        _["is_outlier"] = is_outlier,
        _["threshold"] = threshold,
        _["quantile"] = quantile,
        _["ierr"] = ierr
    );
}

// [[Rcpp::export(.detect_outliers_rcpp)]]
List detect_outliers_rcpp(int n_families, NumericVector distances, IntegerVector gene_to_fam, double percentile) {
    // derived from the inputs, not asked of the caller
    int n_genes = (int) distances.size();

    // outputs and work space
    std::vector<double> tmp_work_array(n_genes);
    std::vector<int> tmp_perm(n_genes);
    std::vector<int> tmp_stack_left(n_genes);
    std::vector<int> tmp_stack_right(n_genes);
    tox::BoolBuffer is_outlier_c(n_genes);
    NumericVector loess_x(n_families);
    NumericVector loess_y(n_families);
    IntegerVector loess_n(n_families);
    NumericVector quantile(n_genes);
    int ierr = 0;

    detect_outliers_c(
        &n_genes,
        &n_families,
        distances.begin(),
        gene_to_fam.begin(),
        tmp_work_array.data(),
        tmp_perm.data(),
        tmp_stack_left.data(),
        tmp_stack_right.data(),
        is_outlier_c.data(),
        loess_x.begin(),
        loess_y.begin(),
        loess_n.begin(),
        quantile.begin(),
        &ierr,
        &percentile
    );

    // convert the outputs back
    LogicalVector is_outlier = is_outlier_c.to_r();

    return List::create(
        _["is_outlier"] = is_outlier,
        _["loess_x"] = loess_x,
        _["loess_y"] = loess_y,
        _["loess_n"] = loess_n,
        _["quantile"] = quantile,
        _["ierr"] = ierr
    );
}
