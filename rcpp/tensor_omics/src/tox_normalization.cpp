// Generated. Do not edit.
#include <Rcpp.h>
#include <numeric>
#include <functional>
#include "tox_marshal.h"

using namespace Rcpp;

extern "C" {
    void normalize_unit_length_c(double*, const int*, int*);
    void normalization_pipeline_c(const int*, const int*, const double*, double*, const int*, const int*, const double*, const int*, const bool*, int*);
    void normalize_by_std_dev_c(const int*, const int*, const double*, double*, const double*, const int*, int*);
    void root_mean_sq_normalization_c(const int*, const int*, const double*, double*, int*);
    void quantile_normalization_c(const int*, const int*, const double*, double*, double*, double*, int*, int*);
    void log2_transformation_c(const int*, const int*, const double*, double*, int*);
    void calc_tiss_avg_c(const int*, const int*, const int*, const double*, double*, int*);
    void calc_fchange_c(const int*, const int*, const int*, const int*, const int*, const double*, double*, int*);
}

// [[Rcpp::export(.normalize_unit_length_rcpp)]]
List normalize_unit_length_rcpp(NumericVector vector) {
    // derived from the inputs, not asked of the caller
    int n_dims = (int) vector.size();

    // copy what is modified in place, so the caller's stays intact
    NumericVector vector_out = clone(vector);

    // outputs and work space
    int ierr = 0;

    normalize_unit_length_c(
        vector_out.begin(),
        &n_dims,
        &ierr
    );

    return List::create(
        _["vector"] = vector_out,
        _["ierr"] = ierr
    );
}

// [[Rcpp::export(.normalization_pipeline_rcpp)]]
List normalization_pipeline_rcpp(NumericVector expr, IntegerVector reps_per_tissue, double span, int degree, bool use_quantile) {
    // derived from the inputs, not asked of the caller
    int n_genes = (int) IntegerVector(expr.attr("dim"))[1];
    int n_replicates = (int) IntegerVector(expr.attr("dim"))[0];
    int n_tissues = (int) reps_per_tissue.size();

    // outputs and work space
    NumericVector log_transformed_expr(n_tissues * n_genes);
    log_transformed_expr.attr("dim") = IntegerVector::create(n_tissues, n_genes);
    int ierr = 0;

    normalization_pipeline_c(
        &n_genes,
        &n_replicates,
        expr.begin(),
        log_transformed_expr.begin(),
        reps_per_tissue.begin(),
        &n_tissues,
        &span,
        &degree,
        &use_quantile,
        &ierr
    );

    return List::create(
        _["log_transformed_expr"] = log_transformed_expr,
        _["ierr"] = ierr
    );
}

// [[Rcpp::export(.normalize_by_std_dev_rcpp)]]
List normalize_by_std_dev_rcpp(NumericVector expr, double span, int degree) {
    // derived from the inputs, not asked of the caller
    int n_genes = (int) IntegerVector(expr.attr("dim"))[1];
    int n_replicates = (int) IntegerVector(expr.attr("dim"))[0];

    // outputs and work space
    NumericVector normalized_expr(n_replicates * n_genes);
    normalized_expr.attr("dim") = IntegerVector::create(n_replicates, n_genes);
    int ierr = 0;

    normalize_by_std_dev_c(
        &n_genes,
        &n_replicates,
        expr.begin(),
        normalized_expr.begin(),
        &span,
        &degree,
        &ierr
    );

    return List::create(
        _["normalized_expr"] = normalized_expr,
        _["ierr"] = ierr
    );
}

// [[Rcpp::export(.root_mean_sq_normalization_rcpp)]]
List root_mean_sq_normalization_rcpp(NumericVector expr) {
    // derived from the inputs, not asked of the caller
    int n_genes = (int) IntegerVector(expr.attr("dim"))[1];
    int n_replicates = (int) IntegerVector(expr.attr("dim"))[0];

    // outputs and work space
    NumericVector normalized_expr(n_replicates * n_genes);
    normalized_expr.attr("dim") = IntegerVector::create(n_replicates, n_genes);
    int ierr = 0;

    root_mean_sq_normalization_c(
        &n_genes,
        &n_replicates,
        expr.begin(),
        normalized_expr.begin(),
        &ierr
    );

    return List::create(
        _["normalized_expr"] = normalized_expr,
        _["ierr"] = ierr
    );
}

// [[Rcpp::export(.quantile_normalization_rcpp)]]
List quantile_normalization_rcpp(NumericVector expr) {
    // derived from the inputs, not asked of the caller
    int n_genes = (int) IntegerVector(expr.attr("dim"))[1];
    int n_replicates = (int) IntegerVector(expr.attr("dim"))[0];

    // outputs and work space
    NumericVector normalized_expr(n_replicates * n_genes);
    normalized_expr.attr("dim") = IntegerVector::create(n_replicates, n_genes);
    NumericVector rank_means(n_genes);
    std::vector<double> tmp_genes_row(n_genes);
    std::vector<int> tmp_perm(n_genes);
    int ierr = 0;

    quantile_normalization_c(
        &n_genes,
        &n_replicates,
        expr.begin(),
        normalized_expr.begin(),
        rank_means.begin(),
        tmp_genes_row.data(),
        tmp_perm.data(),
        &ierr
    );

    return List::create(
        _["normalized_expr"] = normalized_expr,
        _["rank_means"] = rank_means,
        _["ierr"] = ierr
    );
}

// [[Rcpp::export(.log2_transformation_rcpp)]]
List log2_transformation_rcpp(NumericVector expr) {
    // derived from the inputs, not asked of the caller
    int n_genes = (int) IntegerVector(expr.attr("dim"))[1];
    int n_tissues = (int) IntegerVector(expr.attr("dim"))[0];

    // outputs and work space
    NumericVector transformed_expr(n_tissues * n_genes);
    transformed_expr.attr("dim") = IntegerVector::create(n_tissues, n_genes);
    int ierr = 0;

    log2_transformation_c(
        &n_genes,
        &n_tissues,
        expr.begin(),
        transformed_expr.begin(),
        &ierr
    );

    return List::create(
        _["transformed_expr"] = transformed_expr,
        _["ierr"] = ierr
    );
}

// [[Rcpp::export(.calc_tiss_avg_rcpp)]]
List calc_tiss_avg_rcpp(IntegerVector reps_per_tissue, NumericVector expr) {
    // derived from the inputs, not asked of the caller
    int n_genes = (int) IntegerVector(expr.attr("dim"))[1];
    int n_tissues = (int) reps_per_tissue.size();

    // outputs and work space
    NumericVector tissue_averages(n_tissues * n_genes);
    tissue_averages.attr("dim") = IntegerVector::create(n_tissues, n_genes);
    int ierr = 0;

    calc_tiss_avg_c(
        &n_genes,
        &n_tissues,
        reps_per_tissue.begin(),
        expr.begin(),
        tissue_averages.begin(),
        &ierr
    );

    return List::create(
        _["tissue_averages"] = tissue_averages,
        _["ierr"] = ierr
    );
}

// [[Rcpp::export(.calc_fchange_rcpp)]]
List calc_fchange_rcpp(IntegerVector control_tissues, IntegerVector condition_tissues, NumericVector expr) {
    // derived from the inputs, not asked of the caller
    int n_genes = (int) IntegerVector(expr.attr("dim"))[1];
    int n_tissues = (int) IntegerVector(expr.attr("dim"))[0];
    int n_pairs = (int) control_tissues.size();

    // outputs and work space
    NumericVector fold_changes(n_pairs * n_genes);
    fold_changes.attr("dim") = IntegerVector::create(n_pairs, n_genes);
    int ierr = 0;

    calc_fchange_c(
        &n_genes,
        &n_tissues,
        &n_pairs,
        control_tissues.begin(),
        condition_tissues.begin(),
        expr.begin(),
        fold_changes.begin(),
        &ierr
    );

    return List::create(
        _["fold_changes"] = fold_changes,
        _["ierr"] = ierr
    );
}
