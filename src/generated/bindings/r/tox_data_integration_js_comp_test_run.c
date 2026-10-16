// Generated. Do not edit.
#if !defined(NO_R_BINDING) && !defined(NO_C_BINDING)
#include <R.h>
#include <Rinternals.h>
#include "tox_marshal.h"
// tox_marshal.h 0e1e7c507a726932 -- its hash, so that fpm, which only hashes this file, recompiles it when the header changes

// the Fortran C-ABI symbols this module calls
void run_js_comp_test_c(const int*, const int*, const int*, const int*, const int*, const double*, const double*, const double*, int*, int*, int*, double*, double*, int*, unsigned char*, int*, int*, double*, int*, int*, int*, double*, int*, int*, double*, int*, int*, double*, double*, double*, double*, double*, const int*, const int*, const int*, const int*, const int*, const double*, const double*, const double*, const char*, int*);
void run_js_comp_test_expert_c(const int*, const int*, const int*, const int*, const int*, const double*, const int*, const double*, const double*, int*, int*, int*, double*, double*, int*, unsigned char*, int*, int*, double*, int*, int*, int*, double*, int*, int*, double*, int*, int*, double*, double*, double*, double*, double*, double*, double*, int*, int*, double*, int*, int*, int*, int*, int*, int*, int*, int*, double*, double*, double*, double*, const int*, const int*, const int*, const int*, const int*, const double*, const double*, const double*, const char*, int*);
void run_js_comp_test_adaptive_c(const int*, const int*, const int*, const int*, const double*, const double*, const int*, int*, int*, double*, double*, int*, unsigned char*, int*, int*, double*, int*, int*, int*, double*, int*, int*, double*, int*, int*, double*, double*, double*, double*, double*, const int*, const int*, const int*, const int*, const int*, const double*, const double*, const double*, const char*, int*);

SEXP run_js_comp_test_call(SEXP n_neighbors, SEXP gene_means, SEXP residuals, SEXP x_star, SEXP n_permutations, SEXP random_seed, SEXP min_residuals_per_bin, SEXP m_min, SEXP m_max, SEXP gamma_occupancy, SEXP lower_residual_range_quantile, SEXP upper_residual_range_quantile, SEXP bin_count_mode) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_studies = INTEGER(Rf_getAttrib(gene_means, R_DimSymbol))[1];
    int max_n_genes_all_studies = INTEGER(Rf_getAttrib(gene_means, R_DimSymbol))[0];
    int max_n_reps_all_studies = INTEGER(Rf_getAttrib(residuals, R_DimSymbol))[0];
    int n_points = (int) Rf_length(x_star);

    // scalar inputs, pulled from their length-1 vectors
    int n_neighbors_v = Rf_asInteger(n_neighbors);
    int n_permutations_v = Rf_asInteger(n_permutations);
    int random_seed_v = Rf_asInteger(random_seed);
    int min_residuals_per_bin_v = Rf_asInteger(min_residuals_per_bin);
    int m_min_v = Rf_asInteger(m_min);
    int m_max_v = Rf_asInteger(m_max);
    double gamma_occupancy_v = Rf_asReal(gamma_occupancy);
    double lower_residual_range_quantile_v = Rf_asReal(lower_residual_range_quantile);
    double upper_residual_range_quantile_v = Rf_asReal(upper_residual_range_quantile);

    // convert what Fortran cannot take from R directly
    char* bin_count_mode_c = tox_char_in(bin_count_mode, 29);

    // outputs and work space
    SEXP neighborhood_indices = PROTECT(Rf_allocVector(INTSXP, n_neighbors_v * n_points * n_studies)); nprot++;
    { SEXP neighborhood_indices_dim = PROTECT(Rf_allocVector(INTSXP, 3)); INTEGER(neighborhood_indices_dim)[0] = n_neighbors_v; INTEGER(neighborhood_indices_dim)[1] = n_points; INTEGER(neighborhood_indices_dim)[2] = n_studies; Rf_setAttrib(neighborhood_indices, R_DimSymbol, neighborhood_indices_dim); UNPROTECT(1); }
    SEXP neighborhood_range = PROTECT(Rf_allocVector(INTSXP, 2 * n_points * n_studies)); nprot++;
    { SEXP neighborhood_range_dim = PROTECT(Rf_allocVector(INTSXP, 3)); INTEGER(neighborhood_range_dim)[0] = 2; INTEGER(neighborhood_range_dim)[1] = n_points; INTEGER(neighborhood_range_dim)[2] = n_studies; Rf_setAttrib(neighborhood_range, R_DimSymbol, neighborhood_range_dim); UNPROTECT(1); }
    SEXP n_bins_per_point = PROTECT(Rf_allocVector(INTSXP, n_points)); nprot++;
    SEXP shared_residual_range_low = PROTECT(Rf_allocVector(REALSXP, n_points)); nprot++;
    SEXP shared_residual_range_high = PROTECT(Rf_allocVector(REALSXP, n_points)); nprot++;
    int max_n_bins_per_point = 0;
    unsigned char* occupancy_failed_c = tox_bool_alloc(n_points);
    SEXP n_pooled_residuals = PROTECT(Rf_allocVector(INTSXP, n_points)); nprot++;
    SEXP min_bin_occupancy = PROTECT(Rf_allocVector(INTSXP, n_points)); nprot++;
    SEXP mean_bin_occupancy = PROTECT(Rf_allocVector(REALSXP, n_points)); nprot++;
    SEXP max_bin_occupancy = PROTECT(Rf_allocVector(INTSXP, n_points)); nprot++;
    SEXP sturges_bins = PROTECT(Rf_allocVector(INTSXP, n_points)); nprot++;
    SEXP fd_bins = PROTECT(Rf_allocVector(INTSXP, n_points)); nprot++;
    SEXP pmfs = PROTECT(Rf_allocVector(REALSXP, n_points * 256 * n_studies)); nprot++;
    { SEXP pmfs_dim = PROTECT(Rf_allocVector(INTSXP, 3)); INTEGER(pmfs_dim)[0] = n_points; INTEGER(pmfs_dim)[1] = 256; INTEGER(pmfs_dim)[2] = n_studies; Rf_setAttrib(pmfs, R_DimSymbol, pmfs_dim); UNPROTECT(1); }
    SEXP counts = PROTECT(Rf_allocVector(INTSXP, n_points * 256 * n_studies)); nprot++;
    { SEXP counts_dim = PROTECT(Rf_allocVector(INTSXP, 3)); INTEGER(counts_dim)[0] = n_points; INTEGER(counts_dim)[1] = 256; INTEGER(counts_dim)[2] = n_studies; Rf_setAttrib(counts, R_DimSymbol, counts_dim); UNPROTECT(1); }
    SEXP included_n_reps = PROTECT(Rf_allocVector(INTSXP, n_points * n_studies)); nprot++;
    { SEXP included_n_reps_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(included_n_reps_dim)[0] = n_points; INTEGER(included_n_reps_dim)[1] = n_studies; Rf_setAttrib(included_n_reps, R_DimSymbol, included_n_reps_dim); UNPROTECT(1); }
    SEXP mean_pmf = PROTECT(Rf_allocVector(REALSXP, n_points * 256)); nprot++;
    { SEXP mean_pmf_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(mean_pmf_dim)[0] = n_points; INTEGER(mean_pmf_dim)[1] = 256; Rf_setAttrib(mean_pmf, R_DimSymbol, mean_pmf_dim); UNPROTECT(1); }
    SEXP mean_pmf_counts = PROTECT(Rf_allocVector(INTSXP, n_points * 256)); nprot++;
    { SEXP mean_pmf_counts_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(mean_pmf_counts_dim)[0] = n_points; INTEGER(mean_pmf_counts_dim)[1] = 256; Rf_setAttrib(mean_pmf_counts, R_DimSymbol, mean_pmf_counts_dim); UNPROTECT(1); }
    SEXP mean_pmf_included_n_reps = PROTECT(Rf_allocVector(INTSXP, n_points)); nprot++;
    SEXP js_divergences = PROTECT(Rf_allocVector(REALSXP, n_points * n_studies)); nprot++;
    { SEXP js_divergences_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(js_divergences_dim)[0] = n_points; INTEGER(js_divergences_dim)[1] = n_studies; Rf_setAttrib(js_divergences, R_DimSymbol, js_divergences_dim); UNPROTECT(1); }
    SEXP weights = PROTECT(Rf_allocVector(REALSXP, n_points * n_studies)); nprot++;
    { SEXP weights_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(weights_dim)[0] = n_points; INTEGER(weights_dim)[1] = n_studies; Rf_setAttrib(weights, R_DimSymbol, weights_dim); UNPROTECT(1); }
    SEXP global_js_divergence = PROTECT(Rf_allocVector(REALSXP, n_studies)); nprot++;
    SEXP p_values = PROTECT(Rf_allocVector(REALSXP, n_studies)); nprot++;
    SEXP p_values_observed_consensus = PROTECT(Rf_allocVector(REALSXP, n_studies)); nprot++;
    int ierr = 0;

    run_js_comp_test_c(
        &n_studies,
        &max_n_genes_all_studies,
        &max_n_reps_all_studies,
        &n_points,
        &n_neighbors_v,
        REAL(gene_means),
        REAL(residuals),
        REAL(x_star),
        INTEGER(neighborhood_indices),
        INTEGER(neighborhood_range),
        INTEGER(n_bins_per_point),
        REAL(shared_residual_range_low),
        REAL(shared_residual_range_high),
        &max_n_bins_per_point,
        occupancy_failed_c,
        INTEGER(n_pooled_residuals),
        INTEGER(min_bin_occupancy),
        REAL(mean_bin_occupancy),
        INTEGER(max_bin_occupancy),
        INTEGER(sturges_bins),
        INTEGER(fd_bins),
        REAL(pmfs),
        INTEGER(counts),
        INTEGER(included_n_reps),
        REAL(mean_pmf),
        INTEGER(mean_pmf_counts),
        INTEGER(mean_pmf_included_n_reps),
        REAL(js_divergences),
        REAL(weights),
        REAL(global_js_divergence),
        REAL(p_values),
        REAL(p_values_observed_consensus),
        &n_permutations_v,
        &random_seed_v,
        &min_residuals_per_bin_v,
        &m_min_v,
        &m_max_v,
        &gamma_occupancy_v,
        &lower_residual_range_quantile_v,
        &upper_residual_range_quantile_v,
        bin_count_mode_c,
        &ierr
    );

    // convert the outputs back
    SEXP occupancy_failed = PROTECT(tox_bool_out(occupancy_failed_c, n_points)); nprot++;

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 25)); nprot++;
    SET_VECTOR_ELT(_out, 0, neighborhood_indices);
    SET_VECTOR_ELT(_out, 1, neighborhood_range);
    SET_VECTOR_ELT(_out, 2, n_bins_per_point);
    SET_VECTOR_ELT(_out, 3, shared_residual_range_low);
    SET_VECTOR_ELT(_out, 4, shared_residual_range_high);
    SET_VECTOR_ELT(_out, 5, Rf_ScalarInteger(max_n_bins_per_point));
    SET_VECTOR_ELT(_out, 6, occupancy_failed);
    SET_VECTOR_ELT(_out, 7, n_pooled_residuals);
    SET_VECTOR_ELT(_out, 8, min_bin_occupancy);
    SET_VECTOR_ELT(_out, 9, mean_bin_occupancy);
    SET_VECTOR_ELT(_out, 10, max_bin_occupancy);
    SET_VECTOR_ELT(_out, 11, sturges_bins);
    SET_VECTOR_ELT(_out, 12, fd_bins);
    SET_VECTOR_ELT(_out, 13, pmfs);
    SET_VECTOR_ELT(_out, 14, counts);
    SET_VECTOR_ELT(_out, 15, included_n_reps);
    SET_VECTOR_ELT(_out, 16, mean_pmf);
    SET_VECTOR_ELT(_out, 17, mean_pmf_counts);
    SET_VECTOR_ELT(_out, 18, mean_pmf_included_n_reps);
    SET_VECTOR_ELT(_out, 19, js_divergences);
    SET_VECTOR_ELT(_out, 20, weights);
    SET_VECTOR_ELT(_out, 21, global_js_divergence);
    SET_VECTOR_ELT(_out, 22, p_values);
    SET_VECTOR_ELT(_out, 23, p_values_observed_consensus);
    SET_VECTOR_ELT(_out, 24, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 25)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("neighborhood_indices"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("neighborhood_range"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("n_bins_per_point"));
    SET_STRING_ELT(_nms, 3, Rf_mkChar("shared_residual_range_low"));
    SET_STRING_ELT(_nms, 4, Rf_mkChar("shared_residual_range_high"));
    SET_STRING_ELT(_nms, 5, Rf_mkChar("max_n_bins_per_point"));
    SET_STRING_ELT(_nms, 6, Rf_mkChar("occupancy_failed"));
    SET_STRING_ELT(_nms, 7, Rf_mkChar("n_pooled_residuals"));
    SET_STRING_ELT(_nms, 8, Rf_mkChar("min_bin_occupancy"));
    SET_STRING_ELT(_nms, 9, Rf_mkChar("mean_bin_occupancy"));
    SET_STRING_ELT(_nms, 10, Rf_mkChar("max_bin_occupancy"));
    SET_STRING_ELT(_nms, 11, Rf_mkChar("sturges_bins"));
    SET_STRING_ELT(_nms, 12, Rf_mkChar("fd_bins"));
    SET_STRING_ELT(_nms, 13, Rf_mkChar("pmfs"));
    SET_STRING_ELT(_nms, 14, Rf_mkChar("counts"));
    SET_STRING_ELT(_nms, 15, Rf_mkChar("included_n_reps"));
    SET_STRING_ELT(_nms, 16, Rf_mkChar("mean_pmf"));
    SET_STRING_ELT(_nms, 17, Rf_mkChar("mean_pmf_counts"));
    SET_STRING_ELT(_nms, 18, Rf_mkChar("mean_pmf_included_n_reps"));
    SET_STRING_ELT(_nms, 19, Rf_mkChar("js_divergences"));
    SET_STRING_ELT(_nms, 20, Rf_mkChar("weights"));
    SET_STRING_ELT(_nms, 21, Rf_mkChar("global_js_divergence"));
    SET_STRING_ELT(_nms, 22, Rf_mkChar("p_values"));
    SET_STRING_ELT(_nms, 23, Rf_mkChar("p_values_observed_consensus"));
    SET_STRING_ELT(_nms, 24, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP run_js_comp_test_expert_call(SEXP n_neighbors, SEXP gene_means, SEXP gene_means_perms, SEXP residuals, SEXP x_star, SEXP n_permutations, SEXP random_seed, SEXP min_residuals_per_bin, SEXP m_min, SEXP m_max, SEXP gamma_occupancy, SEXP lower_residual_range_quantile, SEXP upper_residual_range_quantile, SEXP bin_count_mode) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_studies = INTEGER(Rf_getAttrib(gene_means, R_DimSymbol))[1];
    int max_n_genes_all_studies = INTEGER(Rf_getAttrib(gene_means, R_DimSymbol))[0];
    int max_n_reps_all_studies = INTEGER(Rf_getAttrib(residuals, R_DimSymbol))[0];
    int n_points = (int) Rf_length(x_star);

    // scalar inputs, pulled from their length-1 vectors
    int n_neighbors_v = Rf_asInteger(n_neighbors);
    int n_permutations_v = Rf_asInteger(n_permutations);
    int random_seed_v = Rf_asInteger(random_seed);
    int min_residuals_per_bin_v = Rf_asInteger(min_residuals_per_bin);
    int m_min_v = Rf_asInteger(m_min);
    int m_max_v = Rf_asInteger(m_max);
    double gamma_occupancy_v = Rf_asReal(gamma_occupancy);
    double lower_residual_range_quantile_v = Rf_asReal(lower_residual_range_quantile);
    double upper_residual_range_quantile_v = Rf_asReal(upper_residual_range_quantile);

    // convert what Fortran cannot take from R directly
    char* bin_count_mode_c = tox_char_in(bin_count_mode, 29);

    // outputs and work space
    SEXP neighborhood_indices = PROTECT(Rf_allocVector(INTSXP, n_neighbors_v * n_points * n_studies)); nprot++;
    { SEXP neighborhood_indices_dim = PROTECT(Rf_allocVector(INTSXP, 3)); INTEGER(neighborhood_indices_dim)[0] = n_neighbors_v; INTEGER(neighborhood_indices_dim)[1] = n_points; INTEGER(neighborhood_indices_dim)[2] = n_studies; Rf_setAttrib(neighborhood_indices, R_DimSymbol, neighborhood_indices_dim); UNPROTECT(1); }
    SEXP neighborhood_range = PROTECT(Rf_allocVector(INTSXP, 2 * n_points * n_studies)); nprot++;
    { SEXP neighborhood_range_dim = PROTECT(Rf_allocVector(INTSXP, 3)); INTEGER(neighborhood_range_dim)[0] = 2; INTEGER(neighborhood_range_dim)[1] = n_points; INTEGER(neighborhood_range_dim)[2] = n_studies; Rf_setAttrib(neighborhood_range, R_DimSymbol, neighborhood_range_dim); UNPROTECT(1); }
    SEXP n_bins_per_point = PROTECT(Rf_allocVector(INTSXP, n_points)); nprot++;
    SEXP shared_residual_range_low = PROTECT(Rf_allocVector(REALSXP, n_points)); nprot++;
    SEXP shared_residual_range_high = PROTECT(Rf_allocVector(REALSXP, n_points)); nprot++;
    int max_n_bins_per_point = 0;
    unsigned char* occupancy_failed_c = tox_bool_alloc(n_points);
    SEXP n_pooled_residuals = PROTECT(Rf_allocVector(INTSXP, n_points)); nprot++;
    SEXP min_bin_occupancy = PROTECT(Rf_allocVector(INTSXP, n_points)); nprot++;
    SEXP mean_bin_occupancy = PROTECT(Rf_allocVector(REALSXP, n_points)); nprot++;
    SEXP max_bin_occupancy = PROTECT(Rf_allocVector(INTSXP, n_points)); nprot++;
    SEXP sturges_bins = PROTECT(Rf_allocVector(INTSXP, n_points)); nprot++;
    SEXP fd_bins = PROTECT(Rf_allocVector(INTSXP, n_points)); nprot++;
    SEXP pmfs = PROTECT(Rf_allocVector(REALSXP, n_points * 256 * n_studies)); nprot++;
    { SEXP pmfs_dim = PROTECT(Rf_allocVector(INTSXP, 3)); INTEGER(pmfs_dim)[0] = n_points; INTEGER(pmfs_dim)[1] = 256; INTEGER(pmfs_dim)[2] = n_studies; Rf_setAttrib(pmfs, R_DimSymbol, pmfs_dim); UNPROTECT(1); }
    SEXP counts = PROTECT(Rf_allocVector(INTSXP, n_points * 256 * n_studies)); nprot++;
    { SEXP counts_dim = PROTECT(Rf_allocVector(INTSXP, 3)); INTEGER(counts_dim)[0] = n_points; INTEGER(counts_dim)[1] = 256; INTEGER(counts_dim)[2] = n_studies; Rf_setAttrib(counts, R_DimSymbol, counts_dim); UNPROTECT(1); }
    SEXP included_n_reps = PROTECT(Rf_allocVector(INTSXP, n_points * n_studies)); nprot++;
    { SEXP included_n_reps_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(included_n_reps_dim)[0] = n_points; INTEGER(included_n_reps_dim)[1] = n_studies; Rf_setAttrib(included_n_reps, R_DimSymbol, included_n_reps_dim); UNPROTECT(1); }
    SEXP mean_pmf = PROTECT(Rf_allocVector(REALSXP, n_points * 256)); nprot++;
    { SEXP mean_pmf_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(mean_pmf_dim)[0] = n_points; INTEGER(mean_pmf_dim)[1] = 256; Rf_setAttrib(mean_pmf, R_DimSymbol, mean_pmf_dim); UNPROTECT(1); }
    SEXP mean_pmf_counts = PROTECT(Rf_allocVector(INTSXP, n_points * 256)); nprot++;
    { SEXP mean_pmf_counts_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(mean_pmf_counts_dim)[0] = n_points; INTEGER(mean_pmf_counts_dim)[1] = 256; Rf_setAttrib(mean_pmf_counts, R_DimSymbol, mean_pmf_counts_dim); UNPROTECT(1); }
    SEXP mean_pmf_included_n_reps = PROTECT(Rf_allocVector(INTSXP, n_points)); nprot++;
    SEXP js_divergences = PROTECT(Rf_allocVector(REALSXP, n_points * n_studies)); nprot++;
    { SEXP js_divergences_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(js_divergences_dim)[0] = n_points; INTEGER(js_divergences_dim)[1] = n_studies; Rf_setAttrib(js_divergences, R_DimSymbol, js_divergences_dim); UNPROTECT(1); }
    SEXP weights = PROTECT(Rf_allocVector(REALSXP, n_points * n_studies)); nprot++;
    { SEXP weights_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(weights_dim)[0] = n_points; INTEGER(weights_dim)[1] = n_studies; Rf_setAttrib(weights, R_DimSymbol, weights_dim); UNPROTECT(1); }
    SEXP global_js_divergence = PROTECT(Rf_allocVector(REALSXP, n_studies)); nprot++;
    SEXP p_values = PROTECT(Rf_allocVector(REALSXP, n_studies)); nprot++;
    SEXP p_values_observed_consensus = PROTECT(Rf_allocVector(REALSXP, n_studies)); nprot++;
    double* tmp_neighborhood_residuals_gathered = (double*) R_alloc(max_n_reps_all_studies * n_neighbors_v * n_points, sizeof(double));
    double* tmp_permutation_mean_pmf = (double*) R_alloc(n_points * 256, sizeof(double));
    int* tmp_permutation_point_pool = (int*) R_alloc(256, sizeof(int));
    int* tmp_permutation_point_draw = (int*) R_alloc(256, sizeof(int));
    double* tmp_pooled_residuals = (double*) R_alloc((max_n_reps_all_studies*n_neighbors_v*n_studies), sizeof(double));
    int* tmp_pooled_residuals_perm = (int*) R_alloc((max_n_reps_all_studies*n_neighbors_v*n_studies), sizeof(int));
    int* tmp_bin_counts_search = (int*) R_alloc(256, sizeof(int));
    int* tmp_residual_to_study_perm = (int*) R_alloc((max_n_reps_all_studies*n_neighbors_v*n_studies), sizeof(int));
    int* tmp_study_block_ends = (int*) R_alloc(n_studies, sizeof(int));
    int* tmp_study_block_sizes = (int*) R_alloc(n_studies, sizeof(int));
    int* tmp_point_n_neighbors = (int*) R_alloc(n_studies, sizeof(int));
    int* tmp_permutation_mean_pmf_counts = (int*) R_alloc(n_points * 256, sizeof(int));
    int* tmp_permutation_counts = (int*) R_alloc(n_points * 256, sizeof(int));
    double* tmp_permutation_pmfs = (double*) R_alloc(n_points * 256 * n_studies, sizeof(double));
    double* tmp_permutation_js_divergences = (double*) R_alloc(n_points * n_studies, sizeof(double));
    double* tmp_permutation_weights = (double*) R_alloc(n_points * n_studies, sizeof(double));
    double* tmp_permutation_global_js_divergence = (double*) R_alloc(n_studies, sizeof(double));
    int ierr = 0;

    run_js_comp_test_expert_c(
        &n_studies,
        &max_n_genes_all_studies,
        &max_n_reps_all_studies,
        &n_points,
        &n_neighbors_v,
        REAL(gene_means),
        INTEGER(gene_means_perms),
        REAL(residuals),
        REAL(x_star),
        INTEGER(neighborhood_indices),
        INTEGER(neighborhood_range),
        INTEGER(n_bins_per_point),
        REAL(shared_residual_range_low),
        REAL(shared_residual_range_high),
        &max_n_bins_per_point,
        occupancy_failed_c,
        INTEGER(n_pooled_residuals),
        INTEGER(min_bin_occupancy),
        REAL(mean_bin_occupancy),
        INTEGER(max_bin_occupancy),
        INTEGER(sturges_bins),
        INTEGER(fd_bins),
        REAL(pmfs),
        INTEGER(counts),
        INTEGER(included_n_reps),
        REAL(mean_pmf),
        INTEGER(mean_pmf_counts),
        INTEGER(mean_pmf_included_n_reps),
        REAL(js_divergences),
        REAL(weights),
        REAL(global_js_divergence),
        REAL(p_values),
        REAL(p_values_observed_consensus),
        tmp_neighborhood_residuals_gathered,
        tmp_permutation_mean_pmf,
        tmp_permutation_point_pool,
        tmp_permutation_point_draw,
        tmp_pooled_residuals,
        tmp_pooled_residuals_perm,
        tmp_bin_counts_search,
        tmp_residual_to_study_perm,
        tmp_study_block_ends,
        tmp_study_block_sizes,
        tmp_point_n_neighbors,
        tmp_permutation_mean_pmf_counts,
        tmp_permutation_counts,
        tmp_permutation_pmfs,
        tmp_permutation_js_divergences,
        tmp_permutation_weights,
        tmp_permutation_global_js_divergence,
        &n_permutations_v,
        &random_seed_v,
        &min_residuals_per_bin_v,
        &m_min_v,
        &m_max_v,
        &gamma_occupancy_v,
        &lower_residual_range_quantile_v,
        &upper_residual_range_quantile_v,
        bin_count_mode_c,
        &ierr
    );

    // convert the outputs back
    SEXP occupancy_failed = PROTECT(tox_bool_out(occupancy_failed_c, n_points)); nprot++;

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 25)); nprot++;
    SET_VECTOR_ELT(_out, 0, neighborhood_indices);
    SET_VECTOR_ELT(_out, 1, neighborhood_range);
    SET_VECTOR_ELT(_out, 2, n_bins_per_point);
    SET_VECTOR_ELT(_out, 3, shared_residual_range_low);
    SET_VECTOR_ELT(_out, 4, shared_residual_range_high);
    SET_VECTOR_ELT(_out, 5, Rf_ScalarInteger(max_n_bins_per_point));
    SET_VECTOR_ELT(_out, 6, occupancy_failed);
    SET_VECTOR_ELT(_out, 7, n_pooled_residuals);
    SET_VECTOR_ELT(_out, 8, min_bin_occupancy);
    SET_VECTOR_ELT(_out, 9, mean_bin_occupancy);
    SET_VECTOR_ELT(_out, 10, max_bin_occupancy);
    SET_VECTOR_ELT(_out, 11, sturges_bins);
    SET_VECTOR_ELT(_out, 12, fd_bins);
    SET_VECTOR_ELT(_out, 13, pmfs);
    SET_VECTOR_ELT(_out, 14, counts);
    SET_VECTOR_ELT(_out, 15, included_n_reps);
    SET_VECTOR_ELT(_out, 16, mean_pmf);
    SET_VECTOR_ELT(_out, 17, mean_pmf_counts);
    SET_VECTOR_ELT(_out, 18, mean_pmf_included_n_reps);
    SET_VECTOR_ELT(_out, 19, js_divergences);
    SET_VECTOR_ELT(_out, 20, weights);
    SET_VECTOR_ELT(_out, 21, global_js_divergence);
    SET_VECTOR_ELT(_out, 22, p_values);
    SET_VECTOR_ELT(_out, 23, p_values_observed_consensus);
    SET_VECTOR_ELT(_out, 24, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 25)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("neighborhood_indices"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("neighborhood_range"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("n_bins_per_point"));
    SET_STRING_ELT(_nms, 3, Rf_mkChar("shared_residual_range_low"));
    SET_STRING_ELT(_nms, 4, Rf_mkChar("shared_residual_range_high"));
    SET_STRING_ELT(_nms, 5, Rf_mkChar("max_n_bins_per_point"));
    SET_STRING_ELT(_nms, 6, Rf_mkChar("occupancy_failed"));
    SET_STRING_ELT(_nms, 7, Rf_mkChar("n_pooled_residuals"));
    SET_STRING_ELT(_nms, 8, Rf_mkChar("min_bin_occupancy"));
    SET_STRING_ELT(_nms, 9, Rf_mkChar("mean_bin_occupancy"));
    SET_STRING_ELT(_nms, 10, Rf_mkChar("max_bin_occupancy"));
    SET_STRING_ELT(_nms, 11, Rf_mkChar("sturges_bins"));
    SET_STRING_ELT(_nms, 12, Rf_mkChar("fd_bins"));
    SET_STRING_ELT(_nms, 13, Rf_mkChar("pmfs"));
    SET_STRING_ELT(_nms, 14, Rf_mkChar("counts"));
    SET_STRING_ELT(_nms, 15, Rf_mkChar("included_n_reps"));
    SET_STRING_ELT(_nms, 16, Rf_mkChar("mean_pmf"));
    SET_STRING_ELT(_nms, 17, Rf_mkChar("mean_pmf_counts"));
    SET_STRING_ELT(_nms, 18, Rf_mkChar("mean_pmf_included_n_reps"));
    SET_STRING_ELT(_nms, 19, Rf_mkChar("js_divergences"));
    SET_STRING_ELT(_nms, 20, Rf_mkChar("weights"));
    SET_STRING_ELT(_nms, 21, Rf_mkChar("global_js_divergence"));
    SET_STRING_ELT(_nms, 22, Rf_mkChar("p_values"));
    SET_STRING_ELT(_nms, 23, Rf_mkChar("p_values_observed_consensus"));
    SET_STRING_ELT(_nms, 24, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP run_js_comp_test_adaptive_call(SEXP gene_means, SEXP residuals, SEXP pooled_neighborhood_range, SEXP n_permutations, SEXP random_seed, SEXP min_residuals_per_bin, SEXP m_min, SEXP m_max, SEXP gamma_occupancy, SEXP lower_residual_range_quantile, SEXP upper_residual_range_quantile, SEXP bin_count_mode) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_studies = INTEGER(Rf_getAttrib(gene_means, R_DimSymbol))[1];
    int max_n_genes_all_studies = INTEGER(Rf_getAttrib(gene_means, R_DimSymbol))[0];
    int max_n_reps_all_studies = INTEGER(Rf_getAttrib(residuals, R_DimSymbol))[0];
    int n_points = INTEGER(Rf_getAttrib(pooled_neighborhood_range, R_DimSymbol))[1];

    // scalar inputs, pulled from their length-1 vectors
    int n_permutations_v = Rf_asInteger(n_permutations);
    int random_seed_v = Rf_asInteger(random_seed);
    int min_residuals_per_bin_v = Rf_asInteger(min_residuals_per_bin);
    int m_min_v = Rf_asInteger(m_min);
    int m_max_v = Rf_asInteger(m_max);
    double gamma_occupancy_v = Rf_asReal(gamma_occupancy);
    double lower_residual_range_quantile_v = Rf_asReal(lower_residual_range_quantile);
    double upper_residual_range_quantile_v = Rf_asReal(upper_residual_range_quantile);

    // convert what Fortran cannot take from R directly
    char* bin_count_mode_c = tox_char_in(bin_count_mode, 29);

    // outputs and work space
    SEXP n_neighbors_per_point = PROTECT(Rf_allocVector(INTSXP, n_studies * n_points)); nprot++;
    { SEXP n_neighbors_per_point_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(n_neighbors_per_point_dim)[0] = n_studies; INTEGER(n_neighbors_per_point_dim)[1] = n_points; Rf_setAttrib(n_neighbors_per_point, R_DimSymbol, n_neighbors_per_point_dim); UNPROTECT(1); }
    SEXP n_bins_per_point = PROTECT(Rf_allocVector(INTSXP, n_points)); nprot++;
    SEXP shared_residual_range_low = PROTECT(Rf_allocVector(REALSXP, n_points)); nprot++;
    SEXP shared_residual_range_high = PROTECT(Rf_allocVector(REALSXP, n_points)); nprot++;
    int max_n_bins_per_point = 0;
    unsigned char* occupancy_failed_c = tox_bool_alloc(n_points);
    SEXP n_pooled_residuals = PROTECT(Rf_allocVector(INTSXP, n_points)); nprot++;
    SEXP min_bin_occupancy = PROTECT(Rf_allocVector(INTSXP, n_points)); nprot++;
    SEXP mean_bin_occupancy = PROTECT(Rf_allocVector(REALSXP, n_points)); nprot++;
    SEXP max_bin_occupancy = PROTECT(Rf_allocVector(INTSXP, n_points)); nprot++;
    SEXP sturges_bins = PROTECT(Rf_allocVector(INTSXP, n_points)); nprot++;
    SEXP fd_bins = PROTECT(Rf_allocVector(INTSXP, n_points)); nprot++;
    SEXP pmfs = PROTECT(Rf_allocVector(REALSXP, n_points * 256 * n_studies)); nprot++;
    { SEXP pmfs_dim = PROTECT(Rf_allocVector(INTSXP, 3)); INTEGER(pmfs_dim)[0] = n_points; INTEGER(pmfs_dim)[1] = 256; INTEGER(pmfs_dim)[2] = n_studies; Rf_setAttrib(pmfs, R_DimSymbol, pmfs_dim); UNPROTECT(1); }
    SEXP counts = PROTECT(Rf_allocVector(INTSXP, n_points * 256 * n_studies)); nprot++;
    { SEXP counts_dim = PROTECT(Rf_allocVector(INTSXP, 3)); INTEGER(counts_dim)[0] = n_points; INTEGER(counts_dim)[1] = 256; INTEGER(counts_dim)[2] = n_studies; Rf_setAttrib(counts, R_DimSymbol, counts_dim); UNPROTECT(1); }
    SEXP included_n_reps = PROTECT(Rf_allocVector(INTSXP, n_points * n_studies)); nprot++;
    { SEXP included_n_reps_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(included_n_reps_dim)[0] = n_points; INTEGER(included_n_reps_dim)[1] = n_studies; Rf_setAttrib(included_n_reps, R_DimSymbol, included_n_reps_dim); UNPROTECT(1); }
    SEXP mean_pmf = PROTECT(Rf_allocVector(REALSXP, n_points * 256)); nprot++;
    { SEXP mean_pmf_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(mean_pmf_dim)[0] = n_points; INTEGER(mean_pmf_dim)[1] = 256; Rf_setAttrib(mean_pmf, R_DimSymbol, mean_pmf_dim); UNPROTECT(1); }
    SEXP mean_pmf_counts = PROTECT(Rf_allocVector(INTSXP, n_points * 256)); nprot++;
    { SEXP mean_pmf_counts_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(mean_pmf_counts_dim)[0] = n_points; INTEGER(mean_pmf_counts_dim)[1] = 256; Rf_setAttrib(mean_pmf_counts, R_DimSymbol, mean_pmf_counts_dim); UNPROTECT(1); }
    SEXP mean_pmf_included_n_reps = PROTECT(Rf_allocVector(INTSXP, n_points)); nprot++;
    SEXP js_divergences = PROTECT(Rf_allocVector(REALSXP, n_points * n_studies)); nprot++;
    { SEXP js_divergences_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(js_divergences_dim)[0] = n_points; INTEGER(js_divergences_dim)[1] = n_studies; Rf_setAttrib(js_divergences, R_DimSymbol, js_divergences_dim); UNPROTECT(1); }
    SEXP weights = PROTECT(Rf_allocVector(REALSXP, n_points * n_studies)); nprot++;
    { SEXP weights_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(weights_dim)[0] = n_points; INTEGER(weights_dim)[1] = n_studies; Rf_setAttrib(weights, R_DimSymbol, weights_dim); UNPROTECT(1); }
    SEXP global_js_divergence = PROTECT(Rf_allocVector(REALSXP, n_studies)); nprot++;
    SEXP p_values = PROTECT(Rf_allocVector(REALSXP, n_studies)); nprot++;
    SEXP p_values_observed_consensus = PROTECT(Rf_allocVector(REALSXP, n_studies)); nprot++;
    int ierr = 0;

    run_js_comp_test_adaptive_c(
        &n_studies,
        &max_n_genes_all_studies,
        &max_n_reps_all_studies,
        &n_points,
        REAL(gene_means),
        REAL(residuals),
        INTEGER(pooled_neighborhood_range),
        INTEGER(n_neighbors_per_point),
        INTEGER(n_bins_per_point),
        REAL(shared_residual_range_low),
        REAL(shared_residual_range_high),
        &max_n_bins_per_point,
        occupancy_failed_c,
        INTEGER(n_pooled_residuals),
        INTEGER(min_bin_occupancy),
        REAL(mean_bin_occupancy),
        INTEGER(max_bin_occupancy),
        INTEGER(sturges_bins),
        INTEGER(fd_bins),
        REAL(pmfs),
        INTEGER(counts),
        INTEGER(included_n_reps),
        REAL(mean_pmf),
        INTEGER(mean_pmf_counts),
        INTEGER(mean_pmf_included_n_reps),
        REAL(js_divergences),
        REAL(weights),
        REAL(global_js_divergence),
        REAL(p_values),
        REAL(p_values_observed_consensus),
        &n_permutations_v,
        &random_seed_v,
        &min_residuals_per_bin_v,
        &m_min_v,
        &m_max_v,
        &gamma_occupancy_v,
        &lower_residual_range_quantile_v,
        &upper_residual_range_quantile_v,
        bin_count_mode_c,
        &ierr
    );

    // convert the outputs back
    SEXP occupancy_failed = PROTECT(tox_bool_out(occupancy_failed_c, n_points)); nprot++;

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 24)); nprot++;
    SET_VECTOR_ELT(_out, 0, n_neighbors_per_point);
    SET_VECTOR_ELT(_out, 1, n_bins_per_point);
    SET_VECTOR_ELT(_out, 2, shared_residual_range_low);
    SET_VECTOR_ELT(_out, 3, shared_residual_range_high);
    SET_VECTOR_ELT(_out, 4, Rf_ScalarInteger(max_n_bins_per_point));
    SET_VECTOR_ELT(_out, 5, occupancy_failed);
    SET_VECTOR_ELT(_out, 6, n_pooled_residuals);
    SET_VECTOR_ELT(_out, 7, min_bin_occupancy);
    SET_VECTOR_ELT(_out, 8, mean_bin_occupancy);
    SET_VECTOR_ELT(_out, 9, max_bin_occupancy);
    SET_VECTOR_ELT(_out, 10, sturges_bins);
    SET_VECTOR_ELT(_out, 11, fd_bins);
    SET_VECTOR_ELT(_out, 12, pmfs);
    SET_VECTOR_ELT(_out, 13, counts);
    SET_VECTOR_ELT(_out, 14, included_n_reps);
    SET_VECTOR_ELT(_out, 15, mean_pmf);
    SET_VECTOR_ELT(_out, 16, mean_pmf_counts);
    SET_VECTOR_ELT(_out, 17, mean_pmf_included_n_reps);
    SET_VECTOR_ELT(_out, 18, js_divergences);
    SET_VECTOR_ELT(_out, 19, weights);
    SET_VECTOR_ELT(_out, 20, global_js_divergence);
    SET_VECTOR_ELT(_out, 21, p_values);
    SET_VECTOR_ELT(_out, 22, p_values_observed_consensus);
    SET_VECTOR_ELT(_out, 23, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 24)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("n_neighbors_per_point"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("n_bins_per_point"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("shared_residual_range_low"));
    SET_STRING_ELT(_nms, 3, Rf_mkChar("shared_residual_range_high"));
    SET_STRING_ELT(_nms, 4, Rf_mkChar("max_n_bins_per_point"));
    SET_STRING_ELT(_nms, 5, Rf_mkChar("occupancy_failed"));
    SET_STRING_ELT(_nms, 6, Rf_mkChar("n_pooled_residuals"));
    SET_STRING_ELT(_nms, 7, Rf_mkChar("min_bin_occupancy"));
    SET_STRING_ELT(_nms, 8, Rf_mkChar("mean_bin_occupancy"));
    SET_STRING_ELT(_nms, 9, Rf_mkChar("max_bin_occupancy"));
    SET_STRING_ELT(_nms, 10, Rf_mkChar("sturges_bins"));
    SET_STRING_ELT(_nms, 11, Rf_mkChar("fd_bins"));
    SET_STRING_ELT(_nms, 12, Rf_mkChar("pmfs"));
    SET_STRING_ELT(_nms, 13, Rf_mkChar("counts"));
    SET_STRING_ELT(_nms, 14, Rf_mkChar("included_n_reps"));
    SET_STRING_ELT(_nms, 15, Rf_mkChar("mean_pmf"));
    SET_STRING_ELT(_nms, 16, Rf_mkChar("mean_pmf_counts"));
    SET_STRING_ELT(_nms, 17, Rf_mkChar("mean_pmf_included_n_reps"));
    SET_STRING_ELT(_nms, 18, Rf_mkChar("js_divergences"));
    SET_STRING_ELT(_nms, 19, Rf_mkChar("weights"));
    SET_STRING_ELT(_nms, 20, Rf_mkChar("global_js_divergence"));
    SET_STRING_ELT(_nms, 21, Rf_mkChar("p_values"));
    SET_STRING_ELT(_nms, 22, Rf_mkChar("p_values_observed_consensus"));
    SET_STRING_ELT(_nms, 23, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

#endif  // R binding
