// Generated. Do not edit.
#if !defined(NO_R_BINDING) && !defined(NO_C_BINDING)
#include <R.h>
#include <Rinternals.h>
#include "tox_marshal.h"

// the Fortran C-ABI symbols this module calls
void serialize_stc_results_as_json_c(const char*, const int*, const int*, const int*, const int*, const int*, const int*, const int*, const double*, const char*, const int*, const unsigned char*, const unsigned char*, const int*, const double*, const double*, const double*, const int*, const double*, const double*, const int*, const unsigned char*, const int*, const unsigned char*, const double*, const int*, const int*, const int*, const int*, const double*, const int*, const double*, const double*, const double*, const int*, const double*, const double*, const char*, const double*, const unsigned char*, const int*, const int*, const double*, const unsigned char*, const unsigned char*, const unsigned char*, const unsigned char*, const int*, const int*, const double*, const double*, const double*, const int*, int*);
void write_stc_interactive_html_report_c(const char*, const int*, const int*, const int*, const int*, const int*, const int*, const int*, const double*, const char*, const int*, const unsigned char*, const unsigned char*, const int*, const double*, const double*, const double*, const int*, const double*, const double*, const int*, const unsigned char*, const int*, const unsigned char*, const double*, const int*, const int*, const int*, const int*, const double*, const int*, const double*, const double*, const double*, const int*, const double*, const double*, const char*, const double*, const unsigned char*, const int*, const int*, const double*, const unsigned char*, const unsigned char*, const unsigned char*, const unsigned char*, const int*, const int*, const double*, const double*, const double*, const int*, int*);

SEXP serialize_stc_results_as_json_call(SEXP filename, SEXP n_super_ensembles, SEXP vectors, SEXP dim_names, SEXP seed_selection_mask, SEXP ensemble_masks, SEXP ensemble_stop_reason, SEXP ensemble_growth_radii, SEXP ensemble_U_history, SEXP ensemble_S_history, SEXP ensemble_d_history, SEXP ensemble_G_history, SEXP ensemble_mu_history, SEXP ensemble_k_history, SEXP ensemble_accepted_history, SEXP ensemble_member_added_at_step, SEXP ensemble_low_confidence_masks, SEXP ensemble_U_first, SEXP ensemble_d_first, SEXP super_ensembles, SEXP k_min, SEXP k_density, SEXP chordal_dist_max_as_prcnt_of_range, SEXP d_max, SEXP G_max, SEXP RMSE_change_max, SEXP f_max, SEXP a, SEXP exclusion_radius_percentile, SEXP bandwidth_percentile, SEXP reconciliation_mode, SEXP min_overlap_coefficient, SEXP allowed_stop_reasons, SEXP filter_d_min, SEXP filter_d_max, SEXP filter_var_explained_min, SEXP ensemble_eligible, SEXP ensemble_eligible_by_stop_condition, SEXP ensemble_eligible_by_dimension, SEXP ensemble_eligible_by_var_explained, SEXP estimated_k_min, SEXP estimated_k_density, SEXP estimated_density_quantile, SEXP estimated_chordal_dist_max_as_prcnt_of_range, SEXP estimated_G_max, SEXP estimated_d_max) {
    int nprot = 0;
    // optionals: a null pointer and size 0 when the caller omits them
    int allowed_stop_reasons_size = 0;
    if (allowed_stop_reasons != R_NilValue) {
        allowed_stop_reasons_size = (int) Rf_length(allowed_stop_reasons);
    }
    const int* filter_d_min_p = NULL;
    int filter_d_min_size = 0;
    if (filter_d_min != R_NilValue) {
        filter_d_min_size = (int) Rf_length(filter_d_min);
        filter_d_min_p = INTEGER(filter_d_min);
    }
    const int* filter_d_max_p = NULL;
    int filter_d_max_size = 0;
    if (filter_d_max != R_NilValue) {
        filter_d_max_size = (int) Rf_length(filter_d_max);
        filter_d_max_p = INTEGER(filter_d_max);
    }
    const double* filter_var_explained_min_p = NULL;
    int filter_var_explained_min_size = 0;
    if (filter_var_explained_min != R_NilValue) {
        filter_var_explained_min_size = (int) Rf_length(filter_var_explained_min);
        filter_var_explained_min_p = REAL(filter_var_explained_min);
    }
    const int* estimated_k_min_p = NULL;
    int estimated_k_min_size = 0;
    if (estimated_k_min != R_NilValue) {
        estimated_k_min_size = (int) Rf_length(estimated_k_min);
        estimated_k_min_p = INTEGER(estimated_k_min);
    }
    const int* estimated_k_density_p = NULL;
    int estimated_k_density_size = 0;
    if (estimated_k_density != R_NilValue) {
        estimated_k_density_size = (int) Rf_length(estimated_k_density);
        estimated_k_density_p = INTEGER(estimated_k_density);
    }
    const double* estimated_density_quantile_p = NULL;
    int estimated_density_quantile_size = 0;
    if (estimated_density_quantile != R_NilValue) {
        estimated_density_quantile_size = (int) Rf_length(estimated_density_quantile);
        estimated_density_quantile_p = REAL(estimated_density_quantile);
    }
    const double* estimated_chordal_dist_max_as_prcnt_of_range_p = NULL;
    int estimated_chordal_dist_max_as_prcnt_of_range_size = 0;
    if (estimated_chordal_dist_max_as_prcnt_of_range != R_NilValue) {
        estimated_chordal_dist_max_as_prcnt_of_range_size = (int) Rf_length(estimated_chordal_dist_max_as_prcnt_of_range);
        estimated_chordal_dist_max_as_prcnt_of_range_p = REAL(estimated_chordal_dist_max_as_prcnt_of_range);
    }
    const double* estimated_G_max_p = NULL;
    int estimated_G_max_size = 0;
    if (estimated_G_max != R_NilValue) {
        estimated_G_max_size = (int) Rf_length(estimated_G_max);
        estimated_G_max_p = REAL(estimated_G_max);
    }
    const int* estimated_d_max_p = NULL;
    int estimated_d_max_size = 0;
    if (estimated_d_max != R_NilValue) {
        estimated_d_max_size = (int) Rf_length(estimated_d_max);
        estimated_d_max_p = INTEGER(estimated_d_max);
    }

    // derived from the inputs, not asked of the caller
    int filename_strlen = tox_max_strlen(filename);
    int n_dimensions = INTEGER(Rf_getAttrib(vectors, R_DimSymbol))[0];
    int n_vectors = INTEGER(Rf_getAttrib(vectors, R_DimSymbol))[1];
    int n_selected_seed = tox_sum_true(seed_selection_mask);
    int o = INTEGER(Rf_getAttrib(ensemble_U_history, R_DimSymbol))[2];
    int max_group_size = INTEGER(Rf_getAttrib(super_ensembles, R_DimSymbol))[0];
    int dim_names_strlen = tox_max_strlen(dim_names);

    // scalar inputs, pulled from their length-1 vectors
    int n_super_ensembles_v = Rf_asInteger(n_super_ensembles);
    int k_min_v = Rf_asInteger(k_min);
    int k_density_v = Rf_asInteger(k_density);
    double chordal_dist_max_as_prcnt_of_range_v = Rf_asReal(chordal_dist_max_as_prcnt_of_range);
    int d_max_v = Rf_asInteger(d_max);
    double G_max_v = Rf_asReal(G_max);
    double RMSE_change_max_v = Rf_asReal(RMSE_change_max);
    double f_max_v = Rf_asReal(f_max);
    int a_v = Rf_asInteger(a);
    double exclusion_radius_percentile_v = Rf_asReal(exclusion_radius_percentile);
    double bandwidth_percentile_v = Rf_asReal(bandwidth_percentile);
    double min_overlap_coefficient_v = Rf_asReal(min_overlap_coefficient);

    // convert what Fortran cannot take from R directly
    char* filename_c = tox_char_in(filename, filename_strlen);
    char* dim_names_c = tox_char_in(dim_names, dim_names_strlen);
    unsigned char* seed_selection_mask_c = tox_bool_in(seed_selection_mask);
    unsigned char* ensemble_masks_c = tox_bool_in(ensemble_masks);
    unsigned char* ensemble_accepted_history_c = tox_bool_in(ensemble_accepted_history);
    unsigned char* ensemble_low_confidence_masks_c = tox_bool_in(ensemble_low_confidence_masks);
    char* reconciliation_mode_c = tox_char_in(reconciliation_mode, 25);
    unsigned char* allowed_stop_reasons_c = tox_bool_in(allowed_stop_reasons);
    unsigned char* ensemble_eligible_c = tox_bool_in(ensemble_eligible);
    unsigned char* ensemble_eligible_by_stop_condition_c = tox_bool_in(ensemble_eligible_by_stop_condition);
    unsigned char* ensemble_eligible_by_dimension_c = tox_bool_in(ensemble_eligible_by_dimension);
    unsigned char* ensemble_eligible_by_var_explained_c = tox_bool_in(ensemble_eligible_by_var_explained);

    // outputs and work space
    int ierr = 0;

    serialize_stc_results_as_json_c(
        filename_c,
        &filename_strlen,
        &n_dimensions,
        &n_vectors,
        &n_selected_seed,
        &o,
        &max_group_size,
        &n_super_ensembles_v,
        REAL(vectors),
        dim_names_c,
        &dim_names_strlen,
        seed_selection_mask_c,
        ensemble_masks_c,
        INTEGER(ensemble_stop_reason),
        REAL(ensemble_growth_radii),
        REAL(ensemble_U_history),
        REAL(ensemble_S_history),
        INTEGER(ensemble_d_history),
        REAL(ensemble_G_history),
        REAL(ensemble_mu_history),
        INTEGER(ensemble_k_history),
        ensemble_accepted_history_c,
        INTEGER(ensemble_member_added_at_step),
        ensemble_low_confidence_masks_c,
        REAL(ensemble_U_first),
        INTEGER(ensemble_d_first),
        INTEGER(super_ensembles),
        &k_min_v,
        &k_density_v,
        &chordal_dist_max_as_prcnt_of_range_v,
        &d_max_v,
        &G_max_v,
        &RMSE_change_max_v,
        &f_max_v,
        &a_v,
        &exclusion_radius_percentile_v,
        &bandwidth_percentile_v,
        reconciliation_mode_c,
        &min_overlap_coefficient_v,
        allowed_stop_reasons != R_NilValue ? allowed_stop_reasons_c : NULL,
        filter_d_min_p,
        filter_d_max_p,
        filter_var_explained_min_p,
        ensemble_eligible_c,
        ensemble_eligible_by_stop_condition_c,
        ensemble_eligible_by_dimension_c,
        ensemble_eligible_by_var_explained_c,
        estimated_k_min_p,
        estimated_k_density_p,
        estimated_density_quantile_p,
        estimated_chordal_dist_max_as_prcnt_of_range_p,
        estimated_G_max_p,
        estimated_d_max_p,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 1)); nprot++;
    SET_VECTOR_ELT(_out, 0, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 1)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP write_stc_interactive_html_report_call(SEXP filename, SEXP n_super_ensembles, SEXP vectors, SEXP dim_names, SEXP seed_selection_mask, SEXP ensemble_masks, SEXP ensemble_stop_reason, SEXP ensemble_growth_radii, SEXP ensemble_U_history, SEXP ensemble_S_history, SEXP ensemble_d_history, SEXP ensemble_G_history, SEXP ensemble_mu_history, SEXP ensemble_k_history, SEXP ensemble_accepted_history, SEXP ensemble_member_added_at_step, SEXP ensemble_low_confidence_masks, SEXP ensemble_U_first, SEXP ensemble_d_first, SEXP super_ensembles, SEXP k_min, SEXP k_density, SEXP chordal_dist_max_as_prcnt_of_range, SEXP d_max, SEXP G_max, SEXP RMSE_change_max, SEXP f_max, SEXP a, SEXP exclusion_radius_percentile, SEXP bandwidth_percentile, SEXP reconciliation_mode, SEXP min_overlap_coefficient, SEXP allowed_stop_reasons, SEXP filter_d_min, SEXP filter_d_max, SEXP filter_var_explained_min, SEXP ensemble_eligible, SEXP ensemble_eligible_by_stop_condition, SEXP ensemble_eligible_by_dimension, SEXP ensemble_eligible_by_var_explained, SEXP estimated_k_min, SEXP estimated_k_density, SEXP estimated_density_quantile, SEXP estimated_chordal_dist_max_as_prcnt_of_range, SEXP estimated_G_max, SEXP estimated_d_max) {
    int nprot = 0;
    // optionals: a null pointer and size 0 when the caller omits them
    int allowed_stop_reasons_size = 0;
    if (allowed_stop_reasons != R_NilValue) {
        allowed_stop_reasons_size = (int) Rf_length(allowed_stop_reasons);
    }
    const int* filter_d_min_p = NULL;
    int filter_d_min_size = 0;
    if (filter_d_min != R_NilValue) {
        filter_d_min_size = (int) Rf_length(filter_d_min);
        filter_d_min_p = INTEGER(filter_d_min);
    }
    const int* filter_d_max_p = NULL;
    int filter_d_max_size = 0;
    if (filter_d_max != R_NilValue) {
        filter_d_max_size = (int) Rf_length(filter_d_max);
        filter_d_max_p = INTEGER(filter_d_max);
    }
    const double* filter_var_explained_min_p = NULL;
    int filter_var_explained_min_size = 0;
    if (filter_var_explained_min != R_NilValue) {
        filter_var_explained_min_size = (int) Rf_length(filter_var_explained_min);
        filter_var_explained_min_p = REAL(filter_var_explained_min);
    }
    const int* estimated_k_min_p = NULL;
    int estimated_k_min_size = 0;
    if (estimated_k_min != R_NilValue) {
        estimated_k_min_size = (int) Rf_length(estimated_k_min);
        estimated_k_min_p = INTEGER(estimated_k_min);
    }
    const int* estimated_k_density_p = NULL;
    int estimated_k_density_size = 0;
    if (estimated_k_density != R_NilValue) {
        estimated_k_density_size = (int) Rf_length(estimated_k_density);
        estimated_k_density_p = INTEGER(estimated_k_density);
    }
    const double* estimated_density_quantile_p = NULL;
    int estimated_density_quantile_size = 0;
    if (estimated_density_quantile != R_NilValue) {
        estimated_density_quantile_size = (int) Rf_length(estimated_density_quantile);
        estimated_density_quantile_p = REAL(estimated_density_quantile);
    }
    const double* estimated_chordal_dist_max_as_prcnt_of_range_p = NULL;
    int estimated_chordal_dist_max_as_prcnt_of_range_size = 0;
    if (estimated_chordal_dist_max_as_prcnt_of_range != R_NilValue) {
        estimated_chordal_dist_max_as_prcnt_of_range_size = (int) Rf_length(estimated_chordal_dist_max_as_prcnt_of_range);
        estimated_chordal_dist_max_as_prcnt_of_range_p = REAL(estimated_chordal_dist_max_as_prcnt_of_range);
    }
    const double* estimated_G_max_p = NULL;
    int estimated_G_max_size = 0;
    if (estimated_G_max != R_NilValue) {
        estimated_G_max_size = (int) Rf_length(estimated_G_max);
        estimated_G_max_p = REAL(estimated_G_max);
    }
    const int* estimated_d_max_p = NULL;
    int estimated_d_max_size = 0;
    if (estimated_d_max != R_NilValue) {
        estimated_d_max_size = (int) Rf_length(estimated_d_max);
        estimated_d_max_p = INTEGER(estimated_d_max);
    }

    // derived from the inputs, not asked of the caller
    int filename_strlen = tox_max_strlen(filename);
    int n_dimensions = INTEGER(Rf_getAttrib(vectors, R_DimSymbol))[0];
    int n_vectors = INTEGER(Rf_getAttrib(vectors, R_DimSymbol))[1];
    int n_selected_seed = tox_sum_true(seed_selection_mask);
    int o = INTEGER(Rf_getAttrib(ensemble_U_history, R_DimSymbol))[2];
    int max_group_size = INTEGER(Rf_getAttrib(super_ensembles, R_DimSymbol))[0];
    int dim_names_strlen = tox_max_strlen(dim_names);

    // scalar inputs, pulled from their length-1 vectors
    int n_super_ensembles_v = Rf_asInteger(n_super_ensembles);
    int k_min_v = Rf_asInteger(k_min);
    int k_density_v = Rf_asInteger(k_density);
    double chordal_dist_max_as_prcnt_of_range_v = Rf_asReal(chordal_dist_max_as_prcnt_of_range);
    int d_max_v = Rf_asInteger(d_max);
    double G_max_v = Rf_asReal(G_max);
    double RMSE_change_max_v = Rf_asReal(RMSE_change_max);
    double f_max_v = Rf_asReal(f_max);
    int a_v = Rf_asInteger(a);
    double exclusion_radius_percentile_v = Rf_asReal(exclusion_radius_percentile);
    double bandwidth_percentile_v = Rf_asReal(bandwidth_percentile);
    double min_overlap_coefficient_v = Rf_asReal(min_overlap_coefficient);

    // convert what Fortran cannot take from R directly
    char* filename_c = tox_char_in(filename, filename_strlen);
    char* dim_names_c = tox_char_in(dim_names, dim_names_strlen);
    unsigned char* seed_selection_mask_c = tox_bool_in(seed_selection_mask);
    unsigned char* ensemble_masks_c = tox_bool_in(ensemble_masks);
    unsigned char* ensemble_accepted_history_c = tox_bool_in(ensemble_accepted_history);
    unsigned char* ensemble_low_confidence_masks_c = tox_bool_in(ensemble_low_confidence_masks);
    char* reconciliation_mode_c = tox_char_in(reconciliation_mode, 25);
    unsigned char* allowed_stop_reasons_c = tox_bool_in(allowed_stop_reasons);
    unsigned char* ensemble_eligible_c = tox_bool_in(ensemble_eligible);
    unsigned char* ensemble_eligible_by_stop_condition_c = tox_bool_in(ensemble_eligible_by_stop_condition);
    unsigned char* ensemble_eligible_by_dimension_c = tox_bool_in(ensemble_eligible_by_dimension);
    unsigned char* ensemble_eligible_by_var_explained_c = tox_bool_in(ensemble_eligible_by_var_explained);

    // outputs and work space
    int ierr = 0;

    write_stc_interactive_html_report_c(
        filename_c,
        &filename_strlen,
        &n_dimensions,
        &n_vectors,
        &n_selected_seed,
        &o,
        &max_group_size,
        &n_super_ensembles_v,
        REAL(vectors),
        dim_names_c,
        &dim_names_strlen,
        seed_selection_mask_c,
        ensemble_masks_c,
        INTEGER(ensemble_stop_reason),
        REAL(ensemble_growth_radii),
        REAL(ensemble_U_history),
        REAL(ensemble_S_history),
        INTEGER(ensemble_d_history),
        REAL(ensemble_G_history),
        REAL(ensemble_mu_history),
        INTEGER(ensemble_k_history),
        ensemble_accepted_history_c,
        INTEGER(ensemble_member_added_at_step),
        ensemble_low_confidence_masks_c,
        REAL(ensemble_U_first),
        INTEGER(ensemble_d_first),
        INTEGER(super_ensembles),
        &k_min_v,
        &k_density_v,
        &chordal_dist_max_as_prcnt_of_range_v,
        &d_max_v,
        &G_max_v,
        &RMSE_change_max_v,
        &f_max_v,
        &a_v,
        &exclusion_radius_percentile_v,
        &bandwidth_percentile_v,
        reconciliation_mode_c,
        &min_overlap_coefficient_v,
        allowed_stop_reasons != R_NilValue ? allowed_stop_reasons_c : NULL,
        filter_d_min_p,
        filter_d_max_p,
        filter_var_explained_min_p,
        ensemble_eligible_c,
        ensemble_eligible_by_stop_condition_c,
        ensemble_eligible_by_dimension_c,
        ensemble_eligible_by_var_explained_c,
        estimated_k_min_p,
        estimated_k_density_p,
        estimated_density_quantile_p,
        estimated_chordal_dist_max_as_prcnt_of_range_p,
        estimated_G_max_p,
        estimated_d_max_p,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 1)); nprot++;
    SET_VECTOR_ELT(_out, 0, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 1)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

#endif  // R binding
