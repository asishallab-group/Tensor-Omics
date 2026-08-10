// Generated. Do not edit.
#if !defined(NO_R_BINDING) && !defined(NO_C_BINDING)
#include <R.h>
#include <Rinternals.h>
#include "tox_marshal.h"

// the Fortran C-ABI symbols this module calls
void ensemble_identification_c(const double*, const int*, const int*, const int*, const int*, const int*, const int*, const double*, const int*, const double*, const double*, const double*, const int*, const int*, unsigned char*, int*, double*, double*, double*, int*, double*, double*, int*, unsigned char*, int*, unsigned char*, double*, int*, int*);
void ensemble_identification_merged_c(const double*, const int*, const int*, const int*, const int*, const unsigned char*, const int*, const int*, const double*, const int*, const double*, const double*, const double*, const int*, const int*, unsigned char*, int*, double*, double*, double*, int*, double*, double*, int*, unsigned char*, int*, unsigned char*, double*, int*, int*);

SEXP ensemble_identification_call(SEXP vectors, SEXP kd_indices, SEXP dimension_order, SEXP seed_index, SEXP k_min, SEXP chordal_dist_max_as_prcnt_of_range, SEXP d_max, SEXP G_max, SEXP RMSE_change_max, SEXP f_max, SEXP a, SEXP o) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_dimensions = INTEGER(Rf_getAttrib(vectors, R_DimSymbol))[0];
    int n_vectors = INTEGER(Rf_getAttrib(vectors, R_DimSymbol))[1];

    // scalar inputs, pulled from their length-1 vectors
    int seed_index_v = Rf_asInteger(seed_index);
    int k_min_v = Rf_asInteger(k_min);
    double chordal_dist_max_as_prcnt_of_range_v = Rf_asReal(chordal_dist_max_as_prcnt_of_range);
    int d_max_v = Rf_asInteger(d_max);
    double G_max_v = Rf_asReal(G_max);
    double RMSE_change_max_v = Rf_asReal(RMSE_change_max);
    double f_max_v = Rf_asReal(f_max);
    int a_v = Rf_asInteger(a);
    int o_v = Rf_asInteger(o);

    // outputs and work space
    unsigned char* final_ensemble_mask_c = tox_bool_alloc(n_vectors);
    int stop_reason = 0;
    double growth_radius = 0;
    SEXP U_history = PROTECT(Rf_allocVector(REALSXP, n_dimensions * n_dimensions * o_v)); nprot++;
    { SEXP U_history_dim = PROTECT(Rf_allocVector(INTSXP, 3)); INTEGER(U_history_dim)[0] = n_dimensions; INTEGER(U_history_dim)[1] = n_dimensions; INTEGER(U_history_dim)[2] = o_v; Rf_setAttrib(U_history, R_DimSymbol, U_history_dim); UNPROTECT(1); }
    SEXP S_history = PROTECT(Rf_allocVector(REALSXP, n_dimensions * o_v)); nprot++;
    { SEXP S_history_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(S_history_dim)[0] = n_dimensions; INTEGER(S_history_dim)[1] = o_v; Rf_setAttrib(S_history, R_DimSymbol, S_history_dim); UNPROTECT(1); }
    SEXP d_history = PROTECT(Rf_allocVector(INTSXP, o_v)); nprot++;
    SEXP G_history = PROTECT(Rf_allocVector(REALSXP, o_v)); nprot++;
    SEXP mu_history = PROTECT(Rf_allocVector(REALSXP, n_dimensions * o_v)); nprot++;
    { SEXP mu_history_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(mu_history_dim)[0] = n_dimensions; INTEGER(mu_history_dim)[1] = o_v; Rf_setAttrib(mu_history, R_DimSymbol, mu_history_dim); UNPROTECT(1); }
    SEXP k_history = PROTECT(Rf_allocVector(INTSXP, o_v)); nprot++;
    unsigned char* accepted_history_c = tox_bool_alloc(o_v);
    SEXP member_added_at_step = PROTECT(Rf_allocVector(INTSXP, n_vectors)); nprot++;
    unsigned char* low_confidence_mask_c = tox_bool_alloc(n_vectors);
    SEXP U_first = PROTECT(Rf_allocVector(REALSXP, n_dimensions * n_dimensions)); nprot++;
    { SEXP U_first_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(U_first_dim)[0] = n_dimensions; INTEGER(U_first_dim)[1] = n_dimensions; Rf_setAttrib(U_first, R_DimSymbol, U_first_dim); UNPROTECT(1); }
    int d_first = 0;
    int ierr = 0;

    ensemble_identification_c(
        REAL(vectors),
        &n_dimensions,
        &n_vectors,
        INTEGER(kd_indices),
        INTEGER(dimension_order),
        &seed_index_v,
        &k_min_v,
        &chordal_dist_max_as_prcnt_of_range_v,
        &d_max_v,
        &G_max_v,
        &RMSE_change_max_v,
        &f_max_v,
        &a_v,
        &o_v,
        final_ensemble_mask_c,
        &stop_reason,
        &growth_radius,
        REAL(U_history),
        REAL(S_history),
        INTEGER(d_history),
        REAL(G_history),
        REAL(mu_history),
        INTEGER(k_history),
        accepted_history_c,
        INTEGER(member_added_at_step),
        low_confidence_mask_c,
        REAL(U_first),
        &d_first,
        &ierr
    );

    // convert the outputs back
    SEXP final_ensemble_mask = PROTECT(tox_bool_out(final_ensemble_mask_c, n_vectors)); nprot++;
    SEXP accepted_history = PROTECT(tox_bool_out(accepted_history_c, o_v)); nprot++;
    SEXP low_confidence_mask = PROTECT(tox_bool_out(low_confidence_mask_c, n_vectors)); nprot++;

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 15)); nprot++;
    SET_VECTOR_ELT(_out, 0, final_ensemble_mask);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(stop_reason));
    SET_VECTOR_ELT(_out, 2, Rf_ScalarReal(growth_radius));
    SET_VECTOR_ELT(_out, 3, U_history);
    SET_VECTOR_ELT(_out, 4, S_history);
    SET_VECTOR_ELT(_out, 5, d_history);
    SET_VECTOR_ELT(_out, 6, G_history);
    SET_VECTOR_ELT(_out, 7, mu_history);
    SET_VECTOR_ELT(_out, 8, k_history);
    SET_VECTOR_ELT(_out, 9, accepted_history);
    SET_VECTOR_ELT(_out, 10, member_added_at_step);
    SET_VECTOR_ELT(_out, 11, low_confidence_mask);
    SET_VECTOR_ELT(_out, 12, U_first);
    SET_VECTOR_ELT(_out, 13, Rf_ScalarInteger(d_first));
    SET_VECTOR_ELT(_out, 14, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 15)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("final_ensemble_mask"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("stop_reason"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("growth_radius"));
    SET_STRING_ELT(_nms, 3, Rf_mkChar("U_history"));
    SET_STRING_ELT(_nms, 4, Rf_mkChar("S_history"));
    SET_STRING_ELT(_nms, 5, Rf_mkChar("d_history"));
    SET_STRING_ELT(_nms, 6, Rf_mkChar("G_history"));
    SET_STRING_ELT(_nms, 7, Rf_mkChar("mu_history"));
    SET_STRING_ELT(_nms, 8, Rf_mkChar("k_history"));
    SET_STRING_ELT(_nms, 9, Rf_mkChar("accepted_history"));
    SET_STRING_ELT(_nms, 10, Rf_mkChar("member_added_at_step"));
    SET_STRING_ELT(_nms, 11, Rf_mkChar("low_confidence_mask"));
    SET_STRING_ELT(_nms, 12, Rf_mkChar("U_first"));
    SET_STRING_ELT(_nms, 13, Rf_mkChar("d_first"));
    SET_STRING_ELT(_nms, 14, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP ensemble_identification_merged_call(SEXP vectors, SEXP kd_indices, SEXP dimension_order, SEXP seed_selection_mask, SEXP k_min, SEXP chordal_dist_max_as_prcnt_of_range, SEXP d_max, SEXP G_max, SEXP RMSE_change_max, SEXP f_max, SEXP a, SEXP o) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_dimensions = INTEGER(Rf_getAttrib(vectors, R_DimSymbol))[0];
    int n_vectors = INTEGER(Rf_getAttrib(vectors, R_DimSymbol))[1];
    int n_selected_seed = tox_sum_true(seed_selection_mask);

    // scalar inputs, pulled from their length-1 vectors
    int k_min_v = Rf_asInteger(k_min);
    double chordal_dist_max_as_prcnt_of_range_v = Rf_asReal(chordal_dist_max_as_prcnt_of_range);
    int d_max_v = Rf_asInteger(d_max);
    double G_max_v = Rf_asReal(G_max);
    double RMSE_change_max_v = Rf_asReal(RMSE_change_max);
    double f_max_v = Rf_asReal(f_max);
    int a_v = Rf_asInteger(a);
    int o_v = Rf_asInteger(o);

    // convert what Fortran cannot take from R directly
    unsigned char* seed_selection_mask_c = tox_bool_in(seed_selection_mask);

    // outputs and work space
    unsigned char* ensemble_masks_c = tox_bool_alloc(n_vectors * n_selected_seed);
    SEXP ensemble_stop_reason = PROTECT(Rf_allocVector(INTSXP, n_selected_seed)); nprot++;
    SEXP ensemble_growth_radii = PROTECT(Rf_allocVector(REALSXP, n_selected_seed)); nprot++;
    SEXP ensemble_U_history = PROTECT(Rf_allocVector(REALSXP, n_dimensions * n_dimensions * o_v * n_selected_seed)); nprot++;
    { SEXP ensemble_U_history_dim = PROTECT(Rf_allocVector(INTSXP, 4)); INTEGER(ensemble_U_history_dim)[0] = n_dimensions; INTEGER(ensemble_U_history_dim)[1] = n_dimensions; INTEGER(ensemble_U_history_dim)[2] = o_v; INTEGER(ensemble_U_history_dim)[3] = n_selected_seed; Rf_setAttrib(ensemble_U_history, R_DimSymbol, ensemble_U_history_dim); UNPROTECT(1); }
    SEXP ensemble_S_history = PROTECT(Rf_allocVector(REALSXP, n_dimensions * o_v * n_selected_seed)); nprot++;
    { SEXP ensemble_S_history_dim = PROTECT(Rf_allocVector(INTSXP, 3)); INTEGER(ensemble_S_history_dim)[0] = n_dimensions; INTEGER(ensemble_S_history_dim)[1] = o_v; INTEGER(ensemble_S_history_dim)[2] = n_selected_seed; Rf_setAttrib(ensemble_S_history, R_DimSymbol, ensemble_S_history_dim); UNPROTECT(1); }
    SEXP ensemble_d_history = PROTECT(Rf_allocVector(INTSXP, o_v * n_selected_seed)); nprot++;
    { SEXP ensemble_d_history_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(ensemble_d_history_dim)[0] = o_v; INTEGER(ensemble_d_history_dim)[1] = n_selected_seed; Rf_setAttrib(ensemble_d_history, R_DimSymbol, ensemble_d_history_dim); UNPROTECT(1); }
    SEXP ensemble_G_history = PROTECT(Rf_allocVector(REALSXP, o_v * n_selected_seed)); nprot++;
    { SEXP ensemble_G_history_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(ensemble_G_history_dim)[0] = o_v; INTEGER(ensemble_G_history_dim)[1] = n_selected_seed; Rf_setAttrib(ensemble_G_history, R_DimSymbol, ensemble_G_history_dim); UNPROTECT(1); }
    SEXP ensemble_mu_history = PROTECT(Rf_allocVector(REALSXP, n_dimensions * o_v * n_selected_seed)); nprot++;
    { SEXP ensemble_mu_history_dim = PROTECT(Rf_allocVector(INTSXP, 3)); INTEGER(ensemble_mu_history_dim)[0] = n_dimensions; INTEGER(ensemble_mu_history_dim)[1] = o_v; INTEGER(ensemble_mu_history_dim)[2] = n_selected_seed; Rf_setAttrib(ensemble_mu_history, R_DimSymbol, ensemble_mu_history_dim); UNPROTECT(1); }
    SEXP ensemble_k_history = PROTECT(Rf_allocVector(INTSXP, o_v * n_selected_seed)); nprot++;
    { SEXP ensemble_k_history_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(ensemble_k_history_dim)[0] = o_v; INTEGER(ensemble_k_history_dim)[1] = n_selected_seed; Rf_setAttrib(ensemble_k_history, R_DimSymbol, ensemble_k_history_dim); UNPROTECT(1); }
    unsigned char* ensemble_accepted_history_c = tox_bool_alloc(o_v * n_selected_seed);
    SEXP ensemble_member_added_at_step = PROTECT(Rf_allocVector(INTSXP, n_vectors * n_selected_seed)); nprot++;
    { SEXP ensemble_member_added_at_step_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(ensemble_member_added_at_step_dim)[0] = n_vectors; INTEGER(ensemble_member_added_at_step_dim)[1] = n_selected_seed; Rf_setAttrib(ensemble_member_added_at_step, R_DimSymbol, ensemble_member_added_at_step_dim); UNPROTECT(1); }
    unsigned char* ensemble_low_confidence_masks_c = tox_bool_alloc(n_vectors * n_selected_seed);
    SEXP ensemble_U_first = PROTECT(Rf_allocVector(REALSXP, n_dimensions * n_dimensions * n_selected_seed)); nprot++;
    { SEXP ensemble_U_first_dim = PROTECT(Rf_allocVector(INTSXP, 3)); INTEGER(ensemble_U_first_dim)[0] = n_dimensions; INTEGER(ensemble_U_first_dim)[1] = n_dimensions; INTEGER(ensemble_U_first_dim)[2] = n_selected_seed; Rf_setAttrib(ensemble_U_first, R_DimSymbol, ensemble_U_first_dim); UNPROTECT(1); }
    SEXP ensemble_d_first = PROTECT(Rf_allocVector(INTSXP, n_selected_seed)); nprot++;
    int ierr = 0;

    ensemble_identification_merged_c(
        REAL(vectors),
        &n_dimensions,
        &n_vectors,
        INTEGER(kd_indices),
        INTEGER(dimension_order),
        seed_selection_mask_c,
        &n_selected_seed,
        &k_min_v,
        &chordal_dist_max_as_prcnt_of_range_v,
        &d_max_v,
        &G_max_v,
        &RMSE_change_max_v,
        &f_max_v,
        &a_v,
        &o_v,
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
        &ierr
    );

    // convert the outputs back
    SEXP ensemble_masks = PROTECT(tox_bool_out(ensemble_masks_c, n_vectors * n_selected_seed)); nprot++;
    { SEXP ensemble_masks_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(ensemble_masks_dim)[0] = n_vectors; INTEGER(ensemble_masks_dim)[1] = n_selected_seed; Rf_setAttrib(ensemble_masks, R_DimSymbol, ensemble_masks_dim); UNPROTECT(1); }
    SEXP ensemble_accepted_history = PROTECT(tox_bool_out(ensemble_accepted_history_c, o_v * n_selected_seed)); nprot++;
    { SEXP ensemble_accepted_history_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(ensemble_accepted_history_dim)[0] = o_v; INTEGER(ensemble_accepted_history_dim)[1] = n_selected_seed; Rf_setAttrib(ensemble_accepted_history, R_DimSymbol, ensemble_accepted_history_dim); UNPROTECT(1); }
    SEXP ensemble_low_confidence_masks = PROTECT(tox_bool_out(ensemble_low_confidence_masks_c, n_vectors * n_selected_seed)); nprot++;
    { SEXP ensemble_low_confidence_masks_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(ensemble_low_confidence_masks_dim)[0] = n_vectors; INTEGER(ensemble_low_confidence_masks_dim)[1] = n_selected_seed; Rf_setAttrib(ensemble_low_confidence_masks, R_DimSymbol, ensemble_low_confidence_masks_dim); UNPROTECT(1); }

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 15)); nprot++;
    SET_VECTOR_ELT(_out, 0, ensemble_masks);
    SET_VECTOR_ELT(_out, 1, ensemble_stop_reason);
    SET_VECTOR_ELT(_out, 2, ensemble_growth_radii);
    SET_VECTOR_ELT(_out, 3, ensemble_U_history);
    SET_VECTOR_ELT(_out, 4, ensemble_S_history);
    SET_VECTOR_ELT(_out, 5, ensemble_d_history);
    SET_VECTOR_ELT(_out, 6, ensemble_G_history);
    SET_VECTOR_ELT(_out, 7, ensemble_mu_history);
    SET_VECTOR_ELT(_out, 8, ensemble_k_history);
    SET_VECTOR_ELT(_out, 9, ensemble_accepted_history);
    SET_VECTOR_ELT(_out, 10, ensemble_member_added_at_step);
    SET_VECTOR_ELT(_out, 11, ensemble_low_confidence_masks);
    SET_VECTOR_ELT(_out, 12, ensemble_U_first);
    SET_VECTOR_ELT(_out, 13, ensemble_d_first);
    SET_VECTOR_ELT(_out, 14, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 15)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("ensemble_masks"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ensemble_stop_reason"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("ensemble_growth_radii"));
    SET_STRING_ELT(_nms, 3, Rf_mkChar("ensemble_U_history"));
    SET_STRING_ELT(_nms, 4, Rf_mkChar("ensemble_S_history"));
    SET_STRING_ELT(_nms, 5, Rf_mkChar("ensemble_d_history"));
    SET_STRING_ELT(_nms, 6, Rf_mkChar("ensemble_G_history"));
    SET_STRING_ELT(_nms, 7, Rf_mkChar("ensemble_mu_history"));
    SET_STRING_ELT(_nms, 8, Rf_mkChar("ensemble_k_history"));
    SET_STRING_ELT(_nms, 9, Rf_mkChar("ensemble_accepted_history"));
    SET_STRING_ELT(_nms, 10, Rf_mkChar("ensemble_member_added_at_step"));
    SET_STRING_ELT(_nms, 11, Rf_mkChar("ensemble_low_confidence_masks"));
    SET_STRING_ELT(_nms, 12, Rf_mkChar("ensemble_U_first"));
    SET_STRING_ELT(_nms, 13, Rf_mkChar("ensemble_d_first"));
    SET_STRING_ELT(_nms, 14, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

#endif  // R binding
