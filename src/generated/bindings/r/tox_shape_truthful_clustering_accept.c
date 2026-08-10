// Generated. Do not edit.
#if !defined(NO_R_BINDING) && !defined(NO_C_BINDING)
#include <R.h>
#include <Rinternals.h>
#include "tox_marshal.h"

// the Fortran C-ABI symbols this module calls
void accept_ensemble_c(const int*, const int*, const double*, const int*, const double*, const int*, const int*, const double*, const double*, const double*, const int*, const double*, const double*, const double*, const int*, const double*, const double*, unsigned char*, int*);

SEXP accept_ensemble_call(SEXP U_first, SEXP d_first, SEXP U_history, SEXP d_history, SEXP history_len, SEXP G_t, SEXP normal_error_t, SEXP U_tp1, SEXP d_tp1, SEXP G_tp1, SEXP normal_error_tp1, SEXP chordal_dist_max_as_prcnt_of_range, SEXP d_max, SEXP G_max, SEXP RMSE_change_max) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_dimensions = INTEGER(Rf_getAttrib(U_first, R_DimSymbol))[0];
    int o = INTEGER(Rf_getAttrib(U_history, R_DimSymbol))[2];

    // scalar inputs, pulled from their length-1 vectors
    int d_first_v = Rf_asInteger(d_first);
    int history_len_v = Rf_asInteger(history_len);
    double G_t_v = Rf_asReal(G_t);
    double normal_error_t_v = Rf_asReal(normal_error_t);
    int d_tp1_v = Rf_asInteger(d_tp1);
    double G_tp1_v = Rf_asReal(G_tp1);
    double normal_error_tp1_v = Rf_asReal(normal_error_tp1);
    double chordal_dist_max_as_prcnt_of_range_v = Rf_asReal(chordal_dist_max_as_prcnt_of_range);
    int d_max_v = Rf_asInteger(d_max);
    double G_max_v = Rf_asReal(G_max);
    double RMSE_change_max_v = Rf_asReal(RMSE_change_max);

    // outputs and work space
    unsigned char is_accepted = 0;
    int ierr = 0;

    accept_ensemble_c(
        &n_dimensions,
        &o,
        REAL(U_first),
        &d_first_v,
        REAL(U_history),
        INTEGER(d_history),
        &history_len_v,
        &G_t_v,
        &normal_error_t_v,
        REAL(U_tp1),
        &d_tp1_v,
        &G_tp1_v,
        &normal_error_tp1_v,
        &chordal_dist_max_as_prcnt_of_range_v,
        &d_max_v,
        &G_max_v,
        &RMSE_change_max_v,
        &is_accepted,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(_out, 0, Rf_ScalarLogical(is_accepted != 0));
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("is_accepted"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

#endif  // R binding
