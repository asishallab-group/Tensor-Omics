// Generated. Do not edit.
#if !defined(NO_R_BINDING) && !defined(NO_C_BINDING)
#include <R.h>
#include <Rinternals.h>
#include "tox_marshal.h"

// the Fortran C-ABI symbols this module calls
void accept_ensemble_c(const int*, const double*, const int*, const double*, const double*, const int*, const double*, const double*, const int*, const double*, unsigned char*, int*);

SEXP accept_ensemble_call(SEXP U_t, SEXP d_t, SEXP G_t, SEXP U_tp1, SEXP d_tp1, SEXP G_tp1, SEXP alpha_max, SEXP d_max, SEXP G_max) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_dimensions = INTEGER(Rf_getAttrib(U_t, R_DimSymbol))[0];

    // scalar inputs, pulled from their length-1 vectors
    int d_t_v = Rf_asInteger(d_t);
    double G_t_v = Rf_asReal(G_t);
    int d_tp1_v = Rf_asInteger(d_tp1);
    double G_tp1_v = Rf_asReal(G_tp1);
    double alpha_max_v = Rf_asReal(alpha_max);
    int d_max_v = Rf_asInteger(d_max);
    double G_max_v = Rf_asReal(G_max);

    // outputs and work space
    unsigned char is_accepted = 0;
    int ierr = 0;

    accept_ensemble_c(
        &n_dimensions,
        REAL(U_t),
        &d_t_v,
        &G_t_v,
        REAL(U_tp1),
        &d_tp1_v,
        &G_tp1_v,
        &alpha_max_v,
        &d_max_v,
        &G_max_v,
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
