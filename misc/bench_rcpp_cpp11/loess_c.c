/* Benchmark shim — pure C .Call path (no C++, no Rcpp, no cpp11).
 *
 * This is the hand-written equivalent of what a C-only emitter would generate:
 * the `SEXP` entry point plus the R_init registration. REAL()/INTEGER() hand
 * back the pointer INTO R's own vector buffer (zero-copy). The output is a
 * freshly Rf_allocVector'd vector that R owns.
 *
 * PROTECT is self-counted (nprot++ / UNPROTECT(nprot)) so the balance is never
 * hard-coded — the exact idiom an emitter would emit. There are no C++
 * destructors here, and R restores the PROTECT stack on any error longjmp, so
 * the cpp11 unwind_protect machinery has nothing to protect against.
 */
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

extern void loess_smooth_2d_c(const int*, const int*, const double*, const double*,
                              const int*, const int*, const double*, const double*,
                              const double*, double*, int*);

SEXP loess_smooth_2d_ccall(SEXP x_ref, SEXP y_ref, SEXP indices_used,
                           SEXP x_query, SEXP kernel_sigma, SEXP kernel_cutoff) {
    /* extents derived from the inputs, not asked of the caller */
    int n_total  = LENGTH(x_ref);
    int n_target = LENGTH(x_query);
    int n_used   = LENGTH(indices_used);
    double ks = Rf_asReal(kernel_sigma);
    double kc = Rf_asReal(kernel_cutoff);
    int ierr = 0, nprot = 0;

    SEXP y_out = PROTECT(Rf_allocVector(REALSXP, n_target)); nprot++;

    loess_smooth_2d_c(&n_total, &n_target,
                      REAL(x_ref), REAL(y_ref),          /* zero-copy: into R's buffer */
                      INTEGER(indices_used), &n_used,
                      REAL(x_query), &ks, &kc,
                      REAL(y_out), &ierr);

    /* return list(y_out=, ierr=) */
    SEXP out = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(out, 0, y_out);
    SET_VECTOR_ELT(out, 1, Rf_ScalarInteger(ierr));
    SEXP nms = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(nms, 0, Rf_mkChar("y_out"));
    SET_STRING_ELT(nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(out, R_NamesSymbol, nms);

    UNPROTECT(nprot);
    return out;
}

static const R_CallMethodDef CallEntries[] = {
    {"loess_smooth_2d_ccall", (DL_FUNC) &loess_smooth_2d_ccall, 6},
    {NULL, NULL, 0}
};

/* Name must match the .so basename: loess_c.so -> R_init_loess_c. */
void R_init_loess_c(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
