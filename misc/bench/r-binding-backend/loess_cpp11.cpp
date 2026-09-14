// Benchmark shim — cpp11 .Call path (the candidate architecture).
//
// This is the hand-written equivalent of what a cpp11 emitter would generate:
// the cpp11-typed core, the `extern "C" SEXP` .Call entry point, and the
// R_init registration. Writing the glue by hand (instead of leaning on
// cpp11::cpp_source) is deliberate — it's exactly what the generator would emit
// into an R package, and it avoids cpp_source's decor/callr/... build-time deps.
//
// The read-only inputs are cpp11::doubles / cpp11::integers; REAL()/INTEGER()
// hand back the pointer INTO R's own vector buffer (zero-copy, like Rcpp's
// .begin()). The output is a freshly-allocated writable::doubles that R owns.
// PROTECT bookkeeping is RAII (cpp11), so a generator emitting this can't miscount.
#include <cpp11.hpp>
#include <cpp11/declarations.hpp>   // BEGIN_CPP11 / END_CPP11

extern "C" {
    void loess_smooth_2d_c(const int*, const int*, const double*, const double*,
                           const int*, const int*, const double*, const double*,
                           const double*, double*, int*);
}

static cpp11::writable::list loess_impl(cpp11::doubles x_ref,
                                        cpp11::doubles y_ref,
                                        cpp11::integers indices_used,
                                        cpp11::doubles x_query,
                                        double kernel_sigma,
                                        double kernel_cutoff) {
    // extents derived from the inputs, not asked of the caller
    int n_total  = (int) x_ref.size();
    int n_target = (int) x_query.size();
    int n_used   = (int) indices_used.size();

    cpp11::writable::doubles y_out(n_target);
    int ierr = 0;

    // Zero-copy: pointers into R-managed memory, no duplication of the inputs.
    SEXP s_x_ref = x_ref;
    SEXP s_y_ref = y_ref;
    SEXP s_idx   = indices_used;
    SEXP s_xq    = x_query;
    SEXP s_yout  = y_out;

    loess_smooth_2d_c(
        &n_total, &n_target,
        REAL(s_x_ref), REAL(s_y_ref),
        INTEGER(s_idx), &n_used,
        REAL(s_xq), &kernel_sigma, &kernel_cutoff,
        REAL(s_yout), &ierr);

    using namespace cpp11::literals;
    return cpp11::writable::list({"y_out"_nm = y_out, "ierr"_nm = ierr});
}

extern "C" SEXP loess_smooth_2d_cpp11_call(SEXP x_ref, SEXP y_ref, SEXP indices_used,
                                           SEXP x_query, SEXP kernel_sigma,
                                           SEXP kernel_cutoff) {
    BEGIN_CPP11
    cpp11::writable::list result = loess_impl(
        cpp11::as_cpp<cpp11::doubles>(x_ref),
        cpp11::as_cpp<cpp11::doubles>(y_ref),
        cpp11::as_cpp<cpp11::integers>(indices_used),
        cpp11::as_cpp<cpp11::doubles>(x_query),
        cpp11::as_cpp<double>(kernel_sigma),
        cpp11::as_cpp<double>(kernel_cutoff));
    return (SEXP) result;
    END_CPP11
}

static const R_CallMethodDef CallEntries[] = {
    {"loess_smooth_2d_cpp11", (DL_FUNC) &loess_smooth_2d_cpp11_call, 6},
    {NULL, NULL, 0}
};

// Name must match the .so basename: loess_cpp11.so -> R_init_loess_cpp11.
extern "C" void R_init_loess_cpp11(DllInfo* dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
