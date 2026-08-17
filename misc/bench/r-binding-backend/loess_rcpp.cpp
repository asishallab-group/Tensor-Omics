// Benchmark shim — Rcpp path (the architecture we currently generate).
//
// Logic is identical to the generated rcpp/tensor_omics/src/f42_utils.cpp entry
// for loess_smooth_2d, trimmed to just this one routine so the comparison is
// self-contained. Inputs are Rcpp vectors that ALIAS R's buffers (zero-copy);
// `.begin()` is the pointer into that buffer.
#include <Rcpp.h>

using namespace Rcpp;

extern "C" {
    void loess_smooth_2d_c(const int*, const int*, const double*, const double*,
                           const int*, const int*, const double*, const double*,
                           const double*, double*, int*);
}

// [[Rcpp::export(.loess_smooth_2d_rcpp)]]
List loess_smooth_2d_rcpp(NumericVector x_ref, NumericVector y_ref,
                          IntegerVector indices_used, NumericVector x_query,
                          double kernel_sigma, double kernel_cutoff) {
    // extents derived from the inputs, not asked of the caller
    int n_total  = (int) x_ref.size();
    int n_target = (int) x_query.size();
    int n_used   = (int) indices_used.size();

    NumericVector y_out(n_target);
    int ierr = 0;

    loess_smooth_2d_c(
        &n_total, &n_target,
        x_ref.begin(), y_ref.begin(),
        indices_used.begin(), &n_used,
        x_query.begin(), &kernel_sigma, &kernel_cutoff,
        y_out.begin(), &ierr);

    return List::create(_["y_out"] = y_out, _["ierr"] = ierr);
}
