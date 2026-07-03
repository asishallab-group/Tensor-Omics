#include <Rcpp.h>

using namespace Rcpp;

static_assert(sizeof(Rcomplex) == sizeof(double _Complex) && alignof(Rcomplex) == alignof(double _Complex), "Rcomplex layout is incompatible with C's 'double _Complex' and thus Fortran's 'c_double_complex'");

extern "C" {
    void tox_loess_required_workspace_c(
        const int* d,
        const int* nvmax,
        int* int_workspace_size,
        int* real_workspace_size,
        const int* setlf,
        int* ierr
    );

    void loess_fit_plain_c(
        const int* n,
        const double* x,
        const double* y,
        const double* w,
        const double* eval_points,
        const double* span,
        const int* degree,
        const int* nvmax,
        const int* infl,
        const int* setlf,
        int* int_workspace,
        const int* int_workspace_size,
        double* real_workspace,
        const int* real_workspace_size,
        double* diagl,
        double* fitted_values,
        int* ierr
    );

    void loess_fit_robust_c(
        const int* n,
        const double* x,
        const double* y,
        const double* w,
        const double* eval_points,
        const double* span,
        const int* degree,
        const int* nvmax,
        const int* infl,
        const int* setlf,
        const int* n_iters,
        int* int_workspace,
        const int* int_workspace_size,
        double* real_workspace,
        const int* real_workspace_size,
        double* diagl,
        double* robust_weights,
        double* combined_weights,
        double* residuals,
        int* permutation_indices,
        double* fitted_values,
        int* ierr
    );

    void loess_c(
        const double* x,
        const int* n_x_elements,
        const double* y,
        const int* n_y_elements,
        const double* span,
        const int* degree,
        double* fitted_values,
        const char* mode,
        const int* n_iters,
        int* ierr
    );
}

List tox_loess_required_workspace_rcpp(
    const int d,
    const int nvmax,
    const LogicalVector setlf
) {



    const int* setlf_p = setlf.begin();
    // Is already 0/1 mapped, thus check only for NA values and return ERR_NAN_INF code
    for (int i = 0; i < setlf.size(); ++i) {
        if (setlf_p[i] == NA_LOGICAL) {
            return List::create(Named("ierr") = 204);
        }
    }

    int int_workspace_size = 0;
    int real_workspace_size = 0;
    int ierr = 0;

    tox_loess_required_workspace_c(
        &d,
        &nvmax,
        &int_workspace_size,
        &real_workspace_size,
        &setlf,
        &ierr
    );

//{format(self.arguments, "type_conversion_outputs")}

    return List::create(
        Named("int_workspace_size") = int_workspace_size,
        Named("real_workspace_size") = real_workspace_size,
        Named("ierr") = ierr
    );
}

List loess_fit_plain_rcpp(
    const NumericVector& x,
    const NumericVector& y,
    const NumericVector& w,
    const NumericMatrix& eval_points,
    const double span,
    const int degree,
    const int nvmax,
    const LogicalVector infl,
    const LogicalVector setlf,
    IntegerVector& int_workspace,
    NumericVector& real_workspace,
    NumericVector& diagl
) {


    IntegerVector x_shape = x.attr("dim");
    n = x_shape[0]
    IntegerVector int_workspace_shape = int_workspace.attr("dim");
    int_workspace_size = int_workspace_shape[0]
    IntegerVector real_workspace_shape = real_workspace.attr("dim");
    real_workspace_size = real_workspace_shape[0]

    const int* infl_p = infl.begin();
    // Is already 0/1 mapped, thus check only for NA values and return ERR_NAN_INF code
    for (int i = 0; i < infl.size(); ++i) {
        if (infl_p[i] == NA_LOGICAL) {
            return List::create(Named("ierr") = 204);
        }
    }
    const int* setlf_p = setlf.begin();
    // Is already 0/1 mapped, thus check only for NA values and return ERR_NAN_INF code
    for (int i = 0; i < setlf.size(); ++i) {
        if (setlf_p[i] == NA_LOGICAL) {
            return List::create(Named("ierr") = 204);
        }
    }

    NumericVector fitted_values(n);
    int ierr = 0;

    loess_fit_plain_c(
        &n,
        x.begin(),
        y.begin(),
        w.begin(),
        eval_points.begin(),
        &span,
        &degree,
        &nvmax,
        &infl,
        &setlf,
        int_workspace.begin(),
        &int_workspace_size,
        real_workspace.begin(),
        &real_workspace_size,
        diagl.begin(),
        fitted_values.begin(),
        &ierr
    );

//{format(self.arguments, "type_conversion_outputs")}

    return List::create(
        Named("fitted_values") = fitted_values,
        Named("ierr") = ierr
    );
}

List loess_fit_robust_rcpp(
    const NumericVector& x,
    const NumericVector& y,
    const NumericVector& w,
    const NumericMatrix& eval_points,
    const double span,
    const int degree,
    const int nvmax,
    const LogicalVector infl,
    const LogicalVector setlf,
    const int n_iters,
    IntegerVector& int_workspace,
    NumericVector& real_workspace,
    NumericVector& diagl,
    NumericVector& robust_weights,
    NumericVector& combined_weights,
    NumericVector& residuals,
    IntegerVector& permutation_indices
) {


    IntegerVector x_shape = x.attr("dim");
    n = x_shape[0]
    IntegerVector int_workspace_shape = int_workspace.attr("dim");
    int_workspace_size = int_workspace_shape[0]
    IntegerVector real_workspace_shape = real_workspace.attr("dim");
    real_workspace_size = real_workspace_shape[0]

    const int* infl_p = infl.begin();
    // Is already 0/1 mapped, thus check only for NA values and return ERR_NAN_INF code
    for (int i = 0; i < infl.size(); ++i) {
        if (infl_p[i] == NA_LOGICAL) {
            return List::create(Named("ierr") = 204);
        }
    }
    const int* setlf_p = setlf.begin();
    // Is already 0/1 mapped, thus check only for NA values and return ERR_NAN_INF code
    for (int i = 0; i < setlf.size(); ++i) {
        if (setlf_p[i] == NA_LOGICAL) {
            return List::create(Named("ierr") = 204);
        }
    }

    NumericVector fitted_values(n);
    int ierr = 0;

    loess_fit_robust_c(
        &n,
        x.begin(),
        y.begin(),
        w.begin(),
        eval_points.begin(),
        &span,
        &degree,
        &nvmax,
        &infl,
        &setlf,
        &n_iters,
        int_workspace.begin(),
        &int_workspace_size,
        real_workspace.begin(),
        &real_workspace_size,
        diagl.begin(),
        robust_weights.begin(),
        combined_weights.begin(),
        residuals.begin(),
        permutation_indices.begin(),
        fitted_values.begin(),
        &ierr
    );

//{format(self.arguments, "type_conversion_outputs")}

    return List::create(
        Named("fitted_values") = fitted_values,
        Named("ierr") = ierr
    );
}

List loess_rcpp(
    const NumericVector& x,
    const NumericVector& y,
    const double span,
    const int degree,
    const String& mode,
    const int n_iters
) {


    IntegerVector x_shape = x.attr("dim");
    n_x_elements = x_shape[0]
    IntegerVector y_shape = y.attr("dim");
    n_y_elements = y_shape[0]
    IntegerVector fitted_values_shape = fitted_values.attr("dim");
    size(y) = fitted_values_shape[0]

    if (mode == NA_STRING) {
        return List::create(Named("ierr") = 204);
    }
    const char* mode_p = mode.get_cstring();

    NumericVector fitted_values(size(y));
    int ierr = 0;

    loess_c(
        x.begin(),
        &n_x_elements,
        y.begin(),
        &n_y_elements,
        &span,
        &degree,
        fitted_values.begin(),
        mode_p,
        &n_iters,
        &ierr
    );

//{format(self.arguments, "type_conversion_outputs")}

    return List::create(
        Named("fitted_values") = fitted_values,
        Named("ierr") = ierr
    );
}