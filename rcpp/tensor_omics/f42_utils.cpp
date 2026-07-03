#include <Rcpp.h>

using namespace Rcpp;

static_assert(sizeof(Rcomplex) == sizeof(double _Complex) && alignof(Rcomplex) == alignof(double _Complex), "Rcomplex layout is incompatible with C's 'double _Complex' and thus Fortran's 'c_double_complex'");

extern "C" {
    void which_c(
        const int* mask,
        const int* n_mask_elements,
        const int* n,
        int* idx_out,
        const int* n_idx_out_elements,
        const int* m_max,
        int* m_out,
        int* ierr
    );

    void loess_smooth_2d_c(
        const int* n_total,
        const int* n_target,
        const double* x_ref,
        const double* y_ref,
        const int* indices_used,
        const int* n_used,
        const double* x_query,
        const double* kernel_sigma,
        const double* kernel_cutoff,
        double* y_out,
        int* ierr
    );

    void compute_edf_c(
        const double* values,
        const int* n_values,
        const int* perm,
        double* unique_values,
        double* cdf_values,
        int* n_unique,
        int* ierr
    );

    void compute_edf_c(
        const double* values,
        const int* n_values,
        double* unique_values,
        double* cdf_values,
        int* n_unique,
        int* ierr
    );

    void compute_empirical_p_values_c(
        const int* n_genes,
        const double* rdi,
        const double* sorted_rdi,
        const int* perm,
        double* p_values,
        const double* c_const,
        int* ierr
    );
}

List which_rcpp(
    const LogicalVector& mask,
    const int n,
    const int m_max
) {


    IntegerVector mask_shape = mask.attr("dim");
    n_mask_elements = mask_shape[0]
    IntegerVector idx_out_shape = idx_out.attr("dim");
    n_idx_out_elements = idx_out_shape[0]

    const int* mask_p = mask.begin();
    // Is already 0/1 mapped, thus check only for NA values and return ERR_NAN_INF code
    for (int i = 0; i < mask.size(); ++i) {
        if (mask_p[i] == NA_LOGICAL) {
            return List::create(Named("ierr") = 204);
        }
    }

    IntegerVector idx_out(n_idx_out_elements);
    int m_out = 0;
    int ierr = 0;

    which_c(
        mask_p,
        &n_mask_elements,
        &n,
        idx_out.begin(),
        &n_idx_out_elements,
        &m_max,
        &m_out,
        &ierr
    );

//{format(self.arguments, "type_conversion_outputs")}

    return List::create(
        Named("idx_out") = idx_out,
        Named("m_out") = m_out,
        Named("ierr") = ierr
    );
}

List loess_smooth_2d_rcpp(
    const NumericVector& x_ref,
    const NumericVector& y_ref,
    const IntegerVector& indices_used,
    const NumericVector& x_query,
    const double kernel_sigma,
    const double kernel_cutoff
) {


    IntegerVector x_ref_shape = x_ref.attr("dim");
    n_total = x_ref_shape[0]
    IntegerVector indices_used_shape = indices_used.attr("dim");
    n_used = indices_used_shape[0]
    IntegerVector x_query_shape = x_query.attr("dim");
    n_target = x_query_shape[0]



    NumericVector y_out(n_target);
    int ierr = 0;

    loess_smooth_2d_c(
        &n_total,
        &n_target,
        x_ref.begin(),
        y_ref.begin(),
        indices_used.begin(),
        &n_used,
        x_query.begin(),
        &kernel_sigma,
        &kernel_cutoff,
        y_out.begin(),
        &ierr
    );

//{format(self.arguments, "type_conversion_outputs")}

    return List::create(
        Named("y_out") = y_out,
        Named("ierr") = ierr
    );
}

List compute_edf_rcpp(
    const NumericVector& values,
    const IntegerVector& perm
) {


    IntegerVector values_shape = values.attr("dim");
    n_values = values_shape[0]



    NumericVector unique_values(n_values);
    NumericVector cdf_values(n_values);
    int n_unique = 0;
    int ierr = 0;

    compute_edf_c(
        values.begin(),
        &n_values,
        perm.begin(),
        unique_values.begin(),
        cdf_values.begin(),
        &n_unique,
        &ierr
    );

//{format(self.arguments, "type_conversion_outputs")}

    return List::create(
        Named("unique_values") = unique_values,
        Named("cdf_values") = cdf_values,
        Named("n_unique") = n_unique,
        Named("ierr") = ierr
    );
}

List compute_edf_rcpp(
    const NumericVector& values
) {


    IntegerVector values_shape = values.attr("dim");
    n_values = values_shape[0]



    NumericVector unique_values(n_values);
    NumericVector cdf_values(n_values);
    int n_unique = 0;
    int ierr = 0;

    compute_edf_c(
        values.begin(),
        &n_values,
        unique_values.begin(),
        cdf_values.begin(),
        &n_unique,
        &ierr
    );

//{format(self.arguments, "type_conversion_outputs")}

    return List::create(
        Named("unique_values") = unique_values,
        Named("cdf_values") = cdf_values,
        Named("n_unique") = n_unique,
        Named("ierr") = ierr
    );
}

List compute_empirical_p_values_rcpp(
    const NumericVector& rdi,
    const NumericVector& sorted_rdi,
    const IntegerVector& perm,
    const double c_const
) {


    IntegerVector rdi_shape = rdi.attr("dim");
    n_genes = rdi_shape[0]



    NumericVector p_values(n_genes);
    int ierr = 0;

    compute_empirical_p_values_c(
        &n_genes,
        rdi.begin(),
        sorted_rdi.begin(),
        perm.begin(),
        p_values.begin(),
        &c_const,
        &ierr
    );

//{format(self.arguments, "type_conversion_outputs")}

    return List::create(
        Named("p_values") = p_values,
        Named("ierr") = ierr
    );
}