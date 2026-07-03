#include <Rcpp.h>

using namespace Rcpp;

static_assert(sizeof(Rcomplex) == sizeof(double _Complex) && alignof(Rcomplex) == alignof(double _Complex), "Rcomplex layout is incompatible with C's 'double _Complex' and thus Fortran's 'c_double_complex'");

extern "C" {
    void normalize_unit_length_c(
        double* vector,
        const int* n_dims,
        int* ierr
    );
}

List normalize_unit_length_rcpp(
    NumericVector& vector
) {


    IntegerVector vector_shape = vector.attr("dim");
    n_dims = vector_shape[0]



    int ierr = 0;

    normalize_unit_length_c(
        vector.begin(),
        &n_dims,
        &ierr
    );

//{format(self.arguments, "type_conversion_outputs")}

    return List::create(
        Named("ierr") = ierr
    );
}