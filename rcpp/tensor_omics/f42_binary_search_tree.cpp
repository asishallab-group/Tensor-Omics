#include <Rcpp.h>

using namespace Rcpp;

static_assert(sizeof(Rcomplex) == sizeof(double _Complex) && alignof(Rcomplex) == alignof(double _Complex), "Rcomplex layout is incompatible with C's 'double _Complex' and thus Fortran's 'c_double_complex'");

extern "C" {
    void get_sorted_value_c(
        const double* values,
        const int* n_values_elements,
        const int* sorted_indices,
        const int* n_sorted_indices_elements,
        const int* position,
        int* ierr,
        double* sorted_value
    );

    void build_bst_index_c(
        const double* values,
        const int* num_values,
        int* sorted_indices,
        int* tmp_left_stack,
        int* tmp_right_stack,
        int* ierr
    );

    void bst_range_query_c(
        const double* values,
        const int* sorted_indices,
        const int* num_values,
        const double* lower_bound,
        const double* upper_bound,
        int* output_indices,
        int* num_matches,
        int* ierr
    );
}

List get_sorted_value_rcpp(
    const NumericVector& values,
    const IntegerVector& sorted_indices,
    const int position
) {


    IntegerVector values_shape = values.attr("dim");
    n_values_elements = values_shape[0]
    IntegerVector sorted_indices_shape = sorted_indices.attr("dim");
    n_sorted_indices_elements = sorted_indices_shape[0]



    int ierr = 0;
    double sorted_value = 0.0;

    get_sorted_value_c(
        values.begin(),
        &n_values_elements,
        sorted_indices.begin(),
        &n_sorted_indices_elements,
        &position,
        &ierr,
        &sorted_value
    );

//{format(self.arguments, "type_conversion_outputs")}

    return List::create(
        Named("ierr") = ierr,
        Named("sorted_value") = sorted_value
    );
}

List build_bst_index_rcpp(
    const NumericVector& values
) {


    IntegerVector values_shape = values.attr("dim");
    num_values = values_shape[0]



    IntegerVector sorted_indices(num_values);
    std::vector<int> tmp_left_stack(num_values);
    std::vector<int> tmp_right_stack(num_values);
    int ierr = 0;

    build_bst_index_c(
        values.begin(),
        &num_values,
        sorted_indices.begin(),
        tmp_left_stack.data(),
        tmp_right_stack.data(),
        &ierr
    );

//{format(self.arguments, "type_conversion_outputs")}

    return List::create(
        Named("sorted_indices") = sorted_indices,
        Named("ierr") = ierr
    );
}

List bst_range_query_rcpp(
    const NumericVector& values,
    const IntegerVector& sorted_indices,
    const double lower_bound,
    const double upper_bound
) {


    IntegerVector values_shape = values.attr("dim");
    num_values = values_shape[0]



    IntegerVector output_indices(num_values);
    int num_matches = 0;
    int ierr = 0;

    bst_range_query_c(
        values.begin(),
        sorted_indices.begin(),
        &num_values,
        &lower_bound,
        &upper_bound,
        output_indices.begin(),
        &num_matches,
        &ierr
    );

//{format(self.arguments, "type_conversion_outputs")}

    return List::create(
        Named("output_indices") = output_indices,
        Named("num_matches") = num_matches,
        Named("ierr") = ierr
    );
}