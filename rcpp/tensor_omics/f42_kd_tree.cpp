#include <Rcpp.h>

using namespace Rcpp;

static_assert(sizeof(Rcomplex) == sizeof(double _Complex) && alignof(Rcomplex) == alignof(double _Complex), "Rcomplex layout is incompatible with C's 'double _Complex' and thus Fortran's 'c_double_complex'");

extern "C" {
    void build_kd_index_c(
        const double* points,
        const int* num_dimensions,
        const int* num_points,
        int* kd_indices,
        const int* dimension_order,
        int* tmp_workspace,
        double* tmp_value_buffer,
        int* tmp_permutation,
        int* tmp_left_stack,
        int* tmp_right_stack,
        int* tmp_recursion_stack,
        int* ierr
    );

    void build_spherical_kd_c(
        const double* vectors,
        const int* num_dimensions,
        const int* num_vectors,
        int* sphere_indices,
        const int* dimension_order,
        int* tmp_workspace,
        double* tmp_value_buffer,
        int* tmp_permutation,
        int* tmp_left_stack,
        int* tmp_right_stack,
        int* tmp_recursion_stack,
        int* ierr
    );

    void get_kd_point_c(
        const double* points,
        const int* n_points_elements_dim_1,
        const int* n_points_elements_dim_2,
        const int* kd_indices,
        const int* n_kd_indices_elements,
        const int* position,
        double* point_values,
        const int* n_point_values_elements,
        int* ierr
    );
}

List build_kd_index_rcpp(
    const NumericMatrix& points,
    const IntegerVector& dimension_order
) {


    IntegerVector points_shape = points.attr("dim");
    num_dimensions = points_shape[0]
    num_points = points_shape[1]



    IntegerVector kd_indices(num_points);
    std::vector<int> tmp_workspace(num_points);
    std::vector<double> tmp_value_buffer(num_points);
    std::vector<int> tmp_permutation(num_points);
    std::vector<int> tmp_left_stack(num_points);
    std::vector<int> tmp_right_stack(num_points);
    std::vector<int> tmp_recursion_stack(3 * num_points);
    int ierr = 0;

    build_kd_index_c(
        points.begin(),
        &num_dimensions,
        &num_points,
        kd_indices.begin(),
        dimension_order.begin(),
        tmp_workspace.data(),
        tmp_value_buffer.data(),
        tmp_permutation.data(),
        tmp_left_stack.data(),
        tmp_right_stack.data(),
        tmp_recursion_stack.data(),
        &ierr
    );

//{format(self.arguments, "type_conversion_outputs")}

    return List::create(
        Named("kd_indices") = kd_indices,
        Named("ierr") = ierr
    );
}

List build_spherical_kd_rcpp(
    const NumericMatrix& vectors,
    const IntegerVector& dimension_order
) {


    IntegerVector vectors_shape = vectors.attr("dim");
    num_dimensions = vectors_shape[0]
    num_vectors = vectors_shape[1]



    IntegerVector sphere_indices(num_vectors);
    std::vector<int> tmp_workspace(num_vectors);
    std::vector<double> tmp_value_buffer(num_vectors);
    std::vector<int> tmp_permutation(num_vectors);
    std::vector<int> tmp_left_stack(num_vectors);
    std::vector<int> tmp_right_stack(num_vectors);
    std::vector<int> tmp_recursion_stack(3 * num_vectors);
    int ierr = 0;

    build_spherical_kd_c(
        vectors.begin(),
        &num_dimensions,
        &num_vectors,
        sphere_indices.begin(),
        dimension_order.begin(),
        tmp_workspace.data(),
        tmp_value_buffer.data(),
        tmp_permutation.data(),
        tmp_left_stack.data(),
        tmp_right_stack.data(),
        tmp_recursion_stack.data(),
        &ierr
    );

//{format(self.arguments, "type_conversion_outputs")}

    return List::create(
        Named("sphere_indices") = sphere_indices,
        Named("ierr") = ierr
    );
}

List get_kd_point_rcpp(
    const NumericMatrix& points,
    const IntegerVector& kd_indices,
    const int position
) {


    IntegerVector points_shape = points.attr("dim");
    n_points_elements_dim_1 = points_shape[0]
    n_points_elements_dim_2 = points_shape[1]
    IntegerVector kd_indices_shape = kd_indices.attr("dim");
    n_kd_indices_elements = kd_indices_shape[0]
    IntegerVector point_values_shape = point_values.attr("dim");
    n_point_values_elements = point_values_shape[0]



    NumericVector point_values(n_point_values_elements);
    int ierr = 0;

    get_kd_point_c(
        points.begin(),
        &n_points_elements_dim_1,
        &n_points_elements_dim_2,
        kd_indices.begin(),
        &n_kd_indices_elements,
        &position,
        point_values.begin(),
        &n_point_values_elements,
        &ierr
    );

//{format(self.arguments, "type_conversion_outputs")}

    return List::create(
        Named("point_values") = point_values,
        Named("ierr") = ierr
    );
}