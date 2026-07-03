#include <Rcpp.h>

using namespace Rcpp;

static_assert(sizeof(Rcomplex) == sizeof(double _Complex) && alignof(Rcomplex) == alignof(double _Complex), "Rcomplex layout is incompatible with C's 'double _Complex' and thus Fortran's 'c_double_complex'");

extern "C" {
    void omics_vector_RAP_projection_c(
        const double* vecs,
        const int* n_axes,
        const int* n_vecs,
        const int* vecs_selection_mask,
        const int* n_selected_vecs,
        const int* axes_selection_mask,
        const int* n_selected_axes,
        double* projections,
        int* ierr
    );

    void omics_field_RAP_projection_c(
        const double* vecs,
        const int* n_axes,
        const int* n_vecs,
        const int* vecs_selection_mask,
        const int* n_selected_vecs,
        const int* axes_selection_mask,
        const int* n_selected_axes,
        double* projections,
        int* ierr
    );

    void clock_hand_angle_between_vectors_c(
        const double* v1,
        const double* v2,
        const int* n_dims,
        double* signed_angle,
        const int* selected_axes_for_signed,
        int* ierr
    );

    void clock_hand_angles_for_shift_vectors_c(
        const double* origins,
        const double* targets,
        const int* n_dims,
        const int* n_vecs,
        const int* vecs_selection_mask,
        const int* n_selected_vecs,
        const int* selected_axes_for_signed,
        double* signed_angles,
        int* ierr
    );

    void compute_relative_axis_contributions_c(
        const double* vec,
        const int* n_axes,
        double* contributions,
        int* ierr
    );

    void relative_axes_changes_from_shift_vector_c(
        const double* vec,
        const int* n_axes,
        double* contributions,
        int* ierr
    );

    void relative_axes_expression_from_expression_vector_c(
        const double* vec,
        const int* n_axes,
        double* contributions,
        int* ierr
    );
}

List omics_vector_RAP_projection_rcpp(
    const NumericMatrix& vecs,
    const LogicalVector& vecs_selection_mask,
    const LogicalVector& axes_selection_mask
) {


    IntegerVector vecs_shape = vecs.attr("dim");
    n_axes = vecs_shape[0]
    n_vecs = vecs_shape[1]
    IntegerVector projections_shape = projections.attr("dim");
    n_selected_axes = projections_shape[0]
    n_selected_vecs = projections_shape[1]

    const int* vecs_selection_mask_p = vecs_selection_mask.begin();
    // Is already 0/1 mapped, thus check only for NA values and return ERR_NAN_INF code
    for (int i = 0; i < vecs_selection_mask.size(); ++i) {
        if (vecs_selection_mask_p[i] == NA_LOGICAL) {
            return List::create(Named("ierr") = 204);
        }
    }
    const int* axes_selection_mask_p = axes_selection_mask.begin();
    // Is already 0/1 mapped, thus check only for NA values and return ERR_NAN_INF code
    for (int i = 0; i < axes_selection_mask.size(); ++i) {
        if (axes_selection_mask_p[i] == NA_LOGICAL) {
            return List::create(Named("ierr") = 204);
        }
    }

    NumericMatrix projections(n_selected_axes * n_selected_vecs);
    int ierr = 0;

    omics_vector_RAP_projection_c(
        vecs.begin(),
        &n_axes,
        &n_vecs,
        vecs_selection_mask_p,
        &n_selected_vecs,
        axes_selection_mask_p,
        &n_selected_axes,
        projections.begin(),
        &ierr
    );

//{format(self.arguments, "type_conversion_outputs")}

    return List::create(
        Named("projections") = projections,
        Named("ierr") = ierr
    );
}

List omics_field_RAP_projection_rcpp(
    const NumericMatrix& vecs,
    const LogicalVector& vecs_selection_mask,
    const LogicalVector& axes_selection_mask
) {


    IntegerVector vecs_shape = vecs.attr("dim");
    n_vecs = vecs_shape[1]
    IntegerVector axes_selection_mask_shape = axes_selection_mask.attr("dim");
    n_axes = axes_selection_mask_shape[0]
    IntegerVector projections_shape = projections.attr("dim");
    n_selected_axes = projections_shape[0]
    n_selected_vecs = projections_shape[1]

    const int* vecs_selection_mask_p = vecs_selection_mask.begin();
    // Is already 0/1 mapped, thus check only for NA values and return ERR_NAN_INF code
    for (int i = 0; i < vecs_selection_mask.size(); ++i) {
        if (vecs_selection_mask_p[i] == NA_LOGICAL) {
            return List::create(Named("ierr") = 204);
        }
    }
    const int* axes_selection_mask_p = axes_selection_mask.begin();
    // Is already 0/1 mapped, thus check only for NA values and return ERR_NAN_INF code
    for (int i = 0; i < axes_selection_mask.size(); ++i) {
        if (axes_selection_mask_p[i] == NA_LOGICAL) {
            return List::create(Named("ierr") = 204);
        }
    }

    NumericMatrix projections(n_selected_axes * n_selected_vecs);
    int ierr = 0;

    omics_field_RAP_projection_c(
        vecs.begin(),
        &n_axes,
        &n_vecs,
        vecs_selection_mask_p,
        &n_selected_vecs,
        axes_selection_mask_p,
        &n_selected_axes,
        projections.begin(),
        &ierr
    );

//{format(self.arguments, "type_conversion_outputs")}

    return List::create(
        Named("projections") = projections,
        Named("ierr") = ierr
    );
}

List clock_hand_angle_between_vectors_rcpp(
    const NumericVector& v1,
    const NumericVector& v2,
    const IntegerVector& selected_axes_for_signed
) {


    IntegerVector v1_shape = v1.attr("dim");
    n_dims = v1_shape[0]



    double signed_angle = 0.0;
    int ierr = 0;

    clock_hand_angle_between_vectors_c(
        v1.begin(),
        v2.begin(),
        &n_dims,
        &signed_angle,
        selected_axes_for_signed.begin(),
        &ierr
    );

//{format(self.arguments, "type_conversion_outputs")}

    return List::create(
        Named("signed_angle") = signed_angle,
        Named("ierr") = ierr
    );
}

List clock_hand_angles_for_shift_vectors_rcpp(
    const NumericMatrix& origins,
    const NumericMatrix& targets,
    const LogicalVector& vecs_selection_mask,
    const IntegerVector& selected_axes_for_signed
) {


    IntegerVector origins_shape = origins.attr("dim");
    n_dims = origins_shape[0]
    n_vecs = origins_shape[1]
    IntegerVector signed_angles_shape = signed_angles.attr("dim");
    n_selected_vecs = signed_angles_shape[0]

    const int* vecs_selection_mask_p = vecs_selection_mask.begin();
    // Is already 0/1 mapped, thus check only for NA values and return ERR_NAN_INF code
    for (int i = 0; i < vecs_selection_mask.size(); ++i) {
        if (vecs_selection_mask_p[i] == NA_LOGICAL) {
            return List::create(Named("ierr") = 204);
        }
    }

    NumericVector signed_angles(n_selected_vecs);
    int ierr = 0;

    clock_hand_angles_for_shift_vectors_c(
        origins.begin(),
        targets.begin(),
        &n_dims,
        &n_vecs,
        vecs_selection_mask_p,
        &n_selected_vecs,
        selected_axes_for_signed.begin(),
        signed_angles.begin(),
        &ierr
    );

//{format(self.arguments, "type_conversion_outputs")}

    return List::create(
        Named("signed_angles") = signed_angles,
        Named("ierr") = ierr
    );
}

List compute_relative_axis_contributions_rcpp(
    const NumericVector& vec
) {


    IntegerVector vec_shape = vec.attr("dim");
    n_axes = vec_shape[0]



    NumericVector contributions(n_axes);
    int ierr = 0;

    compute_relative_axis_contributions_c(
        vec.begin(),
        &n_axes,
        contributions.begin(),
        &ierr
    );

//{format(self.arguments, "type_conversion_outputs")}

    return List::create(
        Named("contributions") = contributions,
        Named("ierr") = ierr
    );
}

List relative_axes_changes_from_shift_vector_rcpp(
    const NumericVector& vec
) {


    IntegerVector vec_shape = vec.attr("dim");
    n_axes = vec_shape[0]



    NumericVector contributions(n_axes);
    int ierr = 0;

    relative_axes_changes_from_shift_vector_c(
        vec.begin(),
        &n_axes,
        contributions.begin(),
        &ierr
    );

//{format(self.arguments, "type_conversion_outputs")}

    return List::create(
        Named("contributions") = contributions,
        Named("ierr") = ierr
    );
}

List relative_axes_expression_from_expression_vector_rcpp(
    const NumericVector& vec
) {


    IntegerVector vec_shape = vec.attr("dim");
    n_axes = vec_shape[0]



    NumericVector contributions(n_axes);
    int ierr = 0;

    relative_axes_expression_from_expression_vector_c(
        vec.begin(),
        &n_axes,
        contributions.begin(),
        &ierr
    );

//{format(self.arguments, "type_conversion_outputs")}

    return List::create(
        Named("contributions") = contributions,
        Named("ierr") = ierr
    );
}