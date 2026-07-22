// Generated. Do not edit.
#include <Rcpp.h>
#include <numeric>
#include <functional>
#include "tox_marshal.h"

using namespace Rcpp;

extern "C" {
    void omics_vector_RAP_projection_c(const double*, const int*, const int*, const bool*, const int*, const bool*, const int*, double*, int*);
    void omics_field_RAP_projection_c(const double*, const int*, const int*, const bool*, const int*, const bool*, const int*, double*, int*);
    void clock_hand_angle_between_vectors_c(const double*, const double*, const int*, double*, const int*, int*);
    void clock_hand_angles_for_shift_vectors_c(const double*, const int*, const int*, const bool*, const int*, const int*, double*, int*);
    void relative_axes_changes_from_shift_vector_c(const double*, const int*, double*, int*);
    void relative_axes_expression_from_expression_vector_c(const double*, const int*, double*, int*);
}

// [[Rcpp::export(.omics_vector_RAP_projection_rcpp)]]
List omics_vector_RAP_projection_rcpp(NumericVector vecs, LogicalVector vecs_selection_mask, LogicalVector axes_selection_mask) {
    // derived from the inputs, not asked of the caller
    int n_axes = (int) IntegerVector(vecs.attr("dim"))[0];
    int n_vecs = (int) IntegerVector(vecs.attr("dim"))[1];
    int n_selected_vecs = (int) sum(vecs_selection_mask);
    int n_selected_axes = (int) sum(axes_selection_mask);

    // convert what C cannot take directly
    tox::BoolBuffer vecs_selection_mask_c(vecs_selection_mask);
    tox::BoolBuffer axes_selection_mask_c(axes_selection_mask);

    // outputs and work space
    NumericVector projections(n_selected_axes * n_selected_vecs);
    projections.attr("dim") = IntegerVector::create(n_selected_axes, n_selected_vecs);
    int ierr = 0;

    omics_vector_RAP_projection_c(
        vecs.begin(),
        &n_axes,
        &n_vecs,
        vecs_selection_mask_c.data(),
        &n_selected_vecs,
        axes_selection_mask_c.data(),
        &n_selected_axes,
        projections.begin(),
        &ierr
    );

    return List::create(
        _["projections"] = projections,
        _["ierr"] = ierr
    );
}

// [[Rcpp::export(.omics_field_RAP_projection_rcpp)]]
List omics_field_RAP_projection_rcpp(NumericVector fields, LogicalVector fields_selection_mask, LogicalVector axes_selection_mask) {
    // derived from the inputs, not asked of the caller
    int n_axes = (int) IntegerVector(fields.attr("dim"))[0];
    int n_fields = (int) IntegerVector(fields.attr("dim"))[2];
    int n_selected_fields = (int) sum(fields_selection_mask);
    int n_selected_axes = (int) sum(axes_selection_mask);

    // convert what C cannot take directly
    tox::BoolBuffer fields_selection_mask_c(fields_selection_mask);
    tox::BoolBuffer axes_selection_mask_c(axes_selection_mask);

    // outputs and work space
    NumericVector projections(n_selected_axes * n_selected_fields);
    projections.attr("dim") = IntegerVector::create(n_selected_axes, n_selected_fields);
    int ierr = 0;

    omics_field_RAP_projection_c(
        fields.begin(),
        &n_axes,
        &n_fields,
        fields_selection_mask_c.data(),
        &n_selected_fields,
        axes_selection_mask_c.data(),
        &n_selected_axes,
        projections.begin(),
        &ierr
    );

    return List::create(
        _["projections"] = projections,
        _["ierr"] = ierr
    );
}

// [[Rcpp::export(.clock_hand_angle_between_vectors_rcpp)]]
List clock_hand_angle_between_vectors_rcpp(NumericVector v1, NumericVector v2, IntegerVector selected_axes_for_signed) {
    // derived from the inputs, not asked of the caller
    int n_dims = (int) v1.size();

    // outputs and work space
    double signed_angle = 0;
    int ierr = 0;

    clock_hand_angle_between_vectors_c(
        v1.begin(),
        v2.begin(),
        &n_dims,
        &signed_angle,
        selected_axes_for_signed.begin(),
        &ierr
    );

    return List::create(
        _["signed_angle"] = signed_angle,
        _["ierr"] = ierr
    );
}

// [[Rcpp::export(.clock_hand_angles_for_shift_vectors_rcpp)]]
List clock_hand_angles_for_shift_vectors_rcpp(NumericVector fields, LogicalVector fields_selection_mask, IntegerVector selected_axes_for_signed) {
    // derived from the inputs, not asked of the caller
    int n_dims = (int) IntegerVector(fields.attr("dim"))[0];
    int n_fields = (int) IntegerVector(fields.attr("dim"))[2];
    int n_selected_fields = (int) sum(fields_selection_mask);

    // convert what C cannot take directly
    tox::BoolBuffer fields_selection_mask_c(fields_selection_mask);

    // outputs and work space
    NumericVector signed_angles(n_selected_fields);
    int ierr = 0;

    clock_hand_angles_for_shift_vectors_c(
        fields.begin(),
        &n_dims,
        &n_fields,
        fields_selection_mask_c.data(),
        &n_selected_fields,
        selected_axes_for_signed.begin(),
        signed_angles.begin(),
        &ierr
    );

    return List::create(
        _["signed_angles"] = signed_angles,
        _["ierr"] = ierr
    );
}

// [[Rcpp::export(.relative_axes_changes_from_shift_vector_rcpp)]]
List relative_axes_changes_from_shift_vector_rcpp(NumericVector vec) {
    // derived from the inputs, not asked of the caller
    int n_axes = (int) vec.size();

    // outputs and work space
    NumericVector contributions(n_axes);
    int ierr = 0;

    relative_axes_changes_from_shift_vector_c(
        vec.begin(),
        &n_axes,
        contributions.begin(),
        &ierr
    );

    return List::create(
        _["contributions"] = contributions,
        _["ierr"] = ierr
    );
}

// [[Rcpp::export(.relative_axes_expression_from_expression_vector_rcpp)]]
List relative_axes_expression_from_expression_vector_rcpp(NumericVector vec) {
    // derived from the inputs, not asked of the caller
    int n_axes = (int) vec.size();

    // outputs and work space
    NumericVector contributions(n_axes);
    int ierr = 0;

    relative_axes_expression_from_expression_vector_c(
        vec.begin(),
        &n_axes,
        contributions.begin(),
        &ierr
    );

    return List::create(
        _["contributions"] = contributions,
        _["ierr"] = ierr
    );
}
