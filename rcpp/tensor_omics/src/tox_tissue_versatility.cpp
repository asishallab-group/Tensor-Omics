// Generated. Do not edit.
#include <Rcpp.h>
#include <numeric>
#include <functional>
#include "tox_marshal.h"

using namespace Rcpp;

extern "C" {
    void compute_tissue_versatility_c(const int*, const int*, const double*, const bool*, const int*, const bool*, const int*, double*, double*, int*);
}

// [[Rcpp::export(.compute_tissue_versatility_rcpp)]]
List compute_tissue_versatility_rcpp(NumericVector expression_vectors, LogicalVector exp_vecs_selection_index, int n_selected_vectors, LogicalVector axes_selection, int n_selected_axes) {
    // derived from the inputs, not asked of the caller
    int n_axes = (int) IntegerVector(expression_vectors.attr("dim"))[0];
    int n_vectors = (int) IntegerVector(expression_vectors.attr("dim"))[1];

    // convert what C cannot take directly
    tox::BoolBuffer exp_vecs_selection_index_c(exp_vecs_selection_index);
    tox::BoolBuffer axes_selection_c(axes_selection);

    // outputs and work space
    NumericVector tissue_versatilities(n_selected_vectors);
    NumericVector tissue_angles_deg(n_selected_vectors);
    int ierr = 0;

    compute_tissue_versatility_c(
        &n_axes,
        &n_vectors,
        expression_vectors.begin(),
        exp_vecs_selection_index_c.data(),
        &n_selected_vectors,
        axes_selection_c.data(),
        &n_selected_axes,
        tissue_versatilities.begin(),
        tissue_angles_deg.begin(),
        &ierr
    );

    return List::create(
        _["tissue_versatilities"] = tissue_versatilities,
        _["tissue_angles_deg"] = tissue_angles_deg,
        _["ierr"] = ierr
    );
}
