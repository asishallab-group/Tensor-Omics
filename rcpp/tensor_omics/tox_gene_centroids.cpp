#include <Rcpp.h>

using namespace Rcpp;

static_assert(sizeof(Rcomplex) == sizeof(double _Complex) && alignof(Rcomplex) == alignof(double _Complex), "Rcomplex layout is incompatible with C's 'double _Complex' and thus Fortran's 'c_double_complex'");

extern "C" {
    void mean_vector_c(
        const double* expression_vectors,
        const int* n_axes,
        const int* n_genes,
        const int* gene_indices,
        const int* n_selected_genes,
        double* centroid,
        int* ierr
    );

    void group_centroid_c(
        const double* expression_vectors,
        const int* n_axes,
        const int* n_genes,
        const int* gene_to_family,
        const int* n_families,
        double* centroid_matrix,
        const char* mode,
        int* tmp_selected_indices,
        int* ierr,
        const int* ortholog_set
    );
}

List mean_vector_rcpp(
    const NumericMatrix& expression_vectors,
    const IntegerVector& gene_indices
) {


    IntegerVector expression_vectors_shape = expression_vectors.attr("dim");
    n_axes = expression_vectors_shape[0]
    n_genes = expression_vectors_shape[1]
    IntegerVector gene_indices_shape = gene_indices.attr("dim");
    n_selected_genes = gene_indices_shape[0]



    NumericVector centroid(n_axes);
    int ierr = 0;

    mean_vector_c(
        expression_vectors.begin(),
        &n_axes,
        &n_genes,
        gene_indices.begin(),
        &n_selected_genes,
        centroid.begin(),
        &ierr
    );

//{format(self.arguments, "type_conversion_outputs")}

    return List::create(
        Named("centroid") = centroid,
        Named("ierr") = ierr
    );
}

List group_centroid_rcpp(
    const NumericMatrix& expression_vectors,
    const IntegerVector& gene_to_family,
    const String& mode,
    const LogicalVector& ortholog_set = 
) {


    IntegerVector expression_vectors_shape = expression_vectors.attr("dim");
    n_axes = expression_vectors_shape[0]
    n_genes = expression_vectors_shape[1]
    IntegerVector centroid_matrix_shape = centroid_matrix.attr("dim");
    n_families = centroid_matrix_shape[1]

    if (mode == NA_STRING) {
        return List::create(Named("ierr") = 204);
    }
    const char* mode_p = mode.get_cstring();
    const int* ortholog_set_p = ortholog_set.begin();
    // Is already 0/1 mapped, thus check only for NA values and return ERR_NAN_INF code
    for (int i = 0; i < ortholog_set.size(); ++i) {
        if (ortholog_set_p[i] == NA_LOGICAL) {
            return List::create(Named("ierr") = 204);
        }
    }

    NumericMatrix centroid_matrix(n_axes * n_families);
    std::vector<int> tmp_selected_indices(n_genes);
    int ierr = 0;

    group_centroid_c(
        expression_vectors.begin(),
        &n_axes,
        &n_genes,
        gene_to_family.begin(),
        &n_families,
        centroid_matrix.begin(),
        mode_p,
        tmp_selected_indices.data(),
        &ierr,
        ortholog_set_p
    );

//{format(self.arguments, "type_conversion_outputs")}

    return List::create(
        Named("centroid_matrix") = centroid_matrix,
        Named("ierr") = ierr
    );
}