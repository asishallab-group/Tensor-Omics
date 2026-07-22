// Generated. Do not edit.
#include <Rcpp.h>
#include <numeric>
#include <functional>
#include "tox_marshal.h"

using namespace Rcpp;

extern "C" {
    void validate_data_structure_c(const int*, const int*, const int*, const char*, const int*, const int*, const char*, const int*, const int*, const int*, const int*, const double*, const int*, const int*, const double*, const int*, const int*, const double*, const int*, const int*, int*);
    void validate_gene_to_family_mapping_c(const int*, const int*, const int*, int*);
    void validate_expression_data_c(const double*, const int*, const int*, const bool*, int*);
    void validate_family_centroids_c(const double*, const int*, const int*, int*);
    void validate_shift_vectors_c(const double*, const int*, const int*, const double*, const int*, const int*, const double*, const int*, const int*, const int*, const int*, const int*, int*);
    void validate_string_array_uniqueness_c(const char*, const int*, const int*, int*);
    void validate_all_data_c(const int*, const int*, const int*, const char*, const int*, const int*, const char*, const int*, const int*, const int*, const int*, const double*, const int*, const int*, const double*, const int*, const int*, const double*, const int*, const int*, int*, const bool*, const bool*);
}

// [[Rcpp::export(.validate_data_structure_rcpp)]]
List validate_data_structure_rcpp(int n_genes, int n_families, int n_samples, CharacterVector gene_ids, CharacterVector gene_family_ids, IntegerVector gene_to_fam, NumericVector expression_vectors, NumericVector family_centroids, NumericVector shift_vectors) {
    // derived from the inputs, not asked of the caller
    int gene_ids_strlen = tox::max_strlen(gene_ids);
    int n_gene_ids_elements = (int) gene_ids.size();
    int gene_family_ids_strlen = tox::max_strlen(gene_family_ids);
    int n_gene_family_ids_elements = (int) gene_family_ids.size();
    int n_gene_to_fam_elements = (int) gene_to_fam.size();
    int n_expression_vectors_elements_dim_1 = (int) IntegerVector(expression_vectors.attr("dim"))[0];
    int n_expression_vectors_elements_dim_2 = (int) IntegerVector(expression_vectors.attr("dim"))[1];
    int n_family_centroids_elements_dim_1 = (int) IntegerVector(family_centroids.attr("dim"))[0];
    int n_family_centroids_elements_dim_2 = (int) IntegerVector(family_centroids.attr("dim"))[1];
    int n_shift_vectors_elements_dim_1 = (int) IntegerVector(shift_vectors.attr("dim"))[0];
    int n_shift_vectors_elements_dim_2 = (int) IntegerVector(shift_vectors.attr("dim"))[1];

    // convert what C cannot take directly
    tox::CharBuffer gene_ids_c(gene_ids, gene_ids_strlen);
    tox::CharBuffer gene_family_ids_c(gene_family_ids, gene_family_ids_strlen);

    // outputs and work space
    int ierr = 0;

    validate_data_structure_c(
        &n_genes,
        &n_families,
        &n_samples,
        gene_ids_c.data(),
        &gene_ids_strlen,
        &n_gene_ids_elements,
        gene_family_ids_c.data(),
        &gene_family_ids_strlen,
        &n_gene_family_ids_elements,
        gene_to_fam.begin(),
        &n_gene_to_fam_elements,
        expression_vectors.begin(),
        &n_expression_vectors_elements_dim_1,
        &n_expression_vectors_elements_dim_2,
        family_centroids.begin(),
        &n_family_centroids_elements_dim_1,
        &n_family_centroids_elements_dim_2,
        shift_vectors.begin(),
        &n_shift_vectors_elements_dim_1,
        &n_shift_vectors_elements_dim_2,
        &ierr
    );

    return List::create(
        _["ierr"] = ierr
    );
}

// [[Rcpp::export(.validate_gene_to_family_mapping_rcpp)]]
List validate_gene_to_family_mapping_rcpp(IntegerVector gene_to_fam, int n_families) {
    // derived from the inputs, not asked of the caller
    int n_gene_to_fam_elements = (int) gene_to_fam.size();

    // outputs and work space
    int ierr = 0;

    validate_gene_to_family_mapping_c(
        gene_to_fam.begin(),
        &n_gene_to_fam_elements,
        &n_families,
        &ierr
    );

    return List::create(
        _["ierr"] = ierr
    );
}

// [[Rcpp::export(.validate_expression_data_rcpp)]]
List validate_expression_data_rcpp(NumericVector expression_vectors, bool check_non_negative) {
    // derived from the inputs, not asked of the caller
    int n_expression_vectors_elements_dim_1 = (int) IntegerVector(expression_vectors.attr("dim"))[0];
    int n_expression_vectors_elements_dim_2 = (int) IntegerVector(expression_vectors.attr("dim"))[1];

    // outputs and work space
    int ierr = 0;

    validate_expression_data_c(
        expression_vectors.begin(),
        &n_expression_vectors_elements_dim_1,
        &n_expression_vectors_elements_dim_2,
        &check_non_negative,
        &ierr
    );

    return List::create(
        _["ierr"] = ierr
    );
}

// [[Rcpp::export(.validate_family_centroids_rcpp)]]
List validate_family_centroids_rcpp(NumericVector family_centroids) {
    // derived from the inputs, not asked of the caller
    int n_family_centroids_elements_dim_1 = (int) IntegerVector(family_centroids.attr("dim"))[0];
    int n_family_centroids_elements_dim_2 = (int) IntegerVector(family_centroids.attr("dim"))[1];

    // outputs and work space
    int ierr = 0;

    validate_family_centroids_c(
        family_centroids.begin(),
        &n_family_centroids_elements_dim_1,
        &n_family_centroids_elements_dim_2,
        &ierr
    );

    return List::create(
        _["ierr"] = ierr
    );
}

// [[Rcpp::export(.validate_shift_vectors_rcpp)]]
List validate_shift_vectors_rcpp(NumericVector shift_vectors, NumericVector expression_vectors, NumericVector family_centroids, IntegerVector gene_to_fam, int n_samples) {
    // derived from the inputs, not asked of the caller
    int n_shift_vectors_elements_dim_1 = (int) IntegerVector(shift_vectors.attr("dim"))[0];
    int n_shift_vectors_elements_dim_2 = (int) IntegerVector(shift_vectors.attr("dim"))[1];
    int n_expression_vectors_elements_dim_1 = (int) IntegerVector(expression_vectors.attr("dim"))[0];
    int n_expression_vectors_elements_dim_2 = (int) IntegerVector(expression_vectors.attr("dim"))[1];
    int n_family_centroids_elements_dim_1 = (int) IntegerVector(family_centroids.attr("dim"))[0];
    int n_family_centroids_elements_dim_2 = (int) IntegerVector(family_centroids.attr("dim"))[1];
    int n_gene_to_fam_elements = (int) gene_to_fam.size();

    // outputs and work space
    int ierr = 0;

    validate_shift_vectors_c(
        shift_vectors.begin(),
        &n_shift_vectors_elements_dim_1,
        &n_shift_vectors_elements_dim_2,
        expression_vectors.begin(),
        &n_expression_vectors_elements_dim_1,
        &n_expression_vectors_elements_dim_2,
        family_centroids.begin(),
        &n_family_centroids_elements_dim_1,
        &n_family_centroids_elements_dim_2,
        gene_to_fam.begin(),
        &n_gene_to_fam_elements,
        &n_samples,
        &ierr
    );

    return List::create(
        _["ierr"] = ierr
    );
}

// [[Rcpp::export(.validate_string_array_uniqueness_rcpp)]]
List validate_string_array_uniqueness_rcpp(CharacterVector str_arr) {
    // derived from the inputs, not asked of the caller
    int str_arr_strlen = tox::max_strlen(str_arr);
    int n_str_arr_elements = (int) str_arr.size();

    // convert what C cannot take directly
    tox::CharBuffer str_arr_c(str_arr, str_arr_strlen);

    // outputs and work space
    int ierr = 0;

    validate_string_array_uniqueness_c(
        str_arr_c.data(),
        &str_arr_strlen,
        &n_str_arr_elements,
        &ierr
    );

    return List::create(
        _["ierr"] = ierr
    );
}

// [[Rcpp::export(.validate_all_data_rcpp)]]
List validate_all_data_rcpp(int n_genes, int n_families, int n_samples, CharacterVector gene_ids, CharacterVector gene_family_ids, IntegerVector gene_to_fam, NumericVector expression_vectors, NumericVector family_centroids, NumericVector shift_vectors, bool check_uniqueness, bool check_shift_consistency) {
    // derived from the inputs, not asked of the caller
    int gene_ids_strlen = tox::max_strlen(gene_ids);
    int n_gene_ids_elements = (int) gene_ids.size();
    int gene_family_ids_strlen = tox::max_strlen(gene_family_ids);
    int n_gene_family_ids_elements = (int) gene_family_ids.size();
    int n_gene_to_fam_elements = (int) gene_to_fam.size();
    int n_expression_vectors_elements_dim_1 = (int) IntegerVector(expression_vectors.attr("dim"))[0];
    int n_expression_vectors_elements_dim_2 = (int) IntegerVector(expression_vectors.attr("dim"))[1];
    int n_family_centroids_elements_dim_1 = (int) IntegerVector(family_centroids.attr("dim"))[0];
    int n_family_centroids_elements_dim_2 = (int) IntegerVector(family_centroids.attr("dim"))[1];
    int n_shift_vectors_elements_dim_1 = (int) IntegerVector(shift_vectors.attr("dim"))[0];
    int n_shift_vectors_elements_dim_2 = (int) IntegerVector(shift_vectors.attr("dim"))[1];

    // convert what C cannot take directly
    tox::CharBuffer gene_ids_c(gene_ids, gene_ids_strlen);
    tox::CharBuffer gene_family_ids_c(gene_family_ids, gene_family_ids_strlen);

    // outputs and work space
    int ierr = 0;

    validate_all_data_c(
        &n_genes,
        &n_families,
        &n_samples,
        gene_ids_c.data(),
        &gene_ids_strlen,
        &n_gene_ids_elements,
        gene_family_ids_c.data(),
        &gene_family_ids_strlen,
        &n_gene_family_ids_elements,
        gene_to_fam.begin(),
        &n_gene_to_fam_elements,
        expression_vectors.begin(),
        &n_expression_vectors_elements_dim_1,
        &n_expression_vectors_elements_dim_2,
        family_centroids.begin(),
        &n_family_centroids_elements_dim_1,
        &n_family_centroids_elements_dim_2,
        shift_vectors.begin(),
        &n_shift_vectors_elements_dim_1,
        &n_shift_vectors_elements_dim_2,
        &ierr,
        &check_uniqueness,
        &check_shift_consistency
    );

    return List::create(
        _["ierr"] = ierr
    );
}
