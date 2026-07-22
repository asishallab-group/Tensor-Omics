// Generated. Do not edit.
#include <Rcpp.h>
#include <numeric>
#include <functional>
#include "tox_marshal.h"

using namespace Rcpp;

extern "C" {
    void create_zip_archive_c(const char*, const int*, const char*, const int*, const int*, const char*, const int*, const int*, int*);
    void save_tox_data_c(const char*, const int*, int*, const char*, const int*, const int*, const char*, const int*, const double*, const int*, const int*, const char*, const int*, const int*, const int*, const char*, const int*, const char*, const int*, const int*, const char*, const int*, const double*, const int*, const int*, const char*, const int*, const double*, const int*, const int*, const char*, const int*);
    void get_tox_data_dims_c(const char*, const int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*);
    void read_tox_data_into_c(const char*, const int*, const int*, const int*, char*, const int*, const int*, double*, const int*, int*, const int*, const int*, char*, const int*, const int*, double*, const int*, const int*, double*, int*);
}

// [[Rcpp::export(.create_zip_archive_rcpp)]]
List create_zip_archive_rcpp(CharacterVector zip_filename, CharacterVector keys, CharacterVector filenames) {
    // derived from the inputs, not asked of the caller
    int zip_filename_strlen = tox::max_strlen(zip_filename);
    int keys_strlen = tox::max_strlen(keys);
    int n_keys_elements = (int) keys.size();
    int filenames_strlen = tox::max_strlen(filenames);
    int n_filenames_elements = (int) filenames.size();

    // convert what C cannot take directly
    tox::CharBuffer zip_filename_c(zip_filename, zip_filename_strlen);
    tox::CharBuffer keys_c(keys, keys_strlen);
    tox::CharBuffer filenames_c(filenames, filenames_strlen);

    // outputs and work space
    int ierr = 0;

    create_zip_archive_c(
        zip_filename_c.data(),
        &zip_filename_strlen,
        keys_c.data(),
        &keys_strlen,
        &n_keys_elements,
        filenames_c.data(),
        &filenames_strlen,
        &n_filenames_elements,
        &ierr
    );

    return List::create(
        _["ierr"] = ierr
    );
}

// [[Rcpp::export(.save_tox_data_rcpp)]]
List save_tox_data_rcpp(CharacterVector zip_filename, Nullable<CharacterVector> gene_ids = R_NilValue, Nullable<CharacterVector> gene_ids_file = R_NilValue, Nullable<NumericVector> expression = R_NilValue, Nullable<CharacterVector> expression_file = R_NilValue, Nullable<IntegerVector> gene_to_family = R_NilValue, Nullable<CharacterVector> gene_to_family_file = R_NilValue, Nullable<CharacterVector> family_ids = R_NilValue, Nullable<CharacterVector> family_ids_file = R_NilValue, Nullable<NumericVector> family_centroids = R_NilValue, Nullable<CharacterVector> family_centroids_file = R_NilValue, Nullable<NumericVector> shift_vectors = R_NilValue, Nullable<CharacterVector> shift_vectors_file = R_NilValue) {
    // optionals: a null pointer and size 0 when the caller omits them
    int gene_ids_size = 0;
    CharacterVector gene_ids_val;
    if (gene_ids.isNotNull()) {
        gene_ids_val = gene_ids.get();
        gene_ids_size = gene_ids_val.size();
    }
    int gene_ids_file_size = 0;
    CharacterVector gene_ids_file_val;
    if (gene_ids_file.isNotNull()) {
        gene_ids_file_val = gene_ids_file.get();
        gene_ids_file_size = gene_ids_file_val.size();
    }
    const double* expression_p = nullptr;
    int expression_size = 0;
    std::vector<int> expression_dim(2, 0);
    NumericVector expression_val;
    if (expression.isNotNull()) {
        expression_val = expression.get();
        expression_p = expression_val.begin();
        expression_size = expression_val.size();
        if (!Rf_isNull(expression_val.attr("dim"))) {
            IntegerVector expression_d(expression_val.attr("dim"));
            for (int i = 0; i < 2 && i < expression_d.size(); ++i) expression_dim[i] = expression_d[i];
        }
    }
    int expression_file_size = 0;
    CharacterVector expression_file_val;
    if (expression_file.isNotNull()) {
        expression_file_val = expression_file.get();
        expression_file_size = expression_file_val.size();
    }
    const int* gene_to_family_p = nullptr;
    int gene_to_family_size = 0;
    IntegerVector gene_to_family_val;
    if (gene_to_family.isNotNull()) {
        gene_to_family_val = gene_to_family.get();
        gene_to_family_p = gene_to_family_val.begin();
        gene_to_family_size = gene_to_family_val.size();
    }
    int gene_to_family_file_size = 0;
    CharacterVector gene_to_family_file_val;
    if (gene_to_family_file.isNotNull()) {
        gene_to_family_file_val = gene_to_family_file.get();
        gene_to_family_file_size = gene_to_family_file_val.size();
    }
    int family_ids_size = 0;
    CharacterVector family_ids_val;
    if (family_ids.isNotNull()) {
        family_ids_val = family_ids.get();
        family_ids_size = family_ids_val.size();
    }
    int family_ids_file_size = 0;
    CharacterVector family_ids_file_val;
    if (family_ids_file.isNotNull()) {
        family_ids_file_val = family_ids_file.get();
        family_ids_file_size = family_ids_file_val.size();
    }
    const double* family_centroids_p = nullptr;
    int family_centroids_size = 0;
    std::vector<int> family_centroids_dim(2, 0);
    NumericVector family_centroids_val;
    if (family_centroids.isNotNull()) {
        family_centroids_val = family_centroids.get();
        family_centroids_p = family_centroids_val.begin();
        family_centroids_size = family_centroids_val.size();
        if (!Rf_isNull(family_centroids_val.attr("dim"))) {
            IntegerVector family_centroids_d(family_centroids_val.attr("dim"));
            for (int i = 0; i < 2 && i < family_centroids_d.size(); ++i) family_centroids_dim[i] = family_centroids_d[i];
        }
    }
    int family_centroids_file_size = 0;
    CharacterVector family_centroids_file_val;
    if (family_centroids_file.isNotNull()) {
        family_centroids_file_val = family_centroids_file.get();
        family_centroids_file_size = family_centroids_file_val.size();
    }
    const double* shift_vectors_p = nullptr;
    int shift_vectors_size = 0;
    std::vector<int> shift_vectors_dim(2, 0);
    NumericVector shift_vectors_val;
    if (shift_vectors.isNotNull()) {
        shift_vectors_val = shift_vectors.get();
        shift_vectors_p = shift_vectors_val.begin();
        shift_vectors_size = shift_vectors_val.size();
        if (!Rf_isNull(shift_vectors_val.attr("dim"))) {
            IntegerVector shift_vectors_d(shift_vectors_val.attr("dim"));
            for (int i = 0; i < 2 && i < shift_vectors_d.size(); ++i) shift_vectors_dim[i] = shift_vectors_d[i];
        }
    }
    int shift_vectors_file_size = 0;
    CharacterVector shift_vectors_file_val;
    if (shift_vectors_file.isNotNull()) {
        shift_vectors_file_val = shift_vectors_file.get();
        shift_vectors_file_size = shift_vectors_file_val.size();
    }

    // derived from the inputs, not asked of the caller
    int zip_filename_strlen = tox::max_strlen(zip_filename);
    int gene_ids_strlen = tox::max_strlen(gene_ids_val);
    int n_gene_ids_elements = gene_ids_size;
    int gene_ids_file_strlen = tox::max_strlen(gene_ids_file_val);
    int n_expression_elements_dim_1 = expression_dim[0];
    int n_expression_elements_dim_2 = expression_dim[1];
    int expression_file_strlen = tox::max_strlen(expression_file_val);
    int n_gene_to_family_elements = gene_to_family_size;
    int gene_to_family_file_strlen = tox::max_strlen(gene_to_family_file_val);
    int family_ids_strlen = tox::max_strlen(family_ids_val);
    int n_family_ids_elements = family_ids_size;
    int family_ids_file_strlen = tox::max_strlen(family_ids_file_val);
    int n_family_centroids_elements_dim_1 = family_centroids_dim[0];
    int n_family_centroids_elements_dim_2 = family_centroids_dim[1];
    int family_centroids_file_strlen = tox::max_strlen(family_centroids_file_val);
    int n_shift_vectors_elements_dim_1 = shift_vectors_dim[0];
    int n_shift_vectors_elements_dim_2 = shift_vectors_dim[1];
    int shift_vectors_file_strlen = tox::max_strlen(shift_vectors_file_val);

    // convert what C cannot take directly
    tox::CharBuffer zip_filename_c(zip_filename, zip_filename_strlen);
    tox::CharBuffer gene_ids_c(gene_ids_val, gene_ids_strlen);
    tox::CharBuffer gene_ids_file_c(gene_ids_file_val, gene_ids_file_strlen);
    tox::CharBuffer expression_file_c(expression_file_val, expression_file_strlen);
    tox::CharBuffer gene_to_family_file_c(gene_to_family_file_val, gene_to_family_file_strlen);
    tox::CharBuffer family_ids_c(family_ids_val, family_ids_strlen);
    tox::CharBuffer family_ids_file_c(family_ids_file_val, family_ids_file_strlen);
    tox::CharBuffer family_centroids_file_c(family_centroids_file_val, family_centroids_file_strlen);
    tox::CharBuffer shift_vectors_file_c(shift_vectors_file_val, shift_vectors_file_strlen);

    // outputs and work space
    int ierr = 0;

    save_tox_data_c(
        zip_filename_c.data(),
        &zip_filename_strlen,
        &ierr,
        gene_ids.isNotNull() ? gene_ids_c.data() : nullptr,
        &gene_ids_strlen,
        &n_gene_ids_elements,
        gene_ids_file.isNotNull() ? gene_ids_file_c.data() : nullptr,
        &gene_ids_file_strlen,
        expression_p,
        &n_expression_elements_dim_1,
        &n_expression_elements_dim_2,
        expression_file.isNotNull() ? expression_file_c.data() : nullptr,
        &expression_file_strlen,
        gene_to_family_p,
        &n_gene_to_family_elements,
        gene_to_family_file.isNotNull() ? gene_to_family_file_c.data() : nullptr,
        &gene_to_family_file_strlen,
        family_ids.isNotNull() ? family_ids_c.data() : nullptr,
        &family_ids_strlen,
        &n_family_ids_elements,
        family_ids_file.isNotNull() ? family_ids_file_c.data() : nullptr,
        &family_ids_file_strlen,
        family_centroids_p,
        &n_family_centroids_elements_dim_1,
        &n_family_centroids_elements_dim_2,
        family_centroids_file.isNotNull() ? family_centroids_file_c.data() : nullptr,
        &family_centroids_file_strlen,
        shift_vectors_p,
        &n_shift_vectors_elements_dim_1,
        &n_shift_vectors_elements_dim_2,
        shift_vectors_file.isNotNull() ? shift_vectors_file_c.data() : nullptr,
        &shift_vectors_file_strlen
    );

    return List::create(
        _["ierr"] = ierr
    );
}

// [[Rcpp::export(.get_tox_data_dims_rcpp)]]
List get_tox_data_dims_rcpp(CharacterVector zip_filename) {
    // derived from the inputs, not asked of the caller
    int zip_filename_strlen = tox::max_strlen(zip_filename);

    // convert what C cannot take directly
    tox::CharBuffer zip_filename_c(zip_filename, zip_filename_strlen);

    // outputs and work space
    int n_gene_ids = 0;
    int gene_id_len = 0;
    int n_expression_rows = 0;
    int n_expression_cols = 0;
    int n_gene_to_family = 0;
    int n_family_ids = 0;
    int family_id_len = 0;
    int n_family_centroids_rows = 0;
    int n_family_centroids_cols = 0;
    int n_shift_vectors_rows = 0;
    int n_shift_vectors_cols = 0;
    int ierr = 0;

    get_tox_data_dims_c(
        zip_filename_c.data(),
        &zip_filename_strlen,
        &n_gene_ids,
        &gene_id_len,
        &n_expression_rows,
        &n_expression_cols,
        &n_gene_to_family,
        &n_family_ids,
        &family_id_len,
        &n_family_centroids_rows,
        &n_family_centroids_cols,
        &n_shift_vectors_rows,
        &n_shift_vectors_cols,
        &ierr
    );

    return List::create(
        _["n_gene_ids"] = n_gene_ids,
        _["gene_id_len"] = gene_id_len,
        _["n_expression_rows"] = n_expression_rows,
        _["n_expression_cols"] = n_expression_cols,
        _["n_gene_to_family"] = n_gene_to_family,
        _["n_family_ids"] = n_family_ids,
        _["family_id_len"] = family_id_len,
        _["n_family_centroids_rows"] = n_family_centroids_rows,
        _["n_family_centroids_cols"] = n_family_centroids_cols,
        _["n_shift_vectors_rows"] = n_shift_vectors_rows,
        _["n_shift_vectors_cols"] = n_shift_vectors_cols,
        _["ierr"] = ierr
    );
}

// [[Rcpp::export(.read_tox_data_into_rcpp)]]
List read_tox_data_into_rcpp(CharacterVector zip_filename, int n_gene_ids, int gene_id_len, int n_expression_rows, int n_expression_cols, int n_gene_to_family, int n_family_ids, int family_id_len, int n_family_centroids_rows, int n_family_centroids_cols, int n_shift_vectors_rows, int n_shift_vectors_cols) {
    // derived from the inputs, not asked of the caller
    int zip_filename_strlen = tox::max_strlen(zip_filename);

    // convert what C cannot take directly
    tox::CharBuffer zip_filename_c(zip_filename, zip_filename_strlen);

    // outputs and work space
    tox::CharBuffer gene_ids_c(gene_id_len, n_gene_ids);
    NumericVector expression(n_expression_rows * n_expression_cols);
    expression.attr("dim") = IntegerVector::create(n_expression_rows, n_expression_cols);
    IntegerVector gene_to_family(n_gene_to_family);
    tox::CharBuffer family_ids_c(family_id_len, n_family_ids);
    NumericVector family_centroids(n_family_centroids_rows * n_family_centroids_cols);
    family_centroids.attr("dim") = IntegerVector::create(n_family_centroids_rows, n_family_centroids_cols);
    NumericVector shift_vectors(n_shift_vectors_rows * n_shift_vectors_cols);
    shift_vectors.attr("dim") = IntegerVector::create(n_shift_vectors_rows, n_shift_vectors_cols);
    int ierr = 0;

    read_tox_data_into_c(
        zip_filename_c.data(),
        &zip_filename_strlen,
        &n_gene_ids,
        &gene_id_len,
        gene_ids_c.data(),
        &n_expression_rows,
        &n_expression_cols,
        expression.begin(),
        &n_gene_to_family,
        gene_to_family.begin(),
        &n_family_ids,
        &family_id_len,
        family_ids_c.data(),
        &n_family_centroids_rows,
        &n_family_centroids_cols,
        family_centroids.begin(),
        &n_shift_vectors_rows,
        &n_shift_vectors_cols,
        shift_vectors.begin(),
        &ierr
    );

    // convert the outputs back
    CharacterVector gene_ids = gene_ids_c.to_r();
    CharacterVector family_ids = family_ids_c.to_r();

    return List::create(
        _["gene_ids"] = gene_ids,
        _["expression"] = expression,
        _["gene_to_family"] = gene_to_family,
        _["family_ids"] = family_ids,
        _["family_centroids"] = family_centroids,
        _["shift_vectors"] = shift_vectors,
        _["ierr"] = ierr
    );
}
