// Generated. Do not edit.
#include <Rcpp.h>
#include <numeric>
#include <functional>
#include "tox_marshal.h"

using namespace Rcpp;

extern "C" {
    void read_expression_vectors_tsv_c(const char*, const int*, const int*, const char*, const int*, const int*, double*, const int*, const int*, const int*, const int*, const int*, const int*, const int*, int*, const char*);
    void read_gene_ids_from_tsv_file_c(const char*, const int*, char*, const int*, const int*, const int*, const int*, int*);
    void read_orthofinder_file_c(const char*, const int*, const char*, const int*, const int*, char*, const int*, const int*, int*, const int*, int*);
    void get_unassigned_mask_c(const int*, const int*, bool*, int*, int*);
}

// [[Rcpp::export(.read_expression_vectors_tsv_rcpp)]]
List read_expression_vectors_tsv_rcpp(CharacterVector file_list, CharacterVector gene_ids, NumericVector expression_vectors, int n_header_rows, int gene_col, IntegerVector value_cols, int start_row, Nullable<CharacterVector> delimiter = R_NilValue) {
    // optionals: a null pointer and size 0 when the caller omits them
    int delimiter_size = 0;
    CharacterVector delimiter_val;
    if (delimiter.isNotNull()) {
        delimiter_val = delimiter.get();
        delimiter_size = delimiter_val.size();
    }

    // derived from the inputs, not asked of the caller
    int file_list_strlen = tox::max_strlen(file_list);
    int n_file_list_elements = (int) file_list.size();
    int gene_ids_strlen = tox::max_strlen(gene_ids);
    int n_gene_ids_elements = (int) gene_ids.size();
    int n_expression_vectors_elements_dim_1 = (int) IntegerVector(expression_vectors.attr("dim"))[0];
    int n_expression_vectors_elements_dim_2 = (int) IntegerVector(expression_vectors.attr("dim"))[1];
    int n_value_cols_elements = (int) value_cols.size();

    // copy what is modified in place, so the caller's stays intact
    NumericVector expression_vectors_out = clone(expression_vectors);

    // convert what C cannot take directly
    tox::CharBuffer file_list_c(file_list, file_list_strlen);
    tox::CharBuffer gene_ids_c(gene_ids, gene_ids_strlen);
    tox::CharBuffer delimiter_c(delimiter_val, 1);

    // outputs and work space
    int ierr = 0;

    read_expression_vectors_tsv_c(
        file_list_c.data(),
        &file_list_strlen,
        &n_file_list_elements,
        gene_ids_c.data(),
        &gene_ids_strlen,
        &n_gene_ids_elements,
        expression_vectors_out.begin(),
        &n_expression_vectors_elements_dim_1,
        &n_expression_vectors_elements_dim_2,
        &n_header_rows,
        &gene_col,
        value_cols.begin(),
        &n_value_cols_elements,
        &start_row,
        &ierr,
        delimiter.isNotNull() ? delimiter_c.data() : nullptr
    );

    return List::create(
        _["expression_vectors"] = expression_vectors_out,
        _["ierr"] = ierr
    );
}

// [[Rcpp::export(.read_gene_ids_from_tsv_file_rcpp)]]
List read_gene_ids_from_tsv_file_rcpp(CharacterVector filename, int gene_ids_strlen, int n_gene_ids_elements, int n_header_rows, int gene_col) {
    // derived from the inputs, not asked of the caller
    int filename_strlen = tox::max_strlen(filename);

    // convert what C cannot take directly
    tox::CharBuffer filename_c(filename, filename_strlen);

    // outputs and work space
    tox::CharBuffer gene_ids_c(gene_ids_strlen, n_gene_ids_elements);
    int ierr = 0;

    read_gene_ids_from_tsv_file_c(
        filename_c.data(),
        &filename_strlen,
        gene_ids_c.data(),
        &gene_ids_strlen,
        &n_gene_ids_elements,
        &n_header_rows,
        &gene_col,
        &ierr
    );

    // convert the outputs back
    CharacterVector gene_ids = gene_ids_c.to_r();

    return List::create(
        _["gene_ids"] = gene_ids,
        _["ierr"] = ierr
    );
}

// [[Rcpp::export(.read_orthofinder_file_rcpp)]]
List read_orthofinder_file_rcpp(CharacterVector filename, CharacterVector gene_ids, int family_ids_strlen, int n_family_ids_elements, int n_gene_to_fam_elements) {
    // derived from the inputs, not asked of the caller
    int filename_strlen = tox::max_strlen(filename);
    int gene_ids_strlen = tox::max_strlen(gene_ids);
    int n_gene_ids_elements = (int) gene_ids.size();

    // convert what C cannot take directly
    tox::CharBuffer filename_c(filename, filename_strlen);
    tox::CharBuffer gene_ids_c(gene_ids, gene_ids_strlen);

    // outputs and work space
    tox::CharBuffer family_ids_c(family_ids_strlen, n_family_ids_elements);
    IntegerVector gene_to_fam(n_gene_to_fam_elements);
    int ierr = 0;

    read_orthofinder_file_c(
        filename_c.data(),
        &filename_strlen,
        gene_ids_c.data(),
        &gene_ids_strlen,
        &n_gene_ids_elements,
        family_ids_c.data(),
        &family_ids_strlen,
        &n_family_ids_elements,
        gene_to_fam.begin(),
        &n_gene_to_fam_elements,
        &ierr
    );

    // convert the outputs back
    CharacterVector family_ids = family_ids_c.to_r();

    return List::create(
        _["family_ids"] = family_ids,
        _["gene_to_fam"] = gene_to_fam,
        _["ierr"] = ierr
    );
}

// [[Rcpp::export(.get_unassigned_mask_rcpp)]]
List get_unassigned_mask_rcpp(IntegerVector gene_to_fam) {
    // derived from the inputs, not asked of the caller
    int n_gene_to_fam_elements = (int) gene_to_fam.size();

    // outputs and work space
    tox::BoolBuffer mask_c(((int) gene_to_fam.size()));
    int n_genes_kept = 0;
    int ierr = 0;

    get_unassigned_mask_c(
        gene_to_fam.begin(),
        &n_gene_to_fam_elements,
        mask_c.data(),
        &n_genes_kept,
        &ierr
    );

    // convert the outputs back
    LogicalVector mask = mask_c.to_r();

    return List::create(
        _["mask"] = mask,
        _["n_genes_kept"] = n_genes_kept,
        _["ierr"] = ierr
    );
}
