// Generated. Do not edit.
#include <Rcpp.h>
#include <numeric>
#include <functional>
#include "tox_marshal.h"

using namespace Rcpp;

extern "C" {
    void mask_check_state_c(const int*, const int*, const int*, bool*, int*);
    void detect_neofunctionalization_c(const double*, const int*, const double*, const int*, const int*, const int*, const double*, bool*, int*);
    void detect_dosage_effect_c(const double*, const double*, const int*, const int*, const int*, const int*, int*, const int*, int*, const int*, int*, double*, int*, const double*, const double*);
    void detect_subfunctionalization_c(const double*, const double*, const int*, const int*, const double*, const int*, const int*, int*, const int*, int*, const int*, int*, double*, const double*, const int*, double*, int*);
    void mask_chunk_count_c(const int*, int*, int*);
    void filter_paralogs_by_pattern_subfunctionalization_c(const double*, const double*, const int*, const int*, const int*, int*, const int*, int*);
    void filter_paralogs_by_pattern_dosage_effect_c(const double*, const double*, const int*, const int*, const int*, int*, const int*, int*);
    void calc_work_arr_paralog_subsets_size_c(int*, const int*, int*, const int*, const int*, int*);
}

// [[Rcpp::export(.mask_check_state_rcpp)]]
List mask_check_state_rcpp(IntegerVector bit_mask, int i_gene) {
    // derived from the inputs, not asked of the caller
    int n_bit_mask_elements = (int) bit_mask.size();

    // outputs and work space
    bool state = 0;
    int ierr = 0;

    mask_check_state_c(
        bit_mask.begin(),
        &n_bit_mask_elements,
        &i_gene,
        &state,
        &ierr
    );

    return List::create(
        _["state"] = state,
        _["ierr"] = ierr
    );
}

// [[Rcpp::export(.detect_neofunctionalization_rcpp)]]
List detect_neofunctionalization_rcpp(NumericVector ancestors, NumericVector genes, IntegerVector gene_to_fam, NumericVector thresholds) {
    // derived from the inputs, not asked of the caller
    int n_families = (int) IntegerVector(ancestors.attr("dim"))[1];
    int n_axes = (int) IntegerVector(ancestors.attr("dim"))[0];
    int n_genes = (int) IntegerVector(genes.attr("dim"))[1];

    // outputs and work space
    tox::BoolBuffer neofunc_c(n_genes * n_axes);
    int ierr = 0;

    detect_neofunctionalization_c(
        ancestors.begin(),
        &n_families,
        genes.begin(),
        &n_axes,
        gene_to_fam.begin(),
        &n_genes,
        thresholds.begin(),
        neofunc_c.data(),
        &ierr
    );

    // convert the outputs back
    LogicalVector neofunc = neofunc_c.to_r();
    neofunc.attr("dim") = IntegerVector::create(n_genes, n_axes);

    return List::create(
        _["neofunc"] = neofunc,
        _["ierr"] = ierr
    );
}

// [[Rcpp::export(.detect_dosage_effect_rcpp)]]
List detect_dosage_effect_rcpp(NumericVector ancestor, NumericVector genes, IntegerVector filtered_paralogs_mask, int max_subset_size, int n_paralog_subsets, Nullable<NumericVector> max_angle = R_NilValue, Nullable<NumericVector> gain_gamma = R_NilValue) {
    // optionals: a null pointer and size 0 when the caller omits them
    const double* max_angle_p = nullptr;
    int max_angle_size = 0;
    NumericVector max_angle_val;
    if (max_angle.isNotNull()) {
        max_angle_val = max_angle.get();
        max_angle_p = max_angle_val.begin();
        max_angle_size = max_angle_val.size();
    }
    const double* gain_gamma_p = nullptr;
    int gain_gamma_size = 0;
    NumericVector gain_gamma_val;
    if (gain_gamma.isNotNull()) {
        gain_gamma_val = gain_gamma.get();
        gain_gamma_p = gain_gamma_val.begin();
        gain_gamma_size = gain_gamma_val.size();
    }

    // derived from the inputs, not asked of the caller
    int n_genes = (int) IntegerVector(genes.attr("dim"))[1];
    int n_dims = (int) ancestor.size();
    int n_mask_chunks = (int) filtered_paralogs_mask.size();

    // outputs and work space
    int n_results = 0;
    IntegerVector work_arr_paralog_subsets(n_mask_chunks * n_paralog_subsets);
    work_arr_paralog_subsets.attr("dim") = IntegerVector::create(n_mask_chunks, n_paralog_subsets);
    std::vector<int> tmp_active_mask(n_mask_chunks);
    std::vector<double> tmp_paralog_vector(n_dims);
    int ierr = 0;

    detect_dosage_effect_c(
        ancestor.begin(),
        genes.begin(),
        &n_genes,
        &n_dims,
        filtered_paralogs_mask.begin(),
        &n_mask_chunks,
        &n_results,
        &max_subset_size,
        work_arr_paralog_subsets.begin(),
        &n_paralog_subsets,
        tmp_active_mask.data(),
        tmp_paralog_vector.data(),
        &ierr,
        max_angle_p,
        gain_gamma_p
    );

    return List::create(
        _["n_results"] = n_results,
        _["work_arr_paralog_subsets"] = work_arr_paralog_subsets,
        _["ierr"] = ierr
    );
}

// [[Rcpp::export(.detect_subfunctionalization_rcpp)]]
List detect_subfunctionalization_rcpp(NumericVector ancestor, NumericVector genes, double rdi_threshold, IntegerVector filtered_paralogs_mask, int max_subset_size, int n_paralog_subsets, NumericVector paralog_norms, IntegerVector sorted_paralog_norms_perm) {
    // derived from the inputs, not asked of the caller
    int n_genes = (int) IntegerVector(genes.attr("dim"))[1];
    int n_dims = (int) ancestor.size();
    int n_mask_chunks = (int) filtered_paralogs_mask.size();

    // outputs and work space
    int n_results = 0;
    IntegerVector work_arr_paralog_subsets(n_mask_chunks * n_paralog_subsets);
    work_arr_paralog_subsets.attr("dim") = IntegerVector::create(n_mask_chunks, n_paralog_subsets);
    std::vector<int> tmp_active_mask(n_mask_chunks);
    std::vector<double> tmp_paralog_vector(n_dims);
    std::vector<double> tmp_work_array(n_genes);
    int ierr = 0;

    detect_subfunctionalization_c(
        ancestor.begin(),
        genes.begin(),
        &n_genes,
        &n_dims,
        &rdi_threshold,
        filtered_paralogs_mask.begin(),
        &n_mask_chunks,
        &n_results,
        &max_subset_size,
        work_arr_paralog_subsets.begin(),
        &n_paralog_subsets,
        tmp_active_mask.data(),
        tmp_paralog_vector.data(),
        paralog_norms.begin(),
        sorted_paralog_norms_perm.begin(),
        tmp_work_array.data(),
        &ierr
    );

    return List::create(
        _["n_results"] = n_results,
        _["work_arr_paralog_subsets"] = work_arr_paralog_subsets,
        _["ierr"] = ierr
    );
}

// [[Rcpp::export(.mask_chunk_count_rcpp)]]
List mask_chunk_count_rcpp(int n_genes) {
    // outputs and work space
    int count = 0;
    int ierr = 0;

    mask_chunk_count_c(
        &n_genes,
        &count,
        &ierr
    );

    return List::create(
        _["count"] = count,
        _["ierr"] = ierr
    );
}

// [[Rcpp::export(.filter_paralogs_by_pattern_subfunctionalization_rcpp)]]
List filter_paralogs_by_pattern_subfunctionalization_rcpp(NumericVector gene_angles, double threshold, int n_families, IntegerVector gene_to_fam, int n_mask_chunks) {
    // derived from the inputs, not asked of the caller
    int n_genes = (int) gene_angles.size();

    // outputs and work space
    IntegerVector masks(n_mask_chunks * n_families);
    masks.attr("dim") = IntegerVector::create(n_mask_chunks, n_families);
    int ierr = 0;

    filter_paralogs_by_pattern_subfunctionalization_c(
        gene_angles.begin(),
        &threshold,
        &n_genes,
        &n_families,
        gene_to_fam.begin(),
        masks.begin(),
        &n_mask_chunks,
        &ierr
    );

    return List::create(
        _["masks"] = masks,
        _["ierr"] = ierr
    );
}

// [[Rcpp::export(.filter_paralogs_by_pattern_dosage_effect_rcpp)]]
List filter_paralogs_by_pattern_dosage_effect_rcpp(NumericVector gene_angles, double threshold, int n_families, IntegerVector gene_to_fam, int n_mask_chunks) {
    // derived from the inputs, not asked of the caller
    int n_genes = (int) gene_angles.size();

    // outputs and work space
    IntegerVector masks(n_mask_chunks * n_families);
    masks.attr("dim") = IntegerVector::create(n_mask_chunks, n_families);
    int ierr = 0;

    filter_paralogs_by_pattern_dosage_effect_c(
        gene_angles.begin(),
        &threshold,
        &n_genes,
        &n_families,
        gene_to_fam.begin(),
        masks.begin(),
        &n_mask_chunks,
        &ierr
    );

    return List::create(
        _["masks"] = masks,
        _["ierr"] = ierr
    );
}

// [[Rcpp::export(.calc_work_arr_paralog_subsets_size_rcpp)]]
List calc_work_arr_paralog_subsets_size_rcpp(int max_subset_size, int n_genes, IntegerVector filtered_paralogs_mask) {
    // derived from the inputs, not asked of the caller
    int n_mask_chunks = (int) filtered_paralogs_mask.size();

    // outputs and work space
    int work_array_size = 0;
    int ierr = 0;

    calc_work_arr_paralog_subsets_size_c(
        &max_subset_size,
        &n_genes,
        &work_array_size,
        filtered_paralogs_mask.begin(),
        &n_mask_chunks,
        &ierr
    );

    return List::create(
        _["max_subset_size"] = max_subset_size,
        _["work_array_size"] = work_array_size,
        _["ierr"] = ierr
    );
}
