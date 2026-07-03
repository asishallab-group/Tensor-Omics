#include <Rcpp.h>

using namespace Rcpp;

static_assert(sizeof(Rcomplex) == sizeof(double _Complex) && alignof(Rcomplex) == alignof(double _Complex), "Rcomplex layout is incompatible with C's 'double _Complex' and thus Fortran's 'c_double_complex'");

extern "C" {
    void mask_get_first_successor_idx_c(
        const int* bit_mask,
        const int* n_bit_mask_elements,
        int* idx,
        int* ierr
    );

    void mask_check_state_c(
        const int* bit_mask,
        const int* n_bit_mask_elements,
        const int* i_gene,
        int* state,
        int* ierr
    );

    void detect_neofunctionalization_c(
        const double* ancestors,
        const int* n_families,
        const double* genes,
        const int* n_axes,
        const int* gene_to_fam,
        const int* n_genes,
        const double* thresholds,
        int* neofunc,
        int* ierr
    );

    void detect_dosage_effect_c(
        const double* ancestor,
        const double* genes,
        const int* n_genes,
        const int* n_dims,
        const int* filtered_paralogs_mask,
        const int* n_mask_chunks,
        int* n_results,
        const int* max_subset_size,
        int* work_arr_paralog_subsets,
        const int* n_paralog_subsets,
        int* active_mask,
        double* temp_paralog_vector,
        int* ierr,
        const double* max_angle,
        const double* gain_gamma
    );

    void detect_subfunctionalization_c(
        const double* ancestor,
        const double* genes,
        const int* n_genes,
        const int* n_dims,
        const double* rdi_threshold,
        const int* filtered_paralogs_mask,
        const int* n_mask_chunks,
        int* n_results,
        const int* max_subset_size,
        int* work_arr_paralog_subsets,
        const int* n_paralog_subsets,
        int* active_mask,
        double* temp_paralog_vector,
        const double* paralog_norms,
        const int* sorted_paralog_norms_perm,
        double* temp_work_array,
        int* ierr
    );

    void detect_patterns_c(
        const double* ancestor,
        const double* genes,
        const int* n_genes,
        const int* n_dims,
        const int* pattern_mode,
        const int* filtered_paralogs_mask,
        const int* n_mask_chunks,
        int* n_results,
        const int* max_subset_size,
        int* work_arr_paralog_subsets,
        const int* n_paralog_subsets,
        int* active_mask,
        double* temp_paralog_vector,
        const double* dosage_max_angle,
        const double* dosage_gain_gamma,
        const double* subfunc_rdi_threshold,
        const double* subfunc_paralog_norms,
        const int* subfunc_sorted_paralog_norms_perm,
        double* subfunc_temp_work_array,
        int* ierr
    );

    void mask_chunk_count_c(
        const int* n_genes,
        int* count,
        int* ierr
    );

    void filter_paralogs_by_pattern_subfunctionalization_c(
        const double* gene_angles,
        const double* threshold,
        const int* n_genes,
        const int* n_families,
        const int* gene_to_fam,
        int* masks,
        const int* n_mask_chunks,
        int* ierr
    );

    void filter_paralogs_by_pattern_dosage_effect_c(
        const double* gene_angles,
        const double* threshold,
        const int* n_genes,
        const int* n_families,
        const int* gene_to_fam,
        int* masks,
        const int* n_mask_chunks,
        int* ierr
    );

    void filter_paralogs_by_pattern_c(
        const int* pattern_mode,
        const double* gene_angles,
        const double* threshold,
        const int* n_genes,
        const int* n_families,
        const int* gene_to_fam,
        int* masks,
        const int* n_mask_chunks,
        int* ierr
    );

    void calc_work_arr_paralog_subsets_size_c(
        int* max_subset_size,
        const int* n_genes,
        int* work_array_size,
        const int* filtered_paralogs_mask,
        const int* n_mask_chunks,
        int* ierr
    );

    void mask_set_state_c(
        int* bit_mask,
        const int* n_bit_mask_elements,
        const int* i_gene,
        const int* state,
        int* ierr
    );
}

List mask_get_first_successor_idx_rcpp(
    const IntegerVector& bit_mask
) {


    IntegerVector bit_mask_shape = bit_mask.attr("dim");
    n_bit_mask_elements = bit_mask_shape[0]



    int idx = 0;
    int ierr = 0;

    mask_get_first_successor_idx_c(
        bit_mask.begin(),
        &n_bit_mask_elements,
        &idx,
        &ierr
    );

//{format(self.arguments, "type_conversion_outputs")}

    return List::create(
        Named("idx") = idx,
        Named("ierr") = ierr
    );
}

List mask_check_state_rcpp(
    const IntegerVector& bit_mask,
    const int i_gene
) {


    IntegerVector bit_mask_shape = bit_mask.attr("dim");
    n_bit_mask_elements = bit_mask_shape[0]



    int state = 0;
    int ierr = 0;

    mask_check_state_c(
        bit_mask.begin(),
        &n_bit_mask_elements,
        &i_gene,
        &state,
        &ierr
    );

//{format(self.arguments, "type_conversion_outputs")}

    return List::create(
        Named("state") = state,
        Named("ierr") = ierr
    );
}

List detect_neofunctionalization_rcpp(
    const NumericMatrix& ancestors,
    const NumericMatrix& genes,
    const IntegerVector& gene_to_fam,
    const NumericVector& thresholds
) {


    IntegerVector ancestors_shape = ancestors.attr("dim");
    n_axes = ancestors_shape[0]
    n_families = ancestors_shape[1]
    IntegerVector genes_shape = genes.attr("dim");
    n_genes = genes_shape[1]



    LogicalVector neofunc(n_genes * n_axes);
    int ierr = 0;

    detect_neofunctionalization_c(
        ancestors.begin(),
        &n_families,
        genes.begin(),
        &n_axes,
        gene_to_fam.begin(),
        &n_genes,
        thresholds.begin(),
        neofunc_p,
        &ierr
    );

//{format(self.arguments, "type_conversion_outputs")}

    return List::create(
        Named("neofunc") = neofunc,
        Named("ierr") = ierr
    );
}

List detect_dosage_effect_rcpp(
    const NumericVector& ancestor,
    const NumericMatrix& genes,
    const IntegerVector& filtered_paralogs_mask,
    const int max_subset_size,
    const double max_angle = 3.141592653589793,
    const double gain_gamma = 0.1
) {


    IntegerVector ancestor_shape = ancestor.attr("dim");
    n_dims = ancestor_shape[0]
    IntegerVector genes_shape = genes.attr("dim");
    n_genes = genes_shape[1]
    IntegerVector filtered_paralogs_mask_shape = filtered_paralogs_mask.attr("dim");
    n_mask_chunks = filtered_paralogs_mask_shape[0]
    IntegerVector work_arr_paralog_subsets_shape = work_arr_paralog_subsets.attr("dim");
    n_paralog_subsets = work_arr_paralog_subsets_shape[1]



    int n_results = 0;
    IntegerMatrix work_arr_paralog_subsets(n_mask_chunks * n_paralog_subsets);
    IntegerVector active_mask(n_mask_chunks);
    NumericVector temp_paralog_vector(n_dims);
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
        active_mask.begin(),
        temp_paralog_vector.begin(),
        &ierr,
        &max_angle,
        &gain_gamma
    );

//{format(self.arguments, "type_conversion_outputs")}

    return List::create(
        Named("n_results") = n_results,
        Named("work_arr_paralog_subsets") = work_arr_paralog_subsets,
        Named("active_mask") = active_mask,
        Named("temp_paralog_vector") = temp_paralog_vector,
        Named("ierr") = ierr
    );
}

List detect_subfunctionalization_rcpp(
    const NumericVector& ancestor,
    const NumericMatrix& genes,
    const double rdi_threshold,
    const IntegerVector& filtered_paralogs_mask,
    const int max_subset_size,
    const NumericVector& paralog_norms,
    const IntegerVector& sorted_paralog_norms_perm
) {


    IntegerVector ancestor_shape = ancestor.attr("dim");
    n_dims = ancestor_shape[0]
    IntegerVector genes_shape = genes.attr("dim");
    n_genes = genes_shape[1]
    IntegerVector filtered_paralogs_mask_shape = filtered_paralogs_mask.attr("dim");
    n_mask_chunks = filtered_paralogs_mask_shape[0]
    IntegerVector work_arr_paralog_subsets_shape = work_arr_paralog_subsets.attr("dim");
    n_paralog_subsets = work_arr_paralog_subsets_shape[1]



    int n_results = 0;
    IntegerMatrix work_arr_paralog_subsets(n_mask_chunks * n_paralog_subsets);
    IntegerVector active_mask(n_mask_chunks);
    NumericVector temp_paralog_vector(n_dims);
    NumericVector temp_work_array(n_genes);
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
        active_mask.begin(),
        temp_paralog_vector.begin(),
        paralog_norms.begin(),
        sorted_paralog_norms_perm.begin(),
        temp_work_array.begin(),
        &ierr
    );

//{format(self.arguments, "type_conversion_outputs")}

    return List::create(
        Named("n_results") = n_results,
        Named("work_arr_paralog_subsets") = work_arr_paralog_subsets,
        Named("active_mask") = active_mask,
        Named("temp_paralog_vector") = temp_paralog_vector,
        Named("temp_work_array") = temp_work_array,
        Named("ierr") = ierr
    );
}

List detect_patterns_rcpp(
    const NumericVector& ancestor,
    const NumericMatrix& genes,
    const int pattern_mode,
    const IntegerVector& filtered_paralogs_mask,
    const int max_subset_size,
    const double dosage_max_angle = 3.141592653589793,
    const double dosage_gain_gamma = 0.1,
    const double subfunc_rdi_threshold = ,
    const NumericVector& subfunc_paralog_norms = ,
    const IntegerVector& subfunc_sorted_paralog_norms_perm = 
) {


    IntegerVector ancestor_shape = ancestor.attr("dim");
    n_dims = ancestor_shape[0]
    IntegerVector genes_shape = genes.attr("dim");
    n_genes = genes_shape[1]
    IntegerVector filtered_paralogs_mask_shape = filtered_paralogs_mask.attr("dim");
    n_mask_chunks = filtered_paralogs_mask_shape[0]
    IntegerVector work_arr_paralog_subsets_shape = work_arr_paralog_subsets.attr("dim");
    n_paralog_subsets = work_arr_paralog_subsets_shape[1]



    int n_results = 0;
    IntegerMatrix work_arr_paralog_subsets(n_mask_chunks * n_paralog_subsets);
    IntegerVector active_mask(n_mask_chunks);
    NumericVector temp_paralog_vector(n_dims);
    NumericVector subfunc_temp_work_array(n_genes);
    int ierr = 0;

    detect_patterns_c(
        ancestor.begin(),
        genes.begin(),
        &n_genes,
        &n_dims,
        &pattern_mode,
        filtered_paralogs_mask.begin(),
        &n_mask_chunks,
        &n_results,
        &max_subset_size,
        work_arr_paralog_subsets.begin(),
        &n_paralog_subsets,
        active_mask.begin(),
        temp_paralog_vector.begin(),
        &dosage_max_angle,
        &dosage_gain_gamma,
        &subfunc_rdi_threshold,
        subfunc_paralog_norms.begin(),
        subfunc_sorted_paralog_norms_perm.begin(),
        subfunc_temp_work_array.begin(),
        &ierr
    );

//{format(self.arguments, "type_conversion_outputs")}

    return List::create(
        Named("n_results") = n_results,
        Named("work_arr_paralog_subsets") = work_arr_paralog_subsets,
        Named("active_mask") = active_mask,
        Named("temp_paralog_vector") = temp_paralog_vector,
        Named("subfunc_temp_work_array") = subfunc_temp_work_array,
        Named("ierr") = ierr
    );
}

List mask_chunk_count_rcpp(
    const int n_genes
) {





    int count = 0;
    int ierr = 0;

    mask_chunk_count_c(
        &n_genes,
        &count,
        &ierr
    );

//{format(self.arguments, "type_conversion_outputs")}

    return List::create(
        Named("count") = count,
        Named("ierr") = ierr
    );
}

List filter_paralogs_by_pattern_subfunctionalization_rcpp(
    const NumericVector& gene_angles,
    const double threshold,
    const IntegerVector& gene_to_fam
) {


    IntegerVector gene_angles_shape = gene_angles.attr("dim");
    n_genes = gene_angles_shape[0]
    IntegerVector masks_shape = masks.attr("dim");
    n_mask_chunks = masks_shape[0]
    n_families = masks_shape[1]



    IntegerMatrix masks(n_mask_chunks * n_families);
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

//{format(self.arguments, "type_conversion_outputs")}

    return List::create(
        Named("masks") = masks,
        Named("ierr") = ierr
    );
}

List filter_paralogs_by_pattern_dosage_effect_rcpp(
    const NumericVector& gene_angles,
    const double threshold,
    const IntegerVector& gene_to_fam
) {


    IntegerVector gene_angles_shape = gene_angles.attr("dim");
    n_genes = gene_angles_shape[0]
    IntegerVector masks_shape = masks.attr("dim");
    n_mask_chunks = masks_shape[0]
    n_families = masks_shape[1]



    IntegerMatrix masks(n_mask_chunks * n_families);
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

//{format(self.arguments, "type_conversion_outputs")}

    return List::create(
        Named("masks") = masks,
        Named("ierr") = ierr
    );
}

List filter_paralogs_by_pattern_rcpp(
    const int pattern_mode,
    const NumericVector& gene_angles,
    const double threshold,
    const IntegerVector& gene_to_fam
) {


    IntegerVector gene_angles_shape = gene_angles.attr("dim");
    n_genes = gene_angles_shape[0]
    IntegerVector masks_shape = masks.attr("dim");
    n_mask_chunks = masks_shape[0]
    n_families = masks_shape[1]



    IntegerMatrix masks(n_mask_chunks * n_families);
    int ierr = 0;

    filter_paralogs_by_pattern_c(
        &pattern_mode,
        gene_angles.begin(),
        &threshold,
        &n_genes,
        &n_families,
        gene_to_fam.begin(),
        masks.begin(),
        &n_mask_chunks,
        &ierr
    );

//{format(self.arguments, "type_conversion_outputs")}

    return List::create(
        Named("masks") = masks,
        Named("ierr") = ierr
    );
}

List calc_work_arr_paralog_subsets_size_rcpp(
    int max_subset_size,
    const int n_genes,
    const IntegerVector& filtered_paralogs_mask
) {


    IntegerVector filtered_paralogs_mask_shape = filtered_paralogs_mask.attr("dim");
    n_mask_chunks = filtered_paralogs_mask_shape[0]



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

//{format(self.arguments, "type_conversion_outputs")}

    return List::create(
        Named("work_array_size") = work_array_size,
        Named("ierr") = ierr
    );
}

List mask_set_state_rcpp(
    const int i_gene,
    const LogicalVector state
) {


    IntegerVector bit_mask_shape = bit_mask.attr("dim");
    n_bit_mask_elements = bit_mask_shape[0]

    const int* state_p = state.begin();
    // Is already 0/1 mapped, thus check only for NA values and return ERR_NAN_INF code
    for (int i = 0; i < state.size(); ++i) {
        if (state_p[i] == NA_LOGICAL) {
            return List::create(Named("ierr") = 204);
        }
    }

    IntegerVector bit_mask(n_bit_mask_elements);
    int ierr = 0;

    mask_set_state_c(
        bit_mask.begin(),
        &n_bit_mask_elements,
        &i_gene,
        &state,
        &ierr
    );

//{format(self.arguments, "type_conversion_outputs")}

    return List::create(
        Named("bit_mask") = bit_mask,
        Named("ierr") = ierr
    );
}