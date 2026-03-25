#include <Rcpp.h>

using namespace Rcpp;


// ===================================================================
// FORTRAN FUNCTIONS
// ===================================================================

extern "C" {

    void compute_gene_means_c( double* expr, int* n_genes, int* n_reps, double* means, int* max_n_genes_all_studies, int* ierr);

    void compute_residuals_c( double* expr, int* n_genes, int* n_reps, double* means, int* max_n_genes_all_studies, int* max_n_reps_all_studies, double* resid, int* ierr);

    void determine_js_comp_test_n_points_n_neighbors_alloc_c(int* n_points, int* n_neighbors, double* residuals, int* max_n_reps_all_studies, int* max_n_genes_all_studies, double* shared_residual_range, int* n_bins, double* gene_means, int* n_studies, int* n_bootstraps, double* best_candidate_pair_confidence_interval, const char* join_method, int* min_count_per_mean_bin, double* min_neighbor_overlap, double* succeeding_ci_overlap, double* two_sided_bootstrapping_significance_level, int* random_seed, double* residual_range_quantile, int* ierr);

    void js_comp_test_alloc_c(double* gene_means, int* max_n_genes_all_studies, int* n_studies, double* residuals, double* shared_residual_range, int* n_bins, int* max_n_reps_all_studies, double* x_star, int* n_pool, int* n_points, int* n_neighbors, int* neighborhood_ranges, int* neighborhood_residuals, double* pmfs, int* counts, int* included_n_reps, double* mean_pmf, int* mean_pmf_counts, int* mean_pmf_included_n_reps, double* js_divergences, double* weights, double* global_js_divergence, double* p_values, int* ierr, int* n_permutations, int* random_seed);

    void normalize_by_std_dev_c(int n_genes, int n_tissues,
                                double *input_matrix, double *output_matrix, int *ierr);
    void quantile_normalization_c(int n_genes, int n_tissues, double *input_matrix, double *output_matrix,
                                  double *temp_col, double *rank_means,
                                  int *perm, int *stack_left, int *stack_right,
                                  int max_stack, int *ierr);
    void log2_transformation_c(int n_genes, int n_tissues,
                               double *input_matrix, double *output_matrix, int *ierr);

    void calc_tiss_avg_c(int n_genes, int n_grps,int *group_s, int *group_c,
                         double *input_matrix, double *output_matrix, int *ierr);
    void calc_fchange_c(int n_genes, int n_cols, int n_pairs,int *control_cols, int *cond_cols,
                        double *input_matrix, double *output_matrix, int *ierr);

    void normalization_pipeline_c(int n_genes, int n_tissues,
                                  double *input_matrix, double *buf_stddev, double *buf_quant,
                                  double *buf_avg, double *buf_log, double *temp_col,
                                  double *rank_means, int *perm, int *stack_left,
                                  int *stack_right, int max_stack,
                                  int *group_s, int *group_c, int n_grps, int *ierr);

      void compute_family_scaling_c(
        int n_genes, int n_families,
        double* distances, int* gene_to_fam,
        double* dscale,
        double* loess_x, double* loess_y, int* indices_used,
        int* ierr
      );

      void compute_family_scaling_expert_c(
        int n_genes, int n_families,
        double* distances, int* gene_to_fam,
        double* dscale,
        double* loess_x, double* loess_y, int* indices_used,
        int* perm_tmp, int* stack_left_tmp, int* stack_right_tmp,
        double* family_distances,
        int* ierr
      );

      void compute_rdi_c(
        int n_genes, int n_families,
        double* distances, int* gene_to_fam,
        double* dscale,
        double* rdi, double* sorted_rdi,
        int* perm, int* stack_left, int* stack_right
      );

      void identify_outliers_c(
        int n_genes,
        double* rdi, double* sorted_rdi,
        int* is_outlier_int,
        double* threshold,
        double percentile
      );

      void detect_outliers_c(
        int n_genes, int n_families,
        double* distances, int* gene_to_fam,
        double* work_array,
        int* perm, int* stack_left, int* stack_right,
        int* is_outlier_int,
        double* loess_x, double* loess_y, int* loess_n,
        int* ierr,
        double percentile
      );


      void euclidean_distance_c(double* vec1, double* vec2, int d, double* result);

      void distance_to_centroid_c(
        int n_genes, int n_families,
        double* genes, double* centroids,
        int* gene_to_fam,
        double* distances,
        int d
      );

      void compute_tissue_versatility_c(
        int n_axes, int n_vectors,
        double* expression_vectors,
        int* exp_vecs_selection_index,
        int n_selected_vectors,
        int* axes_selection,
        int n_selected_axes,
        double* tissue_versatilities,
        double* tissue_angles_deg,
        int* ierr
      );

    void euclidean_distance_c(double* vec1, double* vec2, int d, double* result);
    void distance_to_centroid_c(int n_genes, int n_families, double* genes,
                                double* centroids, int* gene_to_fam,
                                double* distances, int d);
    void compute_tissue_versatility_c(int n_axes, int n_vectors,
                                      double* expression_vectors,
                                      int* exp_vecs_selection_index,
                                      int n_selected_vectors,
                                      int* axes_selection,
                                      int n_selected_axes,
                                      double* tissue_versatilities,
                                      double* tissue_angles_deg,
                                      int* ierr);
    void compute_shift_vector_field_c(int d, int n_genes, int n_families,
                                      double* expression_vectors, double* family_centroids,
                                      int* gene_to_centroid, double* shift_vectors,
                                      int* ierr);

    void mean_vector_c(double* expression_vectors, int n_axes, int n_genes,
                       int* gene_indices, int n_selected_genes,
                       double* centroid_col, int* ierr);

    void group_centroid_c(double* expression_vectors, int n_axes, int n_genes,
                         int* gene_to_family, int n_families,
                         double* centroid_matrix, const char* mode,
                         int* ortholog_set, int* selected_indices, int selected_indices_len, int* ierr);

}

/**
 * Calculate Euclidean distance between two vectors
 */
// [[Rcpp::export]]
double tox_euclidean_distance_rcpp(NumericVector vec1, NumericVector vec2) {
    int d = vec1.length();
    double result = 0.0;

    euclidean_distance_c(vec1.begin(), vec2.begin(), d, &result);
    return result;
}

/**
 * Calculate distances from genes to their family centroids
 */
// [[Rcpp::export]]
NumericVector tox_distance_to_centroid_rcpp(NumericVector genes, NumericVector centroids,
                                       IntegerVector gene_to_fam, int d) {
    int n_genes = genes.length() / d;
    int n_families = centroids.length() / d;

    NumericVector distances(n_genes);

    distance_to_centroid_c(n_genes, n_families, genes.begin(),
                          centroids.begin(), gene_to_fam.begin(),
                          distances.begin(), d);

    return distances;
}

/**
 * Calculate Tissue Versatility
 */
// [[Rcpp::export]]
List tox_calculate_tissue_versatility_rcpp(NumericMatrix expression_vectors,
                                      IntegerVector vector_selection,
                                      IntegerVector axis_selection) {
    int n_axes = expression_vectors.nrow();
    int n_vectors = expression_vectors.ncol();
    int n_selected_vectors = sum(vector_selection);
    int n_selected_axes = sum(axis_selection);

    NumericVector tissue_versatilities(n_selected_vectors);
    NumericVector tissue_angles_deg(n_selected_vectors);
    int ierr = 0;

    compute_tissue_versatility_c(n_axes, n_vectors, expression_vectors.begin(),
                                vector_selection.begin(), n_selected_vectors,
                                axis_selection.begin(), n_selected_axes,
                                tissue_versatilities.begin(),
                                tissue_angles_deg.begin(), &ierr);

    return List::create(
        Named("tissue_versatilities") = tissue_versatilities,
        Named("tissue_angles_deg") = tissue_angles_deg,
        Named("n_selected_vectors") = n_selected_vectors,
        Named("n_selected_axes") = n_selected_axes,
        Named("ierr") = ierr
    );

}

/**
 * Calculate normalization by standard deviation
 */
// [[Rcpp::export]]
List tox_normalize_by_std_dev_rcpp(NumericMatrix input) {
    int n_genes = input.nrow();
    int n_tissues = input.ncol();
    NumericMatrix output(n_genes, n_tissues);
    int ierr = 0;

    normalize_by_std_dev_c(n_genes, n_tissues, input.begin(), output.begin(), &ierr);

    return List::create(
        Named("output_vector") = output,
        Named("ierr") = ierr
    );
}

/**
 * Perform quantile normalization
 */
// [[Rcpp::export]]
List tox_quantile_normalization_rcpp(NumericMatrix input) {
    int n_genes = input.nrow();
    int n_tissues = input.ncol();

    int max_stack = static_cast<int>(std::ceil(std::log2(n_genes)) + 10);

    NumericMatrix output(n_genes, n_tissues);
    NumericVector temp_col(n_genes);
    NumericVector rank_means(n_genes);
    IntegerVector perm(n_genes);
    IntegerVector stack_left(max_stack);
    IntegerVector stack_right(max_stack);
    int ierr = 0;

    quantile_normalization_c(n_genes, n_tissues, input.begin(), output.begin(),
                             temp_col.begin(), rank_means.begin(), perm.begin(),
                             stack_left.begin(), stack_right.begin(), max_stack, &ierr);

    return List::create(
        Named("output_vector")     = output,
        Named("rank_means") = rank_means,
        Named("perm")       = perm,
        Named("ierr")       = ierr
    );
}

/**
 * Perform log2 transformation
 */
// [[Rcpp::export]]
List tox_log2_transformation_rcpp(NumericMatrix input) {
    int n_genes = input.nrow();
    int n_tissues = input.ncol();
    NumericMatrix output(n_genes, n_tissues);
    int ierr = 0;

    log2_transformation_c(n_genes, n_tissues, input.begin(), output.begin(), &ierr);

    return List::create(
        Named("output_vector") = output,
        Named("ierr") = ierr
    );
}

// [[Rcpp::export]]
List tox_compute_shift_vector_field_rcpp(NumericMatrix expression_vectors, NumericMatrix family_centroids, IntegerVector gene_to_centroid) {
    int n_axes_genes = expression_vectors.nrow();
    int n_vectors = expression_vectors.ncol();
    int n_axes_centroids = family_centroids.nrow();
    int n_families = family_centroids.ncol();

    NumericMatrix shift_vectors(2 * n_axes_genes, n_vectors);
    int ierr = 0;

    compute_shift_vector_field_c(n_axes_genes, n_vectors, n_families,
                                 expression_vectors.begin(), family_centroids.begin(),
                                 gene_to_centroid.begin(), shift_vectors.begin(), &ierr);

    NumericVector flat(shift_vectors.begin(), shift_vectors.end());

    return List::create(
        Named("shift_vectors") = flat,
        Named("ierr") = ierr
    );
}

/**
 * Calculate tissue averages
 */
// [[Rcpp::export]]
List tox_calc_tiss_avg_rcpp(NumericMatrix input, IntegerVector group_s, IntegerVector group_c) {
    int n_gene = input.nrow();
    int n_grps = group_s.size();
    NumericMatrix output(n_gene, n_grps);
    int ierr = 0;

    calc_tiss_avg_c(n_gene, n_grps, group_s.begin(), group_c.begin(), input.begin(), output.begin(), &ierr);
    return List::create(
        Named("output_vector") = output,
        Named("ierr") = ierr
    );
}

/**
 * Calculate fold change
 */
// [[Rcpp::export]]
List tox_calc_fchange_rcpp(NumericMatrix input, IntegerVector control_cols, IntegerVector cond_cols) {
    int n_genes = input.nrow();
    int n_cols = input.ncol();
    int n_pairs = control_cols.size();
    NumericMatrix output(n_genes, n_pairs);
    int ierr = 0;

    calc_fchange_c(n_genes, n_cols, n_pairs, control_cols.begin(), cond_cols.begin(), input.begin(), output.begin(), &ierr);
    return List::create(
        Named("output_vector") = output,
        Named("ierr") = ierr
    );
}

// [[Rcpp::export]]
List tox_mean_vector_rcpp(NumericMatrix expression_vectors, IntegerVector gene_indices) {
    int n_axes = expression_vectors.nrow();
    int n_genes = expression_vectors.ncol();
    int n_selected_genes = gene_indices.length();

    NumericVector centroid_col(n_axes);
    int ierr = 0;

    mean_vector_c(expression_vectors.begin(), n_axes, n_genes,
                  gene_indices.begin(), n_selected_genes,
                  centroid_col.begin(), &ierr);

    return List::create(
        Named("centroid_col") = centroid_col,
        Named("ierr") = ierr
    );
}

/**
 * Perform normalization pipeline
 */
// [[Rcpp::export]]
List tox_normalization_pipeline_rcpp(NumericMatrix input, IntegerVector group_s, IntegerVector group_c) {
    int n_genes = input.nrow();
    int n_tissues = input.ncol();
    int n_grps = group_s.size();

   int max_stack = static_cast<int>(std::ceil(std::log2(n_genes)) + 10);

    NumericMatrix buf_stddev(n_genes, n_tissues);
    NumericMatrix buf_quant(n_genes, n_tissues);
    NumericMatrix buf_avg(n_genes, n_grps);
    NumericMatrix buf_log(n_genes, n_grps);
    NumericVector temp_col(n_genes);
    NumericVector rank_means(n_genes);
    IntegerVector perm(n_genes);
    IntegerVector stack_left(max_stack);
    IntegerVector stack_right(max_stack);
    int ierr = 0;


     normalization_pipeline_c(n_genes, n_tissues, input.begin(),
                                     buf_stddev.begin(), buf_quant.begin(),
                                     buf_avg.begin(), buf_log.begin(),
                                     temp_col.begin(), rank_means.begin(), perm.begin(),
                                     stack_left.begin(), stack_right.begin(), max_stack,
                                     group_s.begin(), group_c.begin(), n_grps, &ierr);

    return List::create(
        Named("buf_stddev") = buf_stddev,
        Named("buf_quant") = buf_quant,
        Named("buf_avg") = buf_avg,
        Named("buf_log") = buf_log,
        Named("rank_means") = rank_means,
        Named("perm") = perm,
        Named("ierr") = ierr
    );
}




// ===================================================================
// OUTLIER DETECTION WRAPPERS
// ===================================================================

// [[Rcpp::export]]
List tox_compute_family_scaling_rcpp(NumericVector distances, IntegerVector gene_to_fam, int n_families) {
  int n_genes = distances.size();
  NumericVector dscale(n_families);
  NumericVector loess_x(n_families);
  NumericVector loess_y(n_families);
  IntegerVector indices_used(n_families);
  int ierr = 0;

  compute_family_scaling_c(
    n_genes, n_families,
    distances.begin(),
    gene_to_fam.begin(),
    dscale.begin(),
    loess_x.begin(),
    loess_y.begin(),
    indices_used.begin(),
    &ierr
  );

  return List::create(
    Named("dscale") = dscale,
    Named("loess_x") = loess_x,
    Named("loess_y") = loess_y,
    Named("indices_used") = indices_used,
    Named("ierr") = ierr
  );
}

// [[Rcpp::export]]
List tox_compute_family_scaling_expert_rcpp(int n_families, NumericVector distances, IntegerVector gene_to_fam,
                                            IntegerVector perm_tmp, IntegerVector stack_left_tmp, IntegerVector stack_right_tmp,
                                            NumericVector family_distances) {
  int n_genes = distances.size();
  NumericVector dscale(n_families);
  NumericVector loess_x(n_families);
  NumericVector loess_y(n_families);
  IntegerVector indices_used(n_families);
  int ierr = 0;

  compute_family_scaling_expert_c(
    n_genes, n_families,
    distances.begin(),
    gene_to_fam.begin(),
    dscale.begin(),
    loess_x.begin(), loess_y.begin(), indices_used.begin(),
    perm_tmp.begin(), stack_left_tmp.begin(), stack_right_tmp.begin(),
    family_distances.begin(),
    &ierr
  );

  return List::create(
    Named("dscale") = dscale,
    Named("loess_x") = loess_x,
    Named("loess_y") = loess_y,
    Named("indices_used") = indices_used,
    Named("perm_tmp") = perm_tmp,
    Named("stack_left_tmp") = stack_left_tmp,
    Named("stack_right_tmp") = stack_right_tmp,
    Named("family_distances") = family_distances,
    Named("ierr") = ierr
  );
}

// [[Rcpp::export]]
List tox_compute_rdi_rcpp(NumericVector distances, IntegerVector gene_to_fam, NumericVector dscale) {
  int n_genes = distances.size();
  int n_families = dscale.size();
  NumericVector rdi(n_genes);
  NumericVector sorted_rdi(n_genes);
  IntegerVector perm(n_genes);
  IntegerVector stack_left(n_genes);
  IntegerVector stack_right(n_genes);

  compute_rdi_c(
    n_genes, n_families,
    distances.begin(),
    gene_to_fam.begin(),
    dscale.begin(),
    rdi.begin(), sorted_rdi.begin(),
    perm.begin(), stack_left.begin(), stack_right.begin()
  );
  return List::create(
    Named("rdi") = rdi,
    Named("sorted_rdi") = sorted_rdi,
    Named("perm") = perm,
    Named("stack_left") = stack_left,
    Named("stack_right") = stack_right
  );
}

// [[Rcpp::export]]
List tox_identify_outliers_rcpp(NumericVector rdi, double percentile) {
  int n_genes = rdi.size();
  NumericVector sorted_rdi = clone(rdi);
  std::sort(sorted_rdi.begin(), sorted_rdi.end());
  IntegerVector is_outlier_int(n_genes);
  double threshold = 0.0;

  identify_outliers_c(
    n_genes,
    rdi.begin(), sorted_rdi.begin(),
    is_outlier_int.begin(),
    &threshold,
    percentile
  );

  // Convert integer 0/1 flags to logical vector for R
  LogicalVector is_outlier(n_genes);
  for (int i = 0; i < n_genes; ++i) {
    is_outlier[i] = (is_outlier_int[i] != 0);
  }

  return List::create(
    Named("is_outlier") = is_outlier,
    Named("threshold") = threshold
  );
}

// [[Rcpp::export]]
List tox_detect_outliers_rcpp(NumericVector distances, IntegerVector gene_to_fam, int n_families, double percentile) {
  int n_genes = distances.size();
  NumericVector work_array(n_genes);
  IntegerVector perm(n_genes);
  IntegerVector stack_left(n_genes);
  IntegerVector stack_right(n_genes);
  IntegerVector is_outlier_int(n_genes);
  NumericVector loess_x(n_families);
  NumericVector loess_y(n_families);
  IntegerVector loess_n(n_families);
  int ierr = 0;

  detect_outliers_c(
    n_genes, n_families,
    distances.begin(), gene_to_fam.begin(),
    work_array.begin(),
    perm.begin(), stack_left.begin(), stack_right.begin(),
    is_outlier_int.begin(),
    loess_x.begin(), loess_y.begin(), loess_n.begin(),
    &ierr,
    percentile
  );

  // Convert integer flags to logical vector for R
  LogicalVector is_outlier(n_genes);
  for (int i = 0; i < n_genes; ++i) {
    is_outlier[i] = (is_outlier_int[i] != 0);
  }

  return List::create(
    Named("is_outlier") = is_outlier,
    Named("loess_x") = loess_x,
    Named("loess_y") = loess_y,
    Named("loess_n") = loess_n,
    Named("ierr") = ierr
  );
}

// [[Rcpp::export]]
List tox_group_centroid_rcpp(NumericMatrix expression_vectors, IntegerVector gene_to_family, int n_families, IntegerVector ortholog_set, String mode) {
    int n_axes = expression_vectors.nrow();
    int n_genes = expression_vectors.ncol();

    NumericMatrix centroid_matrix(n_axes, n_families);
    IntegerVector selected_indices(n_genes);
    int selected_indices_len = n_genes;
    int ierr = 0;

    group_centroid_c(expression_vectors.begin(), n_axes, n_genes,
                     gene_to_family.begin(), n_families,
                     centroid_matrix.begin(), mode.get_cstring(),
                     ortholog_set.begin(), selected_indices.begin(),
                     selected_indices_len, &ierr);

    return List::create(
        Named("centroid_matrix") = centroid_matrix,
        Named("ierr") = ierr
    );
}

// [[Rcpp::export]]
Rcpp::List compute_gene_means_rcpp(Rcpp::List expr_list) {

    int n_studies = expr_list.size();

    // Determine max dimensions across all studies
    int max_n_genes_all_studies = 0;
    int max_n_reps_all_studies  = 0;

    for (int s = 0; s < n_studies; ++s) {
        Rcpp::NumericMatrix mat = expr_list[s];
        max_n_reps_all_studies  = std::max(max_n_reps_all_studies,  mat.nrow());
        max_n_genes_all_studies = std::max(max_n_genes_all_studies, mat.ncol());
    }

    // Output matrix: (max_n_genes_all_studies x n_studies)
    Rcpp::NumericMatrix means(max_n_genes_all_studies, n_studies);

    int ierr = 0;

    for (int s = 0; s < n_studies; ++s) {

        Rcpp::NumericMatrix expr = expr_list[s];
        int n_reps  = expr.nrow();
        int n_genes = expr.ncol();

        // Pointer to the start of column s in the output matrix
        double* means_col = means.begin() + s * max_n_genes_all_studies;

        compute_gene_means_c(
            expr.begin(),
            &n_genes,
            &n_reps,
            means_col,                     // direct slice
            &max_n_genes_all_studies,
            &ierr
        );
    }

    return Rcpp::List::create(
        Rcpp::Named("means") = means,
        Rcpp::Named("max_n_genes_all_studies") = max_n_genes_all_studies,
        Rcpp::Named("max_n_reps_all_studies") = max_n_reps_all_studies,
        Rcpp::Named("ierr") = ierr
    );
}

// [[Rcpp::export]]
Rcpp::List compute_residuals_rcpp(
    Rcpp::List expr_list,
    Rcpp::NumericMatrix means,   // from compute_gene_means_rcpp
    int max_n_reps_all_studies
) {
    int n_studies = expr_list.size();

    // Determine max dimensions across all studies
    int max_n_genes_all_studies = means.nrow();

    // Allocate 3D output: (max_n_reps, max_n_genes, n_studies)
    Rcpp::NumericVector resid(
        max_n_reps_all_studies *
        max_n_genes_all_studies *
        n_studies
    );

    resid.attr("dim") = Rcpp::IntegerVector::create(
        max_n_reps_all_studies,
        max_n_genes_all_studies,
        n_studies
    );

    int ierr = 0;

    // Process each study
    for (int s = 0; s < n_studies; ++s) {

        Rcpp::NumericMatrix expr = expr_list[s];
        int n_reps  = expr.nrow();
        int n_genes = expr.ncol();

        // Means slice for this study (column s)
        double* means_col = means.begin() + s * max_n_genes_all_studies;

        // Residuals slice for this study (3D block)
        double* resid_slice =
            resid.begin() +
            s * (max_n_reps_all_studies * max_n_genes_all_studies);

        compute_residuals_c(
            expr.begin(),
            &n_genes,
            &n_reps,
            means_col,
            &max_n_genes_all_studies,
            &max_n_reps_all_studies,
            resid_slice,
            &ierr
        );
    }

    return Rcpp::List::create(
        Rcpp::Named("residuals") = resid,
        Rcpp::Named("max_n_genes_all_studies") = max_n_genes_all_studies,
        Rcpp::Named("max_n_reps_all_studies") = max_n_reps_all_studies,
        Rcpp::Named("ierr") = ierr
    );
}

// [[Rcpp::export]]
Rcpp::List determine_js_comp_test_n_points_n_neighbors_alloc_rcpp(
    Rcpp::NumericVector residuals,   // 3D: [max_n_reps, max_n_genes, n_studies]
    Rcpp::NumericMatrix gene_means,  // [max_n_genes, n_studies]
    int n_bootstraps,
    String join_method,         // "min", "max", "median"
    int min_count_per_mean_bin = 5,
    double min_neighbor_overlap = 0.1,
    double succeeding_ci_overlap = 0.9,
    int random_seed = 42,
    double two_sided_bootstrapping_significance_level = 2.5,
    double residual_range_quantile = 95.0
) {
    // dims(residuals) = [max_n_reps_all_studies, max_n_genes_all_studies, n_studies]
    Rcpp::IntegerVector dims = residuals.attr("dim");
    int max_n_reps_all_studies  = dims[0];
    int max_n_genes_all_studies = dims[1];
    int n_studies               = dims[2];

    int n_points  = 0;
    int n_neighbors = 0;
    double shared_residual_range = 0.0;
    int n_bins = 0;
    int ierr = 0;

    Rcpp::NumericMatrix best_ci(2, n_studies);

    determine_js_comp_test_n_points_n_neighbors_alloc_c(
        &n_points,
        &n_neighbors,
        residuals.begin(),
        &max_n_reps_all_studies,
        &max_n_genes_all_studies,
        &shared_residual_range,
        &n_bins,
        gene_means.begin(),
        &n_studies,
        &n_bootstraps,
        best_ci.begin(),
        join_method.get_cstring(),
        &min_count_per_mean_bin,
        &min_neighbor_overlap,
        &succeeding_ci_overlap,
        &two_sided_bootstrapping_significance_level,
        &random_seed,
        &residual_range_quantile,
        &ierr
    );

    return Rcpp::List::create(
        Rcpp::Named("n_points") = n_points,
        Rcpp::Named("n_neighbors") = n_neighbors,
        Rcpp::Named("shared_residual_range") = shared_residual_range,
        Rcpp::Named("n_bins") = n_bins,
        Rcpp::Named("best_candidate_pair_confidence_interval") = best_ci,
        Rcpp::Named("ierr") = ierr
    );
}

// [[Rcpp::export]]
Rcpp::List js_comp_test_alloc_rcpp(
    Rcpp::NumericVector residuals,   // 3D: [max_n_reps, max_n_genes, n_studies]
    Rcpp::NumericMatrix gene_means,  // [max_n_genes, n_studies]
    double shared_residual_range,
    int n_bins,
    int n_points,
    int n_neighbors,
    int n_permutations,
    int random_seed
) {
    // dims(residuals) = [max_n_reps_all_studies, max_n_genes_all_studies, n_studies]
    Rcpp::IntegerVector dims = residuals.attr("dim");
    int max_n_reps_all_studies  = dims[0];
    int max_n_genes_all_studies = dims[1];
    int n_studies               = dims[2];

    // Outputs
    Rcpp::NumericVector x_star(n_points);
    int n_pool = 0;

    Rcpp::IntegerVector neighborhood_ranges(
        2 * n_points * n_studies
    );
    neighborhood_ranges.attr("dim") =
        Rcpp::IntegerVector::create(2, n_points, n_studies);

    Rcpp::IntegerVector neighborhood_residuals(
        n_neighbors * n_points * n_studies
    );
    neighborhood_residuals.attr("dim") =
        Rcpp::IntegerVector::create(n_neighbors, n_points, n_studies);

    Rcpp::NumericVector pmfs(
        n_bins * n_points * n_studies
    );
    pmfs.attr("dim") =
        Rcpp::IntegerVector::create(n_bins, n_points, n_studies);

    Rcpp::IntegerVector counts(
        n_bins * n_points * n_studies
    );
    counts.attr("dim") =
        Rcpp::IntegerVector::create(n_bins, n_points, n_studies);

    Rcpp::IntegerVector included_n_reps(
        n_points * n_studies
    );
    included_n_reps.attr("dim") =
        Rcpp::IntegerVector::create(n_points, n_studies);

    Rcpp::NumericVector mean_pmf(
        n_bins * n_points
    );
    mean_pmf.attr("dim") =
        Rcpp::IntegerVector::create(n_bins, n_points);

    Rcpp::IntegerVector mean_pmf_counts(
        n_bins * n_points
    );
    mean_pmf_counts.attr("dim") =
        Rcpp::IntegerVector::create(n_bins, n_points);

    Rcpp::IntegerVector mean_pmf_included_n_reps(
        n_points
    );

    Rcpp::NumericVector js_divergences(
        n_points * n_studies
    );
    js_divergences.attr("dim") =
        Rcpp::IntegerVector::create(n_points, n_studies);

    Rcpp::NumericVector weights(
        n_points * n_studies
    );
    weights.attr("dim") =
        Rcpp::IntegerVector::create(n_points, n_studies);

    Rcpp::NumericVector global_js_divergence(n_studies);
    Rcpp::NumericVector p_values(n_studies);

    int ierr = 0;

    js_comp_test_alloc_c(
        gene_means.begin(),
        &max_n_genes_all_studies,
        &n_studies,
        residuals.begin(),
        &shared_residual_range,
        &n_bins,
        &max_n_reps_all_studies,
        x_star.begin(),
        &n_pool,
        &n_points,
        &n_neighbors,
        neighborhood_ranges.begin(),
        neighborhood_residuals.begin(),
        pmfs.begin(),
        counts.begin(),
        included_n_reps.begin(),
        mean_pmf.begin(),
        mean_pmf_counts.begin(),
        mean_pmf_included_n_reps.begin(),
        js_divergences.begin(),
        weights.begin(),
        global_js_divergence.begin(),
        p_values.begin(),
        &ierr,
        &n_permutations,
        &random_seed
    );

    return Rcpp::List::create(
        Rcpp::Named("x_star") = x_star,
        Rcpp::Named("n_pool") = n_pool,
        Rcpp::Named("neighborhood_ranges") = neighborhood_ranges,
        Rcpp::Named("neighborhood_residuals") = neighborhood_residuals,
        Rcpp::Named("pmfs") = pmfs,
        Rcpp::Named("counts") = counts,
        Rcpp::Named("included_n_reps") = included_n_reps,
        Rcpp::Named("mean_pmf") = mean_pmf,
        Rcpp::Named("mean_pmf_counts") = mean_pmf_counts,
        Rcpp::Named("mean_pmf_included_n_reps") = mean_pmf_included_n_reps,
        Rcpp::Named("js_divergences") = js_divergences,
        Rcpp::Named("weights") = weights,
        Rcpp::Named("global_js_divergence") = global_js_divergence,
        Rcpp::Named("p_values") = p_values,
        Rcpp::Named("ierr") = ierr
    );
}
