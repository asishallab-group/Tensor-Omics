#include <Rcpp.h>

using namespace Rcpp;


// ===================================================================
// FORTRAN FUNCTIONS
// ===================================================================

extern "C" {

    void fjct_compute_jsd_c( int* family_idx, int* gene_to_family_S1, int* gene_to_family_S2, int* n_genes_S1, int* n_genes_S2, double* neighborhood_residuals_S1, double* neighborhood_residuals_S2, int* neighborhood_genes_S1, int* neighborhood_genes_S2, int* n_reps_S1, int* n_reps_S2, int* n_neighbors, int* n_points, int* n_bins, double* shared_residual_range, double* js_divergences, int* included_n_reps_S1, int* included_n_reps_S2, int* total_included_n_reps, double* global_js_divergence, double* weights, int* ierr );

    void fjct_compute_jsd_expert_c( double* neighborhood_residuals_S1, double* neighborhood_residuals_S2, int* n_reps_S1, int* n_reps_S2, int* n_neighbors, int* n_points, int* neighbor_mask_S1, int* neighbor_mask_S2, int* n_bins, double* shared_residual_range, double* js_divergences, int* included_n_reps_S1, int* included_n_reps_S2, int* total_included_n_reps, double* global_js_divergence, double* weights, double* pmf_S1, double* pmf_S2, int* tmp_counts, int* ierr );

    void fjct_compute_contribution_scores_c( double* global_js_divergences, int* total_included_n_reps_per_f, int* k_families, double* support_weights, double* contribution_scores, int* ierr);

    void gjct_permutation_test_c( double* neighborhood_residuals_S1, double* neighborhood_residuals_S2, int* n_reps_S1, int* n_reps_S2, int* n_neighbors, int* n_points, double* global_jsd_observed, int* n_bins, double* shared_residual_range, int* n_permutations, double* jsd_null, double* p_value, int* ierr, int* random_seed );

    void compute_weighted_global_divergence_c( double* js_divergences, int* n_points, int* included_n_residuals_S1, int* included_n_residuals_S2, double* global_js_divergence, double* weights, int* ierr );

    void compute_divergence_per_reference_point_c( double* pmf_S1, double* pmf_S2, int* n_points, int* n_bins, double* js_divergences, int* ierr );

    void build_residual_histograms_c( double* neighborhood_residuals, int* n_reps, int* n_neighbors, int* n_points, double* shared_residual_range, int* n_bins, int* counts, double* pmf, int* included_n_residuals, int* ierr );

    void determine_shared_residual_range_c( double* neighborhood_residuals_S1, double* neighborhood_residuals_S2, int* n_reps_S1, int* n_reps_S2, int* n_neighbors, int* n_points, double* residual_range_quantile, double* shared_residual_range, int* ierr );

    void determine_shared_residual_range_expert_c( double* residual_pool, int* residual_pool_perm, int* n_pool, double* residual_range_quantile, double* shared_R, int* ierr );

    void compute_gene_means_c( int* n_genes, int* n_reps, double* expr, double* means, int* ierr);

    void compute_residuals_c( int* n_genes, int* n_reps, double* expr, double* means, double* resid, int* ierr);

    void pool_means_c( int* n_genes_S1, double* mean_S1, int* n_genes_S2, double* mean_S2, int* n_points, int* n_pool, double* x_star, int* ierr);

    void pool_means_expert_c( double* pooled_means, int* pooled_means_perm, int* pool_size, int* n_points, int* n_pool, double* x_star, int* ierr);

    void calc_neighborhood_size_c( int* n_pool, int* n_points, int* n_genes_S, double* mean_S, int* desired_size, int* n_neighbors, int* ierr);

    void construct_neighborhoods_c( int* n_points, double* x_star, int* n_genes_S, double* mean_S, int* n_reps_S, double* resid_S, double* neighborhood_residuals, int* neighborhood_indices, int* n_neighbors, int* ierr);

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
Rcpp::List tox_determine_shared_residual_range_expert_rcpp(
    Rcpp::NumericVector residual_pool,
    Rcpp::IntegerVector residual_pool_perm,
    double residual_range_quantile = 95.0
) {
    int pool_size = residual_pool.size();
    double shared_R = 0.0;
    int ierr = 0;

    determine_shared_residual_range_expert_c(
        residual_pool.begin(),
        residual_pool_perm.begin(),
        &pool_size,
        &residual_range_quantile,
        &shared_R,
        &ierr
    );

    return Rcpp::List::create(
        Rcpp::Named("shared_R") = shared_R,
        Rcpp::Named("ierr") = ierr
    );
}

// [[Rcpp::export]]
Rcpp::List tox_determine_shared_residual_range_rcpp(
    Rcpp::NumericVector neighborhood_residuals_S1,
    Rcpp::NumericVector neighborhood_residuals_S2,
    double residual_range_quantile = 95.0
) {
    Rcpp::IntegerVector dims = neighborhood_residuals_S1.attr("dim");
    int n_reps_S1 = dims[0];
    int n_neighbors = dims[1];
    int n_points = dims[2];
    dims = neighborhood_residuals_S2.attr("dim");
    int n_reps_S2 = dims[0];

    double shared_R = 0.0;
    int ierr = 0;

    determine_shared_residual_range_c(
        neighborhood_residuals_S1.begin(),
        neighborhood_residuals_S2.begin(),
        &n_reps_S1,
        &n_reps_S2,
        &n_neighbors,
        &n_points,
        &residual_range_quantile,
        &shared_R,
        &ierr
    );

    return Rcpp::List::create(
        Rcpp::Named("shared_R") = shared_R,
        Rcpp::Named("ierr") = ierr
    );
}

// [[Rcpp::export]]
Rcpp::List tox_build_residual_histograms_rcpp(
    Rcpp::NumericVector neighborhood_residuals,
    double shared_residual_range,
    int n_bins
) {
    Rcpp::IntegerVector dims = neighborhood_residuals.attr("dim");
    int n_reps = dims[0];
    int n_neighbors = dims[1];
    int n_points = dims[2];

    Rcpp::IntegerMatrix counts(n_points, n_bins);
    Rcpp::NumericMatrix pmf(n_points, n_bins);
    Rcpp::IntegerVector included_n_residuals(n_points);

    int ierr = 0;

    build_residual_histograms_c(
        neighborhood_residuals.begin(),
        &n_reps,
        &n_neighbors,
        &n_points,
        &shared_residual_range,
        &n_bins,
        counts.begin(),
        pmf.begin(),
        included_n_residuals.begin(),
        &ierr
    );

    return Rcpp::List::create(
        Rcpp::Named("counts") = counts,
        Rcpp::Named("pmf") = pmf,
        Rcpp::Named("included_n_residuals") = included_n_residuals,
        Rcpp::Named("ierr") = ierr
    );
}

// [[Rcpp::export]]
Rcpp::List tox_compute_divergence_per_reference_point_rcpp(
    Rcpp::NumericMatrix pmf_S1,
    Rcpp::NumericMatrix pmf_S2
) {
    int n_points = pmf_S1.nrow();
    int n_bins = pmf_S1.ncol();

    Rcpp::NumericVector js_divergences(n_points);
    int ierr = 0;

    compute_divergence_per_reference_point_c(
        pmf_S1.begin(),
        pmf_S2.begin(),
        &n_points,
        &n_bins,
        js_divergences.begin(),
        &ierr
    );

    return Rcpp::List::create(
        Rcpp::Named("js_divergences") = js_divergences,
        Rcpp::Named("ierr") = ierr
    );
}

// [[Rcpp::export]]
Rcpp::List tox_compute_weighted_global_divergence_rcpp(
    Rcpp::NumericVector js_divergences,
    Rcpp::IntegerVector included_n_residuals_S1,
    Rcpp::IntegerVector included_n_residuals_S2
) {
    int n_points = js_divergences.size();

    double global_jsd = 0.0;
    Rcpp::NumericVector weights(n_points);
    int ierr = 0;

    compute_weighted_global_divergence_c(
        js_divergences.begin(),
        &n_points,
        included_n_residuals_S1.begin(),
        included_n_residuals_S2.begin(),
        &global_jsd,
        weights.begin(),
        &ierr
    );

    return Rcpp::List::create(
        Rcpp::Named("global_js_divergence") = global_jsd,
        Rcpp::Named("weights") = weights,
        Rcpp::Named("ierr") = ierr
    );
}

// [[Rcpp::export]]
Rcpp::List tox_gjct_permutation_test_rcpp(
    Rcpp::NumericVector neighborhood_residuals_S1,
    Rcpp::NumericVector neighborhood_residuals_S2,
    double global_jsd_observed,
    int n_bins,
    double shared_residual_range,
    int n_permutations,
    int random_seed
) {
    Rcpp::IntegerVector dims = neighborhood_residuals_S1.attr("dim");
    int n_reps_S1 = dims[0];
    int n_neighbors = dims[1];
    int n_points = dims[2];
    dims = neighborhood_residuals_S2.attr("dim");
    int n_reps_S2 = dims[0];

    Rcpp::NumericVector jsd_null(n_permutations);
    double p_value = 0.0;
    int ierr = 0;

    gjct_permutation_test_c(
        neighborhood_residuals_S1.begin(),
        neighborhood_residuals_S2.begin(),
        &n_reps_S1,
        &n_reps_S2,
        &n_neighbors,
        &n_points,
        &global_jsd_observed,
        &n_bins,
        &shared_residual_range,
        &n_permutations,
        jsd_null.begin(),
        &p_value,
        &ierr,
        &random_seed
    );

    return Rcpp::List::create(
        Rcpp::Named("jsd_null") = jsd_null,
        Rcpp::Named("p_value")  = p_value,
        Rcpp::Named("ierr")     = ierr
    );
}


// [[Rcpp::export]]
Rcpp::List tox_calc_neighborhood_size_rcpp(int n_pool,
                                           int n_points,
                                           int n_genes_S,
                                           Rcpp::NumericVector mean_S,
                                           int desired_size = 0) {
    int n_neighbors = 0;
    int ierr = 0;

    calc_neighborhood_size_c(
        &n_pool,
        &n_points,
        &n_genes_S,
        mean_S.begin(),
        &desired_size,
        &n_neighbors,
        &ierr
    );

    return Rcpp::List::create(
        Rcpp::Named("n_neighbors") = n_neighbors,
        Rcpp::Named("ierr")        = ierr
    );
}

// [[Rcpp::export]]
Rcpp::List tox_construct_neighborhoods_rcpp(
        Rcpp::NumericVector x_star,
        int n_pool,
        Rcpp::NumericVector mean_S,
        Rcpp::NumericMatrix resid_S,
        int desired_n_neighbors = 0) {

    int n_points  = x_star.size();
    int n_genes_S = mean_S.size();
    int n_reps_S  = resid_S.nrow();
    int n_neighbors = 0;
    int ierr = 0;

    calc_neighborhood_size_c( &n_pool, &n_points, &n_genes_S, mean_S.begin(), &desired_n_neighbors, &n_neighbors, &ierr);

    if (ierr != 0)
    {
        Rcpp::NumericVector neigh_res(0);
        Rcpp::IntegerMatrix neigh_idx(0, n_points);

        return Rcpp::List::create(
            Rcpp::Named("neighborhood_residuals") = neigh_res,
            Rcpp::Named("neighborhood_indices")   = neigh_idx,
            Rcpp::Named("ierr")                   = ierr
        );
    }

    // Flat buffer
    Rcpp::NumericVector neigh_res(n_reps_S * n_neighbors * n_points);
    Rcpp::IntegerMatrix neigh_idx(n_neighbors, n_points);

    construct_neighborhoods_c(
        &n_points, x_star.begin(),
        &n_genes_S, mean_S.begin(),
        &n_reps_S, resid_S.begin(),
        neigh_res.begin(),
        neigh_idx.begin(),
        &n_neighbors,
        &ierr
    );

    // Convert to 3D array
    neigh_res.attr("dim") = Rcpp::IntegerVector::create(
        n_reps_S,
        n_neighbors,
        n_points
    );

    return Rcpp::List::create(
        Rcpp::Named("neighborhood_residuals") = neigh_res,
        Rcpp::Named("neighborhood_indices")   = neigh_idx,
        Rcpp::Named("ierr")                   = ierr
    );
}

// [[Rcpp::export]]
Rcpp::List tox_compute_gene_means_rcpp(Rcpp::NumericMatrix expr) {
    int n_reps  = expr.nrow();
    int n_genes = expr.ncol();

    Rcpp::NumericVector means(n_genes);
    int ierr = 0;

    compute_gene_means_c(
        &n_genes, &n_reps, expr.begin(),
        means.begin(),
        &ierr
    );

    return Rcpp::List::create(
        Rcpp::Named("means") = means,
        Rcpp::Named("ierr")  = ierr
    );
}

// [[Rcpp::export]]
Rcpp::List tox_compute_residuals_rcpp(Rcpp::NumericMatrix expr,
                                      Rcpp::NumericVector means) {
    int n_reps  = expr.nrow();
    int n_genes = expr.ncol();

    Rcpp::NumericMatrix resid(n_reps, n_genes);
    int ierr = 0;

    compute_residuals_c(
        &n_genes, &n_reps, expr.begin(), 
        means.begin(),
        resid.begin(),
        &ierr
    );

    return Rcpp::List::create(
        Rcpp::Named("resid") = resid,
        Rcpp::Named("ierr")  = ierr
    );
}

// [[Rcpp::export]]
Rcpp::List tox_pool_means_rcpp(Rcpp::NumericVector mean_S1,
                                     Rcpp::NumericVector mean_S2,
                                     int n_points) {
    int n_genes_S1 = mean_S1.size();
    int n_genes_S2 = mean_S2.size();

    Rcpp::NumericVector x_star(n_points);
    int n_pool = 0;
    int ierr = 0;

    pool_means_c(
        &n_genes_S1, mean_S1.begin(),
        &n_genes_S2, mean_S2.begin(),
        &n_points,
        &n_pool,
        x_star.begin(),
        &ierr
    );

    return Rcpp::List::create(
        Rcpp::Named("n_pool") = n_pool,
        Rcpp::Named("x_star") = x_star,
        Rcpp::Named("ierr")   = ierr
    );
}

// [[Rcpp::export]]
Rcpp::List tox_pool_means_expert_rcpp(Rcpp::NumericVector pooled_means,
                               Rcpp::IntegerVector pooled_perm,
                               int n_points) {

    Rcpp::NumericVector x_star(n_points);
    int pool_size = pooled_means.size();
    int n_pool = 0;
    int ierr = 0;

    pool_means_expert_c(
        pooled_means.begin(),
        pooled_perm.begin(),
        &pool_size,
        &n_points,
        &n_pool,
        x_star.begin(),
        &ierr
    );

    return Rcpp::List::create(
        Rcpp::Named("n_pool") = n_pool,
        Rcpp::Named("x_star") = x_star,
        Rcpp::Named("ierr")   = ierr
    );
}

// [[Rcpp::export]]
Rcpp::List tox_fjct_compute_jsd_alloc_rcpp(
    int family_idx,
    Rcpp::IntegerVector gene_to_family_S1,
    Rcpp::IntegerVector gene_to_family_S2,
    Rcpp::NumericVector neighborhood_residuals_S1,
    Rcpp::NumericVector neighborhood_residuals_S2,
    Rcpp::IntegerMatrix neighborhood_genes_S1,
    Rcpp::IntegerMatrix neighborhood_genes_S2,
    int n_bins,
    double shared_residual_range
) {
    Rcpp::IntegerVector dims = neighborhood_residuals_S1.attr("dim");
    int n_reps_S1 = dims[0];
    int n_neighbors = dims[1];
    int n_points = dims[2];
    dims = neighborhood_residuals_S2.attr("dim");
    int n_reps_S2 = dims[0];
    Rcpp::NumericVector jsd(n_points);
    Rcpp::IntegerVector inc1(n_points);
    Rcpp::IntegerVector inc2(n_points);
    int n_genes_S1 = gene_to_family_S1.size();
    int n_genes_S2 = gene_to_family_S2.size();
    int total_included = 0;
    double global_jsd = 0.0;
    Rcpp::NumericVector weights(n_points);
    int ierr = 0;

    fjct_compute_jsd_c(
        &family_idx,
        gene_to_family_S1.begin(),
        gene_to_family_S2.begin(),
        &n_genes_S1,
        &n_genes_S2,
        neighborhood_residuals_S1.begin(),
        neighborhood_residuals_S2.begin(),
        neighborhood_genes_S1.begin(),
        neighborhood_genes_S2.begin(),
        &n_reps_S1,
        &n_reps_S2,
        &n_neighbors,
        &n_points,
        &n_bins,
        &shared_residual_range,
        jsd.begin(),
        inc1.begin(),
        inc2.begin(),
        &total_included,
        &global_jsd,
        weights.begin(),
        &ierr
    );

    return Rcpp::List::create(
        Rcpp::Named("js_divergences") = jsd,
        Rcpp::Named("included_n_reps_S1") = inc1,
        Rcpp::Named("included_n_reps_S2") = inc2,
        Rcpp::Named("total_included_n_reps") = total_included,
        Rcpp::Named("global_js_divergence") = global_jsd,
        Rcpp::Named("weights") = weights,
        Rcpp::Named("ierr") = ierr
    );
}

// [[Rcpp::export]]
Rcpp::List tox_fjct_compute_jsd_expert_rcpp(
    Rcpp::NumericVector neighborhood_residuals_S1,
    Rcpp::NumericVector neighborhood_residuals_S2,
    Rcpp::IntegerMatrix neighbor_mask_S1,
    Rcpp::IntegerMatrix neighbor_mask_S2,
    int n_bins,
    double shared_residual_range
) {
    Rcpp::IntegerVector dims = neighborhood_residuals_S1.attr("dim");
    int n_reps_S1 = dims[0];
    int n_neighbors = dims[1];
    int n_points = dims[2];
    dims = neighborhood_residuals_S2.attr("dim");
    int n_reps_S2 = dims[0];
    Rcpp::NumericVector jsd(n_points);
    Rcpp::IntegerVector inc1(n_points);
    Rcpp::IntegerVector inc2(n_points);
    int total_included = 0;
    double global_jsd = 0.0;
    Rcpp::NumericVector weights(n_points);

    Rcpp::NumericMatrix pmf_S1(n_points, n_bins);
    Rcpp::NumericMatrix pmf_S2(n_points, n_bins);
    Rcpp::IntegerMatrix tmp_counts(n_points, n_bins);

    int ierr = 0;

    fjct_compute_jsd_expert_c(
        neighborhood_residuals_S1.begin(),
        neighborhood_residuals_S2.begin(),
        &n_reps_S1,
        &n_reps_S2,
        &n_neighbors,
        &n_points,
        neighbor_mask_S1.begin(),
        neighbor_mask_S2.begin(),
        &n_bins,
        &shared_residual_range,
        jsd.begin(),
        inc1.begin(),
        inc2.begin(),
        &total_included,
        &global_jsd,
        weights.begin(),
        pmf_S1.begin(),
        pmf_S2.begin(),
        tmp_counts.begin(),
        &ierr
    );

    return Rcpp::List::create(
        Rcpp::Named("js_divergences") = jsd,
        Rcpp::Named("included_n_reps_S1") = inc1,
        Rcpp::Named("included_n_reps_S2") = inc2,
        Rcpp::Named("total_included_n_reps") = total_included,
        Rcpp::Named("global_js_divergence") = global_jsd,
        Rcpp::Named("weights") = weights,
        Rcpp::Named("pmf_S1") = pmf_S1,
        Rcpp::Named("pmf_S2") = pmf_S2,
        Rcpp::Named("tmp_counts") = tmp_counts,
        Rcpp::Named("ierr") = ierr
    );
}

// [[Rcpp::export]]
Rcpp::List tox_fjct_compute_contribution_scores_rcpp(
    Rcpp::NumericVector global_js_divergences,
    Rcpp::IntegerVector total_included_n_reps_per_f
) {
    int k_families = global_js_divergences.size();
    Rcpp::NumericVector support_weights(k_families);
    Rcpp::NumericVector contribution_scores(k_families);
    int ierr = 0;

    fjct_compute_contribution_scores_c(
        global_js_divergences.begin(),
        total_included_n_reps_per_f.begin(),
        &k_families,
        support_weights.begin(),
        contribution_scores.begin(),
        &ierr
    );

    return Rcpp::List::create(
        Rcpp::Named("support_weights") = support_weights,
        Rcpp::Named("contribution_scores") = contribution_scores,
        Rcpp::Named("ierr") = ierr
    );
}
