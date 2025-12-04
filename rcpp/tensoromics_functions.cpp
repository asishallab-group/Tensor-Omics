#include <Rcpp.h>

using namespace Rcpp;


// ===================================================================
// FORTRAN FUNCTIONS
// ===================================================================

extern "C" {

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
