#include <Rcpp.h>
using namespace Rcpp;


 

// ===================================================================
// FORTRAN FUNCTIONS
// ===================================================================

extern "C" {
    void cluster_factor_trajectories_k_means_c(
      int *n_clusters,
      double *trajectories,
      int *n_factors,
      int *n_samples,
      int *n_timepoints,
      double *centroids,
      int *labels,
      int *label_counts,
      int *ierr,
      int *max_iterations
    );
// [[Rcpp::export]]
List tox_cluster_factor_trajectories_k_means_rcpp(int n_clusters,
                                 NumericVector trajectories,
                                 int n_factors,
                                 int n_samples,
                                 int n_timepoints,
                                 NumericMatrix centroids,
                                 int max_iterations = 300) {
    int n_points = n_samples * n_timepoints;
    IntegerVector labels(n_points);
    IntegerVector label_counts(n_clusters);
    int ierr = 0;
    // Fortran expects column-major, R matrices are column-major
    cluster_factor_trajectories_k_means_c(&n_clusters,
                        trajectories.begin(),
                        &n_factors,
                        &n_samples,
                        &n_timepoints,
                        centroids.begin(),
                        labels.begin(),
                        label_counts.begin(),
                        &ierr,
                        &max_iterations);
    centroids.attr("dim") = Dimension(n_factors, n_clusters);
    return List::create(Named("centroids") = centroids,
              Named("labels") = labels,
              Named("label_counts") = label_counts,
              Named("ierr") = ierr);
}

    void k_means_clustering_c(
    int *n_clusters,
    double *data_points,
    int *n_points,
    int *n_dims,
    double *centroids,
    int *labels,
    int *label_counts,
    int *ierr,
    int *max_iterations
  );



    void cluster_factor_trajectories_k_means_c(
    int *n_clusters,
    double *trajectories,
    int *n_factors,
    int *n_samples,
    int *n_timepoints,
    double *centroids,
    int *labels,
    int *label_counts,
    int *ierr,
    int *max_iterations
  );
  void k_means_clustering_c(
    int *n_clusters,
    double *data_points,
    int *n_points,
    int *n_dims,
    double *centroids,
    int *labels,
    int *label_counts,
    int *ierr,
    int *max_iterations
  );

    void relative_axes_changes_from_shift_vector_c(
      double* vec, int* n_axes, double* contributions, int* ierr
    );
    void relative_axes_expression_from_expression_vector_c(
      double* vec, int* n_axes, double* contributions, int* ierr
    );
    void clock_hand_angle_between_vectors_c(
      double* v1, double* v2, int* n_dims, double* signed_angle, int* selected_axes_for_signed, int* ierr
    );
    void clock_hand_angles_for_shift_vectors_c(
      double* origins, double* targets, int* n_dims, int* n_vecs,
      int* vecs_selection_mask, int* n_selected_vecs, int* selected_axes_for_signed,
      double* signed_angles, int* ierr
    );
   
 void omics_vector_RAP_projection_c(
  const double* vecs,
  int n_axes,
  int n_vecs,
  const int* vecs_selection_mask,
  int n_selected_vecs,
  const int* axes_selection_mask,
  int n_selected_axes,
  double* projections,
  int* ierr
);

 
void omics_field_RAP_projection_c(
  const double* vecs,
  int n_axes,
  int n_vecs,
  const int* vecs_selection_mask,
  int n_selected_vecs,
  const int* axes_selection_mask,
  int n_selected_axes,
  double* projections,
  int* ierr
);


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

void build_kd_index_C(
      double* points,
      int num_dimensions,
      int num_points,
      int* kd_indices,
      int* dimension_order,
      int* workspace,
      double* value_buffer,
      int* permutation,
      int* left_stack,
      int* right_stack,
      int* ierr
    );   
void which_c(
            int* mask,
            int n,
            int* idx_out,
            int m_max,
            int* m_out,
            int* ierr
    );

void loess_smooth_2d_c(
            int n_total,
            int n_target,
            double* x_ref,
            double* y_ref,
            int* indices_used,
            int n_used,
            double* x_query,
            double kernel_sigma,
            double kernel_cutoff,
            double* y_out,
            int* ierr
    );
         
void deserialize_int_nd_C(
        int* arr, int arr_size, 
        int* filename_ascii, int fn_len, int* ierr);

void deserialize_real_nd_C(
        double* arr, int arr_size, 
        int* filename_ascii, int fn_len, int* ierr);

void deserialize_char_nd_C(
        int* ascii_arr, int clen, 
        int total_array_size, int* filename_ascii, int fn_len, int* ierr);

void serialize_int_nd_C(
        void* arr, int* dims, 
        int ndim, int* filename_ascii, int fn_len, int* ierr);
void serialize_real_nd_C(
        void* arr, int* dims, 
        int ndim, int* filename_ascii, int fn_len, int* ierr);

void serialize_char_nd_C(
        int* ascii_arr, int* dims, 
        int ndim, int clen, int* filename_ascii, int fn_len, int* ierr);
void get_array_metadata_C(
      const int* filename_ascii,
      int fn_len, int* dims_out, int* dims_out_capacity,
      int* ndims, int* ierr, int* clen);

void build_bst_index_C(
        const double* values, 
        int num_values, int* sorted_indices, int* left_stack, 
        int* right_stack, int* ierr);

void bst_range_query_C(
        const double* values, 
        const int* sorted_indices, int num_values, 
        double lower_bound, double upper_bound, int* output_indices, 
        int* num_matches, int* ierr);

void build_spherical_kd_C(
        const double* vectors, int num_dimensions, 
        int num_vectors, int* sphere_indices, int* dimension_order, 
        int* workspace, double* value_buffer, int* permutation, 
        int* left_stack, int* right_stack, int* ierr);
                                 

void read_expression_vectors_tsv_C(const char* file_list_raw,
                                   int file_list_len,
                                   int n_files,
                                   const char* gene_ids_raw,
                                   int gene_ids_len,
                                   int n_genes,
                                   double* expression_vectors,
                                   int n_samples,
                                   int n_header_rows,
                                   int gene_col,
                                   const int* value_cols,
                                   int n_value_cols,
                                   int* ierr,
                                   const char* delimiter_raw);

void read_gene_ids_from_tsv_file_C(const char* filename_raw,
                                   int fn_len,
                                   char* gene_ids_raw,
                                   int gene_ids_len,
                                   int n_genes,
                                   int n_header_rows,
                                   int gene_col,
                                   int* ierr);


void read_orthofinder_file_C(const char* filename_raw,
                             int fn_len,
                             const char* gene_ids_raw,
                             int gene_ids_len,
                             int n_genes,
                             char* family_ids_raw,
                             int family_ids_len,
                             int n_families,
                             int* gene_to_fam,
                             int* ierr);
        
void filter_unassigned_genes_C(const char* gene_ids_raw,
                               int gene_ids_len,
                               int n_genes,
                               const int* gene_to_fam,
                               int* mask,
                               int* n_genes_kept,
                               int* ierr);

void validate_data_structure_C(int n_genes,
                               int n_families,
                               int n_samples,
                               const char* gene_ids_raw,
                               int gene_ids_len,
                               const char* gene_family_ids_raw,
                               int fam_len,
                               const int* gene_to_fam,
                               const double* expression_vectors,
                               const double* family_centroids,
                               const double* shift_vectors,
                               int* ierr);

void validate_gene_to_family_mapping_C(const int* gene_to_fam,
                                       int n_genes,
                                       int n_families,
                                       int* ierr);

void validate_expression_data_C(const double* expression_vectors,
                                int n_genes,
                                int n_samples,
                                int check_non_negative,
                                int* ierr);

void validate_family_centroids_C(const double* family_centroids,
                                 int n_families,
                                 int n_samples,
                                 int* ierr);

void validate_shift_vectors_C(const double* shift_vectors,
                              const double* expression_vectors,
                              const double* family_centroids,
                              const int* gene_to_fam,
                              int n_genes,
                              int n_samples,
                              int n_families,
                              int* ierr);

void validate_string_array_uniqueness_C(const char* str_arr,
                                        int str_len,
                                        int n_strings,
                                        int* ierr);

void validate_all_data_C(int n_genes,
                         int n_families,
                         int n_samples,
                         const char* gene_ids_raw,
                         int gene_len,
                         const char* gene_family_ids_raw,
                         int fam_len,
                         const int* gene_to_fam,
                         const double* expression_vectors,
                         const double* family_centroids,
                         const double* shift_vectors,
                         int* ierr);

void create_zip_archive_c(const char* zip_filename,
                                  int zip_len,
                                  const char* keys,
                                  int keys_len,
                                  int keys_count,
                                  const char* filenames,
                                  int filenames_len,
                                  int filenames_count,
                                  int* ierr);

void extract_zip_archive_c(const char* zip_filename,
                                   int filename_len,
                                   int* ierr);

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


// [[Rcpp::export]]
List tox_k_means_clustering_rcpp(int n_clusters,
                                 NumericMatrix data_points,
                                 int n_points,
                                 int n_dims,
                                 NumericMatrix centroids,
                                 int max_iterations = 300) {
    IntegerVector labels(n_points);
    IntegerVector label_counts(n_clusters);
    int ierr = 0;
    // Fortran expects column-major, R matrices are column-major
    k_means_clustering_c(&n_clusters,
                        data_points.begin(),
                        &n_points,
                        &n_dims,
                        centroids.begin(),
                        labels.begin(),
                        label_counts.begin(),
                        &ierr,
                        &max_iterations);
    centroids.attr("dim") = Dimension(n_dims, n_clusters);
    return List::create(Named("centroids") = centroids,
              Named("labels") = labels,
              Named("label_counts") = label_counts,
              Named("ierr") = ierr);
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
 * Calculate mean vector
 */
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

// ===================================================================
// [[Rcpp::export]]
List tox_build_kd_index_rcpp(NumericMatrix X, IntegerVector dim_order) {
  int d = X.nrow();
    int n = X.ncol();
    
    IntegerVector kd_ix(n);
    IntegerVector work(n);
    NumericVector subarray(n);
    IntegerVector perm(n);
    IntegerVector stack_left(n);
    IntegerVector stack_right(n);
    int ierr = 0;
    
  build_kd_index_C(X.begin(), d, n, kd_ix.begin(), dim_order.begin(), work.begin(), subarray.begin(), perm.begin(), stack_left.begin(), stack_right.begin(), &ierr);
        return List::create(
        Named("kd_ix") = kd_ix,
        Named("ierr") = ierr
    );
}

// [[Rcpp::export]]
List tox_which_rcpp(IntegerVector mask,int n,int m_max) {

        IntegerVector idx_out(m_max);
        int m_out = 0;
        int ierr = 0;
        which_c(mask.begin(), n, idx_out.begin(), m_max, &m_out, &ierr);
        return List::create(
        Named("idx_out") = idx_out,
        Named("m_out") = m_out,
        Named("ierr") = ierr
        );
}

// [[Rcpp::export]]
List tox_loess_smooth_2d_rcpp(int n_total,
int n_target,
NumericVector x_ref,
NumericVector y_ref,
IntegerVector indices_used,
int n_used,
NumericVector x_query,
double kernel_sigma,
double kernel_cutoff) {
        NumericVector y_out(n_target);
        int ierr = 0;
        loess_smooth_2d_c(n_total, n_target, x_ref.begin(), y_ref.begin(), indices_used.begin(), n_used, x_query.begin(), kernel_sigma, kernel_cutoff, y_out.begin(), &ierr);
        return List::create(
        Named("y_out") = y_out,
        Named("ierr") = ierr
        );
}


static std::vector<int> filename_to_ascii(const std::string &filename) {
  std::vector<int> v;
  for (unsigned char c : filename) v.push_back((int)c);
  return v;
}

// --- Serialize wrappers (match R signatures: arr, filename) ---

// [[Rcpp::export]]
List tox_serialize_int_array_rcpp(IntegerVector arr, std::string filename) {
  IntegerVector dim = arr.hasAttribute("dim") ? as<IntegerVector>(arr.attr("dim")) : IntegerVector::create((int)arr.size());
  int ndim = dim.size();
  std::vector<int> dims(ndim);
  for (int i = 0; i < ndim; ++i) dims[i] = dim[i];

  auto fname = filename_to_ascii(filename);
  int fn_len = static_cast<int>(fname.size());
  int ierr = 0;
  serialize_int_nd_C((void*)arr.begin(), dims.data(), ndim, fname.data(), fn_len, &ierr);
  return List::create(Named("ierr") = ierr);
}

// [[Rcpp::export]]
List tox_serialize_real_array_rcpp(NumericVector arr, std::string filename) {
  IntegerVector dim = arr.hasAttribute("dim") ? as<IntegerVector>(arr.attr("dim")) : IntegerVector::create((int)arr.size());
  int ndim = dim.size();
  std::vector<int> dims(ndim);
  for (int i = 0; i < ndim; ++i) dims[i] = dim[i];

  auto fname = filename_to_ascii(filename);
  int fn_len = static_cast<int>(fname.size());
  int ierr = 0;
  serialize_real_nd_C((void*)arr.begin(), dims.data(), ndim, fname.data(), fn_len, &ierr);
  return List::create(Named("ierr") = ierr);
}

// [[Rcpp::export]]
List tox_serialize_char_array_rcpp(CharacterVector carr, std::string filename) {
  IntegerVector dim = carr.hasAttribute("dim") ? as<IntegerVector>(carr.attr("dim")) : IntegerVector::create((int)carr.size());
  int ndim = dim.size();
  std::vector<int> dims(ndim);
  for (int i = 0; i < ndim; ++i) dims[i] = dim[i];

  int total = 1;
  for (int i = 0; i < ndim; ++i) total *= dims[i];

  int clen = 0;
  for (int i = 0; i < total; ++i) {
    std::string s = as<std::string>(carr[i]);
    if ((int)s.size() > clen) clen = (int)s.size();
  }
  if (clen == 0) clen = 1;

  std::vector<int> ascii_flat(total * clen, 0);
  for (int i = 0; i < total; ++i) {
    std::string s = as<std::string>(carr[i]);
    for (int j = 0; j < (int)s.size() && j < clen; ++j) ascii_flat[i*clen + j] = (unsigned char)s[j];
  }

  auto fname = filename_to_ascii(filename);
  int fn_len = static_cast<int>(fname.size());
  int ierr = 0;
  serialize_char_nd_C(ascii_flat.data(), dims.data(), ndim, clen, fname.data(), fn_len, &ierr);
  return List::create(Named("ierr") = ierr);
}

// --- Deserialize wrappers (match R signatures: filename, max_dims=5) ---

// [[Rcpp::export]]
List tox_deserialize_int_array_rcpp(std::string filename, int max_dims = 5) {
  auto fname = filename_to_ascii(filename);
  int fn_len = static_cast<int>(fname.size());
  std::vector<int> dims_out(max_dims);
  int ndims = 0;
  int ierr = 0;
  int clen = 0;
  get_array_metadata_C(fname.data(), fn_len, dims_out.data(), &max_dims, &ndims, &ierr, &clen);
  if (ierr != 0) return List::create(Named("ierr") = ierr);

  int total = 1;
  for (int i = 0; i < ndims; ++i) total *= dims_out[i];

  IntegerVector out(total);
  deserialize_int_nd_C(out.begin(), total, fname.data(), fn_len, &ierr);
  return List::create(Named("values") = out, Named("dims") = IntegerVector(dims_out.begin(), dims_out.begin()+ndims), Named("ndim") = ndims, Named("ierr") = ierr);
}

// [[Rcpp::export]]
List tox_deserialize_real_array_rcpp(std::string filename, int max_dims = 5) {
  auto fname = filename_to_ascii(filename);
  int fn_len = static_cast<int>(fname.size());
  std::vector<int> dims_out(max_dims);
  int ndims = 0;
  int ierr = 0;
  int clen = 0;
  get_array_metadata_C(fname.data(), fn_len, dims_out.data(), &max_dims, &ndims, &ierr, &clen);
  if (ierr != 0) return List::create(Named("ierr") = ierr);

  int total = 1;
  for (int i = 0; i < ndims; ++i) total *= dims_out[i];

  NumericVector out(total);
  deserialize_real_nd_C(out.begin(), total, fname.data(), fn_len, &ierr);
  return List::create(Named("values") = out, Named("dims") = IntegerVector(dims_out.begin(), dims_out.begin()+ndims), Named("ndim") = ndims, Named("ierr") = ierr);
}

// [[Rcpp::export]]
List tox_deserialize_char_array_rcpp(std::string filename, int max_dims = 5) {
  auto fname = filename_to_ascii(filename);
  int fn_len = static_cast<int>(fname.size());
  std::vector<int> dims_out(max_dims);
  int ndims = 0;
  int ierr = 0;
  int clen = 0;
  get_array_metadata_C(fname.data(), fn_len, dims_out.data(), &max_dims, &ndims, &ierr, &clen);
  if (ierr != 0) return List::create(Named("ierr") = ierr);

  int total = 1;
  for (int i = 0; i < ndims; ++i) total *= dims_out[i];

  std::vector<int> ascii_out(total * clen);
  deserialize_char_nd_C(ascii_out.data(), clen, total, fname.data(), fn_len, &ierr);

  CharacterVector out(total);
  for (int i = 0; i < total; ++i) {
    std::string s;
    for (int j = 0; j < clen; ++j) {
      int ch = ascii_out[i*clen + j];
      if (ch == 0) break;
      s.push_back((char)ch);
    }
    out[i] = s;
  }
  return List::create(Named("values") = out, Named("dims") = IntegerVector(dims_out.begin(), dims_out.begin()+ndims), Named("ndim") = ndims, Named("ierr") = ierr);
}

// [[Rcpp::export]]
List tox_get_array_metadata_rcpp(std::string filename, int dims_out_capacity = 5, bool with_clen = false) {
  auto fname = filename_to_ascii(filename);
  int fn_len = static_cast<int>(fname.size());
  IntegerVector dims_res(dims_out_capacity);
  int ndims = 0;
  int ierr = 0;
  int clen = 0;
  get_array_metadata_C(fname.data(), fn_len, dims_res.begin(), &dims_out_capacity, &ndims, &ierr, &clen);
  if (with_clen) {
    return List::create(Named("dims") = dims_res, Named("ndim") = ndims, Named("clen") = clen, Named("ierr") = ierr);
  }
  return List::create(Named("dims") = dims_res, Named("ndim") = ndims, Named("ierr") = ierr);
}



// [[Rcpp::export]]
IntegerVector build_bst_index_rcpp(NumericVector values) {
    int num_values = static_cast<int>(values.size());
    IntegerVector sorted_indices(num_values);
    IntegerVector left_stack(num_values);
    IntegerVector right_stack(num_values);
    int ierr = 0;

    build_bst_index_C(values.begin(), num_values, sorted_indices.begin(), left_stack.begin(), right_stack.begin(), &ierr);

    // Return sorted indices (1-based Fortran indices preserved)
    return sorted_indices;
}


// [[Rcpp::export]]
List bst_range_query_rcpp(NumericVector values, IntegerVector sorted_indices, double lower_bound, double upper_bound) {
    int num_values = static_cast<int>(values.size());
    IntegerVector output_indices(num_values);
    int num_matches = 0;
    int ierr = 0;

    bst_range_query_C(values.begin(), sorted_indices.begin(), num_values, lower_bound, upper_bound, output_indices.begin(), &num_matches, &ierr);

    return List::create(Named("output_indices") = output_indices,
                                            Named("num_matches") = num_matches,
                                            Named("ierr") = ierr);
}

// [[Rcpp::export]]
List tox_build_spherical_kd_rcpp(NumericMatrix V,IntegerVector dim_order) {
    int d = V.nrow();
    int n = V.ncol();
    
    IntegerVector sphere_ix(n);
    IntegerVector work(n);
    NumericVector subarray(n);
    IntegerVector perm(n);
    IntegerVector stack_left(n);
    IntegerVector stack_right(n);
    int ierr = 0;
    

    build_spherical_kd_C(V.begin(), d, n, sphere_ix.begin(), dim_order.begin(), work.begin(), subarray.begin(), perm.begin(), stack_left.begin(), stack_right.begin(), &ierr);

    return List::create(Named("sphere_ix") = sphere_ix,
                        Named("ierr") = ierr
    );}

// ===================================================================
// TOX DATA RCPP FORWARDERS
// ===================================================================


// [[Rcpp::export]]
List tox_read_expression_vectors_tsv_rcpp(RawMatrix file_list_raw,
                                          RawMatrix gene_ids_raw,
                                          IntegerVector value_cols,
                                          RawVector delimiter_raw,
                                          int n_samples,
                                          int n_header_rows,
                                          int gene_col) {
    int file_list_len = file_list_raw.nrow();
    int n_files = file_list_raw.ncol();
    int gene_ids_len = gene_ids_raw.nrow();
    int n_genes = gene_ids_raw.ncol();
    int n_value_cols = value_cols.size();

    NumericMatrix expression_vectors(n_samples, n_genes);
    int ierr = 0;

    read_expression_vectors_tsv_C(reinterpret_cast<const char*>(file_list_raw.begin()),
                                  file_list_len,
                                  n_files,
                                  reinterpret_cast<const char*>(gene_ids_raw.begin()),
                                  gene_ids_len,
                                  n_genes,
                                  expression_vectors.begin(),
                                  n_samples,
                                  n_header_rows,
                                  gene_col,
                                  value_cols.begin(),
                                  n_value_cols,
                                  &ierr,
                                  reinterpret_cast<const char*>(delimiter_raw.begin()));

    return List::create(Named("expression_vectors") = expression_vectors,
                        Named("ierr") = ierr);
}

// [[Rcpp::export]]
List tox_read_gene_ids_from_tsv_file_rcpp(RawVector filename_raw,
                                          int n_genes,
                                          int gene_ids_len,
                                          int n_header_rows,
                                          int gene_col) {
    RawMatrix gene_ids_raw(gene_ids_len, n_genes);
    std::fill(gene_ids_raw.begin(), gene_ids_raw.end(), static_cast<Rbyte>(0));

    int ierr = 0;

    read_gene_ids_from_tsv_file_C(reinterpret_cast<const char*>(filename_raw.begin()),
                                  static_cast<int>(filename_raw.size()),
                                  reinterpret_cast<char*>(gene_ids_raw.begin()),
                                  gene_ids_len,
                                  n_genes,
                                  n_header_rows,
                                  gene_col,
                                  &ierr);

    return List::create(Named("gene_ids_raw") = gene_ids_raw,
                        Named("ierr") = ierr);
}

// [[Rcpp::export]]
List tox_read_orthofinder_file_rcpp(RawVector filename_raw,
                                    RawMatrix gene_ids_raw,
                                    int n_families,
                                    int family_ids_len) {
    int gene_ids_len = gene_ids_raw.nrow();
    int n_genes = gene_ids_raw.ncol();

    RawMatrix family_ids_raw(family_ids_len, n_families);
    std::fill(family_ids_raw.begin(), family_ids_raw.end(), static_cast<Rbyte>(0));
    IntegerVector gene_to_fam(n_genes);
    int ierr = 0;

    read_orthofinder_file_C(reinterpret_cast<const char*>(filename_raw.begin()),
                            static_cast<int>(filename_raw.size()),
                            reinterpret_cast<const char*>(gene_ids_raw.begin()),
                            gene_ids_len,
                            n_genes,
                            reinterpret_cast<char*>(family_ids_raw.begin()),
                            family_ids_len,
                            n_families,
                            gene_to_fam.begin(),
                            &ierr);

    return List::create(Named("family_ids_raw") = family_ids_raw,
                        Named("gene_to_fam") = gene_to_fam,
                        Named("ierr") = ierr);
}

// [[Rcpp::export]]
List tox_validate_data_structure_rcpp(RawMatrix gene_ids_raw,
                                      RawMatrix gene_family_ids_raw,
                                      IntegerVector gene_to_fam,
                                      NumericMatrix expression_vectors,
                                      NumericMatrix family_centroids,
                                      NumericMatrix shift_vectors) {
    int n_genes = gene_ids_raw.ncol();
    int gene_ids_len = gene_ids_raw.nrow();
    int n_families = gene_family_ids_raw.ncol();
    int fam_len = gene_family_ids_raw.nrow();
    int n_samples = expression_vectors.nrow();
    int ierr = 0;

    validate_data_structure_C(n_genes,
                              n_families,
                              n_samples,
                              reinterpret_cast<const char*>(gene_ids_raw.begin()),
                              gene_ids_len,
                              reinterpret_cast<const char*>(gene_family_ids_raw.begin()),
                              fam_len,
                              gene_to_fam.begin(),
                              expression_vectors.begin(),
                              family_centroids.begin(),
                              shift_vectors.begin(),
                              &ierr);

    return List::create(Named("ierr") = ierr);
}

// [[Rcpp::export]]
List tox_filter_unassigned_genes_rcpp(RawMatrix gene_ids_raw,
                                      IntegerVector gene_to_fam) {
    int gene_ids_len = gene_ids_raw.nrow();
    int n_genes = gene_ids_raw.ncol();
    IntegerVector mask(n_genes);
    int n_genes_kept = 0;
    int ierr = 0;

    filter_unassigned_genes_C(reinterpret_cast<const char*>(gene_ids_raw.begin()),
                              gene_ids_len,
                              n_genes,
                              gene_to_fam.begin(),
                              mask.begin(),
                              &n_genes_kept,
                              &ierr);

    return List::create(Named("mask") = mask,
                        Named("n_genes_kept") = n_genes_kept,
                        Named("ierr") = ierr);
}

// [[Rcpp::export]]
List tox_validate_gene_to_family_mapping_rcpp(IntegerVector gene_to_fam,
                                              int n_families) {
    int n_genes = gene_to_fam.size();
    int ierr = 0;

    validate_gene_to_family_mapping_C(gene_to_fam.begin(),
                                      n_genes,
                                      n_families,
                                      &ierr);

    return List::create(Named("ierr") = ierr);
}

// [[Rcpp::export]]
List tox_validate_expression_data_rcpp(NumericMatrix expression_vectors,
                                       bool check_non_negative) {
    int n_samples = expression_vectors.nrow();
    int n_genes = expression_vectors.ncol();
    int flag = check_non_negative ? 1 : 0;
    int ierr = 0;

    validate_expression_data_C(expression_vectors.begin(),
                               n_genes,
                               n_samples,
                               flag,
                               &ierr);

    return List::create(Named("ierr") = ierr);
}

// [[Rcpp::export]]
List tox_validate_family_centroids_rcpp(NumericMatrix family_centroids) {
    int n_samples = family_centroids.nrow();
    int n_families = family_centroids.ncol();
    int ierr = 0;

    validate_family_centroids_C(family_centroids.begin(),
                                 n_families,
                                 n_samples,
                                 &ierr);

    return List::create(Named("ierr") = ierr);
}

// [[Rcpp::export]]
List tox_validate_shift_vectors_rcpp(NumericMatrix shift_vectors,
                                     NumericMatrix expression_vectors,
                                     NumericMatrix family_centroids,
                                     IntegerVector gene_to_fam) {
    int n_samples = expression_vectors.nrow();
    int n_genes = expression_vectors.ncol();
    int n_families = family_centroids.ncol();
    int ierr = 0;

    validate_shift_vectors_C(shift_vectors.begin(),
                             expression_vectors.begin(),
                             family_centroids.begin(),
                             gene_to_fam.begin(),
                             n_genes,
                             n_samples,
                             n_families,
                             &ierr);

    return List::create(Named("ierr") = ierr);
}

// [[Rcpp::export]]
List tox_validate_string_array_uniqueness_rcpp(RawMatrix string_arr_raw) {
    int str_len = string_arr_raw.nrow();
    int n_strings = string_arr_raw.ncol();
    int ierr = 0;

    validate_string_array_uniqueness_C(reinterpret_cast<const char*>(string_arr_raw.begin()),
                                       str_len,
                                       n_strings,
                                       &ierr);

    return List::create(Named("ierr") = ierr);
}

// [[Rcpp::export]]
List tox_validate_all_data_rcpp(RawMatrix gene_ids_raw,
                                RawMatrix gene_family_ids_raw,
                                IntegerVector gene_to_fam,
                                NumericMatrix expression_vectors,
                                NumericMatrix family_centroids,
                                NumericMatrix shift_vectors) {
    int n_genes = gene_ids_raw.ncol();
    int gene_len = gene_ids_raw.nrow();
    int n_families = gene_family_ids_raw.ncol();
    int fam_len = gene_family_ids_raw.nrow();
    int n_samples = expression_vectors.nrow();
    int ierr = 0;

    validate_all_data_C(n_genes,
                        n_families,
                        n_samples,
                        reinterpret_cast<const char*>(gene_ids_raw.begin()),
                        gene_len,
                        reinterpret_cast<const char*>(gene_family_ids_raw.begin()),
                        fam_len,
                        gene_to_fam.begin(),
                        expression_vectors.begin(),
                        family_centroids.begin(),
                        shift_vectors.begin(),
                        &ierr);

    return List::create(Named("ierr") = ierr);
}

// [[Rcpp::export]]
List tox_create_zip_archive_generic_rcpp(RawVector zip_filename_raw,
                                         RawMatrix keys_raw,
                                         RawMatrix filenames_raw) {
    int zip_len = zip_filename_raw.size();
    int keys_len = keys_raw.nrow();
    int keys_count = keys_raw.ncol();
    int filenames_len = filenames_raw.nrow();
    int filenames_count = filenames_raw.ncol();
    int ierr = 0;

    create_zip_archive_c(reinterpret_cast<const char*>(zip_filename_raw.begin()),
                                 zip_len,
                                 reinterpret_cast<const char*>(keys_raw.begin()),
                                 keys_len,
                                 keys_count,
                                 reinterpret_cast<const char*>(filenames_raw.begin()),
                                 filenames_len,
                                 filenames_count,
                                 &ierr);

    return List::create(Named("ierr") = ierr);
}

// [[Rcpp::export]]
List tox_extract_zip_archive_generic_rcpp(RawVector zip_filename_raw) {
    int zip_len = zip_filename_raw.size();
    int ierr = 0;

    extract_zip_archive_c(reinterpret_cast<const char*>(zip_filename_raw.begin()),
                                  zip_len,
                                  &ierr);

    return List::create(Named("ierr") = ierr);
}
// ===================================================================

// [[Rcpp::export]]
Rcpp::List tox_omics_vector_RAP_projection_rcpp(Rcpp::NumericMatrix vecs,
                                               Rcpp::IntegerVector vecs_selection_mask,
                                               Rcpp::IntegerVector axes_selection_mask) {
  int n_axes = vecs.nrow();
  int n_vecs = vecs.ncol();

  // ---- Check mask lengths ----
  if (vecs_selection_mask.size() != n_vecs) {
    Rcpp::stop("vecs_selection_mask length must be ncol(vecs)=%d, got %d",
               n_vecs, vecs_selection_mask.size());
  }
  if (axes_selection_mask.size() != n_axes) {
    Rcpp::stop("axes_selection_mask length must be nrow(vecs)=%d, got %d",
               n_axes, axes_selection_mask.size());
  }

  // ---- Check mask contents and count selected ----
  int n_selected_vecs = 0;
  for (int j = 0; j < n_vecs; ++j) {
    int v = vecs_selection_mask[j];
    if (v == NA_INTEGER) {
      Rcpp::stop("vecs_selection_mask contains NA at position %d", j + 1);
    }
    if (v != 0 && v != 1) {
      Rcpp::stop("vecs_selection_mask must be 0/1; got %d at position %d", v, j + 1);
    }
    if (v == 1) ++n_selected_vecs;
  }

  int n_selected_axes = 0;
  for (int i = 0; i < n_axes; ++i) {
    int a = axes_selection_mask[i];
    if (a == NA_INTEGER) {
      Rcpp::stop("axes_selection_mask contains NA at position %d", i + 1);
    }
    if (a != 0 && a != 1) {
      Rcpp::stop("axes_selection_mask must be 0/1; got %d at position %d", a, i + 1);
    }
    if (a == 1) ++n_selected_axes;
  }

  if (n_selected_vecs <= 0) Rcpp::stop("No vectors selected (vecs_selection_mask has no 1s).");
  if (n_selected_axes <= 0) Rcpp::stop("No axes selected (axes_selection_mask has no 1s).");

  // ---- Allocate output ----
  Rcpp::NumericMatrix projections(n_selected_axes, n_selected_vecs);
  int ierr = 0;

  // Copy masks for C/Fortran
  std::vector<int> vecs_mask(vecs_selection_mask.begin(), vecs_selection_mask.end());
  std::vector<int> axes_mask(axes_selection_mask.begin(), axes_selection_mask.end());

  // ---- Call native code ----
omics_vector_RAP_projection_c(
  vecs.begin(),
  n_axes,
  n_vecs,
  vecs_mask.data(),
  n_selected_vecs,
  axes_mask.data(),
  n_selected_axes,
  projections.begin(),
  &ierr
);
  return Rcpp::List::create(
    Rcpp::Named("projections") = projections,
    Rcpp::Named("ierr") = ierr
  );
}

// [[Rcpp::export]]
Rcpp::List tox_omics_field_RAP_projection_rcpp(Rcpp::NumericMatrix vecs,
                                              Rcpp::IntegerVector vecs_selection_mask,
                                              Rcpp::IntegerVector axes_selection_mask) {
  int n_rows = vecs.nrow();
  int n_vecs = vecs.ncol();

  // Field input must be (2*n_axes) x n_vecs
  if (n_rows % 2 != 0) {
    Rcpp::stop("Field vecs must have 2*n_axes rows; nrow(vecs)=%d is not even.", n_rows);
  }
  int n_axes = n_rows / 2;

  // ---- Check mask lengths ----
  if (vecs_selection_mask.size() != n_vecs) {
    Rcpp::stop("vecs_selection_mask length must be ncol(vecs)=%d, got %d",
               n_vecs, vecs_selection_mask.size());
  }
  if (axes_selection_mask.size() != n_axes) {
    Rcpp::stop("axes_selection_mask length must be n_axes=nrow(vecs)/2=%d, got %d",
               n_axes, axes_selection_mask.size());
  }

  // ---- Check mask contents and count selected ----
  int n_selected_vecs = 0;
  for (int j = 0; j < n_vecs; ++j) {
    int v = vecs_selection_mask[j];
    if (v == NA_INTEGER) {
      Rcpp::stop("vecs_selection_mask contains NA at position %d", j + 1);
    }
    if (v != 0 && v != 1) {
      Rcpp::stop("vecs_selection_mask must be 0/1; got %d at position %d", v, j + 1);
    }
    if (v == 1) ++n_selected_vecs;
  }

  int n_selected_axes = 0;
  for (int i = 0; i < n_axes; ++i) {
    int a = axes_selection_mask[i];
    if (a == NA_INTEGER) {
      Rcpp::stop("axes_selection_mask contains NA at position %d", i + 1);
    }
    if (a != 0 && a != 1) {
      Rcpp::stop("axes_selection_mask must be 0/1; got %d at position %d", a, i + 1);
    }
    if (a == 1) ++n_selected_axes;
  }

  if (n_selected_vecs <= 0) Rcpp::stop("No vectors selected (vecs_selection_mask has no 1s).");
  if (n_selected_axes <= 0) Rcpp::stop("No axes selected (axes_selection_mask has no 1s).");

  // ---- Allocate output ----
  Rcpp::NumericMatrix projections(n_selected_axes, n_selected_vecs);
  int ierr = 0;

  std::vector<int> vecs_mask(vecs_selection_mask.begin(), vecs_selection_mask.end());
  std::vector<int> axes_mask(axes_selection_mask.begin(), axes_selection_mask.end());

  // ---- Call native code ----
 omics_field_RAP_projection_c(
  vecs.begin(),
  n_axes,
  n_vecs,
  vecs_mask.data(),
  n_selected_vecs,
  axes_mask.data(),
  n_selected_axes,
  projections.begin(),
  &ierr
);
  return Rcpp::List::create(
    Rcpp::Named("projections") = projections,
    Rcpp::Named("ierr") = ierr
  );
}
