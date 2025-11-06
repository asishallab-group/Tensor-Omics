
#include <Rcpp.h>

using namespace Rcpp;


// ===================================================================
// FORTRAN FUNCTIONS
// ===================================================================

extern "C" {
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
     void normalize_by_std_dev_c(int n_genes, int n_tissues,
                                double *input_matrix, double *output_matrix, int *ierr);

    void quantile_normalization_c(int n_genes, int n_tissues,
                                  double *input_matrix, double *output_matrix,
                                double *temp_col, double *rank_means,
                                  int *perm, int *stack_left, int *stack_right,
                                  int max_stack, int *ierr);

    void normalization_pipeline_c(int n_genes, int n_tissues,
                                  double *input_matrix, double *buf_stddev, double *buf_quant,
                                  double *buf_avg, double *buf_log, double *temp_col,
                                  double *rank_means, int *perm, int *stack_left,
                                  int *stack_right, int max_stack,
                                  int *group_s, int *group_c, int n_grps, int *ierr);

    void log2_transformation_c(int n_genes, int n_tissues,
                               double *input_matrix, double *output_matrix, int *ierr);

    void calc_tiss_avg_c(int n_genes, int n_grps,
                         int *group_s, int *group_c,
                         double *input_matrix, double *output_matrix, int *ierr);

    void calc_fchange_c(int n_genes, int n_cols, int n_pairs,
                        int *control_cols, int *cond_cols,
                        double *input_matrix, double *output_matrix, int *ierr);
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

// ===================================================================
// RCPP WRAPPERS
// ===================================================================

// [[Rcpp::export]]
List normalize_by_std_dev_rcpp(NumericMatrix input) {
    int n_genes = input.nrow();
    int n_tissues = input.ncol();

    NumericMatrix output(n_genes, n_tissues);
    int ierr = 0;

    // Call Fortran routine
    normalize_by_std_dev_c(n_genes, n_tissues,
                           input.begin(),
                           output.begin(),
                           &ierr);

    // Just return output + error code (no stop(), no validation)
    return List::create(
        Named("output") = output,
        Named("ierr") = ierr
    );
}


// [[Rcpp::export]]
List quantile_normalization_rcpp(NumericMatrix input, int max_stack) {
    int n_genes = input.nrow();
    int n_tissues = input.ncol();
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
        Named("output") = output,
        Named("rank_means") = rank_means,
        Named("perm") = perm,
        Named("ierr") = ierr
    );
}

// [[Rcpp::export]]
List tox_normalize_data_rcpp(NumericMatrix input,
                             IntegerVector group_s,
                             IntegerVector group_c,
                             int max_stack) {

    int n_genes = input.nrow();
    int n_tissues = input.ncol();
    int n_grps = group_s.size();

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

    normalization_pipeline_c(
        n_genes,
        n_tissues,
        input.begin(),
        buf_stddev.begin(),
        buf_quant.begin(),
        buf_avg.begin(),
        buf_log.begin(),
        temp_col.begin(),
        rank_means.begin(),
        perm.begin(),
        stack_left.begin(),
        stack_right.begin(),
        max_stack,
        group_s.begin(),
        group_c.begin(),
        n_grps,
        &ierr
    );

    return List::create(
        Named("buf_stddev") = buf_stddev,
        Named("buf_quant") = buf_quant,
        Named("buf_avg") = buf_avg,
        Named("buf_log")  = buf_log,
        Named("rank_means") = rank_means,
        Named("perm") = perm,
        Named("ierr") = ierr
    );
}


// [[Rcpp::export]]
List log2_transformation_rcpp(NumericMatrix input) {
    int n_genes = input.nrow();
    int n_tissues = input.ncol();
    NumericMatrix output(n_genes, n_tissues);
    int ierr = 0;

    
    return List::create(
        Named("output") = output,
        Named("ierr") = ierr
    );
}

// [[Rcpp::export]]
List calc_tiss_avg_rcpp(NumericMatrix input, IntegerVector group_s, IntegerVector group_c) {
    int n_gene = input.nrow();
    int n_grps = group_s.size();
    NumericMatrix output(n_gene, n_grps);
    int ierr = 0;

    calc_tiss_avg_c(
        n_gene,
        n_grps,
        group_s.begin(),
        group_c.begin(),
        input.begin(),
        output.begin(),
        &ierr
    );

    return List::create(
        Named("output") = output,
        Named("ierr") = ierr
    );
}

// [[Rcpp::export]]
List calc_fchange_rcpp(NumericMatrix input,
                       IntegerVector control_cols,
                       IntegerVector cond_cols) {

    int n_genes = input.nrow();
    int n_cols  = input.ncol();
    int n_pairs = control_cols.size();

    NumericMatrix output(n_genes, n_pairs);
    int ierr = 0;

    calc_fchange_c(
        n_genes,
        n_cols,
        n_pairs,
        control_cols.begin(),
        cond_cols.begin(),
        input.begin(),
        output.begin(),
        &ierr
    );

    return List::create(
        Named("output") = output,
        Named("ierr") = ierr
    );
}
