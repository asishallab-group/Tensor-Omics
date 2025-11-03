# Load all TensorOmics helper functions
source("r/tensoromics_functions.R")
source("r/experiments/helper_functions.R")


experiment <- "cardamine_original"
outliers_or_complete <- "outliers"  # "outliers" or "complete"


results_dir <- file.path("results", experiment, outliers_or_complete)  
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}
material_dir <- file.path("material", experiment)  

modes = c("all","orthologs") # all or ortholog mode
stds = c("knn") # "std" or "knn"
source_copy_analysis <- TRUE
calculate_fc <- FALSE
control_pattern <- "dietM"
condition_patterns <- c("dietP")

# === Full normalization pipeline ===

# Load raw expression data
input_file <- file.path(material_dir, "values.tsv")  
filtered_fam <- fread(file.path(material_dir, "filtered_families.tsv"))
global_pairs <- fread(file.path(material_dir, "filtered_orthologs.tsv"))
tandems <- fread(file.path(material_dir, "tandem.tsv"), sep = "\t", header = TRUE)
annotations <- read.table(file.path(material_dir, "prot-scriber_result.txt"), header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)

# Get unique genes from genids column (can be comma-separated)
tandem_protein_set <- unique(trimws(unlist(strsplit(tandems$genids, ","))))

if(source_copy_analysis) {
    paralog_pairs <- fread(file.path(material_dir, "filtered_paralogs.tsv"), header=TRUE)

    # Initialize distance vector
    source_copy_distances <- rep(NA_real_, nrow(paralog_pairs))
}

# Extract species names (all columns except the first one, which is Orthogroup)
species_names <- colnames(filtered_fam)[-1]

# Create a dictionary: species -> gene vector
species_genes_dict <- lapply(species_names, function(species) {
  unique(trimws(unlist(strsplit(filtered_fam[[species]], ", "))))
})

names(species_genes_dict) <- species_names



for (m in modes) {
  for (s in stds) {
    mode = m
    std = s

    print(paste0("Running normalization and analysis for mode=", mode, " and std=", std))

    df <- read.table(input_file, header = TRUE, sep = "\t")
    all_zero_genes <- apply(df[,-1], 1, function(x) all(x == 0))
    df <- df[!all_zero_genes, ]
    all_na_genes <- apply(df[,-1], 1, function(x) all(is.na(x)))
    df <- df[!all_na_genes, ]

    # Prepare matrix for processing (removing the gene ID column)
    gene_ids <- df[[1]]              # Save gene identifiers
    col_names <- colnames(df)[-1]    
    df_matrix <- as.matrix(df[,-1])   # Convert the expression values into a matrix

    # === Diagnose data quality before normalization ===
    diagnostics <- tox_diagnose_data_quality(df_matrix)

    # Clean data if there are problems
    if (diagnostics$problems$na_count > 0 || diagnostics$problems$inf_count > 0 || diagnostics$problems$nan_count > 0) {
      
      # Clean the data using our cleaning function
      df_matrix_clean <- tox_clean_data_for_normalization(
        df_matrix,
        remove_all_zero_genes = TRUE,
        na_strategy = "impute_mean",  # Impute NA with gene means instead of removing genes
        min_expression_threshold = 0.0,  # Don't set a threshold for TPM data
        convert_small_to_zero = FALSE    # Preserve small TPM values
      )
      
      # Update gene_ids to match cleaned matrix
      gene_ids <- gene_ids[1:nrow(df_matrix_clean)]
      
    } else {
      df_matrix_clean <- df_matrix
    }

    # === Apply normalization steps sequentially ===
    if(std=="knn") {
      normalized_matrix_std <- tox_normalize_by_knn_smoothed_std(df_matrix_clean,1)    # Normalize by standard deviation
    } else {
      normalized_matrix_std <- tox_normalize_by_std_dev(df_matrix_clean)    # Normalize by standard deviation
    }

    normalized_matrix_qtl <- tox_quantile_normalization(normalized_matrix_std)  # Apply quantile normalization
    norm_mat <- normalized_matrix_qtl

    if(!grepl("cardamine", tolower(experiment))) {
        norm_mat <- tox_calculate_tissue_averages(normalized_matrix_qtl)         # Average replicates by tissue
    }

    normalized_matrix_log <- tox_log2_transformation(norm_mat)     # Log2(x+1) transformation

    if (calculate_fc) {
      # === Calculate fold changes between specified groups ===
      fc_df <- tox_calculate_fc_by_patterns(
        df = normalized_matrix_log,
        control_pattern = control_pattern,           # Define control pattern
        condition_patterns = c(condition_patterns)       # Define condition pattern(s)
      )
      normalized_matrix_log <- fc_df
    }


    # Create final normalized dataframe with averages
    normalized_df <- data.frame(gene_id = gene_ids, normalized_matrix_log)

    # Use real tissue names
    tissue_names <- colnames(normalized_df)[-1]
    unique_fams <- sort(unique(filtered_fam$Orthogroup))
    fam_indices <- setNames(seq_along(unique_fams), unique_fams)


    fam_map <- list()
    for (i in 1:nrow(filtered_fam)) {
      fam_id <- filtered_fam$Orthogroup[i]
      fam_num <- fam_indices[[fam_id]]
      for (col in colnames(filtered_fam)[-1]) {
        prots <- unlist(strsplit(filtered_fam[[col]][i], ", "))
        for (p in prots) {
          fam_map[[p]] <- fam_num
        }
      }
    }
    gene_to_family <- sapply(gene_ids, function(g) if (!is.null(fam_map[[g]])) fam_map[[g]] else NA)

    n_families <- length(unique_fams)

    # === Build ortholog_set ===
    all_orthologs <- unique(c(global_pairs$Gene1, global_pairs$Gene2))
    ortholog_set <- gene_ids %in% all_orthologs

    # Classify genes: if not ortholog, then paralog
    paralog_set <- !ortholog_set

    # === Calculate centroids by family ===
    expression_vectors <- as.matrix(normalized_df[,-1])

    # Filter genes without family (without NA in gene_to_family)
    valid_idx <- which(!is.na(gene_to_family))
    expression_vectors <- as.matrix(normalized_df[valid_idx, -1]) # genes x tissues

    # Transpose to get axes x genes
    expression_vectors <- t(expression_vectors)

    gene_to_family <- gene_to_family[valid_idx]
    ortholog_set <- ortholog_set[valid_idx]
    gene_ids <- gene_ids[valid_idx]
    # n_families is NOT recalculated here, remains the same
    n_axes <- nrow(expression_vectors)

    # Calculate centroids using all arguments
    centroids <- tox_group_centroid(expression_vectors, as.integer(gene_to_family), n_families, ortholog_set, mode = mode)

    ###########################################################################################
    ## DISTANCE ANALYSIS
    ###########################################################################################

    results_dist <- expression_analysis(expression_vectors, tissue_names, n_axes, n_families, "euclidean_distance", annotations, species_genes_dict, outliers_or_complete, gene_ids, centroids, gene_to_family, ortholog_set, tandem_protein_set)

    data_dist <- results_dist$data
    axis_contribs_df_dist <- results_dist$axis_contribs
    histogram_dist <- results_dist$histogram
    outliers_dist <- results_dist$outliers
    word_freq_dist <- results_dist$term_freq
    outlier_idx <- results_dist$outlier_idx

    plots_dist <- make_all_analysis_plots(
        data = data_dist,
        axis_contribs_df = axis_contribs_df_dist,
        histogram = histogram_dist,
        outliers = outliers_dist,
        experiment = experiment
        )

    

    ###########################################################################################
    ## ANGLE ANALYSIS
    ###########################################################################################

    results_angle <- expression_analysis(expression_vectors, tissue_names, n_axes, n_families, "angle", annotations, species_genes_dict, outliers_or_complete, gene_ids, centroids, gene_to_family, ortholog_set, tandem_protein_set)

    data_angle <- results_angle$data
    axis_contribs_df_angle <- results_angle$axis_contribs
    histogram_angle <- results_angle$histogram
    outliers_angle <- results_angle$outliers
    word_freq_angle <- results_angle$term_freq
    outlier_idx_angle <- results_angle$outlier_idx

    plots_angle <- make_all_analysis_plots(
        data = data_angle,
        axis_contribs_df = axis_contribs_df_angle,
        histogram = histogram_angle,
        outliers = outliers_angle,
        experiment = experiment
        )

    ###########################################################################################
    ## SOURCE COPY ANALYSIS
    ###########################################################################################

    if(source_copy_analysis) {

        sc_results_dist <- analyze_source_copy_data(expression_vectors, paralog_pairs, fam_indices, n_families, experiment, annotations, species_genes_dict, gene_ids, "euclidean_distance", outliers_or_complete, tandems)

        sc_data_dist <- sc_results_dist$data
        sc_axis_contribs_df_dist <- sc_results_dist$axis_contribs
        sc_histogram_dist <- sc_results_dist$histogram
        sc_outliers_dist <- sc_results_dist$outliers
        sc_word_freq_dist <- sc_results_dist$term_freq
        sc_outlier_idx <- sc_results_dist$outlier_idx

        sc_plots_dist <- make_all_analysis_plots(
            data = sc_data_dist,
            axis_contribs_df = sc_axis_contribs_df_dist,
            histogram = sc_histogram_dist,
            outliers = sc_outliers_dist,
            experiment = experiment
            )
    

        ###########################################################################################
        ## ANGLE ANALYSIS
        ###########################################################################################
        sc_results_angle <- analyze_source_copy_data(expression_vectors, paralog_pairs, fam_indices, n_families, experiment, annotations, species_genes_dict, gene_ids, "angle", outliers_or_complete, tandems)

        sc_data_angle <- sc_results_angle$data
        sc_axis_contribs_df_angle <- sc_results_angle$axis_contribs
        sc_histogram_angle <- sc_results_angle$histogram
        sc_outliers_angle <- sc_results_angle$outliers
        sc_word_freq_angle <- sc_results_angle$term_freq
        sc_outlier_idx_angle <- sc_results_angle$outlier_idx

        sc_plots_angle <- make_all_analysis_plots(
            data = sc_data_angle,
            axis_contribs_df = sc_axis_contribs_df_angle,
            histogram = sc_histogram_angle,
            outliers = sc_outliers_angle,
            experiment = experiment
            )


    }


    outlier_ids_dist <- gene_ids[outlier_idx]
    outlier_ids_angle <- gene_ids[outlier_idx_angle]
    jaccard <- length(intersect(outlier_ids_dist, outlier_ids_angle)) / length(union(outlier_ids_dist, outlier_ids_angle))
        
    title_plot <- paste0(outliers_or_complete, " analysis - ", experiment, " - ", mode)
    pdf(file.path(results_dir, paste0(experiment, "_", outliers_or_complete, "_", mode, ".pdf")), width = 10, height = 7)

    grid.newpage()
    grid.text(title_plot, x = 0.5, y = 0.5, gp = gpar(fontsize = 24, fontface = "bold"))
    grid.newpage()
    grid.text("Distance analysis", x = 0.5, y = 0.5, gp = gpar(fontsize = 24, fontface = "bold"))

    print(plots_dist$rdi_histogram)
    print(plots_dist$distance_general)
    print(plots_dist$distance_species)
    print(plots_dist$tissue_versatility_general)
    print(plots_dist$tissue_versatility_species)
    print(plots_dist$change_in_versatility_general)
    print(plots_dist$change_in_versatility_species)
    print(plots_dist$change_in_preference_general)
    print(plots_dist$change_in_preference_species)
    print(plots_dist$axis_contribs)

    wordcloud(words = names(word_freq_dist),
              freq = word_freq_dist,
              scale = c(3, 0.5),
              max.words = 70,
              random.order = FALSE,
              colors = brewer.pal(8, "Dark2"))
    
    grid.newpage()
    grid.text("Angle Analysis", x = 0.5, y = 0.5, gp = gpar(fontsize = 24, fontface = "bold"))
    if(outliers_or_complete == "outliers") {
      grid.newpage()
      venn.plot <- venn.diagram(
      x = list(
          Distance = outlier_ids_dist,
          Angle = outlier_ids_angle
      ),
      filename = NULL,
      fill = c("cornflowerblue", "darkorange"),
      alpha = 0.5,
      cex = 2,
      cat.cex = 2,
      cat.pos = c(-20, 20)
      )
      grid.draw(venn.plot)
      grid.text(sprintf("Jaccard similarity: %.3f", jaccard), x = 0.5, y = 0.1, gp = gpar(fontsize = 18, col = "black"))
    } 
    print(plots_angle$rdi_histogram)
    print(plots_angle$distance_general)
    print(plots_angle$distance_species)
    print(plots_angle$tissue_versatility_general)
    print(plots_angle$tissue_versatility_species)
    print(plots_angle$change_in_versatility_general)
    print(plots_angle$change_in_versatility_species)
    print(plots_angle$change_in_preference_general)
    print(plots_angle$change_in_preference_species)
    print(plots_angle$axis_contribs)

    wordcloud(words = names(word_freq_angle),
              freq = word_freq_angle,
              scale = c(3, 0.5),
              max.words = 70,
              random.order = FALSE,
              colors = brewer.pal(8, "Dark2"))


    if (source_copy_analysis) {
      # --- Function to create bidirectional pairs ---
make_bidirectional_pairs <- function(df) {
  # Direct pairs
  direct <- paste(df$gene_id_source, df$gene_id_copy, sep = "_")
  # Inverse pairs
  inverse <- paste(df$gene_id_copy, df$gene_id_source, sep = "_")
  # Combine both (duplicating data for both directions)
  c(direct, inverse)
}

# --- Generate pairs for both analyses ---
pairs_dist  <- make_bidirectional_pairs(sc_data_dist)
pairs_angle <- make_bidirectional_pairs(sc_data_angle)

# --- Make them unique (regardless of internal gene order) ---
canonicalize <- function(pairs) {
  sapply(strsplit(pairs, "_"), function(x) paste(sort(x), collapse = "_"))
}

unique_pairs_dist  <- unique(canonicalize(pairs_dist))
unique_pairs_angle <- unique(canonicalize(pairs_angle))

# --- Calculate Jaccard index ---
sc_jaccard <- length(intersect(unique_pairs_dist, unique_pairs_angle)) /
              length(union(unique_pairs_dist, unique_pairs_angle))

cat("Number of unique pairs (dist):", length(unique_pairs_dist), "\n")
cat("Number of unique pairs (angle):", length(unique_pairs_angle), "\n")
cat("Jaccard similarity =", round(sc_jaccard, 3), "\n")


      grid.newpage()
      grid.text("Source Copy Analysis - Euclidean distance", x = 0.5, y = 0.5, gp = gpar(fontsize = 24, fontface = "bold"))

      print(sc_plots_dist$rdi_histogram)
      print(sc_plots_dist$distance_general)
      print(sc_plots_dist$distance_species)
      print(sc_plots_dist$change_in_versatility_general)
      print(sc_plots_dist$change_in_versatility_species)
      print(sc_plots_dist$change_in_preference_general)
      print(sc_plots_dist$change_in_preference_species)
      print(sc_plots_dist$axis_contribs)

        wordcloud(words = names(sc_word_freq_dist),
                freq = sc_word_freq_dist,
                scale = c(3, 0.5),
                max.words = 70,
                random.order = FALSE,
                colors = brewer.pal(8, "Dark2"))

      grid.newpage()
      grid.text("Source Copy Analysis - Angle", x = 0.5, y = 0.5, gp = gpar(fontsize = 24, fontface = "bold"))
      
      if(outliers_or_complete == "outliers") {
        grid.newpage()
        venn.plot <- venn.diagram(
        x = list(
            Distance = unique_pairs_dist,
            Angle = unique_pairs_angle
        ),
        filename = NULL,
        fill = c("cornflowerblue", "darkorange"),
        alpha = 0.5,
        cex = 2,
        cat.cex = 2,
        cat.pos = c(-20, 20)
        )
        grid.draw(venn.plot)
        grid.text(sprintf("Jaccard similarity: %.3f", sc_jaccard), x = 0.5, y = 0.1, gp = gpar(fontsize = 18, col = "black"))
      } 

      print(sc_plots_angle$rdi_histogram)
      print(sc_plots_angle$distance_general)
      print(sc_plots_angle$distance_species)
      print(sc_plots_angle$change_in_versatility_general)
      print(sc_plots_angle$change_in_versatility_species)
      print(sc_plots_angle$change_in_preference_general)
      print(sc_plots_angle$change_in_preference_species)
      print(sc_plots_angle$axis_contribs)

        wordcloud(words = names(sc_word_freq_angle),
                freq = sc_word_freq_angle,
                scale = c(3, 0.5),
                max.words = 70,
                random.order = FALSE,
                colors = brewer.pal(8, "Dark2"))
    } 


    dev.off()


  }

}