# ==========================================================
# TensorOmics main analysis pipeline
# Performs normalization, distance and angle analyses, and
# optional source-copy comparisons across gene families.
# ==========================================================

# Load all TensorOmics helper functions
source("r/tensoromics_functions.R")
source("r/experiments/helper_functions.R")

# --- General experiment configuration ---
experiment <- "mammalian_sum"
outliers_or_complete <- "outliers"  # options: "outliers" or "complete"
type_test <- "t_test" # options: "t_test" or "wilcoxon"

# Create output directories
results_dir <- file.path("results")  
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}
material_dir <- file.path("material", experiment)  

# Define modes and normalization methods
modes = c("all", "orthologs")      # how centroids should be calculated
stds = c("knn")                    # "std" or "knn" for normalization
source_copy_analysis <- TRUE       # whether to run source-copy comparisons
calculate_fc <- FALSE              # whether to compute fold changes
control_pattern <- "dietM"         # only needed when calculate_fc is TRUE
condition_patterns <- c("dietP")   # only needed when calculate_fc is TRUE

# ==========================================================
# Full normalization and analysis pipeline
# ==========================================================

# --- Load input data ---
input_file <- file.path(material_dir, "values.tsv")  
filtered_fam <- fread(file.path(material_dir, "filtered_families.tsv"))
global_pairs <- fread(file.path(material_dir, "filtered_orthologs.tsv"))
tandems <- fread(file.path(material_dir, "tandem.tsv"), sep = "\t", header = TRUE)
annotations <- read.table(
  file.path(material_dir, "prot-scriber_result.txt"),
  header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE
)

# Extract all protein IDs from the tandem table (can be comma-separated)
tandem_protein_set <- unique(trimws(unlist(strsplit(tandems$genids, ","))))

# --- Optional: load source-copy pairs for paralog analysis ---
if (source_copy_analysis) {
  paralog_pairs <- fread(file.path(material_dir, "filtered_paralogs.tsv"), header = TRUE)
  source_copy_distances <- rep(NA_real_, nrow(paralog_pairs))
}

# Extract species names (all columns except Orthogroup)
species_names <- colnames(filtered_fam)[-1]

# Create dictionary: species -> list of genes
species_genes_dict <- lapply(species_names, function(species) {
  unique(trimws(unlist(strsplit(filtered_fam[[species]], ", "))))
})
names(species_genes_dict) <- species_names

# ==========================================================
# Main loop: iterate over modes and standardization methods
# ==========================================================
for (m in modes) {
  for (s in stds) {
    mode = m
    std = s

    print(paste0("Running normalization and analysis for mode=", mode, " and std=", std))

    # --- Load and preprocess expression matrix ---
    df <- read.table(input_file, header = TRUE, sep = "\t")
    all_zero_genes <- apply(df[, -1], 1, function(x) all(x == 0))
    df <- df[!all_zero_genes, ]
    all_na_genes <- apply(df[, -1], 1, function(x) all(is.na(x)))
    df <- df[!all_na_genes, ]

    gene_ids <- df[[1]]
    col_names <- colnames(df)[-1]
    df_matrix <- as.matrix(df[, -1])

    # --- Diagnose data quality before normalization ---
    diagnostics <- tox_diagnose_data_quality(df_matrix)

    # Clean or impute missing values if needed
    if (diagnostics$problems$na_count > 0 ||
        diagnostics$problems$inf_count > 0 ||
        diagnostics$problems$nan_count > 0) {

      df_matrix_clean <- tox_clean_data_for_normalization(
        df_matrix,
        remove_all_zero_genes = TRUE,
        na_strategy = "impute_mean",
        min_expression_threshold = 0.0,
        convert_small_to_zero = FALSE
      )

      gene_ids <- gene_ids[1:nrow(df_matrix_clean)]

    } else {
      df_matrix_clean <- df_matrix
    }

    # ==========================================================
    # Normalization steps
    # ==========================================================

    if (std == "knn") {
      normalized_matrix_std <- tox_normalize_by_knn_smoothed_std(df_matrix_clean, 1)
    } else {
      normalized_matrix_std <- tox_normalize_by_std_dev(df_matrix_clean)
    }

    normalized_matrix_qtl <- tox_quantile_normalization(normalized_matrix_std)
    norm_mat <- normalized_matrix_qtl

    # Average replicates by tissue if experiment is not Cardamine
    if (!grepl("cardamine", tolower(experiment))) {
      norm_mat <- tox_calculate_tissue_averages(normalized_matrix_qtl)
    }

    normalized_matrix_log <- tox_log2_transformation(norm_mat)

    # Optional: compute fold-change values if requested
    if (calculate_fc) {
      fc_df <- tox_calculate_fc_by_patterns(
        df = normalized_matrix_log,
        control_pattern = control_pattern,
        condition_patterns = c(condition_patterns)
      )
      normalized_matrix_log <- fc_df
    }

    normalized_df <- data.frame(gene_id = gene_ids, normalized_matrix_log)
    tissue_names <- colnames(normalized_df)[-1]
    unique_fams <- sort(unique(filtered_fam$Orthogroup))
    fam_indices <- setNames(seq_along(unique_fams), unique_fams)

    # --- Build family mapping: protein -> family index ---
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

    # --- Build ortholog set ---
    all_orthologs <- unique(c(global_pairs$Gene1, global_pairs$Gene2))
    ortholog_set <- gene_ids %in% all_orthologs
    paralog_set <- !ortholog_set

    # ==========================================================
    # Compute centroids per family
    # ==========================================================
    expression_vectors <- as.matrix(normalized_df[, -1])
    valid_idx <- which(!is.na(gene_to_family))
    expression_vectors <- as.matrix(normalized_df[valid_idx, -1])
    expression_vectors <- t(expression_vectors)  # transpose to axes x genes

    gene_to_family <- gene_to_family[valid_idx]
    ortholog_set <- ortholog_set[valid_idx]
    gene_ids <- gene_ids[valid_idx]
    n_axes <- nrow(expression_vectors)

    centroids <- tox_group_centroid(expression_vectors, as.integer(gene_to_family),
                                    n_families, ortholog_set, mode = mode)

    # ==========================================================
    # DISTANCE ANALYSIS
    # ==========================================================
    results_dist <- expression_analysis(
      expression_vectors, tissue_names, n_axes, n_families,
      "euclidean_distance", annotations, species_genes_dict,
      outliers_or_complete, gene_ids, centroids, gene_to_family,
      ortholog_set, tandem_protein_set
    )

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
      experiment = experiment,
      type_test = type_test
      
    )

    # ==========================================================
    # ANGLE ANALYSIS
    # ==========================================================
    results_angle <- expression_analysis(
      expression_vectors, tissue_names, n_axes, n_families,
      "angle", annotations, species_genes_dict,
      outliers_or_complete, gene_ids, centroids, gene_to_family,
      ortholog_set, tandem_protein_set
    )

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
      experiment = experiment,
      type_test = type_test
    )

    # ==========================================================
    # SOURCE-COPY ANALYSIS
    # ==========================================================
    if (source_copy_analysis) {
      sc_results_dist <- analyze_source_copy_data(
        expression_vectors, paralog_pairs, fam_indices, n_families,
        experiment, annotations, species_genes_dict, gene_ids,
        "euclidean_distance", outliers_or_complete, tandems
      )

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
        experiment = experiment,
        type_test = type_test,
        source_copy = TRUE
      )

      sc_results_angle <- analyze_source_copy_data(
        expression_vectors, paralog_pairs, fam_indices, n_families,
        experiment, annotations, species_genes_dict, gene_ids,
        "angle", outliers_or_complete, tandems
      )

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
        experiment = experiment,
        type_test = type_test,
        source_copy = TRUE
      )
    }

    # ==========================================================
    # OVERLAP AND VISUALIZATION (DISTANCE vs ANGLE)
    # ==========================================================
    outlier_ids_dist <- gene_ids[outlier_idx]
    outlier_ids_angle <- gene_ids[outlier_idx_angle]
    jaccard <- length(intersect(outlier_ids_dist, outlier_ids_angle)) /
               length(union(outlier_ids_dist, outlier_ids_angle))

    title_plot <- paste0(outliers_or_complete, " analysis - ", experiment, " - ", mode)
    pdf(file.path(results_dir, paste0(experiment, "_", type_test, "_", outliers_or_complete, "_", mode, ".pdf")),
        width = 10, height = 7)

    grid.newpage()
    grid.text(title_plot, x = 0.5, y = 0.5, gp = gpar(fontsize = 24, fontface = "bold"))
    grid.newpage()
    grid.text("Distance analysis", x = 0.5, y = 0.5, gp = gpar(fontsize = 24, fontface = "bold"))

    # Print all distance-based plots
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

    if (outliers_or_complete == "outliers") {
      grid.newpage()
      venn.plot <- venn.diagram(
        x = list(Distance = outlier_ids_dist, Angle = outlier_ids_angle),
        filename = NULL,
        fill = c("cornflowerblue", "darkorange"),
        alpha = 0.5,
        cex = 2, cat.cex = 2, cat.pos = c(-20, 20)
      )
      grid.draw(venn.plot)
      grid.text(sprintf("Jaccard similarity: %.3f", jaccard),
                x = 0.5, y = 0.1, gp = gpar(fontsize = 18, col = "black"))
    }

    # Print all angle-based plots
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

    # ==========================================================
    # SOURCE-COPY ANALYSIS (PAIR-LEVEL COMPARISON)
    # ==========================================================
    if (source_copy_analysis) {
      # Build bidirectional pair sets for comparison
      make_bidirectional_pairs <- function(df) {
        direct <- paste(df$gene_id_source, df$gene_id_copy, sep = "_")
        inverse <- paste(df$gene_id_copy, df$gene_id_source, sep = "_")
        c(direct, inverse)
      }

      pairs_dist  <- make_bidirectional_pairs(sc_data_dist)
      pairs_angle <- make_bidirectional_pairs(sc_data_angle)

      canonicalize <- function(pairs) {
        sapply(strsplit(pairs, "_"), function(x) paste(sort(x), collapse = "_"))
      }

      unique_pairs_dist  <- unique(canonicalize(pairs_dist))
      unique_pairs_angle <- unique(canonicalize(pairs_angle))

      sc_jaccard <- length(intersect(unique_pairs_dist, unique_pairs_angle)) /
                    length(union(unique_pairs_dist, unique_pairs_angle))

      # --- Visualization section ---
      grid.newpage()
      grid.text("Source Copy Analysis - Euclidean distance", x = 0.5, y = 0.5,
                gp = gpar(fontsize = 24, fontface = "bold"))

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
      grid.text("Source Copy Analysis - Angle", x = 0.5, y = 0.5,
                gp = gpar(fontsize = 24, fontface = "bold"))

      if (outliers_or_complete == "outliers") {
        grid.newpage()
        venn.plot <- venn.diagram(
          x = list(Distance = unique_pairs_dist, Angle = unique_pairs_angle),
          filename = NULL,
          fill = c("cornflowerblue", "darkorange"),
          alpha = 0.5,
          cex = 2, cat.cex = 2, cat.pos = c(-20, 20)
        )
        grid.draw(venn.plot)
        grid.text(sprintf("Jaccard similarity: %.3f", sc_jaccard),
                  x = 0.5, y = 0.1, gp = gpar(fontsize = 18, col = "black"))
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
