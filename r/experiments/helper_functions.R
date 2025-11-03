library(ggplot2)
library(ggpubr)
library(data.table)
library(wordcloud)
library(tm)
library(dplyr)
library(grid)
library(VennDiagram)

map_species <- function(ids) {
    vapply(ids, function(gid) {
        sp <- names(species_genes_dict)[vapply(
        species_genes_dict, function(genes) gid %in% genes, logical(1)
        )]
        if (length(sp) > 0) sp[1] else NA_character_
    }, character(1))
}


compute_word_frequencies <- function(annotations_df) {

  if (is.null(annotations_df) || nrow(annotations_df) == 0)
    return(list(
      word_freq = numeric(0),
      word_freq_df = data.frame(word = character(), frequency = integer())
    ))
  
  if (!"Human_Readable_Description" %in% names(annotations_df))
    return(list(
      word_freq = numeric(0),
      word_freq_df = data.frame(word = character(), frequency = integer())
    ))

  texts <- annotations_df$Human_Readable_Description
  texts <- texts[!is.na(texts) & nzchar(texts)]
  if (length(texts) == 0)
    return(list(
      word_freq = numeric(0),
      word_freq_df = data.frame(word = character(), frequency = integer())
    ))

  corpus <- tm::Corpus(tm::VectorSource(texts))
  non.inf.words <- c(
    "protein", "containing", "family", "domain", "isoform", "subunit",
    "associated", "duf", "factor", "dependent", "binding", "of", "and",
    "or", "unknown"
  )

  corpus <- tm::tm_map(corpus, tm::content_transformer(tolower))
  corpus <- tm::tm_map(corpus, tm::removePunctuation)
  corpus <- tm::tm_map(corpus, tm::removeNumbers)
  corpus <- tm::tm_map(corpus, tm::stripWhitespace)
  corpus <- tm::tm_map(corpus, tm::removeWords, non.inf.words)

  tdm <- tm::TermDocumentMatrix(corpus)
  tdm_matrix <- as.matrix(tdm)

  if (nrow(tdm_matrix) == 0 || ncol(tdm_matrix) == 0)
    return(list(
      word_freq = numeric(0),
      word_freq_df = data.frame(word = character(), frequency = integer())
    ))

  num_docs <- ncol(tdm_matrix)
  term_doc_freq <- rowSums(tdm_matrix > 0)
  filtered_terms <- term_doc_freq < (num_docs * 0.8)
  filtered_matrix <- if (any(filtered_terms))
    tdm_matrix[filtered_terms, , drop = FALSE] else tdm_matrix

  wf <- sort(rowSums(filtered_matrix), decreasing = TRUE)
  wf_df <- data.frame(
    word = names(wf),
    frequency = as.integer(wf),
    row.names = NULL
  )

  list(word_freq = wf, word_freq_df = wf_df)
}

expression_analysis <- function(
  expression_vectors,
  tissue_names,
  n_axes,
  n_families,
  analysis_type,                # "euclidean_distance" o "angle"
  annotations,
  species_genes_dict,
  outliers_or_complete,
  gene_ids,
  centroids,
  gene_to_family,
  ortholog_set,
  tandem_protein_set
) {


  # =====================================================
  # =====================================================
  comp1   <- expression_vectors
  comp2   <- centroids
  fam_map <- as.integer(gene_to_family)
  genes   <- gene_ids

  # =====================================================
  # =====================================================
  if (analysis_type == "euclidean_distance") {
    distances <- tox_distance_to_centroid(comp1, comp2, fam_map, n_axes)
    scaling   <- tox_compute_family_scaling(distances, fam_map, n_families)
    histogram <- tox_compute_rdi(distances, fam_map, scaling$dscale)
    outliers  <- tox_identify_outliers(histogram$rdi, 95)
  } else if (analysis_type == "angle") {
    distances <- tox_angle_to_centroid(comp1, comp2, fam_map) * (180/pi)
    res_sd    <- tox_circular_sd_per_family(comp1, comp2, fam_map)
    histogram <- tox_scaled_angles(distances * (pi/180), res_sd$s_c, fam_map)
    outliers  <- tox_angle_outliers(histogram, p = 0.95)
  } else {
    stop("analysis_type must be 'euclidean_distance' or 'angle'")
  }

  # =====================================================
  # =====================================================
  selected_idx <- if (outliers_or_complete == "outliers")
                    which(outliers$is_outlier) else seq_along(genes)
  sel <- seq_along(genes) %in% selected_idx

  # =====================================================
  # =====================================================
  is_tandem <- genes %in% tandem_protein_set
  group <- ifelse(ortholog_set[sel], "Ortholog",
                  ifelse(is_tandem[sel], "Tandem", "Non-Tandem"))
  species_vec <- map_species(genes[sel])

  # =====================================================
  # =====================================================
  axis_sel <- rep(TRUE, n_axes)
  vres <- tox_calculate_tissue_versatility(comp1, sel, axis_sel)
  comp2_angles <- tox_calculate_tissue_versatility(comp2, rep(TRUE, ncol(comp2)), axis_sel)$tissue_angles_deg
  delta_angle <- vres$tissue_angles_deg - comp2_angles[as.integer(fam_map[sel])]

  # =====================================================
  # =====================================================
  proj1 <- omics_vector_RAP_projection(comp1, sel, axis_sel)
  proj2 <- omics_vector_RAP_projection(comp2, rep(1, ncol(comp2)), axis_sel)
  delta_pref <- vapply(seq_along(selected_idx), function(i) {
    fam <- as.integer(fam_map[sel][i])
    if (!is.na(fam))
      tox_clock_hand_angle_between_vectors(proj1[, i], proj2[, fam]) * 180/pi
    else NA_real_
  }, numeric(1))

  shift <- tox_compute_shift_vector_field(comp1, comp2, fam_map)$shift_vectors
  if (is.null(dim(shift))) shift <- matrix(shift, nrow = 2*n_axes)
  V  <- shift[(n_axes+1):(2*n_axes), ] - shift[1:n_axes, ]
  Vn <- apply(V, 2, function(x) { n <- sqrt(sum(x^2)); if (n > 0) x/n else x })
  axis_contribs <- apply(Vn[, sel, drop = FALSE], 2, relative_axes_changes_from_shift_vector)

  # =====================================================
  # =====================================================
  df <- data.frame(
    gene_id = genes[sel],
    species = species_vec,
    group = group,
    distance = distances[sel],
    tissue_angle_deg = vres$tissue_angles_deg,
    change_in_tissue_versatility = delta_angle,
    change_in_tissue_preference  = delta_pref
  )

  axis_df <- data.frame(
    gene_id = rep(df$gene_id, each = n_axes),
    species = rep(df$species, each = n_axes),
    tissue  = rep(tissue_names, times = nrow(df)),
    contribution = as.vector(axis_contribs),
    group = rep(df$group, each = n_axes)
  )

  # =====================================================
  # =====================================================
  outlier_ids <- df$gene_id
  annots_sel <- tryCatch(
    annotations[annotations$Annotee_Identifier %in% outlier_ids, ],
    error = function(e) data.frame()
  )
  wf <- compute_word_frequencies(annots_sel)

  # =====================================================
  # =====================================================
  list(
    data = df,
    axis_contribs = axis_df,
    histogram = histogram,
    outliers = outliers,
    term_freq = wf$word_freq,
    outlier_idx = selected_idx
  )
}

analyze_source_copy_data <- function(
  expression_vectors,
  paralog_pairs,
  fam_indices,
  n_families,
  experiment,
  annotations,
  species_genes_dict,
  gene_ids,
  analysis_type,
  outliers_or_complete,
  tandems 
) {


  # =====================================================
  # =====================================================
  n_pairs <- nrow(paralog_pairs)
  source_vecs <- matrix(NA_real_, nrow = nrow(expression_vectors), ncol = n_pairs)
  copy_vecs   <- matrix(NA_real_, nrow = nrow(expression_vectors), ncol = n_pairs)
  distances   <- rep(NA_real_, n_pairs)

  for (i in seq_len(n_pairs)) {
    idx1 <- match(paralog_pairs$Gene1[i], gene_ids)
    idx2 <- match(paralog_pairs$Gene2[i], gene_ids)
    if (!is.na(idx1) && !is.na(idx2)) {
      v1 <- expression_vectors[, idx1]
      v2 <- expression_vectors[, idx2]
      distances[i] <- if (analysis_type == "angle")
        tox_angle_between_vectors(v1, v2)
      else
        tox_euclidean_distance(v1, v2)
      source_vecs[, i] <- v1
      copy_vecs[, i]   <- v2
    }
  }

  valid_idx <- which(!is.na(distances))
  distances <- distances[valid_idx]
  paralog_pairs <- paralog_pairs[valid_idx, ]
  source_vecs <- source_vecs[, valid_idx, drop = FALSE]
  copy_vecs <- copy_vecs[, valid_idx, drop = FALSE]
  fam_map <- fam_indices[paralog_pairs$Orthogroup]
  source_gene_ids <- paralog_pairs$Gene1
  copy_gene_ids   <- paralog_pairs$Gene2

  # =====================================================
  # =====================================================
  if (analysis_type == "euclidean_distance") {
    scaling   <- tox_compute_family_scaling(distances, fam_map, n_families)
    histogram <- tox_compute_rdi(distances, fam_map, scaling$dscale)
    outliers  <- tox_identify_outliers(histogram$rdi, 95)
  } else if (analysis_type == "angle") {
    res_sd <- tox_circular_sd_per_family(source_vecs, copy_vecs, fam_map)
    histogram <- tox_scaled_angles(distances, res_sd$s_c, fam_map)
    outliers  <- tox_angle_outliers(histogram, p = 0.95)
    distances <- distances * (180/pi)
  } else {
    stop("analysis_type must be 'euclidean_distance' or 'angle'")
  }

  selected_idx <- if (outliers_or_complete == "outliers")
                    which(outliers$is_outlier) else seq_along(distances)
  sel <- seq_along(distances) %in% selected_idx

  # =====================================================
  # =====================================================

  parent_to_tandem_row <- setNames(rep(seq_len(nrow(tandems)), sapply(strsplit(tandems$genids, ","), length)),
                                  trimws(unlist(strsplit(paste(tandems$genids, collapse=","), ","))))

  row1 <- parent_to_tandem_row[paralog_pairs$Gene1]
  row2 <- parent_to_tandem_row[paralog_pairs$Gene2]
  is_tandem <- !is.na(row1) & !is.na(row2) & (row1 == row2)
  group <- ifelse(is_tandem[sel], "Tandem", "Non-Tandem")

  # =====================================================
  # =====================================================
  axis_sel <- rep(TRUE, nrow(expression_vectors))
  vres <- tox_calculate_tissue_versatility(source_vecs, sel, axis_sel)
  copy_angles <- tox_calculate_tissue_versatility(copy_vecs, rep(TRUE, ncol(copy_vecs)), axis_sel)$tissue_angles_deg
  delta_angle <- vres$tissue_angles_deg - copy_angles[as.integer(fam_map[sel])]

  proj_source <- omics_vector_RAP_projection(source_vecs, sel, axis_sel)
  proj_copy   <- omics_vector_RAP_projection(copy_vecs, rep(1, ncol(copy_vecs)), axis_sel)
  delta_pref <- vapply(seq_along(selected_idx), function(i) {
    fam <- as.integer(fam_map[sel][i])
    if (!is.na(fam))
      tox_clock_hand_angle_between_vectors(proj_source[, i], proj_copy[, fam]) * 180/pi
    else NA_real_
  }, numeric(1))

  # =====================================================
  # =====================================================
  shift <- tox_compute_shift_vector_field(source_vecs, copy_vecs, fam_map)$shift_vectors
  if (is.null(dim(shift))) shift <- matrix(shift, nrow = 2 * nrow(expression_vectors))
  V  <- shift[(nrow(expression_vectors) + 1):(2 * nrow(expression_vectors)), ] - shift[1:nrow(expression_vectors), ]
  Vn <- apply(V, 2, function(x) { n <- sqrt(sum(x^2)); if (n > 0) x/n else x })
  axis_contribs <- apply(Vn[, sel, drop = FALSE], 2, relative_axes_changes_from_shift_vector)

  # =====================================================
  # =====================================================
  species_src  <- map_species(source_gene_ids[sel])

  df <- data.frame(
    gene_id_source = source_gene_ids[sel],
    gene_id_copy   = copy_gene_ids[sel],
    species = species_src,
    group = group,
    distance = distances[sel],
    tissue_angle_deg = vres$tissue_angles_deg,
    change_in_tissue_versatility = delta_angle,
    change_in_tissue_preference  = delta_pref
  )

  axis_df <- data.frame(
    gene_id_source = rep(df$gene_id_source, each = nrow(expression_vectors)),
    species = rep(df$species, each = nrow(expression_vectors)),
    tissue  = rep(tissue_names, times = nrow(df)),
    contribution = as.vector(axis_contribs),
    group = rep(df$group, each = nrow(expression_vectors))
  )

  # =====================================================
  # =====================================================
  outlier_ids <- unique(c(df$gene_id_source, df$gene_id_copy))
  annots_sel <- tryCatch(
    annotations[annotations$Annotee_Identifier %in% outlier_ids, ],
    error = function(e) data.frame()
  )
  wf <- compute_word_frequencies(annots_sel)

  # =====================================================
  # =====================================================
  list(
    data = df,
    axis_contribs = axis_df,
    histogram = histogram,
    outliers = outliers,
    term_freq = wf$word_freq,
    outlier_idx = selected_idx
  )
}


make_histogram_plot <- function(histogram, outliers, analysis_type) {
  if (analysis_type == "euclidean_distance") {
    df <- data.frame(value = histogram$sorted_rdi)
    x_label <- "RDI value"
    title <- "Histogram of Relative Distance Index (RDI)"
  } else if (analysis_type == "angle") {
    z_values <- if (is.list(histogram)) histogram$z else histogram
    df <- data.frame(value = z_values)
    x_label <- "Scaled angular deviation (z)"
    title <- "Histogram of Scaled Angular Deviations"
  } else {
    stop("analysis_type must be either 'euclidean_distance' or 'angle'")
  }

  ggplot(df, aes(x = value)) +
    geom_histogram(binwidth = 0.05, fill = "#1f77b4", color = "#1f77b4", alpha = 0.7) +
    geom_vline(xintercept = outliers$threshold, color = "red", linetype = "dashed", size = 1.2) +
    labs(
      title = title,
      x = x_label,
      y = "Frequency"
    ) +
    theme_minimal(base_size = 16) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
}


make_all_analysis_plots <- function(data, axis_contribs_df, histogram, outliers, experiment) {
  desired_order <- c("Ortholog", "Non-Tandem", "Tandem")
  present_groups <- desired_order[desired_order %in% unique(data$group)]
  color_palette <- c("Ortholog" = "#1f77b4", "Non-Tandem" = "#ff7f0e", "Tandem" = "#2ca02c")
  fill_colors <- color_palette[present_groups]
  data$group <- factor(data$group, levels = present_groups)
  comparisons <- list(c("Ortholog", "Non-Tandem"),
                      c("Ortholog", "Tandem"),
                      c("Non-Tandem", "Tandem"))
  group_labels <- sapply(present_groups, function(g) paste0(g, " (n=", sum(data$group == g), ")"))

  # --- 1. Distance to centroid ---
  p_distance_general <- make_violin_plot(
    data, x="group", y="distance", fill="group", color="group",
    labels=group_labels, title="Distance to centroid",
    ylab="Distance to centroid", comparisons=comparisons,
    palette=fill_colors
  )
  p_distance_species <- make_violin_plot(
    data, x="group", y="distance", fill="group", color="group",
    labels=group_labels, title="Distance to centroid",
    ylab="Distance to centroid", comparisons=comparisons,
    palette=fill_colors, facet_by="species"
  )

  # --- 2. Histogram (RDI or Angle depending on data) ---
  hist_plot <- make_histogram_plot(histogram, outliers, analysis_type = if ("sorted_rdi" %in% names(histogram)) "euclidean_distance" else "angle")

  # --- 3. Tissue versatility ---
  present_types <- desired_order[desired_order %in% unique(data$group)]
  fill_colors_v <- color_palette[present_types]
  group_labels_v <- sapply(present_types, function(g) paste0(g, " (n=", sum(data$group == g), ")"))

  p_versatility_general <- make_violin_plot(
    data, x="group", y="tissue_angle_deg", fill="group", color="group",
    labels=group_labels_v, title="Tissue versatility",
    ylab="Tissue Angle (Degrees)", comparisons=comparisons,
    palette=fill_colors_v, abs=TRUE
  )
  p_versatility_species <- make_violin_plot(
    data, x="group", y="tissue_angle_deg", fill="group", color="group",
    labels=group_labels_v, title="Tissue versatility",
    ylab="Tissue Angle (Degrees)", comparisons=comparisons,
    palette=fill_colors_v, facet_by="species", abs=TRUE
  )

  # --- 4. Change in tissue versatility ---
  p_change_vers_general <- make_violin_plot(
    data, x="group", y="change_in_tissue_versatility",
    fill="group", color="group", labels=group_labels_v,
    title="Change in tissue versatility",
    ylab="Change in Tissue Angle (Degrees)",
    comparisons=comparisons, palette=fill_colors_v, abs=TRUE
  )
  p_change_vers_species <- make_violin_plot(
    data, x="group", y="change_in_tissue_versatility",
    fill="group", color="group", labels=group_labels_v,
    title="Change in tissue versatility",
    ylab="Change in Tissue Angle (Degrees)",
    comparisons=comparisons, palette=fill_colors_v, facet_by="species", abs=TRUE
  )

  # --- 5. Change in tissue preference ---
  p_change_pref_general <- make_violin_plot(
    data, x="group", y="change_in_tissue_preference",
    fill="group", color="group", labels=group_labels_v,
    title="Change in tissue preference",
    ylab="Change in Tissue Preference (Degrees)",
    comparisons=comparisons, palette=fill_colors_v, abs=TRUE
  )
  p_change_pref_species <- make_violin_plot(
    data, x="group", y="change_in_tissue_preference",
    fill="group", color="group", labels=group_labels_v,
    title="Change in tissue preference",
    ylab="Change in Tissue Preference (Degrees)",
    comparisons=comparisons, palette=fill_colors_v, facet_by="species", abs=TRUE
  )

  # --- 6. Relative axis change from shift vectors ---
  if (grepl("cardamine", tolower(experiment))) {
    tissue_order <- c("seedling", "cotyledon", "developing_leaf", "flower_stage_9", "flower_stage_16")
  } else {
    tissue_order <- unique(axis_contribs_df$tissue)
  }
  axis_contribs_df$tissue <- factor(axis_contribs_df$tissue, levels = tissue_order)
  axis_contribs_df$group <- factor(axis_contribs_df$group, levels = desired_order)

  p_axis_contribs <- make_tissue_boxplot_with_significance(
    df = axis_contribs_df,
    fill_colors = fill_colors_v,
    title = "Relative axis change from shift vectors"
  )

  return(list(
    distance_general = p_distance_general,
    distance_species = p_distance_species,
    rdi_histogram = hist_plot,
    tissue_versatility_general = p_versatility_general,
    tissue_versatility_species = p_versatility_species,
    change_in_versatility_general = p_change_vers_general,
    change_in_versatility_species = p_change_vers_species,
    change_in_preference_general = p_change_pref_general,
    change_in_preference_species = p_change_pref_species,
    axis_contribs = p_axis_contribs
  ))
}



calculate_wilcox_pvalues <- function(data, group_col, value_col, comparisons) {
  group_col <- rlang::sym(group_col)
  value_col <- rlang::sym(value_col)

  df <- data.frame(
    group = data[[rlang::as_string(group_col)]],
    value = data[[rlang::as_string(value_col)]]
  )

  if (length(unique(df$group)) < 2) return(data.frame())

  results <- tryCatch({
    lapply(comparisons, function(pair) {
      x <- df$value[df$group == pair[1]]
      y <- df$value[df$group == pair[2]]

      if (length(x) > 0 && length(y) > 0) {
        test <- wilcox.test(x, y, alternative = "less")
        data.frame(group1 = pair[1], group2 = pair[2], p.value = test$p.value)
      } else {
        NULL
      }
    }) |> dplyr::bind_rows()
  }, error = function(e) data.frame())

  if (nrow(results) == 0) return(data.frame())

  results$signif <- cut(
    results$p.value,
    breaks = c(-Inf, 0.001, 0.01, 0.05, 1),
    labels = c("***", "**", "*", "ns")
  )

  results
}

calculate_wilcox_intra_tissue <- function(data, tissue_col, group_col, value_col) {
  tissue_col <- rlang::sym(tissue_col)

  comparisons <- list(
    c("Ortholog", "Non-Tandem"),
    c("Ortholog", "Tandem"),
    c("Non-Tandem", "Tandem")
  )

  df <- data.frame(
    tissue = data[[rlang::as_string(tissue_col)]],
    group  = data[[group_col]],
    value  = data[[value_col]]
  )

  if (length(unique(df$group)) < 2) return(data.frame())

  results <- df %>%
    dplyr::group_by(tissue) %>%
    dplyr::group_modify(~ calculate_wilcox_pvalues(.x, "group", "value", comparisons)) %>%
    dplyr::ungroup()

  results
}


make_violin_plot <- function(
  data,
  x,
  y,
  fill,
  color,
  labels,
  title,
  ylab,
  comparisons,
  palette,
  facet_by = NULL,
  reverse_y = FALSE,
  abs = FALSE
) {
  if (!is.null(facet_by)) {
    counts <- data %>%
      group_by(!!sym(x), !!sym(facet_by)) %>%
      summarise(n = n(), .groups = "drop")
  } else {
    counts <- data %>%
      group_by(!!sym(x)) %>%
      summarise(n = n(), .groups = "drop")
  }

  counts <- counts %>%
    mutate(
      offset = rep(c(-0.15, 0.15), length.out = n()),
      n_label = paste0("n = ", n)
    )

  if (abs) data[[y]] <- abs(data[[y]])

  if (!is.null(facet_by)) {
    facet_col <- rlang::sym(facet_by)

    pvals <- data %>%
      dplyr::group_by(!!facet_col) %>%
      dplyr::group_modify(~{
        pv <- tryCatch({
          calculate_wilcox_pvalues(.x, x, y, comparisons)
        }, error = function(e) data.frame())
        if (nrow(pv) > 0) pv[[facet_by]] <- unique(.x[[facet_by]])[1]
        pv
      }) %>%
      dplyr::ungroup()

    y_pos_df <- data %>%
      dplyr::group_by(!!facet_col) %>%
      dplyr::summarise(
        y_max = max(!!sym(y), na.rm = TRUE),
        y_min = min(!!sym(y), na.rm = TRUE),
        y_range = y_max - y_min,
        y_label = y_max + 0.55 * y_range,             
        y_n_label = y_min - 0.20 * y_range,          
        .groups = "drop"
      )

    pvals <- dplyr::left_join(pvals, y_pos_df, by = facet_by)
    counts <- dplyr::left_join(counts, y_pos_df, by = facet_by)

    } else {
    pvals <- tryCatch({
      calculate_wilcox_pvalues(data, x, y, comparisons)
    }, error = function(e) data.frame())

    if (ncol(pvals) == 0) {
      pvals <- data.frame(
        group1 = character(),
        group2 = character(),
        p.value = numeric(),
        signif = character(),
        y_label = numeric()
      )
    }

    y_max <- max(data[[y]], na.rm = TRUE)
    y_min <- min(data[[y]], na.rm = TRUE)
    y_range <- y_max - y_min

    if (nrow(pvals) > 0) {
      pvals$y_label <- y_max + 0.50 * y_range
    }

    counts$y_n_label <- y_min - 0.20 * y_range
  }


  p <- ggplot(data, aes_string(x = x, y = y, fill = fill)) +
    geom_jitter(aes_string(color = color), width = 0.12, alpha = 0.3, size = 2) +
    geom_violin(alpha = 0.4, trim = FALSE) +
    geom_boxplot(
      fill = "white",
      outlier.shape = NA,
      alpha = 0.4,
      width = 0.5,
      position = position_dodge(width = 0.8)
    ) +
    scale_fill_manual(values = palette) +
    scale_color_manual(values = palette) +
    scale_x_discrete(labels = labels) +
    labs(title = title, x = NULL, y = ylab) +
    theme_minimal(base_size = 16) +
    theme(
      legend.position = "top",
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )

  if (nrow(pvals) > 0) {
    p <- p + geom_text(
      data = pvals,
      aes(
        x = (match(group1, levels(factor(data[[x]]))) +
             match(group2, levels(factor(data[[x]])))) / 2,
        y = y_label,
        label = signif
      ),
      vjust = 0,
      fontface = "bold",
      size = 4,
      inherit.aes = FALSE
    )
  }

  p <- p + geom_text(
    data = counts,
    aes_string(
      x = paste0("as.numeric(", x, ") + offset"),
      y = "y_n_label",  
      label = "n_label"
    ),
    inherit.aes = FALSE,
    size = 4,
    color = "gray20",
    fontface = "bold",
    vjust = 1
  )

  if (!is.null(facet_by)) {
    p <- p + facet_wrap(as.formula(paste("~", facet_by)), scales = "free_y", nrow = 2)
  }

  if (reverse_y) p <- p + scale_y_reverse()

  p <- p +
    scale_y_continuous(expand = expansion(mult = c(0.15, 0.25))) +
    coord_cartesian(clip = "off") +
    theme(plot.margin = margin(10, 10, 30, 10))

  return(p)
}

make_tissue_boxplot_with_significance <- function(
  df,
  value_col = "contribution",
  tissue_col = "tissue",
  group_col = "group",
  fill_colors,
  title = "Relative axis change from shift vectors",
  ylab = "Normalized vector change"
) {

    pvals_tissue <- calculate_wilcox_intra_tissue(df, tissue_col, group_col, value_col)

    if (ncol(pvals_tissue) == 0) {
    pvals_tissue <- data.frame(
        tissue = character(),
        group1 = character(),
        group2 = character(),
        p.value = numeric(),
        signif = character()
    )
    }

  y_max_df <- df %>%
    dplyr::group_by(.data[[tissue_col]]) %>%
    dplyr::summarise(y_max = max(.data[[value_col]], na.rm = TRUE), .groups = "drop")

  pvals_tissue <- dplyr::left_join(pvals_tissue, y_max_df, by = setNames(tissue_col, "tissue"))

  unique_groups <- sort(unique(df[[group_col]]))
  n_groups <- length(unique_groups)

  if (n_groups == 3) {
    group_positions <- c("Ortholog" = 1, "Non-Tandem" = 2, "Tandem" = 3)
    x_offset_factor <- 0.25
  } else if (n_groups == 2) {
    group_positions <- setNames(1:2, unique_groups)
    x_offset_factor <- 0.0  
  } else {
    group_positions <- setNames(seq_along(unique_groups), unique_groups)
    x_offset_factor <- 0.25
  }

    y_max_global <- max(df[[value_col]], na.rm = TRUE)
    const_sep <- 0.05 * y_max_global  

  pvals_tissue <- pvals_tissue %>%
    dplyr::mutate(
      x_pos = purrr::map2_dbl(group1, group2,
                              ~ mean(c(group_positions[.x], group_positions[.y]))),
      comp_id = paste(group1, group2, sep = "_")
    ) %>%
    dplyr::group_by(.data[[tissue_col]]) %>%
    dplyr::mutate(
    y_pos = y_max + const_sep + (dplyr::row_number() - 2) * const_sep
    ) %>%
    dplyr::ungroup()

  tissue_levels <- levels(factor(df[[tissue_col]]))
  tissue_numeric <- setNames(seq_along(tissue_levels), tissue_levels)
  pvals_tissue$tissue_num <- tissue_numeric[pvals_tissue$tissue]

  n_df <- df %>%
    dplyr::group_by(.data[[tissue_col]], .data[[group_col]]) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
    dplyr::mutate(
      n_label = paste0("n=", n),
      group_pos = group_positions[.data[[group_col]]],
      tissue_num = tissue_numeric[.data[[tissue_col]]]
    )

  y_min_df <- df %>%
    dplyr::group_by(.data[[tissue_col]]) %>%
    dplyr::summarise(y_min = min(.data[[value_col]], na.rm = TRUE), .groups = "drop")

  n_df <- dplyr::left_join(n_df, y_max_df, by = setNames(tissue_col, "tissue")) %>%
          dplyr::left_join(y_min_df, by = setNames(tissue_col, "tissue"))

  const_sep_n <- 0.1 * y_max_global
  n_df <- n_df %>%
    dplyr::mutate(
      y_pos = y_min - const_sep_n -
              0.5 * (group_pos - mean(group_positions)) * const_sep_n
    )


  p <- ggplot(df, aes_string(x = tissue_col, y = value_col, fill = group_col)) +
    geom_jitter(
      aes_string(color = group_col),
      position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5),
      alpha = 0.3, size = 2
    ) +
    geom_boxplot(outlier.shape = NA, alpha = 0.4, width = 0.5) +

    geom_text(
      data = pvals_tissue,
      aes(
        x = tissue_num + (x_pos - mean(group_positions)) * x_offset_factor,
        y = y_pos,
        label = signif
      ),
      vjust = 0,
      fontface = "bold",
      size = 3,
      color = "black",
      inherit.aes = FALSE
    ) +

    geom_text(
      data = n_df,
      aes(
        x = tissue_num + (group_pos - mean(group_positions)) * 0.25, # ← fijo, no depende de x_offset_factor
        y = y_pos,
        label = n_label
      ),
      vjust = 1,
      fontface = "bold",
      size = 3,
      color = "gray30",
      inherit.aes = FALSE
    ) +

    scale_fill_manual(values = fill_colors) +
    scale_color_manual(values = fill_colors) +
    labs(
      title = title,
      x = "Tissue",
      y = ylab,
      fill = "Group",
      color = "Group"
    ) +
    theme_minimal(base_size = 16) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "top"
    )

  return(p)
}
