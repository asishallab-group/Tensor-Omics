# Protein Function Prediction Pipeline (SVD + Kidera + Lexical Consistency)
# Version 4.3.3
# Author: Aaron Schroeder
# Date: 2025-08-18
#
# Overview
# - Build a chemical structure space from SVD embeddings (variable dims) + 10 Kidera factors.
# - Embed query proteins using BLAST hits with *normalized geometric mean* weights.
# - Supports multiple query proteins in a single BLAST file (separate embeddings per query).
# - Find nearest neighbors via Euclidean KD-tree (RANN) with a sequential stop rule.
# - Per corpus (GO/IPR/HRD): enforce high lexical consistency via cosine similarity (>= threshold).
# - Save results per query as CSV + RDS under results/<query_id>/.
#
# Notes
# - ID normalization removes the "UniRef50_" prefix to match other files.
# - The geometric mean is computed directly from DIAMOND (-f6) output using forward/backward overlap and pident (with 0.1 fallbacks).
# - Combination mode:
#     * "sum": conventional weighted combination:  sum_i w_i * v_i     (default)
#     * "mean_after_weight": mean over weighted vectors: (1/n) * sum_i w_i * v_i
# - Designed for very large reference sets, but memory limits apply. Consider chunking if needed.

suppressPackageStartupMessages({
  require(data.table)
  require(Matrix)
  require(RANN)
  require(stringr)
  require(tools)
  require(optparse)
})

# ------------------------
# Utility helpers
# ------------------------
canon_id <- function(x) {
  x <- as.character(x)
  x <- str_trim(x)
  # drop UniRef50_ prefix if present
  x <- sub("^(?i)UniRef50[_:]?", "", x, perl = TRUE)
  x
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

save_both <- function(obj, df, out_dir, base_name) {
  ensure_dir(out_dir)
  rds_path <- file.path(out_dir, sprintf("%s.rds", base_name))
  csv_path <- file.path(out_dir, sprintf("%s.csv", base_name))
  saveRDS(obj, rds_path)
  fwrite(df, csv_path)
  list(rds = rds_path, csv = csv_path)
}

# ------------------------
# Loading reference data
# ------------------------
load_lexical_axes <- function(path) {
  message("Loading lexical axes from: ", path, "...")
  ax <- readRDS(path)
  if (!"doc_vectors" %in% names(ax)) stop("lexical axes RDS missing 'doc_vectors'")
  dv <- ax$doc_vectors
  if (!is.null(rownames(dv))) {
    rownames(dv) <- canon_id(rownames(dv))
  }
  message("Lexical axes loaded: ", nrow(dv), " Entries.")
  list(doc_vectors = dv, explained_variance = ax$explained_variance)
}

# Load mean vector
load_mean_vector <- function(path) {
  message("Loading named mean vector ", path, "...")
  mean_dt <- fread(path)
  # assume first column is protein_id and second one is the vector
  setnames(mean_dt, names(mean_dt)[1:2], c("Protein_ID", "Mean_Vector_Value"))
  mean_dt[, Protein_ID := canon_id(Protein_ID)]
  setkey(mean_dt, Protein_ID)
  message("Mean vector loaded: ", nrow(mean_dt), " Entries.")
  mean_dt
}


compute_norm_gmean_weights <- function(blast_file, mean_dt = NULL) {
  message("Processing blast results from ", blast_file, "...")
  cols <- c("qseqid","sseqid","qstart","qend","qlen",
            "sstart","send","slen","pident","evalue","bitscore")
  dt <- data.table::fread(blast_file, sep = "\t", header = FALSE, col.names = cols)
  dt[, overlap := ((qend - qstart) + (send - sstart)) / (qlen + slen) * 100]
  dt[, qid := canon_id(qseqid)]
  dt[, sid := canon_id(sseqid)]
  dt_back <- data.table::copy(dt)
  data.table::setnames(dt_back,
                       c("qid","sid","overlap","pident"),
                       c("sid","qid","overlap_b","pident_b"))
  dt_back <- dt_back[, .(qid, sid, overlap_b, pident_b)]
  
  
  merged <- dt[dt_back, on = .(qid = sid, sid = qid)]
  
  merged[is.na(overlap_b), overlap_b := 0.1]
  merged[is.na(pident_b),  pident_b  := 0.1]
  
  gm <- function(vals) exp(mean(log(pmax(vals, 1e-12))))
  merged[, gmean := gm(c(overlap, overlap_b, pident, pident_b)), by = 1:nrow(merged)]
  merged <- merged[is.finite(gmean) & gmean > 0]

  # If mean vector is provided (this should be the case) use it
  if (!is.null(mean_dt)) {
    message("Using mean vector...")

    merged[, mean_val := mean_dt[.(sid), Mean_Vector_Value]]
    merged[is.na(mean_val), mean_val := 0] # Falls eine ID nicht im Mean-Vektor ist, 0 verwenden

    # calculate norm factor
    merged[, norm_factor := sum(gmean), by = qid]

    # get final weights
    merged[, w := (gmean - mean_val) / norm_factor]
    merged[norm_factor == 0, w := 0] # Division by zero vermeiden

  } else {
    # If no vector provided, use pseudo weights based on geom_mean
    merged[, w := gmean / sum(gmean), by = qid]
  }

  message("Calculated blast weights")
  merged[, .(id = sid, w, qid)]
}

# ------------------------
# Query embedding in SVD and chemical space
# ------------------------
combine_hit_vectors <- function(weights_dt, svd_mat, query_id, combine_mode = c("sum", "mean_after_weight")) {
  combine_mode <- match.arg(combine_mode)
  ids <- intersect(weights_dt[qid == query_id, id], rownames(svd_mat))
  if (length(ids) == 0) stop(sprintf("None of the BLAST hits for query %s are present in SVD matrix.", query_id))
  w <- weights_dt[qid == query_id][match(ids, id), w]
  V <- svd_mat[ids, , drop = FALSE]
  weighted <- V * w
  vq <- colSums(weighted)
  if (combine_mode == "mean_after_weight") {
    vq <- vq / nrow(weighted)
  }
  as.numeric(vq)
}

append_query_kidera <- function(v_svd_q, query_kidera_row) {
  if (is.null(query_kidera_row) || length(query_kidera_row) != 10) stop("Query Kidera row must have 10 numeric values.")
  c(v_svd_q, as.numeric(query_kidera_row))
}

# ------------------------
# KD-tree search + lexical gating per corpus
# ------------------------
cosine_sim <- function(x, Y) {
  x <- as.numeric(x)
  Y <- as.matrix(Y)
  x_norm <- sqrt(sum(x * x))
  if (x_norm == 0) return(rep(NA_real_, nrow(Y)))
  y_norms <- sqrt(rowSums(Y * Y))
  as.numeric((Y %*% x) / (y_norms * x_norm))
}

collect_neighbors_for_corpus <- function(neighbor_ids, doc_mat_sparse, threshold = 0.9) {
  if (length(neighbor_ids) == 0) return(list(ids = character(), similarities = numeric(), doc_vectors = NULL, mean_vector = NULL))
  
  available <- intersect(neighbor_ids, rownames(doc_mat_sparse))
  if (length(available) == 0) return(list(ids = character(), similarities = numeric(), doc_vectors = NULL, mean_vector = NULL))
  
  first_id <- available[1]
  selected_ids <- c(first_id)
  selected_mat <- as.matrix(doc_mat_sparse[first_id, , drop = FALSE])
  sim_values <- c(1.0) # Cosine similarity of a vector with itself is 1.0
  
  # check cosine for each neighbor to all neighbors that are currently in the set
  if (length(available) > 1) {
    for (cand in available[-1]) {
      cand_vec <- as.matrix(doc_mat_sparse[cand, , drop = FALSE])
      sims <- cosine_sim(cand_vec, selected_mat)
      
      if (all(is.finite(sims)) && all(sims >= threshold)) {
        selected_ids <- c(selected_ids, cand)
        selected_mat <- rbind(selected_mat, cand_vec)
        # Calculate similarity against the mean vector of selected neighbors
        mean_vec_current <- as.numeric(colMeans(selected_mat))
        current_sim <- cosine_sim(cand_vec, matrix(mean_vec_current, nrow = 1))
        sim_values <- c(sim_values, current_sim)
      } else {
        break
      }
    }
  }
  
  mean_vec <- as.numeric(colMeans(selected_mat))
  list(ids = selected_ids, similarities = sim_values, doc_vectors = selected_mat, mean_vector = mean_vec)
}

# ------------------------
# Driver for multiple queries
# ------------------------
run_for_all_queries <- function(blast_file,
                               chemical_space_rds,
                               query_kidera_rds_path,
                               mean_vector_csv = NULL,
                               lexical_paths = list(GO = NULL, IPR = NULL, HRD = NULL),
                               output_base = "results",
                               # cap at 50 to avoid unneccessary computations
                               k_cap = 50,
                               # default = 0.9, this should work fine but might miss some relevant factors, consider using 0.8
                               cosine_threshold = 0.9,
                               combine_mode = c("sum", "mean_after_weight")) {
  
  combine_mode <- match.arg(combine_mode)
  
  # load data
  message("--- Phase 1: Loading references ---")
  message("Loading chemical space from ", chemical_space_rds, "...")
  chem <- readRDS(chemical_space_rds)
  
  # Extract necessary components from the chemical space object
  chem_mat <- chem$mat
  svd_mat <- chem$mat[, colnames(chem$mat) %in% chem$svd$names]
  
  message("Loading query kidera factors...")
  qk <- readRDS(query_kidera_rds_path)
  qk <- as.data.table(qk)
  setnames(qk, names(qk)[1], "id")
  qk[, id := canon_id(id)]
  k_cols <- grep("^KF\\d+$", names(qk), value = TRUE)
  stopifnot(length(k_cols) == 10)

  # NEU: Lade Mean-Vektor, falls angegeben
  mean_dt <- NULL
  if (!is.null(mean_vector_csv)) {
    mean_dt <- load_mean_vector(mean_vector_csv)
  }
  
  message("Loading lexical vectors...")
  corpora <- list()
  for (nm in names(lexical_paths)) {
    p <- lexical_paths[[nm]]
    if (!is.null(p)) {
      corpora[[nm]] <- load_lexical_axes(p)
    }
  }
  
  # Phase 2: BLAST-Ergebnisse verarbeiten
  message("--- Phase 2: Processing blast results ---")
  # Übergabe des Mean-Vektors an die Funktion
  w_dt <- compute_norm_gmean_weights(blast_file, mean_dt = mean_dt)
  query_ids <- unique(w_dt$qid)
  message(length(query_ids), " query proteins found.")
  
  # Phase 3: Suche nach Nachbarn für jede Abfrage
  message("--- Phase 3: Starting neighborhood search ---")
  results <- list()
  for (i in 1:length(query_ids)) {
    query_id <- query_ids[i]
    message(sprintf("Processing query %d/%d: %s", i, length(query_ids), query_id))
    
    query_row <- qk[id == query_id]
    if (nrow(query_row) == 0) {
      message("  -> No kidera factors found; skipping query protein...")
      next
    }
    q_kidera_vec <- as.numeric(query_row[1, ..k_cols])
    
    message("  -> creating vector...")
    vq_svd <- combine_hit_vectors(w_dt, svd_mat, query_id, combine_mode = combine_mode)
    vq_chem <- append_query_kidera(vq_svd, q_kidera_vec)
    
    message("  -> Finding ", k_cap, " nearest neighbors in the chemical space...")
    nn <- RANN::nn2(data = chem_mat, query = matrix(vq_chem, nrow = 1), k = min(k_cap, nrow(chem_mat)))
    nn_idx <- as.integer(nn$nn.idx[1, ])
    neighbor_ids <- rownames(chem_mat)[nn_idx]
    neighbor_ids <- neighbor_ids[neighbor_ids != query_id]
    message("  -> Top 10 neighbors in the chemical space: ", paste(head(neighbor_ids, 10), collapse = ", "))
    
    per_corpus <- list()
    for (nm in names(corpora)) {
      message("  -> Leical analysis for '", nm, "'...")
      dv <- corpora[[nm]]$doc_vectors
      res <- collect_neighbors_for_corpus(neighbor_ids, dv, threshold = cosine_threshold)
      per_corpus[[nm]] <- res
    }
    
    out_dir <- file.path(output_base, tolower(gsub("[^A-Za-z0-9]+", "_", query_id)))
    ensure_dir(out_dir)
    
    rows <- list()
    for (nm in names(per_corpus)) {
      ids <- per_corpus[[nm]]$ids
      sims <- per_corpus[[nm]]$similarities
      if (length(ids)) {
        rows[[nm]] <- data.table(
          query_id = query_id,
          corpus = nm,
          neighbor_id = ids,
          cosine_similarity = sims,
          included = TRUE
        )
      }
    }
    neighbor_df <- if (length(rows)) rbindlist(rows, use.names = TRUE, fill = TRUE) else data.table()
    
    means_dt <- rbindlist(lapply(names(per_corpus), function(nm) {
      mv <- per_corpus[[nm]]$mean_vector
      if (is.null(mv)) return(NULL)
      data.table(query_id = query_id, corpus = nm, dim = paste0("EW", seq_along(mv)), value = mv)
    }), use.names = TRUE, fill = TRUE)
    
    rds_obj <- list(
      query_id = query_id,
      chemical_query_vector = vq_chem,
      chemical_dims = colnames(chem_mat),
      neighbor_order = neighbor_ids,
      per_corpus = per_corpus
    )
    
    message("  -> Saving results in ", out_dir, "...")
    out_files1 <- save_both(rds_obj, neighbor_df, out_dir, sprintf("%s_output", tolower(query_id)))
    means_csv <- file.path(out_dir, sprintf("%s_mean_doc_vectors.csv", tolower(query_id)))
    if (nrow(means_dt)) fwrite(means_dt, means_csv)
    
    results[[query_id]] <- list(
      out_rds = out_files1$rds,
      out_csv = out_files1$csv,
      out_means_csv = means_csv,
      out_dir = out_dir
    )
  }
  
  message("--- All Proteins processed ---")
  results
}

# ------------------------
# Main entry via command line
# ------------------------
option_list <- list(
  make_option(c("--blast"), type = "character", help = "BLAST tabular file (-f6)"),
  make_option(c("--chemical_space"), type = "character", help = "RDS file with pre-built chemical space"),
  make_option(c("--query_kidera"), type = "character", help = "RDS with query Kidera factors"),
  make_option(c("--mean_vector"), type = "character", default = NULL, help = "CSV file with mean vector values for normalization"), # NEUER PARAMETER
  make_option(c("--go"), type = "character", default = NULL, help = "RDS with GO lexical axes"),
  make_option(c("--ipr"), type = "character", default = NULL, help = "RDS with IPR lexical axes"),
  make_option(c("--hrd"), type = "character", default = NULL, help = "RDS with HRD lexical axes"),
  make_option(c("--out"), type = "character", default = "results", help = "Output base directory"),
  make_option(c("--kcap"), type = "integer", default = 50, help = "Number of KD-tree neighbors to consider"),
  make_option(c("--threshold"), type = "double", default = 0.9, help = "Cosine similarity threshold for lexical gating"),
  make_option(c("--combine"), type = "character", default = "sum", help = "Combination mode: sum | mean_after_weight")
)

args <- parse_args(OptionParser(option_list = option_list))

if (!is.null(args$blast) && !is.null(args$chemical_space) && !is.null(args$query_kidera)) {
  lexical_paths <- list(GO = args$go, IPR = args$ipr, HRD = args$hrd)
  run_for_all_queries(
    blast_file = args$blast,
    chemical_space_rds = args$chemical_space,
    query_kidera_rds_path = args$query_kidera,
    mean_vector_csv = args$mean_vector, # Übergabe des neuen Parameters
    lexical_paths = lexical_paths,
    output_base = args$out,
    k_cap = args$kcap,
    cosine_threshold = args$threshold,
    combine_mode = args$combine
  )
}
