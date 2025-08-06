library(dplyr)
library(tibble)
library(readr)

# ---- Helper ----

# Compute cosine similarity
cosine_sim <- function(q_vec, ref_mat) {
  sims <- as.vector(ref_mat %*% q_vec) / (sqrt(rowSums(ref_mat^2)) * sqrt(sum(q_vec^2)))
  return(sims)
}

# Get top N most similar protein IDs to q_vec in ref_mat
get_top_ids <- function(q_vec, ref_mat, top_n = 5) {
  sims <- cosine_sim(q_vec, ref_mat)
  top <- sort(sims, decreasing = TRUE)[1:top_n]
  return(names(top))
}

# Retrieve annotations for top hits
get_predictions <- function(ids, annotation_df) {
  annotation_df %>%
    filter(protein_id %in% ids) %>%
    distinct(protein_id, term) %>%
    pull(term)
}

# Compute geometric mean of percent identity and alignment overlap
compute_gm_sim <- function(pident, align_len, q_len, s_len) {
  overlap <- align_len / pmax(q_len, s_len)
  gm <- sqrt(pident * overlap)
  return(gm)
}

# Determine where similarity drops sharply — STOP condition
get_top_ids_with_stop <- function(q_vec, ref_mat, max_k = 30, window_k = 3) {
  sims <- sort(cosine_sim(q_vec, ref_mat), decreasing = TRUE)
  selected <- names(sims)[1]

  deltas <- numeric()
  for (i in 2:min(length(sims), max_k)) {
    delta <- sims[i] - sims[i - 1]
    deltas <- c(deltas, delta)

    if (i > window_k) {
      mu <- mean(deltas[(i - window_k):(i - 1)])
      sigma <- sd(deltas[(i - window_k):(i - 1)])
      if (delta < (mu - 2 * sigma)) {
        break
      }
    }

    selected <- c(selected, names(sims)[i])
  }

  return(selected)
}

# ---- Load lexical axes and reference spaces ----
axes_GO  <- readRDS("material/rds_files/lexical_axes_GO.rds")
axes_HRD <- readRDS("material/rds_files/lexical_axes_HRD.rds")
axes_IPR <- readRDS("material/rds_files/lexical_axes_IPR.rds")

go_ref  <- axes_GO$doc_vectors
hrd_ref <- axes_HRD$doc_vectors
ipr_ref <- axes_IPR$doc_vectors

# ---- Load annotation data ----
ann_goa <- readRDS("material/rds_files/protein2goa.rds")
ann_hrd <- readRDS("material/rds_files/protein2hrd.rds")
ann_ipr <- readRDS("material/rds_files/protein2ipr.rds")

# ---- Load Kidera ----
kidera <- readRDS("material/rds_files/kidera_factors_uniref50_morethan1_ref_prots.rds")
kidera_ref <- kidera %>%
  select(protein_id, starts_with("K")) %>%
  column_to_rownames("protein_id") %>%
  as.matrix()

# ---- Main function ----
predict_functions <- function(protein_list, blast_results, kidera) {
  all_preds <- list()

  for (q in protein_list) {
    hits <- blast_results %>%
      filter(query == q) %>%
      mutate(
        gm_sim = compute_gm_sim(pident, align_length, q_len, s_len)
      ) %>%
      filter(!is.na(gm_sim) & gm_sim > 0) %>%
      arrange(desc(gm_sim)) %>%
      head(30)

    # Log-transform gm-sim and normalize
    gm_sim_log <- log1p(hits$gm_sim)
    w <- gm_sim_log / sum(gm_sim_log)

    # Filter matching reference vectors
    V_go  <- go_ref[intersect(hits$subject, rownames(go_ref)), , drop = FALSE]
    V_hrd <- hrd_ref[intersect(hits$subject, rownames(hrd_ref)), , drop = FALSE]
    V_ipr <- ipr_ref[intersect(hits$subject, rownames(ipr_ref)), , drop = FALSE]
    V_kid <- kidera_ref[intersect(hits$subject, rownames(kidera_ref)), , drop = FALSE]

    # Match weights to matrix rows
    idx_go  <- match(rownames(V_go), hits$subject)
    idx_hrd <- match(rownames(V_hrd), hits$subject)
    idx_ipr <- match(rownames(V_ipr), hits$subject)
    idx_kid <- match(rownames(V_kid), hits$subject)

    q_go  <- colSums(V_go  * w[idx_go])
    q_hrd <- colSums(V_hrd * w[idx_hrd])
    q_ipr <- colSums(V_ipr * w[idx_ipr])
    q_kid <- colSums(V_kid * w[idx_kid])

    # --- Get neighbors with STOP condition ---
    top_go_lex  <- get_top_ids_with_stop(q_go, go_ref)
    top_go_kid  <- get_top_ids_with_stop(q_kid, go_ref)

    top_hrd_lex <- get_top_ids_with_stop(q_hrd, hrd_ref)
    top_hrd_kid <- get_top_ids_with_stop(q_kid, hrd_ref)

    top_ipr_lex <- get_top_ids_with_stop(q_ipr, ipr_ref)
    top_ipr_kid <- get_top_ids_with_stop(q_kid, ipr_ref)

    # --- Predict terms ---
    pred_go <- union(
      get_predictions(top_go_lex, ann_goa),
      get_predictions(top_go_kid, ann_goa)
    )
    pred_hrd <- union(
      get_predictions(top_hrd_lex, ann_hrd),
      get_predictions(top_hrd_kid, ann_hrd)
    )
    pred_ipr <- union(
      get_predictions(top_ipr_lex, ann_ipr),
      get_predictions(top_ipr_kid, ann_ipr)
    )

    # Collect results
    all_preds[[q]] <- tibble(
      protein_id = q,
      predicted_term = c(pred_go, pred_hrd, pred_ipr),
      source = c(
        rep("GO",  length(pred_go)),
        rep("HRD", length(pred_hrd)),
        rep("IPR", length(pred_ipr))
      )
    )
  }

  bind_rows(all_preds)
}

# ---- Load and prepare BLAST results ----
blast_results <- read_tsv(
  "material/rds_files/diamond_test_1.txt", # <-- BLAST-Output hier angeben
  col_names = c("query", "subject", "q_start", "q_end", "q_len",
                "s_start", "s_end", "s_len", "pident", "evalue", "bitscore"),
  col_types = cols(
    query = col_character(),
    subject = col_character(),
    q_start = col_double(),
    q_end = col_double(),
    q_len = col_double(),
    s_start = col_double(),
    s_end = col_double(),
    s_len = col_double(),
    pident = col_double(),
    evalue = col_double(),
    bitscore = col_double()
  )
) %>%
  mutate(
    align_length = abs(q_end - q_start) + 1
  )

# ---- Automatically detect all query proteins ----
query_proteins <- unique(blast_results$query)

# ---- Run prediction ----
predictions <- predict_functions(query_proteins, blast_results, kidera)

# ---- Save output ----
write_csv(predictions, "predicted_functions.csv")
