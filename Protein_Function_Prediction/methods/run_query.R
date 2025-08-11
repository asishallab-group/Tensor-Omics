# Token-based prediction version of your script
# - Adds token-level calculations for HRD / GO / IPR
# - Keeps previous doc-vector / kidera / svd logic as well
# Requirements: dplyr, tibble, readr, stringr
#Version 2.1.2
library(dplyr)
library(tibble)
library(readr)
library(stringr)

# ---- Helper functions ----
cosine_sim <- function(q_vec, ref_mat) {
  # ref_mat: rows = entities, cols = vector dims
  if (length(q_vec) == 0 || nrow(ref_mat) == 0) return(numeric(0))
  denom_ref <- sqrt(rowSums(ref_mat^2))
  denom_q <- sqrt(sum(q_vec^2))
  # avoid division by zero
  denom_ref[denom_ref == 0] <- .Machine$double.eps
  if (denom_q == 0) denom_q <- .Machine$double.eps
  sims <- as.vector((ref_mat %*% q_vec) / (denom_ref * denom_q))
  names(sims) <- rownames(ref_mat)
  sims
}

get_top_ids <- function(q_vec, ref_mat, top_n = 5) {
  sims <- cosine_sim(q_vec, ref_mat)
  if (length(sims) == 0) return(character(0))
  top <- sort(sims, decreasing = TRUE)[1:min(top_n, length(sims))]
  names(top)
}

get_predictions <- function(ids, annotation_df) {
  if (length(ids) == 0) return(character(0))
  annotation_df %>%
    filter(protein_id %in% ids) %>%
    distinct(protein_id, term) %>%
    pull(term)
}

compute_gm_sim <- function(pident, align_len, q_len, s_len) {
  overlap <- align_len / pmax(q_len, s_len)
  gm <- sqrt(pident * overlap)
  gm
}

get_top_ids_with_stop <- function(q_vec, ref_mat, max_k = 30, window_k = 3) {
  if (length(q_vec) == 0 || nrow(ref_mat) == 0) return(character(0))
  
  sims <- cosine_sim(q_vec, ref_mat)
  if (length(sims) == 0) return(character(0))
  sims_sorted <- sort(sims, decreasing = TRUE)
  
  selected <- names(sims_sorted)[1]
  if (length(sims_sorted) == 1) return(selected)
  
  deltas <- numeric()
  for (i in 2:min(length(sims_sorted), max_k)) {
    delta <- sims_sorted[i] - sims_sorted[i - 1]
    deltas <- c(deltas, delta)
    
    if (i > window_k) {
      mu <- mean(deltas[(i - window_k):(i - 1)])
      sigma <- sd(deltas[(i - window_k):(i - 1)])
      if (is.na(sigma)) sigma <- 0
      if (delta < (mu - 2 * sigma)) {
        break
      }
    }
    
    selected <- c(selected, names(sims_sorted)[i])
  }
  
  selected
}


named_numeric <- function(names_vec) {
  v <- numeric(length(names_vec)); names(v) <- names_vec; v
}

# ---- Load data (unchanged) ----
axes_GO  <- readRDS("/media/BioNAS2/Tensor_Omics/Protein_Function_Prediction/material/lexical_axes_GO.rds")
axes_HRD <- readRDS("/media/BioNAS2/Tensor_Omics/Protein_Function_Prediction/material/lexical_axes_HRD.rds")
axes_IPR <- readRDS("/media/BioNAS2/Tensor_Omics/Protein_Function_Prediction/material/lexical_axes_IPR.rds")

protein_to_tokens_ipr <- readRDS("/media/BioNAS2/Tensor_Omics/Protein_Function_Prediction/material/protein_to_tokens_IPR.rds")
protein_to_tokens_go  <- readRDS("/media/BioNAS2/Tensor_Omics/Protein_Function_Prediction/material/protein_to_tokens_GO.rds")
protein_to_tokens_hrd <- readRDS("/media/BioNAS2/Tensor_Omics/Protein_Function_Prediction/material/protein_to_tokens_HRD.rds")

get_or_null <- function(x, name) {
  if (!is.null(x) && !is.null(x[[name]])) x[[name]] else NULL
}

go_doc_ref  <- get_or_null(axes_GO,  "doc_vectors")
go_token_mat <- get_or_null(axes_GO, "word_axes")
go_token_idf <- get_or_null(axes_GO, "idf")

hrd_doc_ref <- get_or_null(axes_HRD, "doc_vectors")
hrd_token_mat <- get_or_null(axes_HRD, "word_axes")
hrd_token_idf <- get_or_null(axes_HRD, "idf")

ipr_doc_ref  <- get_or_null(axes_IPR,  "doc_vectors")
ipr_token_mat <- get_or_null(axes_IPR, "word_axes")
ipr_token_idf <- get_or_null(axes_IPR, "idf")

go_ref  <- go_doc_ref 
hrd_ref <- hrd_doc_ref 
ipr_ref <- ipr_doc_ref

ann_goa <- readRDS("/media/BioNAS2/Tensor_Omics/Protein_Function_Prediction/material/protein2goa.rds")
ann_hrd <- readRDS("/media/BioNAS2/Tensor_Omics/Protein_Function_Prediction/material/protein2hrd.rds")
ann_ipr <- readRDS("/media/BioNAS2/Tensor_Omics/Protein_Function_Prediction/material/protein2ipr.rds")

kidera <- readRDS("/media/BioNAS2/Tensor_Omics/Protein_Function_Prediction/material/kidera_factors_uniref50_morethan1_ref_prots.rds")
kidera_ref <- kidera %>% select(id, starts_with("K")) %>% column_to_rownames("id") %>% as.matrix()

svd_ref <- read_csv(
  "/media/BioNAS2/Tensor_Omics/Protein_Function_Prediction/results/svd_reduced_matrix.csv",
  col_names = c("protein_id", as.character(0:17)),
  col_types = paste0("c", paste(rep("d", 18), collapse = ""))
) %>%
  column_to_rownames("protein_id") %>%
  as.matrix()

get_protein_tokens <- function(prot_id, type = c("GO", "HRD", "IPR")) {
  type <- match.arg(type)
  if (type == "GO") {
    tokens <- protein_to_tokens_go[[prot_id]]
  } else if (type == "HRD") {
    tokens <- protein_to_tokens_hrd[[prot_id]]
  } else {
    tokens <- protein_to_tokens_ipr[[prot_id]]
  }
  if (is.null(tokens)) character(0) else tokens
}

# Hilfsfunktion zum Entfernen des UniRef50_ Präfix
strip_uniref_prefix <- function(id) {
  sub("^UniRef50_", "", id)
}

# ---- Map tokens to candidate terms (simple lookup) ----
terms_containing_tokens <- function(tokens, annotation_df, max_terms = 200) {
  # Find annotation terms that contain any of tokens (word boundaries) - case insensitive
  if (length(tokens) == 0) return(character(0))
  pattern <- paste0("\\b(", paste0(sprintf("%s", tokens), collapse = "|"), ")\\b")
  # search in distinct terms
  candidate_terms <- annotation_df %>%
    distinct(term) %>%
    filter(str_detect(tolower(term), pattern)) %>%
    pull(term)
  if (length(candidate_terms) > max_terms) candidate_terms <- candidate_terms[1:max_terms]
  candidate_terms
}

# ---- Main function ----
predict_functions_token_based <- function(protein_list, blast_results, kidera,
                                          top_tokens = 50,
                                          debug_dir = NULL) {
  all_preds <- list()
  intermediate <- list()
  
  total <- length(protein_list)
  message("\nStarting token-based function prediction for ", total, " proteins\n")
  
  for (i in seq_along(protein_list)) {
    q <- protein_list[i]
    message("Processing protein ", i, "/", total, ": ", q)
    debug_info <- list(protein_id = q)
    
    hits <- blast_results %>%
      filter(query == q) %>%
      mutate(
        gm_sim = compute_gm_sim(pident, align_length, q_len, s_len)
      ) %>%
      filter(!is.na(gm_sim) & gm_sim > 0) %>%
      arrange(desc(gm_sim)) %>%
      head(30)
    
    debug_info$hits <- hits
    
    if (nrow(hits) == 0) {
      message(" - No valid hits found")
      all_preds[[q]] <- tibble(protein_id = q, predicted_token = character(), predicted_term = character(), source = character())
      intermediate[[q]] <- debug_info
      next
    }
    
    # weights
    gm_sim_log <- log1p(hits$gm_sim)
    weight_sum <- sum(gm_sim_log)
    if (weight_sum < .Machine$double.eps) {
      message(" - Weight sum near zero, using uniform weights")
      w <- rep(1/nrow(hits), nrow(hits))
    } else {
      w <- gm_sim_log / weight_sum
    }
    
    # For each annotation type do token-based aggregation
    annotation_types <- c("GO", "HRD", "IPR")
    token_results <- list()
    term_results <- list()
    
    for (atype in annotation_types) {
      if (atype == "GO") {
        token_mat <- go_token_mat
        token_idf <- go_token_idf
        ann_df <- ann_goa
        prot_to_tokens <- protein_to_tokens_go
      } else if (atype == "HRD") {
        token_mat <- hrd_token_mat
        token_idf <- hrd_token_idf
        ann_df <- ann_hrd
        prot_to_tokens <- protein_to_tokens_hrd
      } else {
        token_mat <- ipr_token_mat
        token_idf <- ipr_token_idf
        ann_df <- ann_ipr
        prot_to_tokens <- protein_to_tokens_ipr
      }
      
      # If no token matrix available, skip token-based for this type
      if (is.null(token_mat) || nrow(token_mat) == 0) {
        message(" - No token matrix for ", atype, "; skipping token-based for this type")
        token_results[[atype]] <- character(0)
        term_results[[atype]] <- character(0)
        next
      }
      
      # Build protein-level token vectors for hits using mapping
      prot_token_vecs <- list()
      hit_ids <- hits$subject
      hit_ids_stripped <- strip_uniref_prefix(hit_ids)
      
      # Debug: Zeige IDs, die nicht im Mapping sind
      missing_in_mapping <- hit_ids_stripped[!hit_ids_stripped %in% names(prot_to_tokens)]
      if (length(missing_in_mapping) > 0) {
        message("DEBUG: Folgende Hit-IDs fehlen im Token-Mapping für ", atype, ": ", paste(missing_in_mapping, collapse = ", "))
      }
      
      # Debug: Zeige IDs, die nicht in der Token-Matrix sind
      missing_in_matrix <- hit_ids_stripped[!hit_ids_stripped %in% rownames(token_mat)]
      if (length(missing_in_matrix) > 0) {
        message("DEBUG: Folgende Hit-IDs fehlen in der Token-Matrix für ", atype, ": ", paste(missing_in_matrix, collapse = ", "))
      }
      
      for (h_idx in seq_along(hit_ids)) {
        hid <- hit_ids_stripped[h_idx]
        tokens <- prot_to_tokens[[hid]]
        if (is.null(tokens) || length(tokens) == 0) {
          prot_token_vecs[[hid]] <- numeric(0)
        } else {
          # Summiere die Vektoren der Tokens aus der Matrix
          valid_tokens <- intersect(tokens, rownames(token_mat))
          if (length(valid_tokens) == 0) {
            prot_token_vecs[[hid]] <- numeric(0)
          } else {
            vec <- colSums(token_mat[valid_tokens, , drop = FALSE])
            prot_token_vecs[[hid]] <- vec
          }
        }
      }
      
      # Some hits may have no token vectors. Filter them and re-normalize weights accordingly
      available_hits <- names(prot_token_vecs)[sapply(prot_token_vecs, function(x) length(x) > 0)]
      # Passe available_hits ebenfalls an (falls nötig)
      idx_in_hits <- match(available_hits, hit_ids_stripped)
      w_sub <- w[idx_in_hits]
      w_sub_sum <- sum(w_sub)
      if (w_sub_sum <= 0) {
        w_sub <- rep(1/length(w_sub), length(w_sub))
      } else {
        w_sub <- w_sub / w_sub_sum
      }
      
      # Ensure all prot token vectors have same length as token_mat columns
      dim_tok <- ncol(token_mat)
      prot_token_matrix <- matrix(0, nrow = length(available_hits), ncol = dim_tok)
      rownames(prot_token_matrix) <- available_hits
      for (j in seq_along(available_hits)) {
        vec <- prot_token_vecs[[available_hits[j]]]
        if (length(vec) != dim_tok) {
          # if dimensions differ (shouldn't), try to pad/truncate
          if (length(vec) == 0) {
            vec <- rep(0, dim_tok)
          } else {
            vec <- head(c(vec, rep(0, dim_tok)), dim_tok)
          }
        }
        prot_token_matrix[j, ] <- vec
      }
      
      # Query token vector: weighted sum of protein token vectors
      q_token_vec <- colSums(prot_token_matrix * matrix(w_sub, nrow = nrow(prot_token_matrix), ncol = dim_tok))
      
      # If q vector is zero, skip similarity
      if (all(q_token_vec == 0)) {
        message(" - Query token vector zero for ", atype)
        token_results[[atype]] <- character(0)
        term_results[[atype]] <- character(0)
        next
      }
      
      # Compute cosine similarity between query token vector and all tokens in token_mat
      denom_q <- sqrt(sum(q_token_vec^2)); if (denom_q == 0) denom_q <- .Machine$double.eps
      denom_tokens <- sqrt(rowSums(token_mat^2)); denom_tokens[denom_tokens == 0] <- .Machine$double.eps
      sims_tokens <- as.vector((token_mat %*% q_token_vec) / (denom_tokens * denom_q))
      names(sims_tokens) <- rownames(token_mat)
      sims_tokens <- sort(sims_tokens, decreasing = TRUE)
      
      top_k <- min(top_tokens, length(sims_tokens))
      top_tokens_names <- names(sims_tokens)[1:top_k]
      top_tokens_scores <- sims_tokens[1:top_k]
      
      token_results[[atype]] <- tibble(token = top_tokens_names, score = as.numeric(top_tokens_scores))
      
      # Map tokens to candidate terms (annotation strings)
      candidates <- terms_containing_tokens(top_tokens_names, ann_df, max_terms = 500)
      # Optionally we can rank candidate terms by frequency among neighbors or by aggregated token scores
      # For simplicity, compute a basic score: how many top tokens appear in the term (weighted)
      if (length(candidates) == 0) {
        term_results[[atype]] <- character(0)
      } else {
        token_score_named <- setNames(as.numeric(top_tokens_scores), top_tokens_names)
        term_scores <- sapply(candidates, function(term) {
          tkns <- unique(tokenize_text(term))
          common <- intersect(tkns, top_tokens_names)
          if (length(common) == 0) return(0)
          sum(token_score_named[common], na.rm = TRUE)
        })
        ranked_terms <- candidates[order(term_scores, decreasing = TRUE)]
        term_results[[atype]] <- ranked_terms
      }
      
      # Save for debugging
      debug_info[[paste0("tokens_", atype)]] <- token_results[[atype]]
      debug_info[[paste0("terms_", atype)]] <- term_results[[atype]]
    } # end for annotation types
    
    # Additionally compute doc-based predictions (original approach) for comparison
    # (keep original neighbors logic)
    # compute q_doc vectors (weighted sums) if doc matrices available
    compute_qdoc <- function(ref_mat) {
      if (is.null(ref_mat)) return(numeric(0))
      if (!is.matrix(ref_mat) && !is.data.frame(ref_mat)) return(numeric(0))
      if (nrow(ref_mat) == 0) return(numeric(0))
      # Entferne Präfix bei den Hit-IDs
      ids <- intersect(strip_uniref_prefix(hits$subject), rownames(ref_mat))
      if (length(ids) == 0) return(numeric(0))
      V <- ref_mat[ids, , drop = FALSE]
      idx <- match(rownames(V), strip_uniref_prefix(hits$subject))
      colSums(V * w[idx])
    }
    q_go  <- compute_qdoc(go_ref)
    q_hrd <- compute_qdoc(hrd_ref)
    q_ipr <- compute_qdoc(ipr_ref)
    q_kid <- compute_qdoc(kidera_ref)
    q_svd <- compute_qdoc(svd_ref)
    
    # find neighbors (doc-based) using existing stop condition
    get_tops <- function(q_vec, ref_mat) {
      if (length(q_vec) == 0 || is.null(ref_mat) || nrow(ref_mat) == 0) return(character(0))
      get_top_ids_with_stop(q_vec, ref_mat)
    }
    tops <- list(
      GO_doc = get_tops(q_go, go_ref),
      HRD_doc = get_tops(q_hrd, hrd_ref),
      IPR_doc = get_tops(q_ipr, ipr_ref),
      KID_doc = get_tops(q_kid, kidera_ref),
      SVD_doc = get_tops(q_svd, svd_ref)
    )
    
    # get predictions from doc neighbors (original mapping)
    preds_doc <- list(
      GO = unique(c(get_predictions(tops$GO_doc, ann_goa))),
      HRD = unique(c(get_predictions(tops$HRD_doc, ann_hrd))),
      IPR = unique(c(get_predictions(tops$IPR_doc, ann_ipr)))
    )
    
    # Format outputs: we will produce both token-based top tokens + mapped terms and doc-term predictions
    # Create tidy tibble rows for the query:
    # Token-based: one row per (atype, token, token_score, mapped_terms...) - for simplicity we return joined mapped terms as a single string
    token_rows <- bind_rows(lapply(names(token_results), function(at) {
      tr <- token_results[[at]]
      # Robustere Prüfung:
      if (is.null(tr) || !is.data.frame(tr) || nrow(tr) == 0) {
        tibble(protein_id = q, annotation_type = at, token = character(0), token_score = numeric(0), mapped_terms = character(0))
      } else {
        mapped_terms_str <- sapply(seq_len(nrow(tr)), function(r) {
          tkn <- tr$token[r]
          # terms containing this single token
          tms <- terms_containing_tokens(tkn, if (at == "GO") ann_goa else if (at == "HRD") ann_hrd else ann_ipr, max_terms = 50)
          paste0(head(tms, 10), collapse = " | ")
        })
        tibble(
          protein_id = q,
          annotation_type = at,
          token = tr$token,
          token_score = tr$score,
          mapped_terms = mapped_terms_str
        )
      }
    }))
    
    # Doc-based rows
    doc_rows <- bind_rows(list(
      GO = tibble(protein_id = q, annotation_type = "GO_doc", predicted_term = preds_doc$GO),
      HRD = tibble(protein_id = q, annotation_type = "HRD_doc", predicted_term = preds_doc$HRD),
      IPR = tibble(protein_id = q, annotation_type = "IPR_doc", predicted_term = preds_doc$IPR)
    ), .id = "tmp") %>% select(-tmp)
    
    # Save results
    # final predictions: tokens table + doc_terms table (we return both as separate outputs glued into longer frame)
    # For convenience, we'll create a 'predicted_term' column in token_rows by using first mapped term if available
    token_rows <- token_rows %>%
      mutate(predicted_term = ifelse(mapped_terms == "", NA_character_, mapped_terms))
    
    combined_rows <- bind_rows(
      token_rows %>% select(protein_id, annotation_type, token, token_score, predicted_term),
      doc_rows %>% rename(annotation_type = annotation_type, token = predicted_term) %>%
        mutate(token_score = NA_real_) %>%
        select(protein_id, annotation_type, token, token_score, predicted_term = token)
    )
    
    all_preds[[q]] <- combined_rows
    intermediate[[q]] <- debug_info
    message(" - Token-based prediction produced ", nrow(combined_rows), " rows (including doc-based rows)")
  } # end loop proteins
  
  # Save intermediate results if requested
  if (!is.null(debug_dir)) {
    dir.create(debug_dir, showWarnings = FALSE)
    saveRDS(intermediate, file.path(debug_dir, "intermediate_results_token_based.rds"))
    message("\nSaved intermediate results to ", file.path(debug_dir, "intermediate_results_token_based.rds"))
  }
  
  # bind rows while preserving columns
  final <- bind_rows(all_preds)
  final
}

# ---- Load and prepare BLAST results ----
blast_results <- read_tsv(
  "/media/BioNAS2/Tensor_Omics/Protein_Function_Prediction/results/diamond_test_1.txt",
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

# ---- Run token-based prediction ----
query_proteins <- unique(blast_results$query)
predictions_token <- predict_functions_token_based(
  query_proteins,
  blast_results,
  kidera,
  top_tokens = 50,
  debug_dir = "debug_results_token"  # set to NULL to disable saving intermediate debug data
)

# ---- Save outputs ----
write_csv(predictions_token, "/media/BioNAS2/Tensor_Omics/Protein_Function_Prediction/results/predicted_functions_token_based.csv")
message("\nSaved token-based predictions to predicted_functions_token_based.csv")

