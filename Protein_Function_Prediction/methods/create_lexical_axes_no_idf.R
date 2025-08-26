# Version 1.4.4
# Author: Aaron Schroeder
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(purrr)
library(Matrix) 
library(tidytext)
library(text2vec)
library(irlba)
library(progressr)
library(pbapply)

# Global progressbar
handlers(global = TRUE)
handlers("txtprogressbar")

# 1. Load annotation tables
ann_hrd <- readRDS("/media/BioNAS2/Tensor_Omics/Protein_Function_Prediction/material/protein2hrd.rds")
ann_goa <- readRDS("/media/BioNAS2/Tensor_Omics/Protein_Function_Prediction/material/protein2goa.rds")
ann_ipr <- readRDS("/media/BioNAS2/Tensor_Omics/Protein_Function_Prediction/material/protein2ipr.rds")

data("stop_words")

# 2. Create corpora from the annotation tables
corpora <- list(
  GO = ann_goa,
  IPR = ann_ipr,
  HRD = ann_hrd
)

clean_doc <- function(text) {
  text %>%
    str_replace_all("[\\[\\]()]", " ") %>%
    str_replace_all("--+", "-") %>%
    str_replace_all(",", " ") %>%
    str_squish()
}

tokenize_custom <- function(protein_ids, texts, pattern, name = NULL) {
  if (name == "GO") {
    tokens_list <- str_split(texts, pattern)
  } else {
    tokens_list <- str_split(texts, "[\\s.;,]+")
  }

  tokens_list <- map(tokens_list, ~ .x[.x != ""])
  
  message("[", name, "] Tokens before filtering: ", sum(lengths(tokens_list)))
  message("[", name, "] Unique Tokens after filtering: ", length(unique(unlist(tokens_list))))

  tibble(
    protein_id = protein_ids,
    tokens_list = tokens_list
  ) %>%
    unnest(tokens_list) %>%
    rename(token = tokens_list) %>%
    mutate(token = str_to_lower(str_trim(token))) %>%
    filter(
      str_length(token) > 2,
      str_detect(token, "[a-z]")
    )
}

splitter_go <- "[\\s.,:]+"

# create axes without tf-idf weighting
build_axes <- function(docs, name) {
  message("Processing corpus: ", name, " (", length(docs[[2]]), " docs)")
  
  protein_ids <- docs[[1]]
  texts <- docs[[2]]

  # Remove NA documents
  valid_docs <- !is.na(texts)
  texts <- texts[valid_docs]
  protein_ids <- protein_ids[valid_docs]

  # Preprocess
  texts <- clean_doc(texts)

  # Tokenize
  p <- progressor(steps = 6)
  p("Tokenizing")
  if (name == "GO") {
    dt <- tokenize_custom(protein_ids, texts, splitter_go, name)
  } else {
    dt <- tokenize_custom(protein_ids, texts, NULL, name)  # Pattern nicht benötigt für HRD/IPR
  }

  message("[", name, "] Tokens nach Filtern: ", nrow(dt))
  message("[", name, "] Unique Tokens nach Filtern: ", n_distinct(dt$token))

  p("Counting term frequencies")
  wf_doc <- dt %>% 
    distinct(protein_id, token) %>% 
    mutate(tf = 1)

  p("Building co-occurrence matrix")
  it <- wf_doc %>% cast_sparse(protein_id, token, value = "tf")
  cooc <- crossprod(it)
  
  # 500 dims needed since it is slightly more general than when tf-idf is used
  nv <- min(500, nrow(cooc) - 1, ncol(cooc) - 1)
  if (nv < 1) stop("Co-occurrence matrix too small for SVD.")
  sv <- irlba(cooc, nv = nv)

  total_variance <- sum(sv$d^2)

  explained_variance <- sv$d^2 / total_variance
  cumulative_variance <- cumsum(explained_variance)

  k <- min(which(cumulative_variance >= 0.85))

  # use slightly more dimensions due to generalization
  if (name == "GO"){
	  k <- max(k, 5)
  }
  else if (name == "IPR") {
	  k <- max(k, 3)
  } else {
	  k <- max(k, 3)
  }

  k <- min(k, nv)

  message("optimal dimensions for ", name, ": ", k, " explain ",
	  round(cumulative_variance[k]*100, 1), "% variance")

  U <- sv$u[, 1:k, drop = FALSE] 
  colnames(U) <- paste0("EW", 1:k)
  rownames(U) <- colnames(cooc)

  p("projecting documents")

  tokens_final <- rownames(U)
  Uf <- U[tokens_final, , drop = FALSE]

  mat_binary <- cast_sparse(wf_doc, protein_id, token, value="tf")
  doc_proj <- mat_binary %*% Uf
  colnames(doc_proj) <- paste0("EW", 1:k)
  
  list(
    word_axes = Uf,
    doc_vectors = doc_proj,
    explained_variance = explained_variance
  )
}

# Run with progress display
with_progress({
  axes <- map2(corpora, names(corpora), build_axes)
})

# Save result
saveRDS(axes$GO,  "/media/BioNAS2/Tensor_Omics/Protein_Function_Prediction/material/lexical_axes_GO_no_idf.rds")
saveRDS(axes$IPR, "/media/BioNAS2/Tensor_Omics/Protein_Function_Prediction/material/lexical_axes_IPR_no_idf.rds")
saveRDS(axes$HRD, "/media/BioNAS2/Tensor_Omics/Protein_Function_Prediction/material/lexical_axes_HRD_no_idf.rds")

message("Done. Lexical axes built for GO, InterPro & HRDs.")
