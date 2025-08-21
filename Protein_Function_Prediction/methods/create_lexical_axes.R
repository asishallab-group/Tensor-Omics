# Version 1.4.3 - Ohne Filterung, mit angepasster Tokenisierung
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

# Aktiviert globale Fortschrittsanzeige
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

# Neue Tokenizer-Funktion für HRD und IPR (an Leerzeichen, Semikolon, Punkten)
tokenize_custom <- function(protein_ids, texts, pattern, name = NULL) {
  if (name == "GO") {
    tokens_list <- str_split(texts, pattern)
  } else {
    # Für HRD und IPR an typischen Trennzeichen splitten
    tokens_list <- str_split(texts, "[\\s.;,]+")
  }

  # Entferne leere Tokens
  tokens_list <- map(tokens_list, ~ .x[.x != ""])
  
  message("[", name, "] Tokens vor Filtern: ", sum(lengths(tokens_list)))
  message("[", name, "] Unique Tokens vor Filtern: ", length(unique(unlist(tokens_list))))

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

# Funktion zur Erstellung der Achsen ohne Filterung
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
  wf_doc <- dt %>% count(protein_id, token, name = "tf")

  p("Building co-occurrence matrix")
  it <- wf_doc %>% cast_sparse(protein_id, token, value = "tf")
  cooc <- crossprod(it)

  p("SVD decomposition")
  nv <- min(50, nrow(cooc) - 1, ncol(cooc) - 1)
  if (nv < 1) stop("Co-occurrence matrix too small for SVD.")
  sv <- irlba(cooc, nv = nv)
  vars <- sv$d^2 / sum(sv$d^2)
  k <- min(which(cumsum(vars) >= 0.90))

  U <- sv$u[, 1:k, drop = FALSE] 
  colnames(U) <- paste0("EW", 1:k)
  rownames(U) <- colnames(cooc)

  p("Computing IDF and projecting documents")
  ndocs <- length(unique(dt$protein_id))
  doc_freq <- dt %>% distinct(protein_id, token) %>%
    count(token, name = "df")
  idf <- doc_freq %>% 
    mutate(idf = log(ndocs / df))

  tokens_final <- rownames(U)
  Uf <- U[tokens_final, , drop = FALSE]
  idf_f <- idf %>% filter(token %in% tokens_final)

  tfidf <- wf_doc %>% 
    left_join(idf_f, by = "token") %>%
    mutate(tf_idf = tf * idf) %>%
    select(protein_id, token, tf_idf)

  mat_tfidf <- cast_sparse(tfidf, protein_id, token, value = "tf_idf")
  doc_proj <- mat_tfidf %*% Uf
  colnames(doc_proj) <- paste0("EW", 1:k)
  
  list(
    word_axes = Uf,
    idf = idf_f,
    doc_vectors = doc_proj,
    explained_variance = vars[1:k]
  )
}

# Run with progress display
with_progress({
  axes <- map2(corpora, names(corpora), build_axes)
})

# Save result
saveRDS(axes$GO,  "/media/BioNAS2/Tensor_Omics/Protein_Function_Prediction/material/lexical_axes_GO.rds")
saveRDS(axes$IPR, "/media/BioNAS2/Tensor_Omics/Protein_Function_Prediction/material/lexical_axes_IPR.rds")
saveRDS(axes$HRD, "/media/BioNAS2/Tensor_Omics/Protein_Function_Prediction/material/lexical_axes_HRD.rds")

message("Done. Lexical axes built for GO, InterPro & HRDs.")
