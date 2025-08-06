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

is_identifier <- function(token) {
  if (length(token) != 1) stop("Token muss Länge 1 haben")
  
  token <- tolower(token)
  has_special_chars <- grepl("[0-9/-]", token)
  too_short <- nchar(token) < 3
  not_stopword <- !token %in% stop_words$word
  
  return(not_stopword && (has_special_chars || too_short))
}

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

# Optimierte Tokenizer-Funktion ohne rep()
tokenize_custom <- function(protein_ids, texts, pattern, name = NULL) {
  if (name == "GO") {
    tokens_list <- str_split(texts, pattern)
  } else {
    tokens_list <- str_extract_all(texts, pattern)
  }

  tibble(
    protein_id = protein_ids,  # nur einmal pro Dokument
    tokens_list = tokens_list
  ) %>%
    unnest(tokens_list) %>%
    rename(token = tokens_list) %>%
    mutate(token = str_to_lower(str_trim(token))) %>%
    filter(
      str_length(token) > 2,
      str_detect(token, "[a-z]"),
      !str_detect(token, "^[0-9]{2,}[a-z]$")
    )
}

splitter <- "[0-9]*'?-?[A-Za-z]+(?:['/-][0-9A-Za-z]+)*"
splitter_go <- "[\\s.,:]+"

# Function to build axes for one corpus
build_axes <- function(docs, name) {
  message("Processing corpus: ", name, " (", length(docs), " docs)")
  
  protein_ids <- docs[[1]]
  texts <- docs[[2]]

  # Remove NA documents
  valid_docs <- !is.na(texts)
  texts <- texts[valid_docs]
  protein_ids <- protein_ids[valid_docs]

  # Preprocess
  texts <- clean_doc(texts)

  # HRD spezifisches Entfernen von IDs
  if (name == "HRD") {
    texts <- pbsapply(texts, function(txt) {
      parts <- str_split(txt, " ")[[1]]
      last <- tail(parts, 1)
      if (is_identifier(last) && length(parts) > 4) {
        str_trim(str_remove(txt, paste0("\\s*", last, "$")))
      } else {
        txt
      }
    })
  }

  # Tokenize
  p <- progressor(steps = 6)
  p("Tokenizing")
  if (name == "GO") {
    dt <- tokenize_custom(protein_ids, texts, splitter_go, name)
  } else {
    dt <- tokenize_custom(protein_ids, texts, splitter, name)
  }

  p("Counting term frequencies")
  wf_doc <- dt %>% count(protein_id, token, name = "tf")
  fw <- wf_doc %>% group_by(token) %>% summarise(fw = sum(tf), .groups = "drop")

  total_fw <- sum(fw$fw)
  
  p("Calculating IC and filtering bottom 10%")
  fw <- fw %>% mutate(IC = -log(fw / total_fw))
  cutoff_ic <- quantile(fw$IC, 0.10, na.rm = TRUE)
  keep_tokens <- fw %>% filter(IC > cutoff_ic) %>% pull(token)

  p("Building co-occurrence matrix")
  dt_f <- wf_doc %>% filter(token %in% keep_tokens)
  it <- dt_f %>% cast_sparse(protein_id, token, value = "tf")
  cooc <- crossprod(it)

  p("SVD decomposition")
  nv <- min(50, nrow(cooc) - 1, ncol(cooc) - 1)
  if (nv < 1) stop("Co-occurrence matrix too small for SVD.")
  sv <- irlba(cooc, nv = nv)
  vars <- sv$d^2 / sum(sv$d^2)
  k <- min(which(cumsum(vars) >= 0.80))

  U <- sv$u[, 1:k, drop = FALSE] 
  colnames(U) <- paste0("EW", 1:k)
  rownames(U) <- colnames(cooc)

  p("Computing IDF and projecting documents")
  ndocs <- length(unique(dt$protein_id))
  doc_freq <- dt %>% distinct(protein_id, token) %>%
    count(token, name = "df")
  idf <- doc_freq %>% 
    mutate(idf = log(ndocs / df)) %>%
    filter(df < 0.95 * ndocs)

  tokens_final <- intersect(rownames(U), idf$token)
  Uf <- U[tokens_final, , drop = FALSE]
  idf_f <- idf %>% filter(token %in% tokens_final)

  tfidf <- wf_doc %>% 
    filter(token %in% tokens_final) %>%
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