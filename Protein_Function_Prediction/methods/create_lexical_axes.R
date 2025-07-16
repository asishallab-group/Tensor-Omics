library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(purrr)
library(Matrix) 
library(tidytext)
library(text2vec)
library(irlba)

# 1. Load annotation tables
ann_hrd <- readRDS("/media/BioNAS2/Tensor_Omics/Protein_Function_Prediction/material/protein2hrd.rds")
ann_goa <- readRDS("/media/BioNAS2/Tensor_Omics/Protein_Function_Prediction/material/protein2goa.rds")
ann_ipr <- readRDS("/media/BioNAS2/Tensor_Omics/Protein_Function_Prediction/material/protein2ipr.rds")

# 2. Create corpora from the annotation tables
corpora <- list(
  GO = ann_goa$GOA_label,
  IPR = ann_ipr$description,
  HRD = ann_hrd$hrd
)

clean_go <- function(text) {
  text %>%
    str_replace_all("[\\[\\]()]", " ") %>%
    str_replace_all("--+", "-") %>%
    str_squish()
}

clean_ipr <- function(text) {
  text %>%
    str_replace_all("[\\[\\]()]", " ") %>%
    str_replace_all("--+", "-") %>%
    str_replace_all(",", " ") %>%
    str_squish()
}

clean_hrd <- function(text) {
  text %>%
    str_replace_all("[\\[\\]()]", " ") %>%
    str_replace_all("--+", "-") %>%
    str_squish()
}

# Prot‑scriber splitter regex (specialized for HRDs)
splitter <- "(?<=\\b)[A-Za-z0-9/-]+(?:'[A-Za-z0-9]+)?(?=\\b)"

# Function to build axes for one corpus
build_axes <- function(docs, name) {
  message("Processing corpus: ", name, " (", length(docs), " docs)")
  
  # Remove NA documents
  valid_docs <- !is.na(docs)
  docs <- docs[valid_docs]
  
  # Corpus-spezifisches Preprocessing
  docs <- switch(name,
    GO = clean_go(docs),
    IPR = clean_ipr(docs),
    HRD = clean_hrd(docs),
    stop("Unknown corpus: ", name)
  )

  # if existent, remove last identifier from HRD descriptions
  if (name == "HRD") {
    docs <- sapply(docs, function(txt) {
      parts <- str_split(txt, " ", simplify = TRUE)
      last <- tail(parts, 1)
      if (is_identifier(last)) {
        str_trim(str_remove(txt, paste0("\\s*", last, "$")))
      } else {
        txt
      }
    })
  }

  # 3. Tokenize & lowercase, then filter
  dt <- tibble(doc_id = seq_along(docs), text = docs) %>%
    unnest_tokens(token, text, token = "regex", pattern = splitter) %>%
    mutate(token = str_to_lower(token)) %>%
    filter(
      str_length(token) > 2,                    # length > 2
      str_detect(token, "[a-z]")                # contains at least one letter ot remove tokens like "[(-"
    )

  # 4. Word frequencies per document & global f_w
  wf_doc <- dt %>% count(doc_id, token, name = "tf")
  fw <- wf_doc %>% group_by(token) %>% summarise(fw = sum(tf))
  
  total_fw <- sum(fw$fw)
  
  # 5. Shannon information content IC(token) = -log(fw/total_fw)
  fw <- fw %>% mutate(IC = -log(fw / total_fw))
  
  # 6. Drop bottom 10% by IC
  cutoff_ic <- quantile(fw$IC, 0.10, na.rm = TRUE)
  keep_tokens <- fw %>% filter(IC > cutoff_ic) %>% pull(token)
  
  # 7. Build co-occurrence matrix (token × token)
  dt_f <- wf_doc %>% filter(token %in% keep_tokens)
  vocab <- unique(dt_f$token)
  
  # Create a token-by-doc matrix
  it <- dt_f %>% cast_sparse(doc_id, token)
  # Co-occurrence = t(it) %*% it
  cooc <- crossprod(it)
  
  # 8. SVD on cooc, keep components explaining >80% variance
  nv <- min(50, nrow(cooc) - 1, ncol(cooc) - 1)
  if (nv < 1) stop("Co-occurrence matrix too small for SVD.")
  sv <- irlba(cooc, nv = nv)
  vars <- sv$d^2 / sum(sv$d^2)
  k <- min(which(cumsum(vars) >= 0.80))

  # Eigenwords
  U <- sv$u[, 1:k, drop = FALSE] 
  # Eigenword vectors
  colnames(U) <- paste0("EW", 1:k)
  rownames(U) <- vocab
  
  # 9. Compute IDF for each token
  ndocs <- length(unique(dt$doc_id))
  doc_freq <- dt %>% distinct(doc_id, token) %>%
    count(token, name = "df")
  idf <- doc_freq %>% 
    mutate(idf = log(ndocs / df)) %>%
    filter(df < 0.95 * ndocs)   # Drop words in ≥95% docs
  
  # Intersect tokens with both U and idf
  tokens_final <- intersect(rownames(U), idf$token)
  Uf <- U[tokens_final, , drop = FALSE]
  idf_f <- idf %>% filter(token %in% tokens_final)
  
  # 10. Project documents into eigenword-space by TF-IDF weighted sum
  tfidf <- wf_doc %>% 
    filter(token %in% tokens_final) %>%
    left_join(idf_f, by = "token") %>%
    mutate(tf_idf = tf * idf) %>%
    select(doc_id, token, tf_idf)
  
  # Build sparse matrix and project
  mat_tfidf <- cast_sparse(tfidf, doc_id, token, value = "tf_idf")
  doc_proj <- mat_tfidf %*% Uf
  
  colnames(doc_proj) <- paste0("EW", 1:k)
  rownames(doc_proj) <- paste0(name, "_doc", rownames(doc_proj))
  
  list(
    word_axes = Uf,
    idf = idf_f,
    doc_vectors = doc_proj,
    explained_variance = vars[1:k]
  )
}

# Run for all three corpora
axes <- map2(corpora, names(corpora), build_axes)

# Save result
saveRDS(axes, "/media/BioNAS2/Tensor_Omics/Protein_Function_Prediction/material/lexical_axes_GO_IPR_HRD.rds")

cat("Done. Lexical axes built for GO, InterPro & HRDs.")