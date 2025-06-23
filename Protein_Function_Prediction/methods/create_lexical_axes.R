library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(purrr)
library(Matrix) 
library(tidytext)
library(text2vec)
library(irlba)

# 1. load annotation table
ann <- readRDS("./material/annotation_table.rds")

# 2. Create each corpus to handle them indipendently
# GOA corpus: one doc per GO term name
corpus_GO <- ann$GOA_label

# IPR corpus
corpus_IPR <- ann$description

# HRD corpus
corpus_HRD <- ann$hrd

corpora <- list(GO = corpus_GO,
                IPR = corpus_IPR,
                HRD = corpus_HRD)

# prot‑scriber splitter regex (specialized for HRDs)
splitter <- "(?<=\\b)[A-Za-z0-9]+(?:'[A-Za-z0-9]+)?(?=\\b)"

# a function to build axes for one corpus
build_axes <- function(docs, name) {
  message("Processing corpus: ", name, " (", length(docs), " docs)")

  # 3. tokenize & lowercase
  dt <- tibble(doc_id = seq_along(docs), text = docs) %>%
    unnest_tokens(token, text, token = "regex", pattern = splitter) %>%
    mutate(token = str_to_lower(token)) %>%
    filter(str_length(token) > 1)

  # 4. word frequencies per document & global f_w
  wf_doc <- dt %>% count(doc_id, token, name = "tf")
  fw <- wf_doc %>% group_by(token) %>% summarise(fw = sum(tf))

  total_fw <- sum(fw$fw)

  # 5. Shannon information content IC(token) = -log(fw/total_fw)
  fw <- fw %>% mutate(IC = -log(fw / total_fw))

  # 6. drop bottom 10% by IC
  cutoff_ic <- quantile(fw$IC, 0.10)
  keep_tokens <- fw %>% filter(IC > cutoff_ic) %>% pull(token)

  # 7. build co‐occurrence matrix (token × token)
  #    context = “document”: co‐occur if in same doc
  dt_f <- wf_doc %>% filter(token %in% keep_tokens)
  vocab <- unique(dt_f$token)
  
  # create a token‐by‐doc matrix
  it <- dt_f %>% cast_sparse(doc_id, token)
  # co‐occurrence = t(it) %*% it
  cooc <- crossprod(it)

  # 8. SVD on cooc, keep components explaining >80% variance
  #    use irlba for speed
  nv <- min(50, nrow(cooc) - 1, ncol(cooc) - 1)
  if (nv < 1) stop("Co-occurrence matrix too small for SVD.")
  sv <- irlba(cooc, nv = nv)
  vars <- sv$d^2 / sum(sv$d^2)
  k <- min(which(cumsum(vars) >= 0.80))
  U  <- sv$u[, 1:k]              # eigenword vectors
  colnames(U) <- paste0("EW",1:k)

  rownames(U) <- vocab

  # 9. compute IDF for each token
  #    IDF = log(N_docs / n_docs(token))
  ndocs <- length(unique(dt$doc_id))
  doc_freq <- dt %>% distinct(doc_id, token) %>%
    count(token, name="df")
  idf <- doc_freq %>% 
    mutate(idf = log(ndocs / df)) %>%
    filter(df < 0.95 * ndocs)         # drop words in ≥95% docs

  # intersect tokens with both U and idf
  tokens_final <- intersect(rownames(U), idf$token)
  Uf <- U[tokens_final, , drop=FALSE]
  idf_f <- idf %>% filter(token %in% tokens_final)

  # 10. project each document into eigenword‐space by TF‐IDF weighted sum
  #     result: a doc×k matrix
  #    ‑ for each doc, collect tf for tokens_final, then tf*idf * Uf
  tfidf <- wf_doc %>% 
    filter(token %in% tokens_final) %>%
    left_join(idf_f, by="token") %>%
    mutate(tf_idf = tf * idf) %>%
    select(doc_id, token, tf_idf)

  # build a sparse df: doc × token with tf_idf, then %*% Uf
  mat_tfidf <- cast_sparse(tfidf, doc_id, token, value = "tf_idf")
  # project: doc×token %*% token×k → doc×k
  doc_proj <- mat_tfidf %*% Uf

  colnames(doc_proj) <- paste0("EW",1:k)
  rownames(doc_proj) <- paste0(name, "_doc", rownames(doc_proj))

  list(
    word_axes = Uf,
    idf        = idf_f,
    doc_vectors = doc_proj,
    explained_variance = vars[1:k]
  )
}

# run for all three corpora
axes <- map2(corpora, names(corpora), build_axes)

# save result
saveRDS(axes, "lexical_axes_GO_IPR_HRD.rds")

message("Done.  Lexical axes built for GO, InterPro & HRDs.")