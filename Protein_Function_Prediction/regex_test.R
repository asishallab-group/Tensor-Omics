library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(purrr)
library(Matrix) 
library(tidytext)
library(text2vec)
library(irlba)

ann_goa <- data.frame(
  GOA_label = c("proteasome complex", "transferase activity, transferring phosphorus-containing groups",
                "transferase activity", "polysaccharide biosynthetic process", "protein transport", "membrane")
)

ann_ipr <- data.frame(
  description = c("Protein kinase domain", "Transmembrane domain", "DNA binding domain", "Zinc finger domain")
)

ann_hrd <- data.frame(
  hrd = c("DNA repair protein", "Cell cycle regulator", "Transcription factor", "Signal transduction protein")
)

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
  GO = ann_goa$GOA_label,
  IPR = ann_ipr$description,
  HRD = ann_hrd$hrd
)

clean_doc <- function(text) {
  text %>%
    str_replace_all("[\\[\\]()]", " ") %>%
    str_replace_all("--+", "-") %>%
    str_replace_all(",", " ") %>%
    str_squish()
}

tokenize_custom <- function(docs, pattern) {
  tibble(doc_id = seq_along(docs), text = docs) %>%
    mutate(tokens = str_extract_all(text, pattern)) %>%
    unnest(tokens) %>%
    rename(token = tokens) %>%
    mutate(token = str_to_lower(token)) %>%
    filter(
      str_length(token) > 2,
      str_detect(token, "[a-z]")  # nur Tokens mit Buchstaben
    )
}

splitter <- "[0-9]*'?-?[A-Za-z]+(?:['/-][0-9A-Za-z]+)*"

# Function to build axes for one corpus
build_axes <- function(docs, name) {
  message("Processing corpus: ", name, " (", length(docs), " docs)")
  
  # Remove NA documents
  valid_docs <- !is.na(docs)
  docs <- docs[valid_docs]
  
  # Corpus-spezifisches Preprocessing
  docs <- clean_doc(docs)

  # if existent, remove last identifier from HRD descriptions
  if (name == "HRD") {
    docs <- sapply(docs, function(txt) {
      parts <- str_split(txt, " ")[[1]]
      last <- tail(parts, 1)
      if (is_identifier(last) && length(parts) > 3) {
        str_trim(str_remove(txt, paste0("\\s*", last, "$")))
      } else {
        txt
      }
    })
  }

  # 3. Tokenize & lowercase, then filter
  dt <- tokenize_custom(docs, splitter)
  message("Number of unique tokens before filtering: ", n_distinct(dt$token), "\nNumber of tokens: ", length(dt$token))
  print(dt)
}

# Run for all three corpora
axes <- map2(corpora, names(corpora), build_axes)