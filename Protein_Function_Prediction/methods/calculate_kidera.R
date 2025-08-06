# Pakete laden
library(Biostrings)     # FASTA effizient lesen
library(Peptides)       # Kidera-Faktor Berechnung
library(future.apply)   # Parallele Verarbeitung
library(tibble)         # Tidy Dataframe
library(dplyr)          # Filter / Mutate

# Anzahl der Worker setzen
plan(multisession, workers = 8)

# Eingabe-FASTA
fasta_path <- "/media/BioNAS2/Tensor_Omics/Protein_Function_Prediction/material/uniref50_morethan1.fasta"

# Sequenzen laden
aa <- readAAStringSet(fasta_path)
ids <- names(aa)
seqs <- as.character(aa)

# IDs bereinigen (nur erster Teil vor Leerzeichen)
ids_clean <- sub(" .*", "", ids)

# Funktion zum Entfernen nicht-Standard-AAs
clean_sequence <- function(seq) {
  gsub("[^ACDEFGHIKLMNPQRSTVWY]", "", seq)
}

# Safe-Kidera mit Reinigung
safe_kidera <- function(seq) {
  seq <- clean_sequence(seq)
  tryCatch(
    as.numeric(kidera(seq)),
    error = function(e) rep(NA_real_, 10)
  )
}

# Kidera-Faktoren in Parallelberechnung
kidera_list <- future_lapply(seqs, safe_kidera)

# Liste in Matrix umwandeln
kidera_matrix <- do.call(rbind, kidera_list)
colnames(kidera_matrix) <- paste0("k", 1:10)

# In Tibble umwandeln + IDs hinzufügen
kidera_df <- as_tibble(kidera_matrix) %>%
  mutate(id = ids_clean, .before = 1)

# Entferne Reihen mit NA (falls Sequenz nach Reinigung leer war)
kidera_df <- kidera_df %>%
  filter(if_all(starts_with("k"), ~ !is.na(.x)))

# Als .rds speichern
saveRDS(
  kidera_df,
  "/media/BioNAS2/Tensor_Omics/Protein_Function_Prediction/material/kidera_factors_uniref50_morethan1_ref_prots.rds"
)

message("Fertig! ", nrow(kidera_df), " Sequenzen verarbeitet.")

