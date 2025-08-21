#Version 1.2.7

library(Biostrings)     # Für effizientes Einlesen der FASTA
library(Peptides)       # Für kideraFactors()
library(future.apply)   # Für Parallelisierung
library(tibble)         # Für hübsche Dataframes
library(dplyr)          # Für Datenmanipulation

# Parallel-Plan with 8 Workers
plan(multisession, workers = 8)

# Path to fasta
fasta_path <- "/media/BioNAS2/Tensor_Omics/Protein_Function_Prediction/material/p_coccineus_different_root.faa"

cat("Loading fasta:", fasta_path, "\n")
aa <- readAAStringSet(fasta_path)
ids <- names(aa)
seqs <- as.character(aa)
cat("Numbers of loaded sequences:", length(seqs), "\n")

# IDs säubern: 'UniRef50_' entfernen und alles nach Leerzeichen abschneiden
ids_clean <- sub("^UniRef50_", "", sub(" .*", "", ids))

# Funktion zur sicheren Kidera-Faktoren-Berechnung mit AA-Säuberung
safe_kidera <- function(seq) {
  # Nicht-Standard-AAs entfernen (X, B, Z, etc.)
  seq_clean <- gsub("[^ACDEFGHIKLMNPQRSTVWY]", "", seq)
  if (nchar(seq_clean) == 0) {
    # Falls nach Reinigung nichts übrig bleibt: NA-Vektor zurückgeben
    return(rep(NA_real_, 10))
  }
  tryCatch(
    as.numeric(kideraFactors(seq_clean)[[1]]),
    error = function(e) rep(NA_real_, 10)
  )
}

cat("Starting kidera calculation with 8 workers...\n")
kidera_list <- future_lapply(seqs, safe_kidera)
cat("Done. Creating Matrix...\n")

kidera_matrix <- do.call(rbind, kidera_list)
colnames(kidera_matrix) <- paste0("KF", 1:10)

# Dataframe zusammenbauen mit id + Kidera-Faktoren
kidera_df <- as_tibble(kidera_matrix) %>%
  mutate(id = ids_clean, .before = 1) %>%
  filter(if_all(starts_with("KF"), ~ !is.na(.x)))

cat("Final dataframe with: ", nrow(kidera_df), "valid sequences.\n")

# Speichern
saveRDS(kidera_df, "/media/BioNAS2/Tensor_Omics/Protein_Function_Prediction/results/query_kidera.rds")
cat("Data saved.\n")

