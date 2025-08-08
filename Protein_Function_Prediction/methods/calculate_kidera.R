# Pakete laden
library(Biostrings)
library(Peptides)
library(future.apply)
library(tibble)
library(dplyr)

# Anzahl der Worker automatisch aus SLURM oder manuell setzen
library(future)
workers <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", 8))
plan(multisession, workers = workers)

# Eingabe-FASTA
fasta_path <- "/media/BioNAS2/Tensor_Omics/Protein_Function_Prediction/material/uniref50_morethan1.fasta"

message("Lade FASTA-Datei: ", fasta_path)
aa <- readAAStringSet(fasta_path)
ids <- names(aa)
seqs <- as.character(aa)

message("Anzahl geladener Sequenzen: ", length(seqs))
message("Beispiel-Sequenz (roh): ", substr(seqs[1], 1, 60), "...")

# IDs bereinigen
ids_clean <- sub(" .*", "", ids)

# Nicht-Standard-AAs entfernen
clean_sequence <- function(seq) {
  gsub("[^ACDEFGHIKLMNPQRSTVWY]", "", seq)
}

# Safe-Kidera mit Reinigung
safe_kidera <- function(seq) {
  seq <- clean_sequence(seq)
  if (nchar(seq) == 0) {
    return(rep(NA_real_, 10))
  }
  tryCatch(
    kideraFactors(seq),
    error = function(e) {
      message("Fehler bei Sequenz (nach Cleanup): ", seq)
      rep(NA_real_, 10)
    }
  )
}

# Berechne Kidera-Faktoren
message("Starte Kidera-Berechnung mit ", workers, " Workern...")
kidera_list <- future_lapply(seqs, safe_kidera)

message("Berechnung abgeschlossen. Erstelle Matrix...")

kidera_matrix <- do.call(rbind, kidera_list)
colnames(kidera_matrix) <- paste0("k", 1:10)

# In Tibble mit IDs
kidera_df <- as_tibble(kidera_matrix) %>%
  mutate(id = ids_clean, .before = 1)

message("Anzahl Zeilen vor NA-Filter: ", nrow(kidera_df))

# Zeige Beispielzeilen
message("Beispiel-Kidera-Werte:")
print(head(kidera_df, 3))

# Entferne Zeilen mit NA
kidera_df <- kidera_df %>%
  filter(if_all(starts_with("k"), ~ !is.na(.x)))

message("Anzahl Zeilen nach NA-Filter: ", nrow(kidera_df))

# Ausgabe speichern
save_path <- "/media/BioNAS2/Tensor_Omics/Protein_Function_Prediction/material/kidera_factors_uniref50_morethan1_ref_prots.rds"
saveRDS(kidera_df, save_path)
message("Speichern abgeschlossen: ", save_path)