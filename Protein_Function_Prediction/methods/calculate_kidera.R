#Version 1.2.7
#Author: Aaron Schroeder

library(Biostrings)     # For efficient reading of fasta
library(Peptides)       # for kideraFactors()
library(future.apply)   # for parallel
library(tibble)         # for dataframes
library(dplyr)          # For efficient data modification

# Parallel-Plan with 8 Workers
plan(multisession, workers = 8)

# Path to fasta
fasta_path <- "/media/BioNAS2/Tensor_Omics/Protein_Function_Prediction/material/uniref50_morethan1.fasta"

cat("Loading fasta:", fasta_path, "\n")
aa <- readAAStringSet(fasta_path)
ids <- names(aa)
seqs <- as.character(aa)
cat("Numbers of loaded sequences:", length(seqs), "\n")

# Clean IDs
ids_clean <- sub("^UniRef50_", "", sub(" .*", "", ids))

# calc kideras
safe_kidera <- function(seq) {
  if (nchar(seq) == 0) {
	# if sequence empty return NA
    return(rep(NA_real_, 10))
  }
  tryCatch(
    as.numeric(kideraFactors(seq)[[1]]),
    error = function(e) rep(NA_real_, 10)
  )
}

cat("Starting kidera calculation with 8 workers...\n")
kidera_list <- future_lapply(seqs, safe_kidera)
cat("Done. Creating Matrix...\n")

kidera_matrix <- do.call(rbind, kidera_list)
colnames(kidera_matrix) <- paste0("KF", 1:10)

# Build result data frame
kidera_df <- as_tibble(kidera_matrix) %>%
  mutate(id = ids_clean, .before = 1) %>%
  filter(if_all(starts_with("KF"), ~ !is.na(.x)))

cat("Final dataframe with: ", nrow(kidera_df), "valid sequences.\n")

# save
saveRDS(kidera_df, "/media/BioNAS2/Tensor_Omics/Protein_Function_Prediction/material/kidera_factors_uniref50_morethan1_ref_prots.rds")
cat("Data saved.\n")

