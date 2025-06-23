library(data.table)
library(jsonlite)

# 1. Efficiently read protein2ipr.dat line-by-line and extract protein_id and description
cat("Reading protein2ipr.dat\n")
protein2ipr_path <- "/media/BioNAS2/Tensor_Omics/Protein_Function_Prediction/material/protein2ipr.dat"
con <- file(protein2ipr_path, open = "r")
batch_size <- 1e5
protein_id_list <- list()
description_list <- list()

repeat {
  lines <- readLines(con, n = batch_size)
  if (length(lines) == 0) break
  split_lines <- strsplit(trimws(lines), " +")
  protein_id_list[[length(protein_id_list)+1]] <- vapply(split_lines, function(x) x[1], character(1))
  description_list[[length(description_list)+1]] <- vapply(split_lines, function(x) {
    if (length(x) >= 5) paste(x[3:(length(x)-3)], collapse = " ")
    else if (length(x) >= 3) x[3]
    else NA_character_
  }, character(1))
}
close(con)

protein2ipr <- data.table(
  protein_id = unlist(protein_id_list),
  description = unlist(description_list)
)[!is.na(description)]
rm(protein_id_list, description_list, split_lines, lines)
gc()

# 2. Read GO terms JSON (small file, no optimization needed)
cat("Reading go-basic-filtered.json\n")
go_basic_path <- "/media/BioNAS2/Tensor_Omics/Protein_Function_Prediction/material/go-basic-filtered.json"
go_basic <- fromJSON(go_basic_path, simplifyVector = TRUE)
go_basic <- data.table(id = go_basic$id, GOA_label = go_basic$lbl)
go_basic[, GO_term := sub("^GO:", "GO_", id)][, id := NULL]

# 3. Read GAF file line-by-line and convert to data.table
cat("Reading goa_uniprot_all.gaf\n")
gaf_path <- "/media/BioNAS2/Tensor_Omics/Protein_Function_Prediction/material/goa_uniprot_all.gaf"
con <- file(gaf_path, open = "r")
gaf_lines <- list()
repeat {
  lines <- readLines(con, n = batch_size)
  if (length(lines) == 0) break
  lines <- lines[sapply(strsplit(lines, "\\t"), length) >= 5]  # Keep valid lines only
  gaf_lines[[length(gaf_lines) + 1]] <- lines
}
close(con)
all_lines <- unlist(gaf_lines)
tmp_gaf <- tempfile(tmpdir = "/media/BioNAS2/Tensor_Omics/Protein_Function_Prediction/tmp", fileext = ".gaf")
writeLines(all_lines, tmp_gaf)
gaf_dt <- fread(tmp_gaf, header = FALSE, fill = TRUE, sep = "\t")
unlink(tmp_gaf)
rm(gaf_lines, all_lines)
gc()

# Extract protein_id and GO_term (first column starting with GO:)
GO_idx <- which(apply(gaf_dt, 2, function(col) any(grepl("^GO:", col))))[1]
goa_uniprot <- data.table(
  protein_id = gaf_dt[[2]],
  GO_term = sub("^GO:", "GO_", gaf_dt[[GO_idx]])
)
rm(gaf_dt)
gc()

# 4. Parse uniprot_sprot.fasta headers
cat("Reading uniprot_sprot.fasta\n")
sprot_headers <- fread(
  cmd = "grep '^>' /media/BioNAS2/Tensor_Omics/Protein_Function_Prediction/material/uniprot_sprot.fasta",
  sep = "\n", header = FALSE
)[[1]]
sprot_df <- data.table(
  protein_id = sub("^>sp\\|(\\w+)\\|.*", "\\1", sprot_headers),
  hrd = sub("^>sp\\|\\w+\\|[^ ]+ (.+) OS=.*$", "\\1", sprot_headers)
)
# 5. Parse uniprot_trembl.fasta headers
cat("Reading uniprot_trembl.fasta\n")
trembl_headers <- fread(
  cmd = "grep '^>' /media/BioNAS2/Tensor_Omics/Protein_Function_Prediction/material/uniprot_trembl.fasta",
  sep = "\n", header = FALSE
)[[1]]
trembl_df <- data.table(
  protein_id = sub("^>tr\\|(\\w+)\\|.*", "\\1", trembl_headers),
  hrd = sub("^>tr\\|\\w+\\|[^ ]+ (.+) OS=.*$", "\\1", trembl_headers)
)

# Combine headers
hrds <- rbindlist(list(sprot_df, trembl_df))
rm(sprot_headers, trembl_headers, sprot_df, trembl_df)
gc()

# 6. Join all annotation info
cat("Creating annotation table\n")
setkey(protein2ipr, protein_id)
setkey(hrds, protein_id)
setkey(goa_uniprot, protein_id)
setkey(go_basic, GO_term)

annotation_table <- merge(protein2ipr, hrds, all.x = TRUE)
annotation_table <- merge(annotation_table, goa_uniprot, all.x = TRUE)
annotation_table <- merge(annotation_table, go_basic, by = "GO_term", all.x = TRUE)
annotation_table <- unique(annotation_table[!is.na(GOA_label) | !is.na(hrd)])

# 7. Save the unfiltered annotation table
saveRDS(annotation_table, "/media/BioNAS2/Tensor_Omics/Protein_Function_Prediction/material/annotation_table_unfiltered.rds")

# 8. Filter for UniRef50 IDs
cat("Reading UniRef50 and filtering\n")
uniref50_headers <- fread(
  cmd = "grep '^>' /media/BioNAS2/Tensor_Omics/Protein_Function_Prediction/material/uniref50.fasta",
  sep = "\n", header = FALSE
)[[1]]
uniref50_ids <- sub("^>UniRef50_(\\w+)\\s.*", "\\1", uniref50_headers)
uniref50_dt <- data.table(protein_id = uniref50_ids)
setkey(uniref50_dt, protein_id)

# Perform semi-join (retain only matching protein_ids)
annotation_table <- annotation_table[uniref50_dt, nomatch = 0]

# 9. Save filtered table
saveRDS(annotation_table, "/media/BioNAS2/Tensor_Omics/Protein_Function_Prediction/material/annotation_table_uniref50.rds")

cat("Done!\n")
