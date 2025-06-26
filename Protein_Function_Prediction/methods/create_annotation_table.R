library(data.table)
library(dplyr)
library(readr)
library(jsonlite)
library(tibble)

cat("Starting script...\n")

# === 1. Read UniRef50 IDs ===
cat("Reading UniRef50 IDs...\n")
uniref50_headers <- readLines("./material/uniref50.fasta")
uniref50_headers <- uniref50_headers[grepl("^>", uniref50_headers)]
uniref50_ids <- sub("^>UniRef50_(\\w+)\\s.*", "\\1", uniref50_headers)
uniref50_set <- new.env(hash = TRUE, parent = emptyenv())
for (id in uniref50_ids) uniref50_set[[id]] <- TRUE
rm(uniref50_headers, uniref50_ids)
gc()

# === 2. protein2ipr: protein_id + description (nur UniRef50) ===
cat("Reading protein2ipr.dat in chunks...\n")
con <- file("./material/protein2ipr.dat", open = "r")
protein_ids <- character()
descriptions <- character()
while (length(lines <- readLines(con, n = 100000, warn = FALSE)) > 0) {
  splits <- strsplit(trimws(lines), " +")
  for (x in splits) {
    if (length(x) >= 5) {
      protein_ids <- c(protein_ids, x[[1]])
      descriptions <- c(descriptions, paste(x[3:(length(x)-3)], collapse = " "))
    } else if (length(x) >= 3) {
      protein_ids <- c(protein_ids, x[[1]])
      descriptions <- c(descriptions, x[3])
    }
  }
}
close(con)

protein2ipr <- tibble(protein_id = protein_ids, description = descriptions) %>%
  filter(!is.na(description)) %>%
  distinct() %>%
  filter(sapply(protein_id, function(x) !is.null(uniref50_set[[x]])))
saveRDS(protein2ipr, "./material/protein2ipr.rds")
write.csv(protein2ipr, "./material/protein2ipr.csv")
rm(protein_ids, descriptions, splits, lines)
gc()

# === 3. protein2goa: protein_id + description + GOA_label (nur UniRef50) ===
cat("Reading go-basic-filtered.json...\n")
goa_basic <- fromJSON(
  "./material/go-basic-filtered.json",
  simplifyVector = TRUE
) %>%
  as_tibble() %>%
  select(GO_term = id, GOA_label = lbl) %>%
  mutate(GO_term = sub("^GO:", "GO_", GO_term))
gc()

cat("Reading goa_uniprot_all.gaf in chunks...\n")
gaf_path <- "./material/goa_uniprot_all.gaf"
gaf_con <- file(gaf_path, open = "r")
gaf_proteins <- character()
gaf_terms <- character()
while (length(lines <- readLines(gaf_con, n = 100000, warn = FALSE)) > 0) {
  for (line in lines) {
    cols <- strsplit(line, " +")[[1]]
    if (length(cols) >= 5 && grepl("^GO:", cols[5])) {
      gaf_proteins <- c(gaf_proteins, cols[2])
      gaf_terms <- c(gaf_terms, cols[5])
    }
  }
}
close(gaf_con)
goa_uniprot <- tibble(
  protein_id = gaf_proteins,
  GO_term = sub("^GO:", "GO_", gaf_terms)
)
rm(gaf_proteins, gaf_terms, lines)
gc()
protein2goa <- goa_uniprot %>%
  left_join(goa_basic, by = "GO_term") %>%
  select(protein_id, GOA_label) %>%
  filter(!is.na(GOA_label)) %>%
  distinct()
# Merge mit protein2ipr, damit description enthalten ist
protein2goa <- protein2goa %>%
  left_join(protein2ipr, by = "protein_id") %>%
  select(protein_id, description, GOA_label) %>%
  filter(!is.na(description))
saveRDS(protein2goa, "./material/protein2goa.rds")
# write.csv(protein2goa, "./material/protein2goa.csv")
rm(goa_uniprot, goa_basic)
gc()

# === 4. protein2hrd: protein_id + hrd (nur UniRef50) ===
cat("Reading swissprot headers...\n")
sprot_headers <- readLines("./material/uniprot_sprot.fasta")
sprot_headers <- sprot_headers[grepl("^>", sprot_headers)]
sprot_df <- tibble(
  header = sprot_headers,
  protein_id = sub("^>sp\\|(\\w+)\\|.*", "\\1", sprot_headers),
  hrd = sub("^>sp\\|\\w+\\|[^ ]+ (.+) OS=.*$", "\\1", sprot_headers)
) %>% select(protein_id, hrd)
rm(sprot_headers)
gc()

cat("Reading trembl headers...\n")
trembl_headers <- readLines("./material/uniprot_trembl.fasta")
trembl_headers <- trembl_headers[grepl("^>", trembl_headers)]
trembl_df <- tibble(
  header = trembl_headers,
  protein_id = sub("^>tr\\|(\\w+)\\|.*", "\\1", trembl_headers),
  hrd = sub("^>tr\\|\\w+\\|[^ ]+ (.+) OS=.*$", "\\1", trembl_headers)
) %>% select(protein_id, hrd)
protein2hrd <- bind_rows(sprot_df, trembl_df) %>%
  filter(!is.na(hrd)) %>%
  distinct() %>%
  filter(sapply(protein_id, function(x) !is.null(uniref50_set[[x]])))
saveRDS(protein2hrd, "./material/protein2hrd.rds")
#write.csv(protein2hrd, "./material/protein2hrd.csv")
rm(trembl_headers, sprot_df, trembl_df)
gc()

cat("All annotation tables created. Done.\n")