# Create an annotation table
# Goal Structure:

library(data.table)
library(dplyr)
library(readr)
library(jsonlite)
library(tibble)

# 1. Read protein2ipr.dat, extract protein_id and description robustly
protein2ipr_raw <- fread(
  "./material/protein2ipr.dat",
  header = FALSE,
  sep = " ",
  fill = TRUE,
  strip.white = TRUE
)

# 1. Read protein2ipr.dat, extract protein_id and description robustly (line by line)
protein2ipr_lines <- readLines("./material/protein2ipr.dat")
protein2ipr_split <- strsplit(trimws(protein2ipr_lines), " +")  # trim whitespace and split

protein2ipr <- tibble(
  protein_id = sapply(protein2ipr_split, function(x) x[1]),
  description = sapply(protein2ipr_split, function(x) {
    # Description is from 3rd element to 3rd last element
    if (length(x) >= 5) {  # Need at least 5 elements to have middle fields (1,2,...n-3,n-2,n-1,n)
      paste(x[3:(length(x)-3)], collapse = " ")
    } else if (length(x) >= 3) {  # If only 3-4 elements, just take 3rd
      x[3]
    } else {
      NA_character_  # Handle unexpected cases
    }
  })
) %>% 
  filter(!is.na(description))  # Remove any rows where description couldn't be extracted

# 3. Read go-basic-filtered.json (GO Term ID and label)
goa_basic <- fromJSON(
  "./material/go-basic-filtered.json",
  simplifyVector = TRUE
) %>% as_tibble() %>%
  select(GO_term = id, GOA_label = lbl) %>%
  mutate(GO_term = sub("^GO:", "GO_", GO_term)) # ensure same format as goa_uniprot

# 4. Read goa_uniprot_all.gaf - select protein_id and the column that starts with "GO:"
gaf_lines <- readLines("./material/goa_uniprot_all.gaf")
gaf_lines <- gaf_lines[sapply(strsplit(gaf_lines, " +"), length) >= 5] # changed "\t" to " +"
tmp_gaf <- tempfile(fileext = ".gaf")
writeLines(gaf_lines, tmp_gaf)
gaf_dt <- fread(
  tmp_gaf,
  header = FALSE,
  fill = TRUE,
  strip.white = TRUE
)
unlink(tmp_gaf)

# Only build goa_uniprot if gaf_dt has at least 2 columns and 1 row
if (nrow(gaf_dt) > 0 && ncol(gaf_dt) >= 2) {
  goa_uniprot <- tibble(
    protein_id = gaf_dt[[2]],
    GO_term = apply(gaf_dt, 1, function(row) {
      go_col <- which(grepl("^GO:", row))
      if (length(go_col) > 0) row[go_col[1]] else NA_character_
    })
  )
} else {
  goa_uniprot <- tibble(protein_id = character(), GO_term = character())
}

if (nrow(goa_uniprot) == 0) {
  goa_uniprot <- tibble(protein_id = character(), GO_term = character())
}

# GO-Terms vereinheitlichen: GO: zu GO_
goa_uniprot <- goa_uniprot %>%
  mutate(GO_term = sub("^GO:", "GO_", GO_term))

# 5. Parse uniprot_sprot.fasta headers for ID and human readable description (HRD)
# Read only header lines to save RAM
sprot_headers <- readLines("./material/uniprot_sprot.fasta")
sprot_headers <- sprot_headers[grepl("^>", sprot_headers)]
sprot_df <- tibble(
  header = sprot_headers
) %>%
  mutate(
    protein_id = sub("^>sp\\|(\\w+)\\|.*", "\\1", header),
    # Extract HRD - text after second '|' and before OS=
    hrd = sub("^>sp\\|\\w+\\|[^ ]+ (.+) OS=.*$", "\\1", header)
  ) %>%
  select(protein_id, hrd)

print(sprot_df)
# 6. Parse uniprot_trembl.fasta headers for ID and HRD
trembl_headers <- readLines("./material/uniprot_trembl.fasta")
trembl_headers <- trembl_headers[grepl("^>", trembl_headers)]
trembl_df <- tibble(
  header = trembl_headers
) %>%
  mutate(
    # Extract Uniprot ID e.g. >tr|A0A7C4TNZ0|...
    protein_id = sub("^>tr\\|(\\w+)\\|.*", "\\1", header),
    hrd = sub("^>tr\\|\\w+\\|[^ ]+ (.+) OS=.*$", "\\1", header)
  ) %>%
  select(protein_id, hrd)

hrds <- rbind(sprot_df, trembl_df)
print(trembl_df)
# 7. Combine all data into a single annotation table
annotation_table <- protein2ipr %>%
    left_join(hrds, by = "protein_id") %>% # join HRD first
    left_join(goa_uniprot, by = "protein_id") %>%
    left_join(goa_basic, by = c("GO_term" = "GO_term")) %>%
    select(protein_id, description, GOA_label, hrd) %>% # GO_term entfernt
    distinct() %>%
    filter(!is.na(GOA_label) | !is.na(hrd)) # GOA_label statt GO_term

# 8. Save the annotation table as RDS file
saveRDS(annotation_table, "./material/annotation_table.rds")

# Tabelle als CSV speichern für VSCode-Ansicht
write.csv(annotation_table, "./material/annotation_table.csv", row.names = FALSE)

# Ausgabe des vollständigen Data Frames
# print(annotation_table, n = Inf)
# Alternativ in RStudio:
# View(annotation_table)