library(data.table)
library(doParallel)
library(dplyr)
library(jsonlite)
library(tibble)

# Set up parallel processing
registerDoParallel(cores = min(4, detectCores() - 1))  # Use 4 cores or available-1
cat("Starting parallel processing with", getDoParWorkers(), "cores...\n")

# === 1. Read UniRef50 IDs (required by all processes) ===
cat("Reading UniRef50 IDs...\n")
uniref50_headers <- readLines("./material/uniref50.fasta")
uniref50_ids <- sub("^>UniRef50_(\\w+)\\s.*", "\\1", grep("^>", uniref50_headers, value = TRUE))
uniref50_set <- new.env(hash = TRUE, parent = emptyenv())
for (id in uniref50_ids) uniref50_set[[id]] <- TRUE
rm(uniref50_headers, uniref50_ids)
gc()

# === 2. Parallel Processing of Large Files ===
results <- foreach(file_type = c("protein2ipr", "goa_uniprot"), .combine = list) %dopar% {
  
  if (file_type == "protein2ipr") {
    # === Process protein2ipr.dat (95GB) ===
    cat("Parallel processing of protein2ipr.dat...\n")
    chunk_size <- 1e7  # 10M lines per chunk
    file_length <- as.integer(system("wc -l < ./material/protein2ipr.dat", intern = TRUE))
    chunks <- seq(1, file_length, by = chunk_size)
    
    protein2ipr <- foreach(
      start = chunks,
      .combine = rbind,
      .packages = c("data.table", "dplyr"),
      .final = function(x) distinct(as_tibble(x))
    ) %do% {
      dt <- fread(
        "./material/protein2ipr.dat",
        skip = start - 1,
        nrows = chunk_size,
        header = FALSE,
        select = c(1, 3),
        col.names = c("protein_id", "description")
      )
      dt[protein_id %chin% names(as.list(uniref50_set))]
    }
    
    saveRDS(protein2ipr, "./material/protein2ipr.rds")
    write.csv(protein2ipr, "./material/protein2ipr.csv")
    return(protein2ipr)
    
  } else if (file_type == "goa_uniprot") {
    # === Process goa_uniprot_all.gaf (250GB) ===
    cat("Parallel processing of goa_uniprot_all.gaf...\n")
    
    # First get GO term info (small file)
    goa_basic <- fromJSON("./material/go-basic-filtered.json") %>%
      as_tibble() %>%
      select(GO_term = id, GOA_label = lbl) %>%
      mutate(GO_term = sub("^GO:", "GO_", GO_term))
    
    # Process GAF in chunks
    gaf_chunk_size <- 5e6  # 5M lines per chunk
    gaf_length <- as.integer(system("wc -l < ./material/goa_uniprot_all.gaf", intern = TRUE))
    gaf_chunks <- seq(1, gaf_length, by = gaf_chunk_size)
    
    goa_uniprot <- foreach(
      start = gaf_chunks,
      .combine = rbind,
      .packages = c("data.table", "dplyr"),
      .final = function(x) distinct(as_tibble(x))
    ) %do% {
      dt <- fread(
        "./material/goa_uniprot_all.gaf",
        skip = start - 1,
        nrows = gaf_chunk_size,
        header = FALSE,
        select = c(2, 5),  # protein_id and GO_term
        col.names = c("protein_id", "GO_term")
      )
      dt <- dt[protein_id %chin% names(as.list(uniref50_set)) & grepl("^GO:", GO_term)]
      dt[, GO_term := sub("^GO:", "GO_", GO_term)]
      dt
    }
    
    protein2goa <- goa_uniprot %>%
      left_join(goa_basic, by = "GO_term") %>%
      select(protein_id, GOA_label) %>%
      filter(!is.na(GOA_label)) %>%
      distinct()
    
    # Merge with descriptions (from already processed protein2ipr if available)
    if (file.exists("./material/protein2ipr.rds")) {
      protein2ipr <- readRDS("./material/protein2ipr.rds")
      protein2goa <- protein2goa %>%
        left_join(select(protein2ipr, protein_id, description), by = "protein_id") %>%
        filter(!is.na(description))
    }
    
    saveRDS(protein2goa, "./material/protein2goa.rds")
    return(protein2goa)
  }
}

# === 3. Process SwissProt/Trembl (sequential as it's faster) ===
cat("Processing SwissProt/Trembl files...\n")
process_uniprot_fasta <- function(file_path, prefix) {
  headers <- readLines(file_path)
  headers <- headers[grepl("^>", headers)]
  tibble(
    protein_id = sub(paste0("^", prefix, "\\|(\\w+)\\|.*"), "\\1", headers),
    hrd = sub(paste0("^", prefix, "\\|\\w+\\|[^ ]+ (.+) OS=.*$"), "\\1", headers)
  ) %>% select(protein_id, hrd)
}

sprot_df <- process_uniprot_fasta("./material/uniprot_sprot.fasta", "sp")
trembl_df <- process_uniprot_fasta("./material/uniprot_trembl.fasta", "tr")

sprot_dt <- as.data.table(sprot_df)
trembl_dt <- as.data.table(trembl_df)

# Filter NA
sprot_dt <- sprot_dt[!is.na(hrd)]
trembl_dt <- trembl_dt[!is.na(hrd)]

# Combine data frames
combined_dt <- rbindlist(list(sprot_dt, trembl_dt))

# remove duplicates (if included)
combined_dt <- unique(combined_dt, by = c("protein_id", "hrd")) 

# merge with efficient hash lookup
uniref50_ids <- names(as.list(uniref50_set))
protein2hrd <- combined_dt[protein_id %chin% uniref50_ids]

saveRDS(protein2hrd, "./material/protein2hrd.rds")

cat("All annotation tables created successfully!\n")
stopImplicitCluster()  # Clean up parallel workers