# create_chemical_space.R
# Date: 2025-08-18
# Version 1.3.2
# Author: Aaron Schroeder
#
# Overview
# - This script combines SVD embeddings and Kidera factors into a single,
#   pre-built chemical space RDS file for use in run_query.R.
#
# Notes
# - The final RDS object is a list containing:
#   - mat: the combined matrix of SVD and Kidera factors
#   - svd: a list containing the original SVD data table, matrix, dimensions, and names
#   - kidera: a list containing the original Kidera data table and matrix
# - This separation ensures that run_query.R can correctly identify and use
#   the SVD-only dimensions for query embedding.

suppressPackageStartupMessages({
  require(data.table)
  require(Matrix)
  require(stringr)
  require(tools)
  require(optparse)
})

# ------------------------
# Utility helpers
# ------------------------
canon_id <- function(x) {
  x <- as.character(x)
  x <- str_trim(x)
  x <- sub("^(?i)UniRef50[_:]?", "", x, perl = TRUE)
  x
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

# ------------------------
# Loading and building data
# ------------------------
load_svd_csv <- function(path, id_col = 1L) {
  message("Loading SVD data from ", path, "...")
  dt <- fread(path)
  if (is.numeric(id_col)) id_col <- names(dt)[id_col]
  stopifnot(id_col %in% names(dt))
  setnames(dt, id_col, "id_raw")
  dt[, id := canon_id(id_raw)]
  dt[, id_raw := NULL]
  num_cols <- setdiff(names(dt), "id")
  for (nc in num_cols) dt[[nc]] <- as.numeric(dt[[nc]])
  dt <- dt[!is.na(id) & id != ""]
  svd_mat <- as.matrix(dt[, ..num_cols])
  rownames(svd_mat) <- dt$id
  message("SVD data read: ", nrow(svd_mat), " Entries with ", ncol(svd_mat), " dimensions.")
  list(dt = dt, mat = svd_mat, dims = length(num_cols), names = num_cols)
}

load_kidera_rds <- function(path, id_col = "id") {
  message("Loading kidera factors from ", path, "...")
  kd <- readRDS(path)
  kd <- as.data.table(kd)
  stopifnot(id_col %in% names(kd))
  setnames(kd, id_col, "id")
  kd[, id := canon_id(id)]
  k_cols <- grep("^KF\\d+$", names(kd), value = TRUE)
  stopifnot(length(k_cols) == 10)
  setcolorder(kd, c("id", k_cols))
  message("Kidera factors loaded: ", nrow(kd), " entries.")
  list(dt = kd, mat = as.matrix(kd[, ..k_cols]), names = k_cols)
}

build_chemical_space <- function(svd, kidera) {
  message("Building chemical reference space...")
  svd_dt <- as.data.table(svd$dt)
  kd_dt  <- as.data.table(kidera$dt)
  merged <- merge(svd_dt, kd_dt, by = "id", all.x = TRUE, all.y = FALSE, sort = FALSE)
  missing_kd <- merged[is.na(get(kidera$names[1]))]
  if (nrow(missing_kd) > 0) {
    warning(sprintf("%d SVD entries without Kidera factors. They will be dropped.", nrow(missing_kd)))
    merged <- merged[!is.na(get(kidera$names[1]))]
  }
  chem_cols <- c(svd$names, kidera$names)
  chem_mat <- as.matrix(merged[, ..chem_cols])
  rownames(chem_mat) <- merged$id
  message("Chemical room build: ", nrow(chem_mat), " proteins with ", ncol(chem_mat), " dimensions.")
  list(dt = merged, mat = chem_mat, dims = ncol(chem_mat), svd = svd, kidera = kidera)
}

# ------------------------
# Main entry via command line
# ------------------------
option_list <- list(
  make_option(c("--svd"), type = "character", help = "CSV with SVD embeddings"),
  make_option(c("--kidera"), type = "character", help = "RDS with Kidera factors"),
  make_option(c("--out"), type = "character", default = "results/chemical_space.rds", help = "Output RDS file path")
)

args <- parse_args(OptionParser(option_list = option_list))

if (!is.null(args$svd) && !is.null(args$kidera)) {
  svd <- load_svd_csv(args$svd)
  kidera <- load_kidera_rds(args$kidera)
  chem_space <- build_chemical_space(svd, kidera)
  
  ensure_dir(dirname(args$out))
  message("Saving chemical space to ", args$out, "...")
  saveRDS(chem_space, args$out)
  message("Done.")
} else {
  stop("Both --svd and --kidera must be provided.")
}
