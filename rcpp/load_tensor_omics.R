# Load the generated TensorOmics R interface into the calling environment.
#
# Hand-written, and it lives here rather than inside `rcpp/tensor_omics/` because that
# whole directory is generated and gitignored. It is the R counterpart of the generated
# `python/tensor_omics/library.py`.
#
# Usage from the repository root:
#     source("rcpp/load_tensor_omics.R")
#
# Compilation is cached under `rcpp/rcpp_cache`, so only changed wrappers are rebuilt.

library(Rcpp)

.tox_load <- function(root = "rcpp/tensor_omics", cache = "rcpp/rcpp_cache") {
    build <- shQuote(normalizePath("build", mustWork = TRUE))

    # every wrapper links against the one Fortran shared library
    Sys.setenv(PKG_LIBS = paste0(
        "-Wl,-rpath,", build, " -L", build, " -ltensor-omics -lgfortran"
    ))

    dir.create(cache, showWarnings = FALSE, recursive = TRUE)

    sources <- sort(list.files(file.path(root, "src"), pattern = "\\.cpp$", full.names = TRUE))
    for (source_file in sources) {
        sourceCpp(source_file, env = .GlobalEnv, cacheDir = cache)
    }

    # error_handling and tox_validate first: the wrappers call into them
    scripts <- file.path(root, "R", c("error_handling.R", "tox_validate.R"))
    rest <- sort(list.files(file.path(root, "R"), pattern = "\\.R$", full.names = TRUE))
    scripts <- c(scripts, setdiff(rest, scripts))

    for (script in scripts) {
        source(script, local = FALSE)
    }

    invisible(length(sources))
}

.tox_load()
