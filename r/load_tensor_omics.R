# Load the generated TensorOmics R interface into the calling environment.
#
# Hand-written, and it lives here rather than inside `r/tensor_omics/` because that
# whole directory is generated. It is the R counterpart of the generated
# `python/tensor_omics/library.py`.
#
# Usage from the repository root:
#     source("r/load_tensor_omics.R")
#
# The generated shims are pure C (`.Call`), compiled into one shared object with
# `R CMD SHLIB` and cached under `r/r_cache`; only rebuilt when a source changes.

.tox_load <- function(root = "r/tensor_omics", cache = "r/r_cache") {
    build <- normalizePath("build", mustWork = TRUE)
    src <- normalizePath(file.path(root, "src"), mustWork = TRUE)
    dir.create(cache, showWarnings = FALSE, recursive = TRUE)
    cache <- normalizePath(cache, mustWork = TRUE)
    so <- file.path(cache, paste0("tensoromics", .Platform$dynlib.ext))

    sources <- sort(list.files(src, pattern = "\\.c$", full.names = TRUE))
    headers <- list.files(src, pattern = "\\.h$", full.names = TRUE)
    # rebuild only when a shim (or the shared header) is newer than the compiled object
    newest <- max(file.info(c(sources, headers))$mtime)
    if (!file.exists(so) || file.info(so)$mtime < newest) {
        # compile in the cache, so the generated src/ stays free of .o litter
        file.copy(c(sources, headers), cache, overwrite = TRUE)
        # one shared object, linked against the Fortran library; the library carries its own
        # runtime as NEEDED, so the rpath to build/ is all the shims need
        Sys.setenv(PKG_LIBS = paste0("-L", shQuote(build), " -Wl,-rpath,", shQuote(build),
                                     " -ltensor-omics"))
        old <- setwd(cache)
        on.exit(setwd(old), add = TRUE)
        status <- system2(
            file.path(R.home("bin"), "R"),
            c("CMD", "SHLIB", shQuote(basename(sources)), "-o", "tensoromics.so"),
            stdout = TRUE, stderr = TRUE
        )
        setwd(old)
        if (!is.null(attr(status, "status")) && attr(status, "status") != 0)
            stop("failed to compile the R interface:\n", paste(status, collapse = "\n"))
    }

    if (is.null(getLoadedDLLs()[["tensoromics"]])) dyn.load(so)

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
