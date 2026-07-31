# Load the generated TensorOmics R interface into the calling environment.
#
# Hand-written, and it lives here rather than inside `r/tensor_omics/` because that
# whole directory is generated. It is the R counterpart of the generated
# `python/tensor_omics/library.py`.
#
# Usage from the repository root:
#     source("r/load_tensor_omics.R")
#
# The C `.Call` shims are compiled into `build/libtensor-omics.so` by the ordinary build
# (`./build.sh`), so there is nothing to compile here -- this loads that library and sources
# the generated R wrappers. Build the library first, with the R interface included (the
# default; i.e. not `--directive=NO_R_INTERFACE`).

.tox_load <- function(root = "r/tensor_omics", lib = "build/libtensor-omics.so") {
    so <- normalizePath(lib, mustWork = TRUE)

    # A rebuilt .so is only picked up by a fresh R session; within one session dyn.load is a
    # no-op once loaded (dyn.unload first to reload). Fresh Rscript runs -- the test flow --
    # are unaffected.
    if (is.null(getLoadedDLLs()[["libtensor-omics"]])) dyn.load(so)

    # a wrapper reaches C through .Call by name (resolved by dynamic lookup); if the R shims
    # were left out of the build (--directive=NO_R_INTERFACE, or R absent when it ran) the
    # symbol is not there
    if (!is.loaded("loess_smooth_2d_call"))
        stop("'", so, "' has no R interface -- rebuild it with the R interface included ",
             "(drop --directive=NO_R_INTERFACE, and make sure R is installed).")

    # error_handling and tox_validate first: the wrappers call into them
    scripts <- file.path(root, c("error_handling.R", "tox_validate.R"))
    rest <- sort(list.files(root, pattern = "\\.R$", full.names = TRUE))
    scripts <- c(scripts, setdiff(rest, scripts))

    for (script in scripts) {
        source(script, local = FALSE)
    }

    invisible(so)
}

.tox_load()
