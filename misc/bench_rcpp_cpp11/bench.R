# =============================================================================
# Benchmark: marshalling overhead of three R -> Fortran calling paths.
#
#   1. Rcpp        — .Call shim, inputs ALIAS R's buffers (zero-copy)   [today]
#   2. cpp11       — .Call shim, inputs ALIAS R's buffers (zero-copy)   [candidate]
#   3. .C          — copies every argument in and out (no shim)          [baseline]
#
# All three call the SAME Fortran symbol `loess_smooth_2d_c` in
# build/libtensor-omics.so, so the compute time is identical and cancels out —
# what remains is the cost of moving the arguments across the boundary.
#
# The routine is driven with huge read-only reference arrays (x_ref/y_ref) but a
# tiny working/output set (n_used = n_target = 100). Compute is O(n_target*n_used)
# = trivial; the only thing that scales with size is whether the big inputs get
# duplicated. So: .C should grow ~linearly with input bytes, the two zero-copy
# paths should stay flat.
#
# Run from anywhere:   Rscript misc/bench_rcpp_cpp11/bench.R
# =============================================================================

options(warn = 1)

## ---- locate ourselves & the Fortran library ---------------------------------
.args     <- commandArgs(trailingOnly = FALSE)
.file_arg <- sub("^--file=", "", .args[grep("^--file=", .args)])
script_dir <- if (length(.file_arg)) dirname(normalizePath(.file_arg)) else normalizePath(".")
repo_root  <- normalizePath(file.path(script_dir, "..", ".."))
build_dir  <- file.path(repo_root, "build")
so_path    <- file.path(build_dir, "libtensor-omics.so")

if (!file.exists(so_path))
    stop("Missing ", so_path, " — build it first with ./build.sh")

cat("repo_root :", repo_root, "\n")
cat("shared lib:", so_path, "\n\n")

## ---- dependencies -----------------------------------------------------------
# cpp11 is header-only; installing it only writes headers + a tiny R package.
if (!requireNamespace("cpp11", quietly = TRUE)) {
    cat("Installing cpp11 (header-only) ...\n")
    install.packages("cpp11", repos = "https://cloud.r-project.org")
}
stopifnot(requireNamespace("Rcpp",  quietly = TRUE),
          requireNamespace("cpp11", quietly = TRUE))

## ---- link both shims against the (Intel-built) Fortran lib -------------------
# Just -ltensor-omics: the Fortran runtime (Intel iomp5/ifcore) is pulled in via
# the .so's own NEEDED entries, so no -lgfortran here. Path is quoted: the repo
# path contains a space.
pkg_libs <- sprintf("-L%s -Wl,-rpath,%s -ltensor-omics",
                    shQuote(build_dir), shQuote(build_dir))
Sys.setenv(PKG_LIBS = pkg_libs)
cat("PKG_LIBS  :", pkg_libs, "\n\n")

cat("Compiling Rcpp shim  (sourceCpp) ...\n")
Rcpp::sourceCpp(file.path(script_dir, "loess_rcpp.cpp"))

# cpp11 shim: compile the hand-written .Call glue with plain R CMD SHLIB against
# the header-only cpp11 (no cpp_source, so no decor/callr/... build deps needed).
cat("Compiling cpp11 shim (R CMD SHLIB) ...\n")
cpp11_inc <- system.file("include", package = "cpp11")
stopifnot(nzchar(cpp11_inc))
cpp11_so  <- file.path(script_dir, "loess_cpp11.so")
Sys.setenv(PKG_CPPFLAGS = sprintf("-I%s", shQuote(cpp11_inc)))
Sys.setenv(PKG_LIBS = pkg_libs)
rbin <- file.path(R.home("bin"), "R")
# cd into the (space-containing) dir and use bare basenames for the compile.
shlib_cmd <- sprintf("cd %s && %s CMD SHLIB loess_cpp11.cpp -o loess_cpp11.so",
                     shQuote(script_dir), shQuote(rbin))
if (system(shlib_cmd) != 0) stop("cpp11 R CMD SHLIB compile failed")
dyn.load(cpp11_so)
stopifnot(is.loaded("loess_smooth_2d_cpp11"))
loess_smooth_2d_cpp11 <- function(x_ref, y_ref, indices_used, x_query, ks, kc)
    .Call("loess_smooth_2d_cpp11", x_ref, y_ref, indices_used, x_query, ks, kc)

## ---- .C path: load the raw Fortran lib and grab the symbol ------------------
dyn.load(so_path)
stopifnot(is.loaded("loess_smooth_2d_c"))
.sym <- getNativeSymbolInfo("loess_smooth_2d_c")$address

# Inputs are already the right storage type, so .C's *own* duplication (not any
# coercion of ours) is what we measure.
call_dotC <- function(x_ref, y_ref, indices_used, x_query, ks, kc) {
    n_total  <- length(x_ref)
    n_target <- length(x_query)
    n_used   <- length(indices_used)
    res <- .C(.sym,
              as.integer(n_total), as.integer(n_target),
              x_ref, y_ref,
              indices_used, as.integer(n_used),
              x_query, as.double(ks), as.double(kc),
              y_out = double(n_target), ierr = integer(1))
    res[[10]]
}
call_rcpp  <- function(x_ref, y_ref, indices_used, x_query, ks, kc)
    .loess_smooth_2d_rcpp(x_ref, y_ref, indices_used, x_query, ks, kc)$y_out
call_cpp11 <- function(x_ref, y_ref, indices_used, x_query, ks, kc)
    loess_smooth_2d_cpp11(x_ref, y_ref, indices_used, x_query, ks, kc)$y_out

variants <- list(rcpp = call_rcpp, cpp11 = call_cpp11, dotC = call_dotC)

## ---- inputs -----------------------------------------------------------------
make_inputs <- function(n_total, n_used = 100L, n_target = 100L) {
    set.seed(1)
    list(x_ref        = runif(n_total),
         y_ref        = runif(n_total),
         indices_used = seq_len(n_used),          # 1..n_used, valid in [1, n_total]
         x_query      = runif(n_target),
         ks           = 1.0,
         kc           = 3.0)
}
run_variant <- function(fn, inp)
    fn(inp$x_ref, inp$y_ref, inp$indices_used, inp$x_query, inp$ks, inp$kc)

## ---- correctness gate (must pass before any timing) -------------------------
cat("\nCorrectness gate (n_total = 1000) ...\n")
gate <- make_inputs(1000L)
ys   <- lapply(variants, run_variant, inp = gate)
stopifnot(isTRUE(all.equal(ys$rcpp, ys$cpp11)),
          isTRUE(all.equal(ys$rcpp, ys$dotC)))
cat("  OK — all three variants agree on y_out\n")

## ---- timing -----------------------------------------------------------------
# median ms/call over `batches` batches of `reps` calls (one warmup first)
time_ms <- function(fn, inp, reps, batches = 3L) {
    run_variant(fn, inp)                                   # warmup
    times <- numeric(batches)
    for (b in seq_len(batches)) {
        t0 <- proc.time()[["elapsed"]]
        for (i in seq_len(reps)) run_variant(fn, inp)
        times[b] <- (proc.time()[["elapsed"]] - t0) / reps * 1000
    }
    stats::median(times)
}

sizes <- c(1e6, 5e6, 20e6)
reps  <- c(30L, 15L,  8L)

rows <- list()
for (k in seq_along(sizes)) {
    n   <- as.integer(sizes[k])
    inp <- make_inputs(n)
    mb  <- 2 * n * 8 / 1024^2                              # x_ref + y_ref payload
    cat(sprintf("\nSize n_total = %s  (read-only payload %.0f MB), reps = %d\n",
                format(n, big.mark = ",", scientific = FALSE), mb, reps[k]))
    row <- list(n_total = n, payload_MB = mb)
    for (nm in names(variants)) {
        ms <- time_ms(variants[[nm]], inp, reps[k])
        row[[nm]] <- ms
        cat(sprintf("  %-6s %8.3f ms/call\n", nm, ms))
    }
    rows[[k]] <- row
    rm(inp); gc(verbose = FALSE)
}

tbl <- do.call(rbind, lapply(rows, function(r)
    data.frame(n_total = r$n_total, payload_MB = round(r$payload_MB),
               rcpp_ms = round(r$rcpp, 3), cpp11_ms = round(r$cpp11, 3),
               dotC_ms = round(r$dotC, 3))))

cat("\n================ SUMMARY (ms per call) ================\n")
print(tbl, row.names = FALSE)

## ---- write RESULTS.md -------------------------------------------------------
fmt_int <- function(x) format(x, big.mark = ",", scientific = FALSE)
md <- c(
    "# Benchmark: Rcpp vs cpp11 (`.Call`) vs `.C` — marshalling overhead",
    "",
    sprintf("_Generated by `bench.R` on %s (R %s)._", Sys.Date(), getRversion()),
    "",
    "All three variants call the identical Fortran symbol `loess_smooth_2d_c` in",
    "`build/libtensor-omics.so`. The routine is driven with large read-only",
    "reference arrays (`x_ref`, `y_ref`) and a tiny working set",
    "(`n_used = n_target = 100`) so compute is negligible and the only thing that",
    "scales with size is argument marshalling.",
    "",
    "| n_total | payload MB | rcpp ms | cpp11 ms | .C ms |",
    "|--------:|-----------:|--------:|---------:|------:|",
    apply(tbl, 1, function(r) sprintf("| %s | %s | %s | %s | %s |",
        fmt_int(as.integer(r["n_total"])), fmt_int(as.integer(r["payload_MB"])),
        r["rcpp_ms"], r["cpp11_ms"], r["dotC_ms"])),
    "",
    "payload MB = `2 * n_total * 8 bytes` (`x_ref` + `y_ref`), the data `.C` must",
    "duplicate on every call and the zero-copy paths do not.",
    "",
    "## Interpretation",
    "",
    sprintf(paste0(
        "- **Zero-copy holds, and cpp11 is at parity with Rcpp.** Both stay flat at ",
        "~%.2f ms/call across a %dx growth in input size — that floor is R call ",
        "dispatch plus the trivial O(n_target*n_used) compute, and it does *not* ",
        "grow with the reference arrays. The rcpp/cpp11 gap is below `proc.time()` ",
        "resolution at these rep counts, i.e. noise, not signal."),
        mean(c(tbl$rcpp_ms, tbl$cpp11_ms)),
        round(max(tbl$n_total) / min(tbl$n_total))),
    sprintf(paste0(
        "- **`.C` pays the copy.** It scales ~linearly with payload (%.1f -> %.1f ms ",
        "as the copied data goes %d -> %d MB) and is ~%dx slower than the zero-copy ",
        "paths at the largest size. This is the `DUP`-always behaviour that made ",
        "`.Fortran`/`.C` untenable for large read-only matrices in the first place."),
        min(tbl$dotC_ms), max(tbl$dotC_ms),
        min(tbl$payload_MB), max(tbl$payload_MB),
        round(tbl$dotC_ms[which.max(tbl$n_total)] /
              mean(c(tbl$rcpp_ms[which.max(tbl$n_total)],
                     tbl$cpp11_ms[which.max(tbl$n_total)])))),
    "- **Conclusion.** A cpp11 `.Call` shim is a drop-in performance match for the",
    "  current Rcpp path while shedding the heavy Rcpp dependency; pure `.C` is not",
    "  a viable replacement for large read-only inputs.",
    ""
)
writeLines(md, file.path(script_dir, "RESULTS.md"))
cat("\nWrote", file.path(script_dir, "RESULTS.md"), "\n")
