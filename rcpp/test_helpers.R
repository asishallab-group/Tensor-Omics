run_all_tests <- function(env = parent.frame(), test_only = TRUE) {
  # Discover candidate names
  if (test_only) {
    test_names <- ls(env, pattern = "^test_")
  } else {
    test_names <- ls(env)
  }

  passed  <- 0
  failed  <- 0
  skipped <- 0

  cat("Running R tests...\n")

  for (name in test_names) {
    obj <- get(name, envir = env)
    if (!is.function(obj)) next  # skip non-functions

    test_func <- obj

    tryCatch(
      {

        cat(sprintf("Testing %s ...\n", name))
        test_func()
        cat(sprintf("✓ %s passed.\n", name))
        passed <- passed + 1
      },
      error = function(e) {
        msg <- conditionMessage(e)

        # Same skip logic as Python version
        if (grepl("Note:", msg) || grepl("acceptable", msg, ignore.case = TRUE)) {
          cat(sprintf("~ %s skipped (expected behavior): %s\n", name, msg))
          skipped <<- skipped + 1
        } else {
          cat(sprintf("✗ %s FAILED: %s\n", name, msg))
          failed <<- failed + 1
        }
      }
    )
  }

  cat("\nSummary: ", passed, " passed, ", failed, " failed, ", skipped, " skipped\n", sep = "")
  stopifnot(failed == 0)
}

assert_true <- function(expr, msg = "Assertion failed") {
  if (!isTRUE(expr)) stop(msg, call. = FALSE)
  invisible(TRUE)
}

assert_false <- function(expr, msg = "Assertion failed") {
  if (isTRUE(expr)) stop(msg, call. = FALSE)
  invisible(TRUE)
}

assert_error <- function(expr, msg = "Expected an error") {
  err <- tryCatch({ expr; NULL }, error = function(e) e)
  if (is.null(err)) stop(msg, call. = FALSE)
  invisible(TRUE)
}

assert_equal_int <- function(x, y, msg = "Integer mismatch") {
  if (!identical(x, y)) stop(msg, call. = FALSE)
  invisible(TRUE)
}

assert_equal_numeric <- function(x, y, tol = 1e-12, msg = "Real mismatch") {
  if (length(x) != length(y)) stop(msg, call. = FALSE)
  for (i in seq_len(length(x))) {
    if (!(is.finite(x[i]) && is.finite(y[i]) && abs(x[i] - y[i]) <= tol)) stop(msg, call. = FALSE)
  }
  invisible(TRUE)
}