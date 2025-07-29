# Test script for the Fortran CSV reader library.
# This script loads the shared library and calls the R-compatible
# Fortran subroutines to test the full data reading pipeline.

# --- 1. Library Loading ---

# The build.sh script creates a symlink in the build directory
lib_path <- file.path("build", "libtensor-omics.so")

if (!file.exists(lib_path)) {
  stop(paste("Shared library not found at:", lib_path, 
             "\nPlease build the Fortran project using the './build.sh' script first."))
}
dyn.load(lib_path)

# --- 2. Helper Functions for Data Conversion ---

MAX_FIELD_LEN <- 512 # This MUST match the Fortran parameter

#' Converts an R string to an integer vector of ASCII codes.
string_to_ascii <- function(s) {
  if (nchar(s) == 0) return(integer(0))
  return(as.integer(charToRaw(s)))
}

#' Converts a 3D integer array of ASCII codes from Fortran back to an R character matrix.
ascii_to_string_matrix <- function(ascii_array) {
  dims <- dim(ascii_array)
  num_rows <- dims[2]
  num_cols <- dims[3]
  mat <- matrix("", nrow = num_rows, ncol = num_cols)
  
  for (j in 1:num_cols) {
    for (i in 1:num_rows) {
      raw_vec <- ascii_array[, i, j]
      null_pos <- which(raw_vec == 0)
      end_pos <- if (length(null_pos) > 0) null_pos[1] - 1 else length(raw_vec)
      if (end_pos > 0) {
        mat[i, j] <- rawToChar(as.raw(raw_vec[1:end_pos]))
      }
    }
  }
  return(mat)
}

# --- 3. Test Function ---

test_full_pipeline <- function() {
  cat("=== Testing Full CSV Reading Pipeline via R Interface ===\n")
  
  csv_content <- paste(
    "ID,Status,Score,Name,Coords",
    "101,T,95.5,Alpha,(1.1, 2.2)",
    "102,F,88.0,Beta,(3.3, -4.4)",
    "103,TRUE,72.1,Gamma,(5.5, 6.6)",
    "104,FALSE,-5.0,Delta,(7.7, 8.8)",
    sep = "\n"
  )
  test_filename <- "r_pipeline_test.csv"
  writeLines(csv_content, test_filename)
  
  tryCatch({
    # Step 1: Get dimensions
    fname_ascii <- string_to_ascii(test_filename)
    res_dims <- .Fortran("get_csv_dims_r",
                         filename_ascii = as.integer(fname_ascii),
                         fn_len = as.integer(length(fname_ascii)),
                         has_header = as.logical(TRUE),
                         delimiter_ascii = as.integer(charToRaw(",")),
                         num_rows = integer(1),
                         num_cols = integer(1),
                         status = integer(1))
    
    if (res_dims$status != 0) stop(paste("get_csv_dims_r failed with status", res_dims$status))
    cat("Step 1: get_csv_dims_r PASSED\n")
    
    nr <- res_dims$num_rows
    nc <- res_dims$num_cols
    
    # Step 2: Read all data as strings
    data_out_ascii <- array(0L, dim = c(MAX_FIELD_LEN, nr, nc))
    res_read <- .Fortran("read_csv_to_strings_r",
                         filename_ascii = as.integer(fname_ascii),
                         fn_len = as.integer(length(fname_ascii)),
                         has_header = as.logical(TRUE),
                         delimiter_ascii = as.integer(charToRaw(",")),
                         num_rows = as.integer(nr),
                         num_cols = as.integer(nc),
                         header_out_ascii = array(0L, dim=c(MAX_FIELD_LEN, nc)),
                         data_out_ascii = data_out_ascii,
                         status = integer(1))

    if (res_read$status != 0) stop(paste("read_csv_to_strings_r failed with status", res_read$status))
    data_str_matrix <- ascii_to_string_matrix(res_read$data_out_ascii)
    stopifnot(data_str_matrix[2, 4] == "Beta")
    cat("Step 2: read_csv_to_strings_r PASSED\n")

    # Step 3: Convert to specific types and test
    
    # Integers (Column 1)
    cols_to_read <- 1L
    int_out <- matrix(0L, nrow = nr, ncol = length(cols_to_read))
    res_int <- .Fortran("read_integer_columns_r",
                        data_in_ascii = res_read$data_out_ascii,
                        num_rows = as.integer(nr), num_cols_in = as.integer(nc),
                        cols_to_read = as.integer(cols_to_read), num_cols_to_read = as.integer(length(cols_to_read)),
                        int_data_out = int_out, status = integer(1))
    if (res_int$status != 0) stop(paste("read_integer_columns_r failed with status", res_int$status))
    stopifnot(all(as.vector(res_int$int_data_out) == c(101, 102, 103, 104)))
    cat("Step 3: read_integer_columns_r PASSED\n")

    # Reals (Column 3)
    cols_to_read <- 3L
    real_out <- matrix(0.0, nrow = nr, ncol = length(cols_to_read))
    res_real <- .Fortran("read_real_columns_r",
                         data_in_ascii = res_read$data_out_ascii,
                         num_rows = as.integer(nr), num_cols_in = as.integer(nc),
                         cols_to_read = as.integer(cols_to_read), num_cols_to_read = as.integer(length(cols_to_read)),
                         real_data_out = real_out, status = integer(1))
    if (res_real$status != 0) stop(paste("read_real_columns_r failed with status", res_real$status))
    # CORRECTED: Convert matrix to vector before comparison
    stopifnot(all.equal(as.vector(res_real$real_data_out), c(95.5, 88.0, 72.1, -5.0)))
    cat("Step 4: read_real_columns_r PASSED\n")

    # Logicals (Column 2)
    cols_to_read <- 2L
    logical_out <- matrix(FALSE, nrow = nr, ncol = length(cols_to_read))
    res_logical <- .Fortran("read_logical_columns_r",
                            data_in_ascii = res_read$data_out_ascii,
                            num_rows = as.integer(nr), num_cols_in = as.integer(nc),
                            cols_to_read = as.integer(cols_to_read), num_cols_to_read = as.integer(length(cols_to_read)),
                            logical_data_out = logical_out, status = integer(1))
    if (res_logical$status != 0) stop(paste("read_logical_columns_r failed with status", res_logical$status))
    stopifnot(all(as.vector(res_logical$logical_data_out) == c(TRUE, FALSE, TRUE, FALSE)))
    cat("Step 5: read_logical_columns_r PASSED\n")

    # Characters (Column 4)
    cols_to_read <- 4L
    char_out_ascii <- array(0L, dim = c(MAX_FIELD_LEN, nr, length(cols_to_read)))
    res_char <- .Fortran("read_character_columns_r",
                         data_in_ascii = res_read$data_out_ascii,
                         num_rows = as.integer(nr), num_cols_in = as.integer(nc),
                         cols_to_read = as.integer(cols_to_read), num_cols_to_read = as.integer(length(cols_to_read)),
                         char_data_out_ascii = char_out_ascii, status = integer(1))
    if (res_char$status != 0) stop(paste("read_character_columns_r failed with status", res_char$status))
    char_out_str <- ascii_to_string_matrix(res_char$char_data_out_ascii)
    stopifnot(all(as.vector(char_out_str) == c("Alpha", "Beta", "Gamma", "Delta")))
    cat("Step 6: read_character_columns_r PASSED\n")

    # Complex (Column 5)
    cols_to_read <- 5L
    complex_out <- matrix(complex(real=0, imaginary=0), nrow = nr, ncol = length(cols_to_read))
    res_complex <- .Fortran("read_complex_columns_r",
                            data_in_ascii = res_read$data_out_ascii,
                            num_rows = as.integer(nr), num_cols_in = as.integer(nc),
                            cols_to_read = as.integer(cols_to_read), num_cols_to_read = as.integer(length(cols_to_read)),
                            complex_data_out = complex_out, status = integer(1))
    if (res_complex$status != 0) stop(paste("read_complex_columns_r failed with status", res_complex$status))
    stopifnot(all.equal(as.vector(res_complex$complex_data_out), c(1.1+2.2i, 3.3-4.4i, 5.5+6.6i, 7.7+8.8i)))
    cat("Step 7: read_complex_columns_r PASSED\n")

  }, finally = {
    if (file.exists(test_filename)) file.remove(test_filename)
  })
}

# --- 4. Main Execution Block ---
main <- function() {
  cat("=================================================\n")
  cat("        R WRAPPER TESTS FOR FORTRAN CSV READER\n")
  cat("=================================================\n")
  
  tryCatch({
    test_full_pipeline()
    cat("\n=================================================\n")
    cat("           ALL R TESTS COMPLETED\n")
    cat("=================================================\n")
    cat("If you see this message, all R interface tests passed! ✓\n")
  }, error = function(e) {
    cat("\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")
    cat("                  A TEST FAILED\n")
    cat("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")
    cat("ERROR:", conditionMessage(e), "\n")
  })
}

main()
