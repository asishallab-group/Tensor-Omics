# test_read_csv.R

# Source the new API functions
source("r/tensoromics_functions.R")

test_full_pipeline <- function() {
  cat("=== Testing Full CSV Reading Pipeline via R API ===\n")
  
  csv_content <- paste(
    "ID,Status,Score,Name,Coords",
    "101,T,95.5,Alpha,(1.1, 2.2)",
    "102,F,88.0,Beta,(3.3, -4.4)",
    "103,TRUE,72.1,Gamma,(5.5, 6.6)",
    "104,FALSE,-5.0,Delta,(7.7, 8.8)",
    sep = "\n"
  )
  test_filename <- "r_api_test.csv"
  writeLines(csv_content, test_filename)
  
  tryCatch({
    # Step 1: Read raw ASCII data
    raw_csv <- tox.read_csv_raw(test_filename)
    cat("Step 1: tox.read_csv_raw PASSED\n")
    
    # Step 2: Convert raw data to strings for verification
    header <- tox.convert_ascii_to_strings(raw_csv$header)
    data_str <- tox.convert_ascii_to_strings(raw_csv$data)
    stopifnot(all(header == c("ID", "Status", "Score", "Name", "Coords")))
    stopifnot(data_str[2, 4] == "Beta")
    cat("Step 2: tox.convert_ascii_to_strings PASSED\n")

    # Step 3: Convert raw ASCII data to specific types
    
    # Integers (Column 1)
    int_data <- tox.read_integer_columns(raw_csv$data, columns = 1)
    stopifnot(all(int_data == c(101, 102, 103, 104)))
    cat("Step 3: tox.read_integer_columns PASSED\n")

    # Reals (Column 3)
    real_data <- tox.read_real_columns(raw_csv$data, columns = 3)
    stopifnot(all.equal(real_data, matrix(c(95.5, 88.0, 72.1, -5.0))))
    cat("Step 4: tox.read_real_columns PASSED\n")

    # Logicals (Column 2)
    logical_data <- tox.read_logical_columns(raw_csv$data, columns = 2)
    stopifnot(all(logical_data == c(TRUE, FALSE, TRUE, FALSE)))
    cat("Step 5: tox.read_logical_columns PASSED\n")

    # Characters (Column 4)
    char_data <- tox.read_character_columns(raw_csv$data, columns = 4)
    stopifnot(all(char_data == c("Alpha", "Beta", "Gamma", "Delta")))
    cat("Step 6: tox.read_character_columns PASSED\n")

    # Complex (Column 5)
    complex_data <- tox.read_complex_columns(raw_csv$data, columns = 5)
    stopifnot(all.equal(complex_data, matrix(c(1.1+2.2i, 3.3-4.4i, 5.5+6.6i, 7.7+8.8i))))
    cat("Step 7: tox.read_complex_columns PASSED\n")

  }, finally = {
    if (file.exists(test_filename)) file.remove(test_filename)
  })
}

# Run the tests
test_full_pipeline()
cat("\nAll R API tests passed! ✓\n")