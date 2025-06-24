# R Test Script for F42CsvReader Library (Revised for Robustness)

# 1. Define constants and library path
# Make sure this path is correct for your project structure.
SHARED_LIB_PATH <- "./build/libF42CsvReader.so"
CSV_FILE <- "data.csv"

# --- Debugging Steps ---
# Check R's current working directory. The CSV file must be here.
print(paste("Current R working directory:", getwd()))
if (!file.exists(CSV_FILE)) {
  stop(paste("CRITICAL ERROR: Cannot find '", CSV_FILE, "' in the current directory.",
             "Please run this script from your project's root directory.", sep=""))
} else {
  print(paste("Found '", CSV_FILE, "'. Proceeding.", sep=""))
}
# --- End Debugging ---


# 2. Load the shared library
# On Windows, this would be "libF42CsvReader.dll"
# On macOS, "libF42CsvReader.dylib"
dyn.load(SHARED_LIB_PATH)

# Convert the filename to an absolute path.
ABSOLUTE_CSV_PATH <- normalizePath(CSV_FILE)
print(paste("Passing absolute path to Fortran:", ABSOLUTE_CSV_PATH))

# 3. Call the CSV reader using POSITIONAL arguments
# The warnings about "passing a char vector" are normal and can usually be ignored.
# The order of arguments must exactly match the Fortran subroutine definition:
# SUBROUTINE read_csv_file_R(filename, has_header, delimiter, io_status)
io_status <- 0
result_read <- .Fortran("read_csv_file_R",
                        as.character(ABSOLUTE_CSV_PATH),
                        as.logical(TRUE),
                        as.character(","),
                        as.integer(io_status))

# The returned list is also positional, io_status is the 4th argument
if (result_read[[4]] != 0) {
  stop(paste("Error reading CSV file. IO Status:", result_read[[4]]))
}
print("CSV file read successfully.")

# 4. Get dimensions
num_rows <- .Fortran("get_num_rows_R", as.integer(0))[[1]]
num_cols <- .Fortran("get_num_cols_R", as.integer(0))[[1]]
print(paste("Dimensions:", num_rows, "rows x", num_cols, "cols"))

# 5. Get and print header
header <- ""
for (j in 1:num_cols) {
  # Pre-allocate space for the character return value
  h_name <- paste(rep(" ", 512), collapse = "")
  h_name <- .Fortran("get_header_R", as.integer(j), h_name)[[2]]
  header <- paste(header, trimws(h_name), sep = "\t")
}
print("Header:")
print(header)
print(paste(rep("-", 50), collapse=""))


# 6. Get and print cell data
print("Data:")
for (i in 1:num_rows) {
  row_str <- ""
  for (j in 1:num_cols) {
    # .Fortran needs variables pre-defined for all potential return types
    cell_result <- .Fortran("get_cell_R",
                             as.integer(i),
                             as.integer(j),
                             as.integer(0),
                             as.integer(0),
                             as.double(0.0),
                             paste(rep(" ", 512), collapse = ""))

    # Process the returned value based on its type (using positional indices)
    # 3:data_type, 4:i_val, 5:r_val, 6:c_val
    cell_str <- ""
    if (cell_result[[3]] == 1) { # Integer
      cell_str <- as.character(cell_result[[4]])
    } else if (cell_result[[3]] == 2) { # Real
      cell_str <- as.character(cell_result[[5]])
    } else if (cell_result[[3]] == 3) { # Character
      cell_str <- trimws(cell_result[[6]])
    }
    row_str <- paste(row_str, cell_str, sep = "\t")
  }
  print(trimws(row_str))
}

# 7. IMPORTANT: Call the cleanup routine
.Fortran("cleanup_csv_data_R")
print("Cleanup function called. Memory released.")

# 8. Unload the library (optional)
dyn.unload(SHARED_LIB_PATH)
