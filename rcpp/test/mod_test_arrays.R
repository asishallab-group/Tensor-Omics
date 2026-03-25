source("rcpp/tensoromics_functions.R")

# --- Testing for all functions ---
test_array_wrapper <- function(tmpdir = tempdir()) {
  cat("Start Wrapper-Tests...\n")
  fn <- function(name) file.path(tmpdir, name)

  cat("Integer tests...\n")
  cat("Dim 1: ")
  arr1 <- 1:10
  # Debug: check storage mode and length
  stopifnot(is.integer(arr1))
  stopifnot(length(arr1) == 10)
  tox_serialize_int_array(arr1, fn("int1d.bin"))
  stopifnot(all(tox_deserialize_int_array(fn("int1d.bin")) == arr1))
  cat("okay\n")

  cat("Dim 2: ")
  arr2 <- matrix(1:12, nrow=3, ncol=4)
  stopifnot(is.integer(arr2))
  stopifnot(length(arr2) == 12)
  tox_serialize_int_array(arr2, fn("int2d.bin"))
  stopifnot(all(tox_deserialize_int_array(fn("int2d.bin")) == arr2))
  cat("okay\n")

  cat("Dim 3: ")
  arr3 <- array(1:24, dim = c(2,3,4))
  tox_serialize_int_array(arr3, fn("int3d.bin"))
  stopifnot(all(tox_deserialize_int_array(fn("int3d.bin")) == arr3))
  cat("okay\n")

  cat("Dim 4: ")
  arr4 <- array(1:48, dim = c(2,3,4,2))
  tox_serialize_int_array(arr4, fn("int4d.bin"))
  stopifnot(all(tox_deserialize_int_array(fn("int4d.bin")) == arr4))
  cat("okay\n")

  cat("Dim 5: ")
  arr5 <- array(1:96, dim = c(2,3,4,2,2))
  tox_serialize_int_array(arr5, fn("int5d.bin"))
  stopifnot(all(tox_deserialize_int_array(fn("int5d.bin")) == arr5))
  cat("okay\n")

  # REAL
  cat("Real tests...\n")
  arr1r <- as.numeric(1:10) * 0.5
  tox_serialize_real_array(arr1r, fn("real1d.bin"))
  stopifnot(all(tox_deserialize_real_array(fn("real1d.bin")) == arr1r))
  cat("okay\n")

  arr2r <- matrix(runif(12), nrow=3, ncol=4)
  tox_serialize_real_array(arr2r, fn("real2d.bin"))
  stopifnot(all(tox_deserialize_real_array(fn("real2d.bin")) == arr2r))
  cat("okay\n")

  arr3r <- array(runif(24), dim = c(2,3,4))
  tox_serialize_real_array(arr3r, fn("real3d.bin"))
  stopifnot(all(tox_deserialize_real_array(fn("real3d.bin")) == arr3r))
  cat("okay\n")

  arr4r <- array(runif(48), dim = c(2,3,4,2))
  tox_serialize_real_array(arr4r, fn("real4d.bin"))
  stopifnot(all(tox_deserialize_real_array(fn("real4d.bin")) == arr4r))
  cat("okay\n")

  arr5r <- array(runif(96), dim = c(2,3,4,2,2))
  tox_serialize_real_array(arr5r, fn("real5d.bin"))
  stopifnot(all(tox_deserialize_real_array(fn("real5d.bin")) == arr5r))
  cat("okay\n")
  cat("Real tests successful!\n")

  # LOGICAL
  cat("Logical tests...\n")
  cat("Dim 1: ")
  arr1l <- c(TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE)
  tox_serialize_logical_array(arr1l, fn("logical1d.bin"))
  stopifnot(all(tox_deserialize_logical_array(fn("logical1d.bin")) == arr1l))
  cat("okay\n")

  cat("Dim 2: ")
  arr2l <- matrix(c(TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE), 
                  nrow=3, ncol=4)
  tox_serialize_logical_array(arr2l, fn("logical2d.bin"))
  stopifnot(all(tox_deserialize_logical_array(fn("logical2d.bin")) == arr2l))
  cat("okay\n")

  cat("Dim 3: ")
  arr3l <- array(rep(c(TRUE, FALSE), 12), dim = c(2,3,4))
  tox_serialize_logical_array(arr3l, fn("logical3d.bin"))
  stopifnot(all(tox_deserialize_logical_array(fn("logical3d.bin")) == arr3l))
  cat("okay\n")

  cat("Dim 4: ")
  arr4l <- array(rep(c(TRUE, FALSE), 24), dim = c(2,3,4,2))
  tox_serialize_logical_array(arr4l, fn("logical4d.bin"))
  stopifnot(all(tox_deserialize_logical_array(fn("logical4d.bin")) == arr4l))
  cat("okay\n")

  cat("Dim 5: ")
  arr5l <- array(rep(c(TRUE, FALSE), 48), dim = c(2,3,4,2,2))
  tox_serialize_logical_array(arr5l, fn("logical5d.bin"))
  stopifnot(all(tox_deserialize_logical_array(fn("logical5d.bin")) == arr5l))
  cat("okay\n")
  
  # Test all TRUE and all FALSE cases
  cat("Special logical cases: ")
  all_true <- rep(TRUE, 8)
  all_false <- rep(FALSE, 8)
  tox_serialize_logical_array(all_true, fn("logical_all_true.bin"))
  tox_serialize_logical_array(all_false, fn("logical_all_false.bin"))
  stopifnot(all(tox_deserialize_logical_array(fn("logical_all_true.bin")) == all_true))
  stopifnot(all(tox_deserialize_logical_array(fn("logical_all_false.bin")) == all_false))
  cat("okay\n")
  cat("Logical tests successful!\n")

  # COMPLEX
  cat("Complex tests...\n")
  cat("Dim 1: ")
  arr1c <- complex(real = 1:10, imaginary = (1:10) * 0.5)
  tox_serialize_complex_array(arr1c, fn("complex1d.bin"))
  result1c <- tox_deserialize_complex_array(fn("complex1d.bin"))
  stopifnot(all(Re(result1c) == Re(arr1c)))
  stopifnot(all(Im(result1c) == Im(arr1c)))
  cat("okay\n")

  cat("Dim 2: ")
  arr2c <- matrix(complex(real = 1:12, imaginary = (1:12) * 0.25), nrow=3, ncol=4)
  tox_serialize_complex_array(arr2c, fn("complex2d.bin"))
  result2c <- tox_deserialize_complex_array(fn("complex2d.bin"))
  stopifnot(all(Re(result2c) == Re(arr2c)))
  stopifnot(all(Im(result2c) == Im(arr2c)))
  cat("okay\n")

  cat("Dim 3: ")
  arr3c <- array(complex(real = 1:24, imaginary = (1:24) * 0.1), dim = c(2,3,4))
  tox_serialize_complex_array(arr3c, fn("complex3d.bin"))
  result3c <- tox_deserialize_complex_array(fn("complex3d.bin"))
  stopifnot(all(Re(result3c) == Re(arr3c)))
  stopifnot(all(Im(result3c) == Im(arr3c)))
  cat("okay\n")

  cat("Dim 4: ")
  arr4c <- array(complex(real = 1:48, imaginary = (1:48) * 0.05), dim = c(2,3,4,2))
  tox_serialize_complex_array(arr4c, fn("complex4d.bin"))
  result4c <- tox_deserialize_complex_array(fn("complex4d.bin"))
  stopifnot(all(Re(result4c) == Re(arr4c)))
  stopifnot(all(Im(result4c) == Im(arr4c)))
  cat("okay\n")

  cat("Dim 5: ")
  arr5c <- array(complex(real = 1:96, imaginary = (1:96) * 0.01), dim = c(2,3,4,2,2))
  tox_serialize_complex_array(arr5c, fn("complex5d.bin"))
  result5c <- tox_deserialize_complex_array(fn("complex5d.bin"))
  stopifnot(all(Re(result5c) == Re(arr5c)))
  stopifnot(all(Im(result5c) == Im(arr5c)))
  cat("okay\n")
  
  # Test special complex cases
  cat("Special complex cases: ")
  # Pure real numbers
  pure_real <- complex(real = 1:6, imaginary = 0)
  # Pure imaginary numbers  
  pure_imag <- complex(real = 0, imaginary = 1:6)
  # Mixed complex numbers
  mixed <- complex(real = c(1, -1, 2, -2, 3, -3), imaginary = c(2, -2, 4, -4, 6, -6))
  
  tox_serialize_complex_array(pure_real, fn("complex_pure_real.bin"))
  tox_serialize_complex_array(pure_imag, fn("complex_pure_imag.bin"))
  tox_serialize_complex_array(mixed, fn("complex_mixed.bin"))
  
  stopifnot(all(tox_deserialize_complex_array(fn("complex_pure_real.bin")) == pure_real))
  stopifnot(all(tox_deserialize_complex_array(fn("complex_pure_imag.bin")) == pure_imag))
  stopifnot(all(tox_deserialize_complex_array(fn("complex_mixed.bin")) == mixed))
  cat("okay\n")
  cat("Complex tests successful!\n")

  # CHAR
  cat("Character tests...\n")
  clen <- 8
  arr1c <- sprintf("%0*d", clen, 1:10)
  tox_serialize_char_array(arr1c, fn("char1d.bin"))
  deserialized_arr1c <- tox_deserialize_char_array(fn("char1d.bin"))
  stopifnot(all(deserialized_arr1c == arr1c))

  arr2c <- matrix(sprintf("%0*d", clen, 1:12), nrow=3, ncol=4)
  tox_serialize_char_array(arr2c, fn("char2d.bin"))
  stopifnot(all(tox_deserialize_char_array(fn("char2d.bin")) == arr2c))

  arr3c <- array(sprintf("%0*d", clen, 1:24), dim = c(2,3,4))
  tox_serialize_char_array(arr3c, fn("char3d.bin"))
  stopifnot(all(tox_deserialize_char_array(fn("char3d.bin")) == arr3c))

  arr4c <- array(sprintf("%0*d", clen, 1:48), dim = c(2,3,4,2))
  tox_serialize_char_array(arr4c, fn("char4d.bin"))
  stopifnot(all(tox_deserialize_char_array(fn("char4d.bin")) == arr4c))

  arr5c <- array(sprintf("%0*d", clen, 1:96), dim = c(2,3,4,2,2))
  tox_serialize_char_array(arr5c, fn("char5d.bin"))
  stopifnot(all(tox_deserialize_char_array(fn("char5d.bin")) == arr5c))

  # 1D-Array with different char lengths
  arr_ascii1 <- c("A", "G1", "GENE003", "BRCA1", "XYZ", "", "12345678", "SEQ")
  tox_serialize_char_array(arr_ascii1, fn("char_ascii1d.bin"))
  stopifnot(all(tox_deserialize_char_array(fn("char_ascii1d.bin")) == arr_ascii1))

  # 2D-Matrix with different length
  arr_ascii2 <- matrix(c("GENE1", "GENE22", "GENE333", "", "ID", "SEQ9999"), nrow = 2, byrow = TRUE)
  tox_serialize_char_array(arr_ascii2, fn("char_ascii2d.bin"))
  stopifnot(all(tox_deserialize_char_array(fn("char_ascii2d.bin")) == arr_ascii2))

  # 3D-Array with realistic data
  arr_ascii3 <- array(c("TP53", "BRCA1", "MT-ATP6", "CYTB", "ND1", "", "NRAS", "EGFR"), dim = c(2, 2, 2))
  tox_serialize_char_array(arr_ascii3, fn("char_ascii3d.bin"))
  stopifnot(all(tox_deserialize_char_array(fn("char_ascii3d.bin")) == arr_ascii3))

  cat("All ASCII string variation tests passed!\n")
  
  cat("\n=== ALL TESTS PASSED SUCCESSFULLY! ===\n")
  cat("Tested: Integer, Real, Logical, Complex, and Character arrays\n")
  cat("Dimensions: 1D through 5D\n")
  cat("All data types working correctly!\n")
}

# start tests
test_array_wrapper()