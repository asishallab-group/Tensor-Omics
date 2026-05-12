source("rcpp/tensoromics_functions.R")
source("rcpp/test_helpers.R")

fn <- function(name, tmpdir = tempdir()) file.path(tmpdir, name)

test_integer_serialization <- function() {
  arr1 <- 1:10
  # Debug: check storage mode and length
  assert_true(is.integer(arr1))
  assert_true(length(arr1) == 10)
  tox_serialize_int_array(arr1, fn("int1d.test.bin"))
  assert_true(all(tox_deserialize_int_array(fn("int1d.test.bin")) == arr1))

  arr2 <- matrix(1:12, nrow=3, ncol=4)
  assert_true(is.integer(arr2))
  assert_true(length(arr2) == 12)
  tox_serialize_int_array(arr2, fn("int2d.test.bin"))
  assert_true(all(tox_deserialize_int_array(fn("int2d.test.bin")) == arr2))

  arr3 <- array(1:24, dim = c(2,3,4))
  tox_serialize_int_array(arr3, fn("int3d.test.bin"))
  assert_true(all(tox_deserialize_int_array(fn("int3d.test.bin")) == arr3))

  arr4 <- array(1:48, dim = c(2,3,4,2))
  tox_serialize_int_array(arr4, fn("int4d.test.bin"))
  assert_true(all(tox_deserialize_int_array(fn("int4d.test.bin")) == arr4))

  arr5 <- array(1:96, dim = c(2,3,4,2,2))
  tox_serialize_int_array(arr5, fn("int5d.test.bin"))
  assert_true(all(tox_deserialize_int_array(fn("int5d.test.bin")) == arr5))
}
test_real_serialization <- function() {
  arr1r <- as.numeric(1:10) * 0.5
  tox_serialize_real_array(arr1r, fn("real1d.test.bin"))
  assert_true(all(tox_deserialize_real_array(fn("real1d.test.bin")) == arr1r))

  arr2r <- matrix(runif(12), nrow=3, ncol=4)
  tox_serialize_real_array(arr2r, fn("real2d.test.bin"))
  assert_true(all(tox_deserialize_real_array(fn("real2d.test.bin")) == arr2r))

  arr3r <- array(runif(24), dim = c(2,3,4))
  tox_serialize_real_array(arr3r, fn("real3d.test.bin"))
  assert_true(all(tox_deserialize_real_array(fn("real3d.test.bin")) == arr3r))

  arr4r <- array(runif(48), dim = c(2,3,4,2))
  tox_serialize_real_array(arr4r, fn("real4d.test.bin"))
  assert_true(all(tox_deserialize_real_array(fn("real4d.test.bin")) == arr4r))

  arr5r <- array(runif(96), dim = c(2,3,4,2,2))
  tox_serialize_real_array(arr5r, fn("real5d.test.bin"))
  assert_true(all(tox_deserialize_real_array(fn("real5d.test.bin")) == arr5r))
}
test__serialization <- function() {
  arr1l <- c(TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE)
  tox_serialize_logical_array(arr1l, fn("logical1d.test.bin"))
  assert_true(all(tox_deserialize_logical_array(fn("logical1d.test.bin")) == arr1l))

  arr2l <- matrix(c(TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE), 
                  nrow=3, ncol=4)
  tox_serialize_logical_array(arr2l, fn("logical2d.test.bin"))
  assert_true(all(tox_deserialize_logical_array(fn("logical2d.test.bin")) == arr2l))

  arr3l <- array(rep(c(TRUE, FALSE), 12), dim = c(2,3,4))
  tox_serialize_logical_array(arr3l, fn("logical3d.test.bin"))
  assert_true(all(tox_deserialize_logical_array(fn("logical3d.test.bin")) == arr3l))

  arr4l <- array(rep(c(TRUE, FALSE), 24), dim = c(2,3,4,2))
  tox_serialize_logical_array(arr4l, fn("logical4d.test.bin"))
  assert_true(all(tox_deserialize_logical_array(fn("logical4d.test.bin")) == arr4l))

  arr5l <- array(rep(c(TRUE, FALSE), 48), dim = c(2,3,4,2,2))
  tox_serialize_logical_array(arr5l, fn("logical5d.test.bin"))
  assert_true(all(tox_deserialize_logical_array(fn("logical5d.test.bin")) == arr5l))
  
  # Test all TRUE and all FALSE cases
  all_true <- rep(TRUE, 8)
  all_false <- rep(FALSE, 8)
  tox_serialize_logical_array(all_true, fn("logical_all_true.test.bin"))
  tox_serialize_logical_array(all_false, fn("logical_all_false.test.bin"))
  assert_true(all(tox_deserialize_logical_array(fn("logical_all_true.test.bin")) == all_true))
  assert_true(all(tox_deserialize_logical_array(fn("logical_all_false.test.bin")) == all_false))
}
test_complex_serialization <- function() {
  arr1c <- complex(real = 1:10, imaginary = (1:10) * 0.5)
  tox_serialize_complex_array(arr1c, fn("complex1d.test.bin"))
  result1c <- tox_deserialize_complex_array(fn("complex1d.test.bin"))
  assert_true(all(Re(result1c) == Re(arr1c)))
  assert_true(all(Im(result1c) == Im(arr1c)))

  arr2c <- matrix(complex(real = 1:12, imaginary = (1:12) * 0.25), nrow=3, ncol=4)
  tox_serialize_complex_array(arr2c, fn("complex2d.test.bin"))
  result2c <- tox_deserialize_complex_array(fn("complex2d.test.bin"))
  assert_true(all(Re(result2c) == Re(arr2c)))
  assert_true(all(Im(result2c) == Im(arr2c)))

  arr3c <- array(complex(real = 1:24, imaginary = (1:24) * 0.1), dim = c(2,3,4))
  tox_serialize_complex_array(arr3c, fn("complex3d.test.bin"))
  result3c <- tox_deserialize_complex_array(fn("complex3d.test.bin"))
  assert_true(all(Re(result3c) == Re(arr3c)))
  assert_true(all(Im(result3c) == Im(arr3c)))

  arr4c <- array(complex(real = 1:48, imaginary = (1:48) * 0.05), dim = c(2,3,4,2))
  tox_serialize_complex_array(arr4c, fn("complex4d.test.bin"))
  result4c <- tox_deserialize_complex_array(fn("complex4d.test.bin"))
  assert_true(all(Re(result4c) == Re(arr4c)))
  assert_true(all(Im(result4c) == Im(arr4c)))

  arr5c <- array(complex(real = 1:96, imaginary = (1:96) * 0.01), dim = c(2,3,4,2,2))
  tox_serialize_complex_array(arr5c, fn("complex5d.test.bin"))
  result5c <- tox_deserialize_complex_array(fn("complex5d.test.bin"))
  assert_true(all(Re(result5c) == Re(arr5c)))
  assert_true(all(Im(result5c) == Im(arr5c)))
  
  # Test special complex cases
  # Pure real numbers
  pure_real <- complex(real = 1:6, imaginary = 0)
  # Pure imaginary numbers  
  pure_imag <- complex(real = 0, imaginary = 1:6)
  # Mixed complex numbers
  mixed <- complex(real = c(1, -1, 2, -2, 3, -3), imaginary = c(2, -2, 4, -4, 6, -6))
  
  tox_serialize_complex_array(pure_real, fn("complex_pure_real.test.bin"))
  tox_serialize_complex_array(pure_imag, fn("complex_pure_imag.test.bin"))
  tox_serialize_complex_array(mixed, fn("complex_mixed.test.bin"))
  
  assert_true(all(tox_deserialize_complex_array(fn("complex_pure_real.test.bin")) == pure_real))
  assert_true(all(tox_deserialize_complex_array(fn("complex_pure_imag.test.bin")) == pure_imag))
  assert_true(all(tox_deserialize_complex_array(fn("complex_mixed.test.bin")) == mixed))
}
test_character_serialization <- function() {
  clen <- 8
  arr1c <- sprintf("%0*d", clen, 1:10)
  tox_serialize_char_array(arr1c, fn("char1d.test.bin"))
  deserialized_arr1c <- tox_deserialize_char_array(fn("char1d.test.bin"))
  assert_true(all(deserialized_arr1c == arr1c))
  arr2c <- matrix(sprintf("%0*d", clen, 1:12), nrow=3, ncol=4)
  tox_serialize_char_array(arr2c, fn("char2d.test.bin"))
  assert_true(all(tox_deserialize_char_array(fn("char2d.test.bin")) == arr2c))

  arr3c <- array(sprintf("%0*d", clen, 1:24), dim = c(2,3,4))
  tox_serialize_char_array(arr3c, fn("char3d.test.bin"))
  assert_true(all(tox_deserialize_char_array(fn("char3d.test.bin")) == arr3c))

  arr4c <- array(sprintf("%0*d", clen, 1:48), dim = c(2,3,4,2))
  tox_serialize_char_array(arr4c, fn("char4d.test.bin"))
  assert_true(all(tox_deserialize_char_array(fn("char4d.test.bin")) == arr4c))

  arr5c <- array(sprintf("%0*d", clen, 1:96), dim = c(2,3,4,2,2))
  tox_serialize_char_array(arr5c, fn("char5d.test.bin"))
  assert_true(all(tox_deserialize_char_array(fn("char5d.test.bin")) == arr5c))

  # 1D-Array with different char lengths
  arr_ascii1 <- c("A", "G1", "GENE003", "BRCA1", "XYZ", "", "12345678", "SEQ")
  tox_serialize_char_array(arr_ascii1, fn("char_ascii1d.test.bin"))
  assert_true(all(tox_deserialize_char_array(fn("char_ascii1d.test.bin")) == arr_ascii1))

  # 2D-Matrix with different length
  arr_ascii2 <- matrix(c("GENE1", "GENE22", "GENE333", "", "ID", "SEQ9999"), nrow = 2, byrow = TRUE)
  tox_serialize_char_array(arr_ascii2, fn("char_ascii2d.test.bin"))
  assert_true(all(tox_deserialize_char_array(fn("char_ascii2d.test.bin")) == arr_ascii2))

  # 3D-Array with realistic data
  arr_ascii3 <- array(c("TP53", "BRCA1", "MT-ATP6", "CYTB", "ND1", "", "NRAS", "EGFR"), dim = c(2, 2, 2))
  tox_serialize_char_array(arr_ascii3, fn("char_ascii3d.test.bin"))
  assert_true(all(tox_deserialize_char_array(fn("char_ascii3d.test.bin")) == arr_ascii3))
}

run_all_tests()