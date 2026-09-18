# R binding and tests for clock hand angle calculations
# Uses Fortran wrappers for RAP projections and angle calculations

# Load the compiled Fortran library
source("r/load_tensor_omics.R")
source("r/test_helpers.R")

# The reference that makes a counter-clockwise turn in the (1, 2) plane positive: `v`
# rotated a quarter turn there. In two dimensions this is the familiar determinant
# convention; above two there is no canonical quarter turn, which is exactly why the
# caller has to state one.
ccw_reference <- function(v) {
  reference <- numeric(length(v))
  reference[1] <- -v[2]
  reference[2] <- v[1]
  reference
}

# Constants
PI <- pi
TOL <- 1e-12


# ==================== 2D TESTS ====================

test_identical_vectors_2d <- function() {
  v1 <- c(1.0, 0.0)
  v2 <- c(1.0, 0.0)
  signed_angle <- clock_hand_angle_between_vectors(v1, v2, ccw_reference(v1))
  assert_equal_numeric(signed_angle, 0.0, TOL, "Identical 2D vectors should give 0 angle")
}

test_opposite_vectors_2d <- function() {
  v1 <- c(1.0, 0.0)
  v2 <- c(-1.0, 0.0)
  signed_angle <- clock_hand_angle_between_vectors(v1, v2, ccw_reference(v1))
  assert_equal_numeric(abs(signed_angle), PI, TOL, "Opposite 2D vectors should give ±π")
}

test_perpendicular_vectors_2d <- function() {
  v1 <- c(1.0, 0.0)
  v2 <- c(0.0, 1.0)
  signed_angle <- clock_hand_angle_between_vectors(v1, v2, ccw_reference(v1))
  assert_equal_numeric(abs(signed_angle), PI/2, TOL, "Perpendicular 2D vectors magnitude")
  assert_true(signed_angle > 0, "Counterclockwise rotation should be positive")
}

test_45_degree_rotation_2d <- function() {
  v1 <- c(1.0, 0.0)
  v2 <- c(sqrt(2)/2, sqrt(2)/2)  # 45 degrees
  expected <- PI/4
  signed_angle <- clock_hand_angle_between_vectors(v1, v2, ccw_reference(v1))
  assert_equal_numeric(signed_angle, expected, TOL, "45-degree counterclockwise rotation")
}

test_clockwise_vs_counterclockwise_2d <- function() {
  v1 <- c(1.0, 0.0)
  v2_ccw <- c(0.0, 1.0)   # 90° counterclockwise
  v2_cw <- c(0.0, -1.0)   # 90° clockwise
  angle_ccw <- clock_hand_angle_between_vectors(v1, v2_ccw, ccw_reference(v1))
  angle_cw <- clock_hand_angle_between_vectors(v1, v2_cw, ccw_reference(v1))
  assert_true(angle_ccw > 0, "Counterclockwise should be positive")
  assert_true(angle_cw < 0, "Clockwise should be negative")
  assert_equal_numeric(abs(angle_ccw), abs(angle_cw), TOL, "Magnitudes should be equal")
}

# ==================== 3D TESTS ====================

test_identical_vectors_3d <- function() {
  v1 <- c(1.0, 1.0, 1.0)
  v2 <- c(1.0, 1.0, 1.0)
  signed_angle <- clock_hand_angle_between_vectors(v1, v2, ccw_reference(v1))
  assert_equal_numeric(signed_angle, 0.0, TOL, "Identical 3D vectors should give 0 angle")
}

test_perpendicular_vectors_3d <- function() {
  v1 <- c(1.0, 0.0, 0.0)
  v2 <- c(0.0, 1.0, 0.0)
  signed_angle <- clock_hand_angle_between_vectors(v1, v2, ccw_reference(v1))
  assert_equal_numeric(abs(signed_angle), PI/2, TOL, "Perpendicular 3D vectors")
}

test_arbitrary_3d_rotation <- function() {
  v1 <- c(1.0, 2.0, 3.0)
  v2 <- c(2.0, 1.0, 3.0)
  v1 <- v1 / sqrt(sum(v1^2))
  v2 <- v2 / sqrt(sum(v2^2))
  signed_angle <- clock_hand_angle_between_vectors(v1, v2, ccw_reference(v1))
  dot_product <- sum(v1 * v2)
  expected_magnitude <- acos(max(-1, min(1, dot_product)))
  assert_equal_numeric(abs(signed_angle), expected_magnitude, TOL, "3D arbitrary rotation magnitude")
}

# ==================== HIGH DIMENSIONAL TESTS ====================

test_high_dimensional_basic <- function() {
  v1 <- c(1.0, 0.0, 0.0, 0.0, 0.0)
  v2 <- c(0.0, 1.0, 0.0, 0.0, 0.0)
  signed_angle <- clock_hand_angle_between_vectors(v1, v2, ccw_reference(v1))
  assert_equal_numeric(abs(signed_angle), PI/2, TOL, "High-dimensional perpendicular vectors")
}

test_high_dimensional_reference_orients <- function() {
  # The turn is in the (3, 5) plane, so only the target itself says which way round it goes
  v1 <- c(0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0)
  v2 <- c(0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0)
  signed_angle <- clock_hand_angle_between_vectors(v1, v2, v2)
  assert_equal_numeric(signed_angle, PI/2, TOL, "A reference along the turn signs it positive")
}

test_reference_that_orients_nothing <- function() {
  v1 <- c(1.0, 0.0, 0.0, 0.0, 0.0)
  v2 <- c(0.0, 1.0, 0.0, 0.0, 0.0)
  # the turn is in the (1, 2) plane; axis 3 says nothing about which way round it is
  reference <- c(0.0, 0.0, 1.0, 0.0, 0.0)
  assert_error(clock_hand_angle_between_vectors(v1, v2, reference),
               "A reference orthogonal to the rotation cannot sign it", ERR_INVALID_INPUT)
}

test_reference_picks_the_orientation <- function() {
  v1 <- c(1.0, 0.0, 0.0, 0.0, 0.0)
  v2 <- c(0.0, 1.0, 0.0, 0.0, 0.0)
  reference <- ccw_reference(v1)

  one_way <- clock_hand_angle_between_vectors(v1, v2, reference)
  the_other <- clock_hand_angle_between_vectors(v1, v2, -reference)

  assert_equal_numeric(one_way, PI/2, TOL, "The stated reference signs the turn positive")
  assert_equal_numeric(the_other, -PI/2, TOL, "Reversing it flips the sign")
}

# ==================== EDGE CASES ====================

test_denormalized_vectors <- function() {
  # Large magnitude vectors
  v1 <- c(100.0, 0.0)
  v2 <- c(0.0, 50.0)
  
  signed_angle <- clock_hand_angle_between_vectors(v1, v2, ccw_reference(v1))
  assert_equal_numeric(abs(signed_angle), PI/2, TOL, "Denormalized vectors should work")
}

test_tiny_vectors_precision <- function() {
  tiny <- 1e-14
  v1 <- c(tiny, 0.0)
  v2 <- c(0.0, tiny)
  
  signed_angle <- clock_hand_angle_between_vectors(v1, v2, ccw_reference(v1))
  assert_equal_numeric(abs(signed_angle), PI/2, 1e-10, "Tiny vectors precision")
}

test_huge_vectors_precision <- function() {
  huge_val <- 1e14
  v1 <- c(huge_val, 0.0)
  v2 <- c(0.0, huge_val)
  
  signed_angle <- clock_hand_angle_between_vectors(v1, v2, ccw_reference(v1))
  assert_equal_numeric(abs(signed_angle), PI/2, TOL, "Huge vectors precision")
}

test_nearly_identical_vectors <- function() {
  epsilon <- 1e-15
  v1 <- c(1.0, 0.0)
  v2 <- c(1.0, epsilon)
  
  signed_angle <- clock_hand_angle_between_vectors(v1, v2, ccw_reference(v1))
  assert_true(abs(signed_angle) < 1e-10, "Nearly identical vectors should have tiny angle")
}

test_nearly_opposite_vectors <- function() {
  epsilon <- 1e-15
  v1 <- c(1.0, 0.0)
  v2 <- c(-1.0, epsilon)
  
  signed_angle <- clock_hand_angle_between_vectors(v1, v2, ccw_reference(v1))
  assert_true(abs(abs(signed_angle) - PI) < 1e-10, "Nearly opposite vectors should be close to π")
}

test_mixed_positive_negative <- function() {
  v1 <- c(1.0, -2.0, 3.0)
  v2 <- c(-2.0, 1.0, -3.0)
  
  # Normalize
  v1 <- v1 / sqrt(sum(v1^2))
  v2 <- v2 / sqrt(sum(v2^2))
  
  signed_angle <- clock_hand_angle_between_vectors(v1, v2, ccw_reference(v1))
  assert_true(abs(signed_angle) >= 0 && abs(signed_angle) <= PI, "Mixed sign vectors in valid range")
}

# ==================== SHIFT VECTORS TESTS ====================

test_single_pair_shift_vectors <- function() {
  origins <- matrix(c(1.0, 0.0), nrow = 2, ncol = 1)
  targets <- matrix(c(0.0, 1.0), nrow = 2, ncol = 1)
  fields <- array(dim=c(nrow(origins), 2, ncol(origins)))
  fields[,1,] <- origins
  fields[,2,] <- targets
  signed_angles <- clock_hand_angles_for_shift_vectors(fields, rep(TRUE, dim(fields)[3]),
                                                     ccw_reference(fields[, 1, 1]))
  assert_equal_numeric(abs(signed_angles[1]), PI/2, TOL, "Single pair shift vectors")
}

test_multiple_pairs_shift_vectors <- function() {
  # Three different rotations
  origins <- matrix(c(1.0, 0.0,   # 90° CCW
                     1.0, 0.0,   # 180°
                     1.0, 0.0),  # 90° CW
                   nrow = 2, ncol = 3)
  
  targets <- matrix(c(0.0, 1.0,   # 90° CCW
                     -1.0, 0.0,   # 180°
                     0.0, -1.0),  # 90° CW
                   nrow = 2, ncol = 3)
  fields <- array(dim=c(nrow(origins), 2, ncol(origins)))
  fields[,1,] <- origins
  fields[,2,] <- targets
  
  signed_angles <- clock_hand_angles_for_shift_vectors(fields, rep(TRUE, dim(fields)[3]),
                                                     ccw_reference(fields[, 1, 1]))
  
  assert_equal_numeric(signed_angles[1], PI/2, TOL, "First rotation (90° CCW)")
  assert_equal_numeric(abs(signed_angles[2]), PI, TOL, "Second rotation (180°)")
  assert_equal_numeric(signed_angles[3], -PI/2, TOL, "Third rotation (90° CW)")
}

test_shift_vectors_with_selection_mask <- function() {
  # Four vectors, but only select 2nd and 4th
  origins <- matrix(c(1.0, 0.0,   # Not selected
                     1.0, 0.0,   # Selected (180°)
                     1.0, 0.0,   # Not selected
                     1.0, 0.0),  # Selected (45°)
                   nrow = 2, ncol = 4)
  
  targets <- matrix(c(0.0, 1.0,                              # Not selected
                     -1.0, 0.0,                             # Selected (180°)
                     0.0, -1.0,                             # Not selected
                     sqrt(2)/2, sqrt(2)/2),                 # Selected (45°)
                   nrow = 2, ncol = 4)
  fields <- array(dim=c(nrow(origins), 2, ncol(origins)))
  fields[,1,] <- origins
  fields[,2,] <- targets
  vecs_selection_mask <- c(FALSE, TRUE, FALSE, TRUE)
  
  signed_angles <- clock_hand_angles_for_shift_vectors(fields, vecs_selection_mask, ccw_reference(fields[, 1, 1]))
  
  assert_equal_numeric(abs(signed_angles[1]), PI, TOL, "Second vector (180°)")
  assert_equal_numeric(signed_angles[2], PI/4, TOL, "Fourth vector (45°)")
}

# ==================== CONSISTENCY TESTS ====================

test_consistency_between_functions <- function() {
  v1 <- c(1.0, 0.0)
  v2 <- c(0.0, 1.0)
  
  # Single function
  single_angle <- clock_hand_angle_between_vectors(v1, v2, ccw_reference(v1))
  
  # Batch function
  origins <- matrix(v1, nrow = 2, ncol = 1)
  targets <- matrix(v2, nrow = 2, ncol = 1)
  fields <- array(dim=c(nrow(origins), 2, ncol(origins)))
  fields[,1,] <- origins
  fields[,2,] <- targets
  batch_angles <- clock_hand_angles_for_shift_vectors(fields, rep(TRUE, dim(fields)[3]),
                                                     ccw_reference(fields[, 1, 1]))
  
  assert_equal_numeric(single_angle, batch_angles[1], TOL, "Single vs batch consistency")
}

test_mathematical_properties <- function() {
  v1 <- c(1.0, 2.0)
  v2 <- c(3.0, 1.0)
  
  # Normalize
  v1 <- v1 / sqrt(sum(v1^2))
  v2 <- v2 / sqrt(sum(v2^2))
  
  angle_12 <- clock_hand_angle_between_vectors(v1, v2, ccw_reference(v1))
  angle_21 <- clock_hand_angle_between_vectors(v2, v1, ccw_reference(v2))
  
  assert_equal_numeric(angle_12, -angle_21, TOL, "Anti-commutativity of signed angles")
}

run_all_tests()