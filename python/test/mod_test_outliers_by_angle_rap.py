import numpy as np
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from tensoromics_functions import (
    tox_compute_family_direction_rap,
    tox_compute_family_direction_rap_expert,
    tox_compute_angular_deviations_rap,
    tox_z_scores_by_dispersion_rap,
    tox_angle_outliers_rap,
    tox_detect_angle_outliers_pipeline_rap
)
from test_helpers import run_all_tests, assert_error

# Constants
PI = np.pi
EPS = 1.0e-12
MAX_SC = 10.0
STAT_NO_STABLE_DIRECTION = 401
STAT_NO_ANGULAR_VARIATION = 402

# ============================================================================
# Test 1: Circular family direction computation with multiple families
# ============================================================================
def test_compute_family_direction_rap_basic():
    """Test circular family direction computation with multiple families."""
    n_genes = 6
    n_families = 2

    # Family 1: Three angles clustered around 0°
    # Family 2: Three angles clustered around 90° (π/2)
    rap_angles = np.array([
        0.0,           # Gene 1, Family 1
        0.1,           # Gene 2, Family 1
        0.1,          # Gene 3, Family 1
        PI/2.0,        # Gene 4, Family 2
        PI/2.0 + 0.1,  # Gene 5, Family 2
        PI/2.0 - 0.1   # Gene 6, Family 2
    ], dtype=np.float64)

    gene_to_fam = np.array([1, 1, 1, 2, 2, 2], dtype=np.int32)

    # Call the expert version with permutation
    result = tox_compute_family_direction_rap_expert(
        rap_angles,
        gene_to_fam
    )

    family_mean_angles = result['family_mean_angles']
    family_dispersions = result['family_dispersions']
    status = result['status']
    assert np.all(status == 0), f"Status is {status}"


    # Family 1 mean should be near 0°
    np.testing.assert_almost_equal(family_mean_angles[0], 0.0, decimal=1)

    # Family 2 mean should be near π/2
    np.testing.assert_almost_equal(family_mean_angles[1], PI/2.0, decimal=2)

    # Dispersions should be positive and reasonable
    assert family_dispersions[0] > 0.0 and family_dispersions[0] < PI, f"1: {family_dispersions[0]}: {status}"
    assert family_dispersions[1] > 0.0 and family_dispersions[1] < PI, "2"



# ============================================================================
# Test 2: Circular family direction with single family
# ============================================================================
def test_compute_family_direction_rap_single_family():
    """Test circular family direction with single family."""
    n_genes = 5
    n_families = 1

    # Angles clustered around 45° (π/4)
    rap_angles = np.array([
        PI/4.0,
        PI/4.0 + 0.05,
        PI/4.0 - 0.05,
        PI/4.0 + 0.03,
        PI/4.0 - 0.03
    ], dtype=np.float64)

    gene_to_fam = np.array([1, 1, 1, 1, 1], dtype=np.int32)

    # Call the allocatable version
    result = tox_compute_family_direction_rap(rap_angles, gene_to_fam)

    family_mean_angles = result['family_mean_angles']
    family_dispersions = result['family_dispersions']
    status = result['status']
    assert np.all(status == 0), f"Status is {status}"

    np.testing.assert_almost_equal(family_mean_angles[0], PI/4.0, decimal=2)
    assert family_dispersions[0] > 0.0


# ============================================================================
# Test 3: Circular family direction with no stable direction
# ============================================================================
def test_compute_family_direction_rap_no_stable_direction():
    """Test circular family direction with no stable direction."""
    n_genes = 4
    n_families = 1

    # Angles at 0°, 90°, 180°, 270° - they cancel out
    rap_angles = np.array([
        0.0,            # 0°
        PI/2.0,         # 90°
        PI,             # 180°
        PI/4.0         # 45°
    ], dtype=np.float64)

    gene_to_fam = np.array([1, 1, 1, 1], dtype=np.int32)

    result = tox_compute_family_direction_rap_expert(
        rap_angles,
        gene_to_fam
    )

    family_dispersions = result['family_dispersions']
    status = result['status']
    assert np.all(status == 401), f"Status should be 401, is {status}"

    np.testing.assert_almost_equal(family_dispersions[0], -1.0, decimal=12)
    # Note: The status code checking might need adjustment based on actual implementation
    # Fortran would set ERR_NO_STABLE_DIRECTION in status

# ============================================================================
# Test 4: Circular family direction with minimal angular variation
# ============================================================================
def test_compute_family_direction_rap_minimal_variation():
    """Test circular family direction with minimal angular variation."""
    n_genes = 4
    n_families = 1

    # Nearly identical angles
    rap_angles = np.array([
        0.0,
        0.0001,
        PI-0.0001,
        0.0002
    ], dtype=np.float64)

    gene_to_fam = np.array([1, 1, 1, 1], dtype=np.int32)

    result = tox_compute_family_direction_rap(rap_angles, gene_to_fam)

    family_dispersions = result['family_dispersions']
    family_mean_angles = result['family_mean_angles']
    status = result['status']
    assert np.all(status == 0), f"Status is {status}"

    if family_dispersions[0] < 0.0:
        # Status check for minimal variation
        # Fortran would set ERR_NO_ANGULAR_VARIATION in status
        pass  # Minimal variation case
    else:
        # If dispersion > TAU, mean should be near 0
        np.testing.assert_almost_equal(family_mean_angles[0], 0.0, decimal=3)

# ============================================================================
# Test 5: Circular family direction with uniform distribution
# ============================================================================
def test_compute_family_direction_rap_uniform_distribution():
    """Test circular family direction with uniform distribution."""
    n_genes = 8
    n_families = 1

    # Angles uniformly spaced around the circle
    rap_angles = np.zeros(n_genes, dtype=np.float64)
    for i in range(n_genes):
        if i % 2 == 0:
            angle = PI / n_genes
        else:
            angle = PI - PI / n_genes

    gene_to_fam = np.ones(n_genes, dtype=np.int32)

    result = tox_compute_family_direction_rap(rap_angles, gene_to_fam)

    family_dispersions = result['family_dispersions']

    # Uniform distribution should have no stable direction
    np.testing.assert_almost_equal(family_dispersions[0], -1.0, decimal=12)

# ============================================================================
# Test 6: Circular family direction with all same direction
# ============================================================================
def test_compute_family_direction_rap_all_same_direction():
    """Test circular family direction with all same direction."""
    n_genes = 5
    n_families = 1

    # All angles exactly the same
    rap_angles = np.full(n_genes, PI/3.0, dtype=np.float64)  # All 60°

    gene_to_fam = np.ones(n_genes, dtype=np.int32)

    result = tox_compute_family_direction_rap_expert(
        rap_angles,
        gene_to_fam
    )

    family_mean_angles = result['family_mean_angles']
    family_dispersions = result['family_dispersions']
    status = result['status']
    assert np.all(status == 402), f"Status should be 402, is {status}"

    np.testing.assert_almost_equal(family_mean_angles[0], PI/3.0, decimal=12)
    # With identical angles, resultant length R=1, so σ = sqrt(-2*ln(1)) = 0
    # This is below TAU, so should be marked as minimal variation
    np.testing.assert_almost_equal(family_dispersions[0], -1.0, decimal=12)
    # Fortran would set ERR_NO_ANGULAR_VARIATION in status

# ============================================================================
# Test 7: Angular deviations computation for circular case
# ============================================================================
def test_compute_angular_deviations_rap_basic():
    """Test angular deviations computation for circular case."""
    n_genes = 5
    n_families = 2

    rap_angles = np.array([
        0.0,     # Gene 1, Family 1
        0.2,     # Gene 2, Family 1
        1.5,     # Gene 3, Family 2
        1.7,     # Gene 4, Family 2
        -1.0     # Gene 5, invalid
    ], dtype=np.float64)

    # Family means: Family 1 at 0.1, Family 2 at 1.6
    family_mean_angles = np.array([0.1, 1.6], dtype=np.float64)

    gene_to_fam = np.array([1, 1, 2, 2, 0], dtype=np.int32)  # Gene 5: unassigned

    angular_deviations = tox_compute_angular_deviations_rap(
        rap_angles,
        family_mean_angles,
        gene_to_fam
    )

    # Gene 1: |wrap(0.0 - 0.1)| = |wrap(-0.1)| = 0.1
    np.testing.assert_almost_equal(angular_deviations[0], 0.1, decimal=12)

    # Gene 2: |wrap(0.2 - 0.1)| = |0.1| = 0.1
    np.testing.assert_almost_equal(angular_deviations[1], 0.1, decimal=12)

    # Gene 3: |wrap(1.5 - 1.6)| = |wrap(-0.1)| = 0.1
    np.testing.assert_almost_equal(angular_deviations[2], 0.1, decimal=12)

    # Gene 4: |wrap(1.7 - 1.6)| = |0.1| = 0.1
    np.testing.assert_almost_equal(angular_deviations[3], 0.1, decimal=12)

    # Gene 5: unassigned -> -1.0
    np.testing.assert_almost_equal(angular_deviations[4], -1.0, decimal=12)

# ============================================================================
# Test 8: Angular deviations with invalid family mapping
# ============================================================================
def test_compute_angular_deviations_rap_invalid_family():
    """Test angular deviations with invalid family mapping."""
    n_genes = 4
    n_families = 2

    rap_angles = np.array([0.1, 0.2, 0.3, 0.4], dtype=np.float64)
    family_mean_angles = np.array([0.15, 0.35], dtype=np.float64)

    # Gene 3: invalid family (3 > n_families)
    gene_to_fam = np.array([1, 1, -3, 2], dtype=np.int32)

    assert_error(lambda: tox_compute_angular_deviations_rap(rap_angles, family_mean_angles, gene_to_fam), "Expected Error for invalid family")


# ============================================================================
# Test 9: Angular deviations with wrapping across π boundary
# ============================================================================
def test_compute_angular_deviations_rap_wrapping():
    """Test angular deviations with wrapping across π boundary."""
    n_genes = 3
    n_families = 1

    # Test case where difference crosses the π boundary
    rap_angles = np.array([
        3.0,     # ~171.9° (π = 3.14159)
        PI,
        0.0      # 0°
    ], dtype=np.float64)

    # Family mean near π (180°)
    family_mean_angles = np.array([PI], dtype=np.float64)

    gene_to_fam = np.array([1, 1, 1], dtype=np.int32)

    angular_deviations = tox_compute_angular_deviations_rap(
        rap_angles,
        family_mean_angles,
        gene_to_fam
    )

    # Gene 1: wrap(3.0 - π) = wrap(~-0.14159) = -0.14159, abs = ~0.14159
    np.testing.assert_almost_equal(angular_deviations[0], abs(3.0 - PI), decimal=12)

    np.testing.assert_almost_equal(angular_deviations[1], 0.0, decimal=12)

    # Gene 3: wrap(0.0 - π) = wrap(-π) = -π, abs = π
    np.testing.assert_almost_equal(angular_deviations[2], PI, decimal=12)

# ============================================================================
# Test 10: Scaling circular angles by dispersion (z-scores)
# ============================================================================
def test_z_scores_by_dispersion_rap_basic():
    """Test scaling circular angles by dispersion (z-scores)."""
    n_genes = 5
    n_families = 2

    angular_deviations = np.array([0.2, 0.4, -1.0, 0.3, 0.6], dtype=np.float64)
    family_dispersions = np.array([0.1, 0.2], dtype=np.float64)  # Family 1: σ=0.1, Family 2: σ=0.2
    gene_to_fam = np.array([1, 1, 0, 2, 2], dtype=np.int32)  # Gene 3: unassigned

    z_scores = tox_z_scores_by_dispersion_rap(
        angular_deviations,
        family_dispersions,
        gene_to_fam
    )

    # Gene 1: 0.2 / 0.1 = 2.0
    np.testing.assert_almost_equal(z_scores[0], 2.0, decimal=12)

    # Gene 2: 0.4 / 0.1 = 4.0
    np.testing.assert_almost_equal(z_scores[1], 4, decimal=12)

    # Gene 3: unassigned -> -1.0
    np.testing.assert_almost_equal(z_scores[2], -1.0, decimal=12)

    # Gene 4: 0.3 / 0.2 = 1.5
    np.testing.assert_almost_equal(z_scores[3], 1.5, decimal=12)

    # Gene 5: 0.6 / 0.2 = 3.0
    np.testing.assert_almost_equal(z_scores[4], 3.0, decimal=12)


# ============================================================================
# Test 11: Scaling circular angles with invalid inputs
# ============================================================================
def test_z_scores_by_dispersion_rap_invalid():
    """Test scaling circular angles with invalid inputs."""
    n_genes = 4
    n_families = 2

    angular_deviations = np.array([-1.0, 0.5, 0.3, 0.4], dtype=np.float64)
    family_dispersions = np.array([0.0, 0.1], dtype=np.float64)  # Family 1: invalid
    gene_to_fam = np.array([1, 1, -3, 2], dtype=np.int32)  # Gene 3: invalid family

    assert_error(lambda: tox_z_scores_by_dispersion_rap( angular_deviations, family_dispersions, gene_to_fam), "Expected Error for invalid family")


# ============================================================================
# Test 12: Scaling circular angles with zero dispersion
# ============================================================================
def test_z_scores_by_dispersion_rap_zero_dispersion():
    """Test scaling circular angles with zero dispersion."""
    n_genes = 3
    n_families = 1

    angular_deviations = np.array([0.1, 0.2, 0.3], dtype=np.float64)
    family_dispersions = np.array([0.0], dtype=np.float64)  # Zero dispersion
    gene_to_fam = np.array([1, 1, 1], dtype=np.int32)

    z_scores = tox_z_scores_by_dispersion_rap(
        angular_deviations,
        family_dispersions,
        gene_to_fam
    )

    # With zero dispersion, all should be marked as invalid
    np.testing.assert_almost_equal(z_scores[0], -1.0, decimal=12)
    np.testing.assert_almost_equal(z_scores[1], -1.0, decimal=12)
    np.testing.assert_almost_equal(z_scores[2], -1.0, decimal=12)

# ============================================================================
# Test 13: Outlier detection for circular angles
# ============================================================================
def test_angle_outliers_rap_alloc_basic():
    """Test outlier detection for circular angles."""
    n_genes = 8

    # Create scaled angles: some normal, some outliers
    z_scores = np.array([
        1.0, 1.5, 2.0, 2.5,
        3.0, 2.5, 2.0, -1.0
    ], dtype=np.float64)

    percentile = 0.8  # 80th percentile

    result = tox_angle_outliers_rap(z_scores, percentile)

    threshold = result['threshold']
    is_outlier = result['is_outlier']

    # With 7 valid scores, 80th percentile should be around 3.5
    assert threshold == 2.5, f"threshold {threshold} should be 2.5"

    # Genes with z_scores >= threshold should be outliers
    assert not is_outlier[0], "gene 1 should not be outlier"
    assert not is_outlier[1], "gene 2 should not be outlier"
    assert not is_outlier[2], "gene 3 should not be outlier"
    assert is_outlier[3], "gene 4 should be outlier"
    assert is_outlier[4], "gene 5 should be outlier"
    assert is_outlier[5], "gene 6 should be outlier"
    assert not is_outlier[6], "gene 7 should not be outlier"
    assert not is_outlier[7], "gene 8 (invalid) should not be outlier"

# ============================================================================
# Test 14: Outlier detection with no valid scaled angles
# ============================================================================
def test_angle_outliers_rap_alloc_no_valid():
    """Test outlier detection with no valid scaled angles."""
    n_genes = 3

    z_scores = np.array([-1.0, -1.0, -1.0], dtype=np.float64)  # All invalid
    percentile = 0.9

    result = tox_angle_outliers_rap(z_scores, percentile)
    assert not any(result["is_outlier"])

# ============================================================================
# Test 15: Outlier detection where all valid scores are outliers
# ============================================================================
def test_angle_outliers_rap_alloc_all_outliers():
    """Test outlier detection where all valid scores are outliers."""
    n_genes = 5

    z_scores = np.array([1.0, 1.1, 1.2, 1.3, 1.4], dtype=np.float64)
    percentile = 0.0  # 0th percentile - everything is outlier

    result = tox_angle_outliers_rap(z_scores, percentile)

    is_outlier = result['is_outlier']

    assert np.all(is_outlier), "all should be outliers with 0th percentile"
    percentile = 1.0  # 0th percentile - everything is inlier

    result = tox_angle_outliers_rap(z_scores, percentile)

    is_outlier = result['is_outlier']

    assert sum(is_outlier) == 1, "only last should be inlier with 100th percentile"

# ============================================================================
# Test 16: Complete pipeline for RAP angle-based outlier detection
# ============================================================================
def test_detect_angle_outliers_pipeline_rap_basic():
    """Test complete pipeline for RAP angle-based outlier detection."""
    n_genes = 10
    n_families = 2

    # Create RAP angles where most genes in each family point in similar directions
    # but one gene in each family is an outlier
    rap_angles = np.array([
        # Family 1 genes (1-5): clustered around 0°
        0.0,           # Gene 1: normal
        0.1,           # Gene 2: normal
        0.01,          # Gene 3: normal
        0.05,          # Gene 4: normal
        1.5,           # Gene 5: OUTLIER (~86°, far from 0°)
        # Family 2 genes (6-10): clustered around π/2 (90°)
        PI/2.0,        # Gene 6: normal
        PI/2.0 + 0.1,  # Gene 7: normal
        PI/2.0 - 0.1,  # Gene 8: normal
        PI/2.0 + 0.05, # Gene 9: normal
        PI        # Gene 10: OUTLIER (far from 90°)
    ], dtype=np.float64)

    gene_to_fam = np.array([1, 1, 1, 1, 1, 2, 2, 2, 2, 2], dtype=np.int32)
    percentile_threshold = 0.85  # 85th percentile

    result = tox_detect_angle_outliers_pipeline_rap(
        rap_angles,
        gene_to_fam,
        percentile_threshold
    )

    is_outlier = result['is_outlier']
    z_scores = result['z_scores']
    status = result['status']

    # Gene 5 and Gene 10 should be outliers (they point in very different directions)
    # With 85th percentile and 10 genes (5 per family), we expect 1 outlier per family
    assert is_outlier[4], "gene 5 should be outlier"  # Python indexing: 0-based
    assert is_outlier[9], "gene 10 should be outlier"

    # Most genes should not be outliers
    assert not np.any(is_outlier[0:4]), "genes 1-4 should not be outliers"
    assert not np.any(is_outlier[5:9]), "genes 6-9 should not be outliers"


# ============================================================================
# Test 17: Pipeline with single family
# ============================================================================
def test_detect_angle_outliers_pipeline_rap_single_family():
    """Test pipeline with single family."""
    n_genes = 8
    n_families = 1

    # Create RAP angles: 7 genes clustered around 45°, 1 gene at -90° (outlier)
    rap_angles = np.array([
        PI/4.0,             # Gene 1: normal (45°)
        PI/4.0 + 0.05,      # Gene 2: normal
        PI/4.0 - 0.05,      # Gene 3: normal
        PI/4.0 + 0.1,       # Gene 4: normal
        PI/4.0 - 0.1,       # Gene 5: normal
        PI/4.0 + 0.03,      # Gene 6: normal
        PI/4.0 - 0.03,      # Gene 7: normal
        PI/2.0             # Gene 8: OUTLIER (90°)
    ], dtype=np.float64)

    gene_to_fam = np.ones(n_genes, dtype=np.int32)
    percentile_threshold = 0.9

    result = tox_detect_angle_outliers_pipeline_rap(
        rap_angles,
        gene_to_fam,
        percentile_threshold
    )

    is_outlier = result['is_outlier']

    # With 8 genes and 90th percentile, expect 1 outlier (gene 8)
    assert is_outlier[7], "gene 8 should be outlier"
    assert not np.any(is_outlier[0:7]), "genes 1-7 should not be outliers"

# ============================================================================
# Test 18: Pipeline with no valid families
# ============================================================================
def test_detect_angle_outliers_pipeline_rap_no_valid_families():
    """Test pipeline with no valid families."""
    n_genes = 3
    n_families = 3

    rap_angles = np.array([0.1, 0.2, 0.3], dtype=np.float64)

    # Each gene in its own family -> families have only 1 gene each
    gene_to_fam = np.array([1, 2, 3], dtype=np.int32)
    percentile_threshold = 0.9

    result = tox_detect_angle_outliers_pipeline_rap( rap_angles, gene_to_fam, percentile_threshold)
    assert all(result["status"] == STAT_NO_STABLE_DIRECTION), result["status"]

# ============================================================================
# Test 19: Pipeline with invalid percentile threshold
# ============================================================================
def test_detect_angle_outliers_pipeline_rap_invalid_percentile():
    """Test pipeline with invalid percentile threshold."""
    n_genes = 5
    n_families = 1

    rap_angles = np.array([0.1, 0.2, 0.3, 0.4, 0.5], dtype=np.float64)
    gene_to_fam = np.ones(n_genes, dtype=np.int32)
    percentile_threshold = 1.01  # Invalid: must be <= 100

    assert_error(lambda: tox_detect_angle_outliers_pipeline_rap( rap_angles, gene_to_fam, percentile_threshold), "Expected Error for invalid quantile")


# ============================================================================
# Test 20: Pipeline with nearly parallel genes
# ============================================================================
def test_detect_angle_outliers_pipeline_rap_parallel_genes():
    """Test pipeline with nearly parallel genes (low angular dispersion)."""
    n_genes = 6
    n_families = 1

    # Create nearly identical angles
    rap_angles = np.array([
        0.0,
        0.0001,
        0.0003,
        0.0002,
        0.0004,
        0.0003
    ], dtype=np.float64)

    gene_to_fam = np.ones(n_genes, dtype=np.int32)
    percentile_threshold = 0.9

    result = tox_detect_angle_outliers_pipeline_rap( rap_angles, gene_to_fam, percentile_threshold)
    assert all(result["status"] == STAT_NO_ANGULAR_VARIATION)


# ============================================================================
# Test 21: Pipeline with clear orthogonal outlier
# ============================================================================
def test_detect_angle_outliers_pipeline_rap_orthogonal_outlier():
    """Test pipeline with clear orthogonal outlier."""
    n_genes = 7
    n_families = 1

    # Create 6 genes pointing around 0°, 1 gene at 90° (orthogonal)
    rap_angles = np.array([
        0.05,     # Normal
        0.0,       # Normal
        0.05,      # Normal
        0.1,       # Normal
        0.1,      # Normal
        0.03,      # Normal
        PI/2.0     # OUTLIER: orthogonal (90°)
    ], dtype=np.float64)

    gene_to_fam = np.ones(n_genes, dtype=np.int32)
    percentile_threshold = 0.85

    result = tox_detect_angle_outliers_pipeline_rap(
        rap_angles,
        gene_to_fam,
        percentile_threshold
    )

    is_outlier = result['is_outlier']

    # Gene 7 should be a clear outlier
    assert is_outlier[6], "orthogonal gene should be outlier"
    assert not np.any(is_outlier[0:6]), "normal genes should not be outliers"

# ============================================================================
# Test 22: Pipeline with mixed directions
# ============================================================================
def test_detect_angle_outliers_pipeline_rap_mixed_directions():
    """Test pipeline with mixed directions."""
    n_genes = 12
    n_families = 2

    # Family 1: angles around 0° with one outlier at 120°
    # Family 2: angles around 180° (-π) with one outlier at -60°
    rap_angles = np.array([
        # Family 1
        0.2,                      # Gene 1: normal
        0.0,                       # Gene 2: normal
        0.1,                       # Gene 3: normal
        0.06,                     # Gene 4: normal
        0.05,                      # Gene 5: normal
        2.0 * PI/3.0,              # Gene 6: OUTLIER (120°)
        # Family 2
        PI,                        # Gene 7: normal (180°)
        PI - 0.1,                  # Gene 8: normal
        PI - 0.2,                  # Gene 9: normal
        PI,                       # Gene 10: normal (180°)
        PI - 0.05,                # Gene 11: normal
        PI/3.0                    # Gene 12: OUTLIER (60°)
    ], dtype=np.float64)

    gene_to_fam = np.array([1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2], dtype=np.int32)
    percentile_threshold = 0.9

    result = tox_detect_angle_outliers_pipeline_rap(
        rap_angles,
        gene_to_fam,
        percentile_threshold
    )

    is_outlier = result['is_outlier']

    # With 90th percentile and 12 genes (6 per family), expect 1 outlier per family
    assert is_outlier[5], "gene 6 should be outlier"
    assert is_outlier[11], "gene 12 should be outlier"

    # Most genes should not be outliers
    assert not np.any(is_outlier[0:5]), "genes 1-5 should not be outliers"
    assert not np.any(is_outlier[6:11]), "genes 7-11 should not be outliers"

# ============================================================================
# Test 23: Pipeline edge cases
# ============================================================================
def test_detect_angle_outliers_pipeline_rap_edge_cases():
    """Test pipeline edge cases."""
    n_genes = 4
    n_families = 1

    # Test with angles that wrap across boundaries
    rap_angles = np.array([
        2.9,     # ~166° (close to π)
        PI,
        0.1,     # ~6°
        3.0      # ~172° (close to π)
    ], dtype=np.float64)

    gene_to_fam = np.ones(n_genes, dtype=np.int32)
    percentile_threshold = 0.75

    result = tox_detect_angle_outliers_pipeline_rap(
        rap_angles,
        gene_to_fam,
        percentile_threshold
    )


if __name__ == "__main__":
    run_all_tests(globals().values())
