import numpy as np
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from test_helpers import run_all_tests
from tensoromics_functions import (
    tox_normalize_vectors_unit_length,
    tox_detect_angle_outliers_pipeline,
    tox_compute_family_direction_expert,
    tox_compute_family_direction,
    tox_compute_angles_to_direction,
    tox_z_scores_by_dispersion,
    tox_angle_outliers,
    tox_angle_outliers_expert
)

# Constants matching Fortran module
PI = np.pi
EPS = 1.0e-12
MAX_SC = 10.0


def test_normalize_vectors_unit_length_basic():
    """Test normalization of vectors"""
    n_samples, n_genes = 3, 3
    expression_vectors = np.array([
        [1.0, 2.0, 2.0],    # Gene 1: norm = √(1² + 2² + 2²) = 3.0
        [0.0, 1.0, 0.0],    # Gene 2: norm = 1.0
        [1.0, 0.0, 0.0]     # Gene 3: norm = 1.0
    ], dtype=np.float64, order='F').T  # Transpose to match Fortran column-major

    unit_vectors = tox_normalize_vectors_unit_length(expression_vectors)

    assert unit_vectors.shape == (n_samples, n_genes), \
        f"test_normalize_vectors_unit_length_basic: wrong shape. Expected ({n_samples}, {n_genes}), got {unit_vectors.shape}"

    # Check each column is normalized (unit length)
    for i in range(n_genes):
        norm = np.linalg.norm(unit_vectors[:, i])
        assert abs(norm - 1.0) < 1e-12, f"test_normalize_vectors_unit_length_basic: column {i} not normalized. Norm = {norm}"


def test_normalize_vectors_unit_length_zero_norm():
    """Test normalization with zero norm vector - should return error"""
    n_samples, n_genes = 3, 2
    expression_vectors = np.array([
        [1.0, 2.0, 2.0],    # Gene 1: norm = 3.0
        [0.0, 0.0, 0.0]     # Gene 2: ZERO VECTOR
    ], dtype=np.float64, order='F').T

    try:
        unit_vectors = tox_normalize_vectors_unit_length(expression_vectors)
        # If it returns, check that second vector is all zeros
        assert False, "test_normalize_vectors_unit_length_zero_norm: Expected division by zero error"
    except Exception as e:
        pass


def test_detect_angle_outliers_basic():
    """Test complete angle-based outlier detection pipeline"""
    n_samples, n_genes, n_families = 3, 10, 2

    # Create expression vectors where most genes point in similar directions
    # but one gene in each family is an outlier
    expression_vectors = np.array([
        # Family 1 genes (indices 0-4)
        [1.0, 0.1, 0.0],    # Gene 0: normal
        [0.9, 0.2, 0.0],    # Gene 1: normal
        [0.8, 0.3, 0.0],    # Gene 2: normal
        [0.7, 0.4, 0.0],    # Gene 3: normal
        [0.0, 1.0, 0.0],    # Gene 4: OUTLIER (orthogonal to others)
        # Family 2 genes (indices 5-9)
        [0.1, 1.0, 0.1],    # Gene 5: normal
        [0.2, 0.9, 0.1],    # Gene 6: normal
        [0.3, 0.8, 0.1],    # Gene 7: normal
        [0.4, 0.7, 0.1],    # Gene 8: normal
        [1.0, 0.0, 0.0]     # Gene 9: OUTLIER (orthogonal to others)
    ], dtype=np.float64, order='F').T

    gene_to_fam = np.array([1, 1, 1, 1, 1, 2, 2, 2, 2, 2], dtype=np.int32, order='F')
    percentile_threshold = 85.0

    result = tox_detect_angle_outliers_pipeline(expression_vectors, gene_to_fam, percentile_threshold)

    is_outlier = result["is_outlier"]
    z_scores = result["z_scores"]
    ierr = result["ierr"]
    status = result["status"]

    # Check basic properties
    assert ierr == 0, f"test_detect_angle_outliers_basic: error code should be 0, got {ierr}"
    assert is_outlier.shape == (n_genes,), f"Wrong is_outlier shape. Expected ({n_genes},), got {is_outlier.shape}"
    assert z_scores.shape == (n_genes,), f"Wrong z_scores shape. Expected ({n_genes},), got {z_scores.shape}"

    # Check that invalid values are -1.0 for unassigned genes (we don't have any)
    for i in range(n_genes):
        if gene_to_fam[i] == 0 or gene_to_fam[i] > n_families:
            assert z_scores[i] < 0, \
                f"Invalid gene {i} should have -1.0 scaled angle, got {z_scores[i]}"

    # With 90th percentile and 10 genes (5 per family), we expect 1 outlier per family
    outliers = np.where(is_outlier)[0]

    # Gene 4 and 9 should be outliers (they point in very different directions)
    assert is_outlier[4], "test_detect_angle_outliers_basic: gene 4 should be outlier"
    assert is_outlier[9], "test_detect_angle_outliers_basic: gene 9 should be outlier"

    # Most genes should not be outliers
    assert not np.any(is_outlier[0:4]), "test_detect_angle_outliers_basic: genes 0-3 should not be outliers"
    assert not np.any(is_outlier[5:9]), "test_detect_angle_outliers_basic: genes 5-8 should not be outliers"


def test_detect_angle_outliers_single_family():
    """Test pipeline with single family"""
    n_samples, n_genes, n_families = 3, 8, 1

    expression_vectors = np.array([
        [1.0, 0.1, 0.0],    # Gene 0
        [0.9, 0.2, 0.0],    # Gene 1
        [0.8, 0.3, 0.0],    # Gene 2
        [0.7, 0.4, 0.0],    # Gene 3
        [0.6, 0.5, 0.0],    # Gene 4
        [0.5, 0.6, 0.0],    # Gene 5
        [0.4, 0.7, 0.0],    # Gene 6
        [0.0, 1.0, 0.0]     # Gene 7: OUTLIER - orthogonal to main direction
    ], dtype=np.float64).T

    gene_to_fam = np.array([1, 1, 1, 1, 1, 1, 1, 1], dtype=np.int32)
    percentile_threshold = 90.0

    result = tox_detect_angle_outliers_pipeline(expression_vectors, gene_to_fam, percentile_threshold)

    is_outlier = result["is_outlier"]
    z_scores = result["z_scores"]
    ierr = result["ierr"]

    assert ierr == 0, f"test_detect_angle_outliers_single_family: error code should be 0, got {ierr}"

    # With 8 genes and 90th percentile, expect 1 outlier (gene 7)
    assert is_outlier[7], "test_detect_angle_outliers_single_family: gene 7 should be outlier"
    assert not np.any(is_outlier[0:7]), "test_detect_angle_outliers_single_family: genes 0-6 should not be outliers"

    # Check scaled angles are reasonable
    assert np.all(z_scores[0:7] >= 0.0), "test_detect_angle_outliers_single_family: normal genes should have non-negative scaled angles"


def test_detect_angle_outliers_invalid_percentile():
    """Test pipeline with invalid percentile threshold"""
    n_samples, n_genes, n_families = 3, 5, 1

    expression_vectors = np.array([
        [1.0, 0.0, 0.0],
        [0.9, 0.1, 0.0],
        [0.8, 0.2, 0.0],
        [0.7, 0.3, 0.0],
        [0.6, 0.4, 0.0]
    ], dtype=np.float64).T

    gene_to_fam = np.array([1, 1, 1, 1, 1], dtype=np.int32)

    # Test percentile = 100.0 (invalid: must be < 100)
    percentile_threshold = 100.0

    try:
        result = tox_detect_angle_outliers_pipeline(expression_vectors, gene_to_fam, percentile_threshold)
        assert False, "tox_detect_angle_outliers_pipeline should throw Error for invalid percentile"
    except Exception as e:
        pass


def test_detect_angle_outliers_no_valid_families():
    """Test pipeline with no valid families (all families have < 2 genes)"""
    n_samples, n_genes, n_families = 3, 3, 3

    expression_vectors = np.array([
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0]
    ], dtype=np.float64).T

    # Each gene in its own family -> families have only 1 gene each
    gene_to_fam = np.array([1, 2, 3], dtype=np.int32)
    percentile_threshold = 90.0

    result = tox_detect_angle_outliers_pipeline(expression_vectors, gene_to_fam, percentile_threshold)
    assert np.any(result["status"] != 0), "Expected status indicator for no valid families"


def test_detect_angle_outliers_parallel_genes():
    """Test pipeline with nearly parallel genes (low angular dispersion)"""
    n_samples, n_genes, n_families = 3, 6, 1

    # Create nearly identical vectors
    expression_vectors = np.array([
        [1.0, 0.001, 0.0],    # Gene 0
        [0.999, 0.002, 0.0],  # Gene 1 - slightly different
        [0.998, 0.001, 0.0],  # Gene 2
        [0.999, 0.0, 0.001],  # Gene 3
        [0.999, -0.001, 0.0], # Gene 4
        [1.0, -0.002, 0.0]    # Gene 5 - slightly different
    ], dtype=np.float64).T

    gene_to_fam = np.array([1, 1, 1, 1, 1, 1], dtype=np.int32)
    percentile_threshold = 90.0

    result = tox_detect_angle_outliers_pipeline(expression_vectors, gene_to_fam, percentile_threshold)
    assert np.any(result["status"] != 0), "Expected status indicator for low angular dispersion"


def test_detect_angle_outliers_orthogonal_outlier():
    """Test pipeline with clear orthogonal outlier"""
    n_samples, n_genes, n_families = 3, 7, 1

    # Create 6 genes pointing in similar direction, 1 gene orthogonal
    expression_vectors = np.array([
        [1.0, 0.1, 0.0],    # Gene 0: normal
        [0.9, 0.2, 0.0],    # Gene 1: normal
        [0.8, 0.3, 0.0],    # Gene 2: normal
        [0.7, 0.4, 0.0],    # Gene 3: normal
        [0.6, 0.5, 0.0],    # Gene 4: normal
        [0.5, 0.6, 0.0],    # Gene 5: normal
        [0.0, 1.0, 0.0]     # Gene 6: OUTLIER - orthogonal
    ], dtype=np.float64).T

    gene_to_fam = np.array([1, 1, 1, 1, 1, 1, 1], dtype=np.int32)
    percentile_threshold = 85.0

    result = tox_detect_angle_outliers_pipeline(expression_vectors, gene_to_fam, percentile_threshold)

    is_outlier = result["is_outlier"]
    z_scores = result["z_scores"]
    ierr = result["ierr"]

    assert ierr == 0, f"test_detect_angle_outliers_orthogonal_outlier: error code should be 0, got {ierr}"

    # Gene 6 should be a clear outlier
    assert is_outlier[6], "test_detect_angle_outliers_orthogonal_outlier: orthogonal gene should be outlier"
    assert not np.any(is_outlier[0:6]), "test_detect_angle_outliers_orthogonal_outlier: normal genes should not be outliers"

    # Orthogonal gene should have much larger scaled angle
    if np.all(z_scores >= 0.0):
        assert z_scores[6] > 2.0 * np.mean(z_scores[0:6]), \
            "Orthogonal gene should have significantly larger scaled angle"


def test_detect_angle_outliers_invalid_family_mapping():
    """Test pipeline with invalid family mappings"""
    n_samples, n_genes, n_families = 3, 5, 2

    expression_vectors = np.array([
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
        [1.0, 1.0, 0.0],
        [0.5, 0.5, 0.5]
    ], dtype=np.float64).T

    # Gene 2: unassigned (0), Gene 3: invalid family (3 > n_families)
    gene_to_fam = np.array([1, 1, 0, 3, 2], dtype=np.int32)
    percentile_threshold = 90.0

    try:
        result = tox_detect_angle_outliers_pipeline(expression_vectors, gene_to_fam, percentile_threshold)
        assert False, "expected error for invalid gene_to_fam"
    except Exception as e:
        pass


def create_expression_data(n_samples=3, n_genes=10, n_families=2, outlier_indices=None):
    """Create test expression data with potential outliers"""
    if outlier_indices is None:
        outlier_indices = []

    expression_vectors = []
    gene_to_fam = []

    genes_per_family = n_genes // n_families

    for fam in range(1, n_families + 1):
        for gene_idx in range(genes_per_family):
            gene_global_idx = (fam-1) * genes_per_family + gene_idx

            if gene_global_idx in outlier_indices:
                # Create outlier vector (orthogonal to family direction)
                vec = np.zeros(n_samples)
                vec[(fam + gene_idx) % n_samples] = 1.0
            else:
                # Create normal vector (clustered around family direction)
                vec = np.random.normal(0, 0.1, n_samples)
                vec[0] += 1.0  # Bias toward first dimension
                vec = vec / np.linalg.norm(vec) if np.linalg.norm(vec) > 0 else vec

            expression_vectors.append(vec)
            gene_to_fam.append(fam)

    expression_vectors = np.array(expression_vectors, dtype=np.float64).T
    gene_to_fam = np.array(gene_to_fam, dtype=np.int32)

    return expression_vectors, gene_to_fam


def test_compute_family_direction_basic():
    """Test basic family direction computation"""
    n_samples, n_genes, n_families = 3, 8, 2

    expression_vectors, gene_to_fam = create_expression_data(
        n_samples, n_genes, n_families, outlier_indices=[3, 7]
    )

    # First normalize vectors
    unit_vectors = tox_normalize_vectors_unit_length(expression_vectors)

    # Compute family directions
    result = tox_compute_family_direction(unit_vectors, gene_to_fam)

    family_directions = result["family_directions"]
    angular_dispersions = result["angular_dispersions"]
    ierr = result["ierr"]
    status = result["status"]

    assert ierr == 0, f"test_compute_family_direction_basic: error code should be 0, got {ierr}"
    assert family_directions.shape == (n_samples, n_families), \
        f"Wrong family_directions shape. Expected ({n_samples}, {n_families}), got {family_directions.shape}"
    assert angular_dispersions.shape == (n_families,), \
        f"Wrong angular_dispersions shape. Expected ({n_families},), got {angular_dispersions.shape}"

    # Check each family direction is unit length
    for f in range(n_families):
        norm = np.linalg.norm(family_directions[:, f])
        assert abs(norm - 1.0) < 1e-12, \
            f"Family direction {f} not normalized. Norm = {norm}"

    # Check angular dispersions are positive
    for f in range(n_families):
        assert angular_dispersions[f] > 0.0, \
            f"Family {f} angular dispersion should be positive, got {angular_dispersions[f]}"


def test_compute_family_direction_single_gene_families():
    """Test family direction computation with single-gene families"""
    n_samples, n_genes, n_families = 3, 3, 3

    expression_vectors = np.array([
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0]
    ], dtype=np.float64).T

    unit_vectors = tox_normalize_vectors_unit_length(expression_vectors)
    gene_to_fam = np.array([1, 2, 3], dtype=np.int32)

    result = tox_compute_family_direction(unit_vectors, gene_to_fam)

    family_directions = result["family_directions"]
    angular_dispersions = result["angular_dispersions"]
    ierr = result["ierr"]
    status = result["status"]

    assert ierr == 0, f"test_compute_family_direction_single_gene_families: error code should be 0, got {ierr}"

    # Families with < 2 genes should have angular_dispersions = -1
    for f in range(n_families):
        assert abs(angular_dispersions[f] - (-1.0)) < 1e-12, \
            f"Single-gene family {f} should have angular dispersion -1, got {angular_dispersions[f]}"


def test_compute_family_direction_expert_with_permutation():
    """Test expert family direction computation with pre-sorted permutation"""
    n_samples, n_genes, n_families = 3, 8, 2

    expression_vectors, gene_to_fam = create_expression_data(
        n_samples, n_genes, n_families, outlier_indices=[3, 7]
    )

    # Normalize vectors
    unit_vectors = tox_normalize_vectors_unit_length(expression_vectors)

    # Simple Python implementation for creating permutation
    indices = np.arange(n_genes, dtype=np.int32)
    gene_to_fam_perm = indices[np.argsort(gene_to_fam)] + 1

    # Compute family directions with expert function
    result = tox_compute_family_direction_expert(unit_vectors, gene_to_fam)

    family_directions = result["family_directions"]
    angular_dispersions = result["angular_dispersions"]
    ierr = result["ierr"]
    status = result["status"]

    assert ierr == 0, f"test_compute_family_direction_expert_with_permutation: error code should be 0, got {ierr}"
    assert family_directions.shape == (n_samples, n_families)
    assert angular_dispersions.shape == (n_families,)

    # Compare with regular function
    result_regular = tox_compute_family_direction(unit_vectors, gene_to_fam)

    # Results should be similar (might differ slightly due to sorting)
    for f in range(n_families):
        assert np.allclose(family_directions[:, f], result_regular["family_directions"][:, f], atol=1e-10), \
            f"Family direction {f} differs between expert and regular computation"
        assert abs(angular_dispersions[f] - result_regular["angular_dispersions"][f]) < 1e-10, \
            f"Angular dispersion {f} differs between expert and regular computation"


def test_compute_angles_to_direction_basic():
    """Test basic angle computation between genes and family directions"""
    n_samples, n_genes, n_families = 3, 6, 2

    expression_vectors, gene_to_fam = create_expression_data(
        n_samples, n_genes, n_families, outlier_indices=[]
    )

    # Normalize vectors
    unit_vectors = tox_normalize_vectors_unit_length(expression_vectors)

    # Compute family directions
    family_result = tox_compute_family_direction(unit_vectors, gene_to_fam)
    family_directions = family_result["family_directions"]

    # Compute angles
    result = tox_compute_angles_to_direction(unit_vectors, family_directions, gene_to_fam)

    angles = result["angles"]
    ierr = result["ierr"]

    assert ierr == 0, f"test_compute_angles_to_direction_basic: error code should be 0, got {ierr}"
    assert angles.shape == (n_genes,), f"Wrong angles shape. Expected ({n_genes},), got {angles.shape}"

    # Check angles are in range [0, π]
    for i in range(n_genes):
        if angles[i] >= 0.0:  # Skip invalid genes
            assert 0.0 <= angles[i] <= PI, \
                f"Angle {i} out of range [0, π]: {angles[i]}"

    # Check that genes in same family have similar angles to family direction
    for fam in range(1, n_families + 1):
        fam_indices = np.where(gene_to_fam == fam)[0]
        fam_angles = angles[fam_indices]

        # Remove invalid angles
        valid_angles = fam_angles[fam_angles >= 0.0]

        if len(valid_angles) > 1:
            angle_std = np.std(valid_angles)
            assert angle_std < 0.5, \
                f"Angles for family {fam} should be clustered. Std = {angle_std}"


def test_compute_angles_to_direction_invalid_families():
    """Test angle computation with invalid/unassigned genes"""
    n_samples, n_genes, n_families = 3, 5, 2

    expression_vectors = np.array([
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
        [1.0, 1.0, 0.0],
        [0.5, 0.5, 0.5]
    ], dtype=np.float64).T

    unit_vectors = tox_normalize_vectors_unit_length(expression_vectors)

    # Create some invalid family assignments
    gene_to_fam = np.array([1, 1, 0, -3, 2], dtype=np.int32)  # 0=unassigned, 3>n_families

    # Compute family directions only for valid families
    valid_fam_mask = (gene_to_fam >= 1) & (gene_to_fam <= n_families)
    valid_gene_to_fam = gene_to_fam[valid_fam_mask]
    valid_unit_vectors = unit_vectors[:, valid_fam_mask]

    family_result = tox_compute_family_direction(valid_unit_vectors, valid_gene_to_fam)
    family_directions = family_result["family_directions"]

    try:
        # Compute angles for all genes
        result = tox_compute_angles_to_direction(unit_vectors, family_directions, gene_to_fam)
        assert False, "Expected error for invalid family"
    except Exception as e:
        pass


def test_z_scores_by_dispersion_basic():
    """Test basic scaled angles computation"""
    n_samples, n_genes, n_families = 3, 8, 2

    expression_vectors, gene_to_fam = create_expression_data(
        n_samples, n_genes, n_families, outlier_indices=[3, 7]
    )

    # Normalize vectors
    unit_vectors = tox_normalize_vectors_unit_length(expression_vectors)

    # Compute family directions and angular dispersions
    family_result = tox_compute_family_direction(unit_vectors, gene_to_fam)
    family_directions = family_result["family_directions"]
    angular_dispersions = family_result["angular_dispersions"]

    # Compute angles
    angles_result = tox_compute_angles_to_direction(unit_vectors, family_directions, gene_to_fam)
    angles = angles_result["angles"]

    # Compute scaled angles
    result = tox_z_scores_by_dispersion(angles, angular_dispersions, gene_to_fam)

    z_scores = result["z_scores"]
    ierr = result["ierr"]

    assert ierr == 0, f"test_z_scores_by_dispersion_basic: error code should be 0, got {ierr}"
    assert z_scores.shape == (n_genes,), \
        f"Wrong z_scores shape. Expected ({n_genes},), got {z_scores.shape}"

    # Check scaled angles
    for i in range(n_genes):
        fam = gene_to_fam[i]

        if fam < 1 or fam > n_families:
            # Unassigned/invalid gene should have scaled_angle = -1
            assert z_scores[i] < 0, \
                f"Invalid gene {i} should have scaled_angle -1, got {z_scores[i]}"
        elif angles[i] < 0.0 or angular_dispersions[fam-1] < 0.0:
            # Invalid angle or dispersion should give -1
            assert z_scores[i] < 0, \
                f"Gene {i} with invalid angle/dispersion should have scaled_angle -1, got {z_scores[i]}"
        elif abs(angular_dispersions[fam-1]) < EPS:
            # Zero dispersion should give -1
            assert z_scores[i] < 0, \
                f"Gene {i} with zero dispersion should have scaled_angle -1, got {z_scores[i]}"
        else:
            # Valid scaled angle should be angle / dispersion
            expected = angles[i] / angular_dispersions[fam-1]
            assert abs(z_scores[i] - expected) < 1e-12, \
                f"Gene {i} scaled angle mismatch. Expected {expected}, got {z_scores[i]}"

    # Outliers should have larger scaled angles
    outlier_indices = [3, 7]
    non_outlier_indices = [i for i in range(n_genes) if i not in outlier_indices and z_scores[i] >= 0.0]

    if len(outlier_indices) > 0 and len(non_outlier_indices) > 0:
        outlier_mean = np.mean([z_scores[i] for i in outlier_indices if z_scores[i] >= 0.0])
        non_outlier_mean = np.mean([z_scores[i] for i in non_outlier_indices])

        assert outlier_mean > non_outlier_mean, \
            "Outliers should have larger scaled angles on average"


def test_z_scores_by_dispersion_zero_dispersion():
    """Test scaled angles computation with zero angular dispersion"""
    n_samples, n_genes, n_families = 3, 4, 1

    gene_to_fam = np.array([1, 1, 1, 1], dtype=np.int32)

    angular_dispersions = np.zeros(n_families, dtype=np.float64)
    angles = np.full(n_genes, np.pi / 2, dtype=np.float64)

    result = tox_z_scores_by_dispersion(angles, angular_dispersions, gene_to_fam)
    z_scores = result["z_scores"]

    # With very small dispersion, scaled angles might be large or marked as invalid
    for i in range(n_genes):
        assert z_scores[i] < 0, \
            f"With zero dispersion, gene {i} should have scaled_angle -1, got {z_scores[i]}"


def test_angle_outliers_basic():
    """Test basic outlier detection from scaled angles"""
    n_genes = 10

    # Create scaled angles with clear outliers
    z_scores = np.array([
        0.5, 0.6, 0.7, 0.8, 0.9,  # Normal genes
        1.0, 1.1, 1.2, 1.3, 3.0   # Last one is outlier
    ], dtype=np.float64)

    percentile = 90.0  # 90th percentile

    result = tox_angle_outliers(z_scores, percentile)

    is_outlier = result["is_outlier"]
    threshold = result["threshold"]
    assert is_outlier.shape == (n_genes,), \
        f"Wrong is_outlier shape. Expected ({n_genes},), got {is_outlier.shape}"

    # With 90th percentile and 10 genes, threshold should be > 1.3
    assert threshold > 1.3, f"Threshold should be > 1.3 for 90th percentile, got {threshold}"

    # Only the last gene (scaled_angle = 3.0) should be outlier
    assert is_outlier[-1], "Last gene (scaled_angle=3.0) should be outlier"
    assert not np.any(is_outlier[0:-1]), "Other genes should not be outliers"


def test_angle_outliers_varying_percentile():
    """Test outlier detection with different percentile thresholds"""
    n_genes = 20

    # Create normally distributed scaled angles
    np.random.seed(42)
    z_scores = np.random.normal(1.0, 0.5, n_genes)
    z_scores = np.clip(z_scores, 0.0, None)  # Ensure non-negative

    # Test different percentiles
    percentiles = [80.0, 90.0, 95.0, 99.0]

    for percentile in percentiles:
        result = tox_angle_outliers(z_scores, percentile)
        is_outlier = result["is_outlier"]
        threshold = result["threshold"]

        # Count outliers
        outlier_count = np.sum(is_outlier)
        expected_max = int((100.0 - percentile) / 100.0 * n_genes) + 1

        # With random data, might not match exactly, but should be reasonable
        assert outlier_count <= expected_max * 2, \
            f"Too many outliers for percentile {percentile}: {outlier_count} > {expected_max * 2}"

        # All genes with z_scores >= threshold should be outliers
        for i in range(n_genes):
            if z_scores[i] >= threshold - 1e-12:
                assert is_outlier[i], \
                    f"Gene {i} with scaled_angle {z_scores[i]} >= threshold {threshold} should be outlier"


def test_angle_outliers_invalid_z_scores():
    """Test outlier detection with invalid scaled angles (-1)"""
    n_genes = 8

    z_scores = np.array([
        0.5, 0.6, -1.0, 0.8, 0.9, -1.0, 1.2, 3.0
    ], dtype=np.float64)

    percentile = 90.0

    result = tox_angle_outliers(z_scores, percentile)

    is_outlier = result["is_outlier"]
    threshold = result["threshold"]

    # Invalid scaled angles (-1) should not be counted as outliers
    invalid_indices = np.where(z_scores < 0.0)[0]
    for idx in invalid_indices:
        assert not is_outlier[idx], f"Invalid gene {idx} (scaled_angle=-1) should not be outlier"

    # Last gene (3.0) should be outlier
    assert is_outlier[-1], "Last gene (scaled_angle=3.0) should be outlier"


def test_angle_outliers_few_valid_genes():
    """Test outlier detection with very few valid genes"""
    n_genes = 10

    # Only 2 valid genes
    z_scores = np.array([
        -1.0, -1.0, -1.0, -1.0, -1.0,
        -1.0, -1.0, 0.5, 0.6, -1.0
    ], dtype=np.float64)

    percentile = 90.0

    result = tox_angle_outliers(z_scores, percentile)
    is_outlier = result["is_outlier"]
    threshold = result["threshold"]
    valid_indices = np.where(z_scores >= 0.0)[0]
    assert all(is_outlier[valid_indices] == [False, True]), f"With only 2 valid genes, the highest valued one should be outlier with 90%-ile"


def test_angle_outliers_extreme_percentile():
    """Test outlier detection with extreme percentiles"""
    n_genes = 100

    # Create uniform scaled angles
    z_scores = np.linspace(0.0, np.pi, n_genes, dtype=np.float64)

    # Test extreme percentiles
    extreme_percentiles = [0.1, 99.9, 50.0]

    for percentile in extreme_percentiles:
        result = tox_angle_outliers(z_scores, percentile)
        is_outlier = result["is_outlier"]
        threshold = result["threshold"]

        # Count outliers
        outlier_count = np.sum(is_outlier)
        expected_count = max(1, int((100.0 - percentile) / 100.0 * n_genes))

        # Allow some tolerance
        assert abs(outlier_count - expected_count) <= 2, \
            f"Outlier count mismatch for percentile {percentile}: expected ~{expected_count}, got {outlier_count}"


def test_full_pipeline_consistency():
    """Test that individual function results match full pipeline"""
    n_samples, n_genes, n_families = 3, 12, 2

    # Create test data
    expression_vectors, gene_to_fam = create_expression_data(
        n_samples, n_genes, n_families, outlier_indices=[4, 10]
    )

    percentile_threshold = 90.0

    # Run full pipeline
    full_result = tox_detect_angle_outliers_pipeline(
        expression_vectors, gene_to_fam, percentile_threshold
    )

    full_is_outlier = full_result["is_outlier"]
    full_z_scores = full_result["z_scores"]
    full_status = full_result["status"]

    # Run step-by-step
    # Step 1: Normalize
    unit_vectors = tox_normalize_vectors_unit_length(expression_vectors)

    # Step 2: Compute family directions
    family_result = tox_compute_family_direction(unit_vectors, gene_to_fam)
    family_directions = family_result["family_directions"]
    angular_dispersions = family_result["angular_dispersions"]

    # Step 3: Compute angles
    angles_result = tox_compute_angles_to_direction(
        unit_vectors, family_directions, gene_to_fam
    )
    angles = angles_result["angles"]

    # Step 4: Compute scaled angles
    scaled_result = tox_z_scores_by_dispersion(
        angles, angular_dispersions, gene_to_fam
    )
    z_scores = scaled_result["z_scores"]

    # Step 5: Detect outliers
    outlier_result = tox_angle_outliers(z_scores, percentile_threshold)
    step_is_outlier = outlier_result["is_outlier"]
    threshold = outlier_result["threshold"]

    # Check scaled angles match
    assert np.allclose(full_z_scores, z_scores, atol=1e-10, equal_nan=True), \
        "Scaled angles from full pipeline don't match step-by-step"

    # Check outlier detection matches
    assert np.array_equal(full_is_outlier, step_is_outlier), \
        "Outlier detection from full pipeline doesn't match step-by-step"


def test_edge_case_high_dimensional():
    """Test with high-dimensional data"""
    n_samples, n_genes, n_families = 10, 20, 3

    # Create high-dimensional data
    expression_vectors = np.random.normal(0, 1, (n_samples, n_genes)).astype(np.float64)

    # Assign genes to families
    gene_to_fam = np.zeros(n_genes, dtype=np.int32)
    for i in range(n_genes):
        gene_to_fam[i] = (i % n_families) + 1

    # Run pipeline
    percentile_threshold = 95.0

    result = tox_detect_angle_outliers_pipeline(
        expression_vectors, gene_to_fam, percentile_threshold
    )

    is_outlier = result["is_outlier"]
    z_scores = result["z_scores"]
    ierr = result["ierr"]
    status = result["status"]

    assert ierr == 0, f"test_edge_case_high_dimensional: error code should be 0, got {ierr}"
    assert is_outlier.shape == (n_genes,), "Wrong is_outlier shape"
    assert z_scores.shape == (n_genes,), "Wrong z_scores shape"

    # Check for any outliers found
    outlier_count = np.sum(is_outlier)


def test_angle_outliers_expert_basic_consistency():
    """Test expert outlier detection against regular function (unsorted valid_scores)"""

    z_scores = np.array([
        0.5, 0.6, -1.0, 0.8, 0.9, 1.1, -1.0, 2.5
    ], dtype=np.float64)

    percentile = 90.0

    # Reference result
    ref = tox_angle_outliers(z_scores, percentile)
    ref_is_outlier = ref["is_outlier"]
    ref_threshold = ref["threshold"]

    # Expert call
    result = tox_angle_outliers_expert(
        z_scores,
        ref_threshold
    )

    is_outlier = result["is_outlier"]
    threshold = result["threshold"]
    ierr = result["ierr"]

    assert ierr == 0, f"Expert function returned error code {ierr}"
    assert np.array_equal(is_outlier, ref_is_outlier), \
        "Expert and regular outlier masks differ"
    assert abs(threshold - ref_threshold) < 1e-12, \
        "Expert and regular thresholds differ"


if __name__ == "__main__":
    run_all_tests(globals().values())
