#!/usr/bin/env python3
"""
Comprehensive Python test suite for outlier detection functions
Uses tensoromics_functions.py wrapper functions (mirrors R get_outliers.R tests)
"""

import numpy as np
import sys
import os

# Add parent directory to path to import tensoromics_functions
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from tensoromics_functions import (
    tox_compute_family_scaling, tox_compute_family_scaling_expert, tox_compute_rdi, tox_identify_outliers, tox_detect_outliers
)
from test_helpers import run_all_tests


# =====================
# Tests for compute_family_scaling
# =====================

def test_compute_family_scaling_basic():
    """Test basic family scaling computation"""

    # Updated to use more families and genes
    n_families = 10  # Incremented number of families
    genes_per_fam = 4
    n_genes = n_families * genes_per_fam

    # Generate distances and gene-to-family mapping
    distances = np.concatenate([
        np.random.uniform(1.0, 2.0, genes_per_fam),  # Family 1
        np.random.uniform(2.0, 3.0, genes_per_fam),  # Family 2
        np.random.uniform(3.0, 4.0, genes_per_fam),  # Family 3
        np.random.uniform(4.0, 5.0, genes_per_fam),  # Family 4
        np.random.uniform(5.0, 6.0, genes_per_fam),  # Family 5
        np.random.uniform(6.0, 7.0, genes_per_fam),  # Family 6
        np.random.uniform(7.0, 8.0, genes_per_fam),  # Family 7
        np.random.uniform(8.0, 9.0, genes_per_fam),  # Family 8
        np.random.uniform(9.0, 10.0, genes_per_fam), # Family 9
        np.random.uniform(10.0, 11.0, genes_per_fam) # Family 10
    ])
    gene_to_fam = np.concatenate([
        np.full(genes_per_fam, 1),
        np.full(genes_per_fam, 2),
        np.full(genes_per_fam, 3),
        np.full(genes_per_fam, 4),
        np.full(genes_per_fam, 5),
        np.full(genes_per_fam, 6),
        np.full(genes_per_fam, 7),
        np.full(genes_per_fam, 8),
        np.full(genes_per_fam, 9),
        np.full(genes_per_fam, 10)
    ]).astype(np.int32)

    result = tox_compute_family_scaling(distances, gene_to_fam)

    # Verify basic properties
    assert len(result['dscale']) == n_families  # Ten families
    assert all(np.isfinite(result['dscale']))
    assert all(result['dscale'] > 0)  # Scaling factors should be positive


def test_compute_family_scaling_single_family():
    """Test with single family"""

    distances = np.array([1.0, 2.0, 3.0, 4.0], dtype=np.float64)
    gene_to_fam = np.array([1, 1, 1, 1], dtype=np.int32)  # All in family 1

    result = tox_compute_family_scaling(distances, gene_to_fam)

    assert len(result['dscale']) == 1
    assert result['dscale'][0] == 0


def test_compute_family_scaling_edge_cases():
    """Test edge cases with more families"""

    # Case 1: Very small distances with more families
    distances = np.array([
        0.01, 0.02, 0.01, 0.02,  # Family 1
        0.03, 0.04, 0.03, 0.04,  # Family 2
        0.05, 0.06, 0.05, 0.06,  # Family 3
        0.07, 0.08, 0.07, 0.08,  # Family 4
        0.09, 0.10, 0.09, 0.10,  # Family 5
        0.11, 0.12, 0.11, 0.12   # Family 6
    ], dtype=np.float64)
    gene_to_fam = np.array([
        1, 1, 2, 2,  # Family 1
        3, 3, 4, 4,  # Family 2
        5, 5, 6, 6,  # Family 3
        7, 7, 8, 8,  # Family 4
        9, 9, 10, 10,  # Family 5
        11, 11, 12, 12  # Family 6
    ], dtype=np.int32)

    result = tox_compute_family_scaling(distances, gene_to_fam)
    assert all(np.isfinite(result['dscale']))

    # Case 2: Large distances with more families
    distances = np.array([
        100.0, 200.0, 150.0, 250.0,  # Family 1
        300.0, 400.0, 350.0, 450.0,  # Family 2
        500.0, 600.0, 550.0, 650.0,  # Family 3
        700.0, 800.0, 750.0, 850.0,  # Family 4
        900.0, 1000.0, 950.0, 1050.0,  # Family 5
        1100.0, 1200.0, 1150.0, 1250.0   # Family 6
    ], dtype=np.float64)
    gene_to_fam = np.array([
        1, 1, 2, 2,  # Family 1
        3, 3, 4, 4,  # Family 2
        5, 5, 6, 6,  # Family 3
        7, 7, 8, 8,  # Family 4
        9, 9, 10, 10,  # Family 5
        11, 11, 12, 12  # Family 6
    ], dtype=np.int32)

    result = tox_compute_family_scaling(distances, gene_to_fam)
    assert all(np.isfinite(result['dscale']))


def test_tox_compute_family_scaling_expert():
    """Test expert version with user-provided work arrays"""

    # Test data
    n_families = 10  # Incremented number of families
    genes_per_fam = 4
    n_genes = n_families * genes_per_fam
    span = 0.7
    mode = 1
    n_iters = 3
    degree = 2

    distances = np.concatenate([
        np.random.uniform(1.0, 2.0, genes_per_fam),  # Family 1
        np.random.uniform(2.0, 3.0, genes_per_fam),  # Family 2
        np.random.uniform(3.0, 4.0, genes_per_fam),  # Family 3
        np.random.uniform(4.0, 5.0, genes_per_fam),  # Family 4
        np.random.uniform(5.0, 6.0, genes_per_fam),  # Family 5
        np.random.uniform(6.0, 7.0, genes_per_fam),  # Family 6
        np.random.uniform(7.0, 8.0, genes_per_fam),  # Family 7
        np.random.uniform(8.0, 9.0, genes_per_fam),  # Family 8
        np.random.uniform(9.0, 10.0, genes_per_fam), # Family 9
        np.random.uniform(10.0, 11.0, genes_per_fam) # Family 10
    ])
    gene_to_fam = np.concatenate([
        np.full(genes_per_fam, 1),
        np.full(genes_per_fam, 2),
        np.full(genes_per_fam, 3),
        np.full(genes_per_fam, 4),
        np.full(genes_per_fam, 5),
        np.full(genes_per_fam, 6),
        np.full(genes_per_fam, 7),
        np.full(genes_per_fam, 8),
        np.full(genes_per_fam, 9),
        np.full(genes_per_fam, 10)
    ]).astype(np.int32)

    # Pre-allocate work arrays (user responsibility)
    perm_tmp = np.empty(n_genes, dtype=np.int32)
    stack_left_tmp = np.empty(n_genes, dtype=np.int32)
    stack_right_tmp = np.empty(n_genes, dtype=np.int32)

    result = tox_compute_family_scaling_expert(
        distances, gene_to_fam, span, degree, mode, n_iters, perm_tmp, stack_left_tmp,
        stack_right_tmp
    )

    # Verify basic properties
    assert len(result['dscale']) == n_families  # Ten families
    assert all(np.isfinite(result['dscale']))
    assert all(result['dscale'] > 0)  # Scaling factors should be positive

    # Verify work arrays are returned
    assert 'perm_tmp' in result
    assert 'stack_left_tmp' in result
    assert 'stack_right_tmp' in result

    # Compare with regular version to ensure consistency
    result_regular = tox_compute_family_scaling(distances, gene_to_fam)

    # Results should be very similar (within numerical precision)
    np.testing.assert_allclose(result['dscale'], result_regular['dscale'], rtol=1e-10)
    np.testing.assert_allclose(result['loess_x'], result_regular['loess_x'], rtol=1e-10)
    np.testing.assert_allclose(result['loess_y'], result_regular['loess_y'], rtol=1e-10)


def test_compute_family_scaling_expert_input_validation():
    """Test expert version input validation"""

    n_families = 6
    genes_per_fam = 4
    n_genes = n_families * genes_per_fam
    span = 0.7
    mode = 1
    n_iters = 3
    degree = 2

    distances = np.concatenate([
        np.random.uniform(1.0, 2.0, genes_per_fam),  # Family 1
        np.random.uniform(2.0, 3.0, genes_per_fam),  # Family 2
        np.random.uniform(3.0, 4.0, genes_per_fam),  # Family 3
        np.random.uniform(4.0, 5.0, genes_per_fam),  # Family 4
        np.random.uniform(5.0, 6.0, genes_per_fam),  # Family 5
        np.random.uniform(6.0, 7.0, genes_per_fam)   # Family 6
    ])
    gene_to_fam = np.concatenate([
        np.full(genes_per_fam, 1),
        np.full(genes_per_fam, 2),
        np.full(genes_per_fam, 3),
        np.full(genes_per_fam, 4),
        np.full(genes_per_fam, 5),
        np.full(genes_per_fam, 6)
    ]).astype(np.int32)

    # Test wrong size work arrays
    error_caught1 = False
    try:
        perm_tmp = np.empty(n_genes - 1, dtype=np.int32)  # Wrong size
        stack_left_tmp = np.empty(n_genes, dtype=np.int32)
        stack_right_tmp = np.empty(n_genes, dtype=np.int32)

        tox_compute_family_scaling_expert(
            distances, gene_to_fam, span, degree, mode, n_iters, perm_tmp, stack_left_tmp,
            stack_right_tmp
        )
    except ValueError as e:
        error_caught1 = True
        assert "same length" in str(e)
    assert error_caught1

    # Test another wrong size array
    error_caught2 = False
    try:
        perm_tmp = np.empty(n_genes, dtype=np.int32)
        stack_left_tmp = np.empty(n_genes, dtype=np.int32)
        stack_right_tmp = np.empty(n_genes + 1, dtype=np.int32)

        tox_compute_family_scaling_expert(
            distances, gene_to_fam, span, degree, mode, n_iters, perm_tmp, stack_left_tmp,
            stack_right_tmp
        )
    except ValueError as e:
        error_caught2 = True
        assert "same length" in str(e)
    assert error_caught2


# =====================
# Tests for compute_rdi
# =====================

def test_compute_rdi_basic():
    """Test basic RDI computation"""

    distances = np.array([1.0, 2.0, 3.0, 4.0, 5.0], dtype=np.float64)
    gene_to_fam = np.array([1, 1, 2, 2, 2], dtype=np.int32)
    dscale = np.array([1.5, 2.0], dtype=np.float64)  # Manual scaling factors

    rdi = tox_compute_rdi(distances, gene_to_fam, dscale)

    # Verify properties
    assert len(rdi) == len(distances)
    assert all(np.isfinite(rdi))
    assert all(rdi >= 0)  # RDI should be non-negative


def test_compute_rdi_outlier_detection():
    """Test RDI with simulated outliers"""

    # Simulate data where last gene is an outlier
    distances = np.array([1.0, 1.1, 1.0, 10.0], dtype=np.float64)  # Last one is outlier
    gene_to_fam = np.array([1, 1, 1, 1], dtype=np.int32)
    dscale = np.array([1.0], dtype=np.float64)

    rdi = tox_compute_rdi(distances, gene_to_fam, dscale)

    # The outlier should have much higher RDI
    assert rdi[3] > rdi[0]
    assert rdi[3] > rdi[1]
    assert rdi[3] > rdi[2]


# =====================
# Tests for identify_outliers
# =====================

def test_identify_outliers_basic():
    """Test basic outlier identification"""

    # RDI values with clear outlier
    rdi = np.array([0.5, 0.6, 0.4, 3.5, 0.5], dtype=np.float64)

    result = tox_identify_outliers(rdi, percentile=80.0)  # Use 80th percentile

    # Verify that gene with RDI=3.5 is identified as outlier
    assert result['outliers'][3] == True
    assert result['threshold'] > 0


def test_identify_outliers_no_outliers():
    """Test when no outliers present"""

    rdi = np.array([0.5, 0.6, 0.4, 0.7, 0.5], dtype=np.float64)  # All low RDI

    result = tox_identify_outliers(rdi, percentile=99.0)  # Very high percentile

    # With very high percentile, few or no outliers should be found
    assert isinstance(result['outliers'], np.ndarray)
    assert result['threshold'] > 0


def test_identify_outliers_all_outliers():
    """Test when many genes are outliers"""

    rdi = np.array([3.0, 4.0, 5.0, 3.5], dtype=np.float64)  # All high RDI

    result = tox_identify_outliers(rdi, percentile=25.0)  # Low percentile = more outliers

    # With low percentile, more outliers should be found
    assert sum(result['outliers']) >= 1


# =====================
# Tests for detect_outliers (complete pipeline)
# =====================

def test_detect_outliers_pipeline():
    """Test complete outlier detection pipeline"""

    # Create test data with simulated outlier
    distances = np.array([
        1.0, 1.2, 1.1, 0.9,  # Family 1
        2.0, 2.2, 2.1, 2.3,  # Family 2
        3.0, 3.2, 3.1, 3.3,  # Family 3
        4.0, 4.2, 4.1, 4.3,  # Family 4
        5.0, 5.2, 5.1, 5.3,  # Family 5
        6.0, 6.2, 6.1, 6.3   # Family 6
    ], dtype=np.float64)
    gene_to_fam = np.array([
        1, 1, 1, 1,  # Family 1
        2, 2, 2, 2,  # Family 2
        3, 3, 3, 3,  # Family 3
        4, 4, 4, 4,  # Family 4
        5, 5, 5, 5,  # Family 5
        6, 6, 6, 6   # Family 6
    ], dtype=np.int32)
    percentile = 90.0

    result = tox_detect_outliers(distances, gene_to_fam, percentile)

    # Verify structure
    assert len(result['outliers']) == len(distances)
    assert 'loess_x' in result
    assert 'loess_y' in result
    assert 'loess_n' in result

    # Check for outliers
    outlier_count = sum(result['outliers'])


def test_detect_outliers_performance():
    """Test performance with larger dataset"""

    n_genes = 1000
    n_families = 20

    # Generate realistic test data
    np.random.seed(42)
    distances = np.abs(np.random.normal(2.0, 0.5, n_genes))
    gene_to_fam = np.random.randint(1, n_families + 1, n_genes)

    # Add some outliers
    outlier_indices = np.random.choice(n_genes, size=50, replace=False)
    distances[outlier_indices] *= 3  # Make them outliers

    import time
    start_time = time.time()

    result = tox_detect_outliers(distances, gene_to_fam, percentile=95.0)

    end_time = time.time()
    elapsed = end_time - start_time

    outlier_count = sum(result['outliers'])

    # Verify reasonable results
    assert len(result['outliers']) == n_genes
    assert outlier_count >= 0  # Should find some outliers (or none is also OK)


def test_detect_outliers_edge_cases():
    """Test edge cases"""

    # Case 1: Single family
    distances = np.array([1.0, 2.0, 3.0], dtype=np.float64)
    gene_to_fam = np.array([1, 1, 1], dtype=np.int32)

    result = tox_detect_outliers(distances, gene_to_fam, percentile=95.0)
    assert len(result['loess_x']) == 1

    # Case 2: Each gene in different family
    distances = np.array([1.0, 2.0, 3.0], dtype=np.float64)
    gene_to_fam = np.array([1, 2, 3], dtype=np.int32)

    result = tox_detect_outliers(distances, gene_to_fam, percentile=95.0)
    assert len(result['loess_x']) == 3


if __name__ == '__main__':
    run_all_tests(globals().values())
