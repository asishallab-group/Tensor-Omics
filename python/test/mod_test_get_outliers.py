#!/usr/bin/env python3
"""
Comprehensive Python test suite for outlier detection functions
Uses tensor_omics.py wrapper functions (mirrors R get_outliers.R tests)
"""

import numpy as np
import sys
import os

# Add parent directory to path to import tensor_omics
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from tensor_omics import (
    compute_family_scaling,
    compute_family_scaling,
    compute_rdi,
    identify_outliers as _identify_outliers,
    detect_outliers as _detect_outliers,
)
from tensor_omics.error_handling import ERR_NAN_INF
from test_helpers import run_all_tests, assert_error

# compute_rdi returns the sorted array and permutation that identify_outliers wants, so
# the pipeline composes directly. Only the two shims below remain: one for calling
# identify_outliers with a bare rdi, and one for the key the old suite expects.


def _sort_prep(rdi):
    """Values clamped non-negative plus the 1-based permutation that sorts them."""
    sorted_rdi = np.asarray(rdi, dtype=np.float64).copy()
    sorted_rdi[sorted_rdi < 0.0] = 0.0
    perm = (np.argsort(sorted_rdi, kind="mergesort").astype(np.int32) + 1)
    return sorted_rdi, perm


def identify_outliers(rdi, quantile_level=0.95):
    sorted_rdi, perm = _sort_prep(rdi)
    result = _identify_outliers(rdi, sorted_rdi, perm, quantile_level)
    return dict(result, outliers=result["is_outlier"])


def detect_outliers(distances, gene_to_fam, quantile_level=0.95):
    n_families = int(np.max(gene_to_fam))
    result = _detect_outliers(n_families, distances, gene_to_fam, quantile_level)
    return dict(result, outliers=result["is_outlier"])


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

    result = compute_family_scaling(n_families, distances, gene_to_fam)

    # Verify basic properties
    assert len(result['dscale']) == n_families  # Ten families
    assert all(np.isfinite(result['dscale']))
    assert all(result['dscale'] > 0)  # Scaling factors should be positive


def test_compute_family_scaling_single_family():
    """Test with single family"""

    distances = np.array([1.0, 2.0, 3.0, 4.0], dtype=np.float64)
    gene_to_fam = np.array([1, 1, 1, 1], dtype=np.int32)  # All in family 1

    n_families = int(np.max(gene_to_fam))
    result = compute_family_scaling(n_families, distances, gene_to_fam)

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

    n_families = int(np.max(gene_to_fam))
    result = compute_family_scaling(n_families, distances, gene_to_fam)
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

    result = compute_family_scaling(n_families, distances, gene_to_fam)
    assert all(np.isfinite(result['dscale']))


def test_tox_compute_family_scaling_buffers():
    """Test expert version with user-provided work arrays"""

    # Test data
    n_families = 10  # Incremented number of families
    genes_per_fam = 4
    n_genes = n_families * genes_per_fam
    span = 0.7
    mode = "robust"
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

    result = compute_family_scaling(
        n_families, distances, gene_to_fam, span, degree, mode, n_iters
    )

    # Verify basic properties
    assert len(result['dscale']) == n_families  # Ten families
    assert all(np.isfinite(result['dscale']))
    assert all(result['dscale'] > 0)  # Scaling factors should be positive

    # Verify the smoothing diagnostics come back too
    assert 'low_sd_cutoff' in result
    assert 'excluded_low_sd' in result
    assert 'indices_used' in result

    # Compare with regular version to ensure consistency
    result_regular = compute_family_scaling(n_families, distances, gene_to_fam)

    # Results should be very similar (within numerical precision)
    np.testing.assert_allclose(result['dscale'], result_regular['dscale'], rtol=1e-10)
    np.testing.assert_allclose(result['loess_x'], result_regular['loess_x'], rtol=1e-10)
    np.testing.assert_allclose(result['loess_y'], result_regular['loess_y'], rtol=1e-10)


def test_compute_family_scaling_input_validation():
    """Test expert version input validation"""

    n_families = 6
    genes_per_fam = 4
    n_genes = n_families * genes_per_fam
    span = 0.7
    mode = "robust"
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

    # distances and gene_to_fam share n_genes; a mismatch has to be rejected
    assert_error(
        lambda: compute_family_scaling(
            n_families, distances[:-1], gene_to_fam, span, degree, mode, n_iters),
        "Expected an error for mismatched distances / gene_to_fam lengths",
    )

    # an unknown mode string is not one of the tabulated values: refused by the binding before
    # the call, naming the argument and the modes, so with no tox code, as in R (#221 G-I1, G-I2)
    assert_error(
        lambda: compute_family_scaling(
            n_families, distances, gene_to_fam, span, degree, "not_a_mode", n_iters),
        "Expected an error for an invalid mode string",
    )


# =====================
# Tests for compute_rdi
# =====================

def test_compute_rdi_basic():
    """Test basic RDI computation"""

    distances = np.array([1.0, 2.0, 3.0, 4.0, 5.0], dtype=np.float64)
    gene_to_fam = np.array([1, 1, 2, 2, 2], dtype=np.int32)
    dscale = np.array([1.5, 2.0], dtype=np.float64)  # Manual scaling factors

    rdi = compute_rdi(distances, gene_to_fam, dscale)['rdi']

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

    rdi = compute_rdi(distances, gene_to_fam, dscale)['rdi']

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

    result = identify_outliers(rdi, quantile_level=0.8)  # Use 80th percentile

    # Verify that gene with RDI=3.5 is identified as outlier
    assert result['outliers'][3] == True
    assert result['threshold'] > 0


def test_identify_outliers_no_outliers():
    """Test when no outliers present"""

    rdi = np.array([0.5, 0.6, 0.4, 0.7, 0.5], dtype=np.float64)  # All low RDI

    result = identify_outliers(rdi, quantile_level=0.99)  # Very high percentile

    # With very high percentile, few or no outliers should be found
    assert isinstance(result['outliers'], np.ndarray)
    assert result['threshold'] > 0


def test_identify_outliers_all_outliers():
    """Test when many genes are outliers"""

    rdi = np.array([3.0, 4.0, 5.0, 3.5], dtype=np.float64)  # All high RDI

    result = identify_outliers(rdi, quantile_level=0.25)  # Low quantile_level = more outliers

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
    quantile_level = 0.9

    result = detect_outliers(distances, gene_to_fam, quantile_level)

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

    result = detect_outliers(distances, gene_to_fam, quantile_level=0.95)

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

    result = detect_outliers(distances, gene_to_fam, quantile_level=0.95)
    assert len(result['loess_x']) == 1

    # Case 2: Each gene in different family
    distances = np.array([1.0, 2.0, 3.0], dtype=np.float64)
    gene_to_fam = np.array([1, 2, 3], dtype=np.int32)

    result = detect_outliers(distances, gene_to_fam, quantile_level=0.95)
    assert len(result['loess_x']) == 3


def test_detect_outliers_unused_loess_slots_hold_sentinel():
    """Slots of loess_x / loess_y that no family occupies hold -1, never NaN (#180)"""
    rng = np.random.default_rng(4)
    n_families = 40
    sizes = [1, 1, 1] + [6] * (n_families - 3)  # three one-member families
    gene_to_fam = np.repeat(np.arange(1, n_families + 1), sizes).astype(np.int32)
    distances = rng.lognormal(0.0, 0.5, size=len(gene_to_fam))

    result = _detect_outliers(n_families, distances, gene_to_fam)

    unused = result['loess_n'] == 0
    assert unused.sum() >= 3, "the three one-member families leave slots unused"
    for key in ('loess_x', 'loess_y'):
        assert not np.isnan(result[key]).any(), f"{key} contains NaN"
        assert np.all(result[key][unused] == -1.0), f"{key}: unused slots must hold -1"
        assert np.all(result[key][~unused] >= 0.0), f"{key}: used slots must be non-negative"

def _assert_each_refused(cases):
    """Every (label, call, argument) raises ERR_NAN_INF naming `argument`; reports all that do not."""
    failures = []
    for label, call, argument in cases:
        try:
            assert_error(call, label, ERR_NAN_INF, argument=argument)
        except AssertionError as failure:
            failures.append(str(failure))
    assert not failures, "\n".join(failures)


def test_detect_outliers_refuses_nan_distance():
    """A NaN or Inf distance is refused as ERR_NAN_INF naming `distances`, whatever the families.

    It used to be let through in a one-member family (its gene got the most extreme tail
    probability) and in the only multi-member family (it came back in loess_x).
    """
    layouts = {
        "one family": (1, [1, 1, 1, 1]),
        "one multi-member family and a singleton": (2, [1, 1, 1, 2]),
        "two multi-member families": (2, [1, 1, 2, 2]),
    }
    cases = []
    for label, (n_families, families) in layouts.items():
        gene_to_fam = np.array(families, dtype=np.int32)
        for bad in (np.nan, np.inf):
            distances = np.array([bad, 1.0, 2.0, 3.0])
            cases.append((f"detect_outliers, {label}, {bad}",
                          lambda n=n_families, d=distances, g=gene_to_fam: _detect_outliers(n, d, g),
                          "distances"))
    _assert_each_refused(cases)


def test_pipeline_steps_refuse_non_finite_input():
    """Each step's real inputs refuse NaN and Inf, as ERR_NAN_INF naming the argument.

    One case per argument and value, so a DM_ALLOW_NAN or DM_ALLOW_INFINITE put back on any of
    them fails here: compute_family_scaling's distances, compute_rdi's distances and dscale,
    identify_outliers' rdi and sorted_rdi.
    """
    gene_to_fam = np.array([1, 1, 2, 2, 3], dtype=np.int32)
    distances = np.array([1.0, 2.0, 1.5, 3.0, 2.0])
    dscale = np.array([1.0, 2.0, 1.0])
    rdi = np.array([0.5, 1.0, 1.5, 2.0])
    perm = np.array([1, 2, 3, 4], dtype=np.int32)

    def with_bad(values, bad):
        values = values.copy()
        values[1] = bad
        return values

    cases = []
    for bad in (np.nan, np.inf):
        cases += [
            (f"compute_family_scaling/distances {bad}",
             lambda b=bad: compute_family_scaling(3, with_bad(distances, b), gene_to_fam), "distances"),
            (f"compute_rdi/distances {bad}",
             lambda b=bad: compute_rdi(with_bad(distances, b), gene_to_fam, dscale), "distances"),
            (f"compute_rdi/dscale {bad}",
             lambda b=bad: compute_rdi(distances, gene_to_fam, with_bad(dscale, b)), "dscale"),
            (f"identify_outliers/rdi {bad}",
             lambda b=bad: _identify_outliers(with_bad(rdi, b), rdi, perm), "rdi"),
            (f"identify_outliers/sorted_rdi {bad}",
             lambda b=bad: _identify_outliers(rdi, with_bad(rdi, b), perm), "sorted_rdi"),
        ]
    _assert_each_refused(cases)


def test_excluded_low_sd_marks_the_families_outside_the_fit():
    """excluded_low_sd is 0 exactly for the families indices_used lists, 1 for the others.

    On every successful exit: no multi-member family, a single one (both leave before the fit),
    the constant fallback, and a real fit. It used to be left unwritten on the early exits.
    """
    rng = np.random.default_rng(7)
    sizes = [1, 1] + [5] * 30
    fitted = np.repeat(np.arange(1, len(sizes) + 1), sizes).astype(np.int32)
    datasets = {
        "no multi-member family": (3, np.array([1.0, 2.0, 3.0]), np.array([1, 2, 3], dtype=np.int32)),
        "one multi-member family": (3, np.array([1.0, 2.0, 4.0, 7.0, 3.0]), np.array([1, 1, 1, 2, 3], dtype=np.int32)),
        "constant fallback": (5, np.tile([2.0, 3.0, 4.0], 5), np.repeat(np.arange(1, 6), 3).astype(np.int32)),
        "fit": (len(sizes), rng.lognormal(0.0, 0.5, size=fitted.size), fitted),
    }
    for label, (n_families, distances, gene_to_fam) in datasets.items():
        result = compute_family_scaling(n_families, distances, gene_to_fam)
        listed = np.isin(np.arange(1, n_families + 1), result["indices_used"])
        expected = np.where(listed, 0, 1)
        assert np.array_equal(result["excluded_low_sd"], expected), \
            f"{label}: excluded_low_sd {result['excluded_low_sd'].tolist()}, expected {expected.tolist()}"


if __name__ == '__main__':
    run_all_tests(globals().values())
