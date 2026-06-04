"""
Test suite for tox_normalization_pipeline (Python)
Mirrors R/Fortran tests: basic, edge case, pipeline vs manual
"""
import numpy as np
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from tensoromics_functions import tox_normalization_pipeline, tox_normalize_by_std_dev, tox_quantile_normalization, tox_calculate_tissue_averages, tox_log2_transformation
from test_helpers import run_all_tests


def test_basic():
    n_genes, n_tissues, n_grps = 10, 6, 2
    input_matrix = np.empty((n_tissues, n_genes), dtype=np.float64, order="F")
    for i in range(n_genes):
        for j in range(n_tissues):
            input_matrix[j, i] = (i + 1) * 2.0 + (j + 1) * 0.5
    group_s = np.array([1,4], dtype=np.int32)
    group_c = np.array([3,3], dtype=np.int32)
    result = tox_normalization_pipeline(input_matrix, group_c, span=0.75, degree=2, use_quantile=1)
    assert np.all(~np.isnan(result)), "NaN in output"
    assert np.all(result >= 0), "Negative values in output"
    assert result.shape == (n_grps, n_genes), f"Output shape mismatch: {result.shape}"


def test_edge_case():
    n_genes, n_tissues, n_grps = 10, 6, 2
    input_matrix = np.zeros((n_genes, n_tissues), dtype=np.float64)
    group_s = np.array([1,4], dtype=np.int32)
    group_c = np.array([3,3], dtype=np.int32)
    try:
        tox_normalization_pipeline(input_matrix, group_c, span=0.75, degree=2, use_quantile=1)
        assert False, "Expected error for zero-variance input"
    except Exception as err:
        pass


def test_pipeline_vs_manual():
    """Compare normalization pipeline with manual stepwise operations using LOESS"""
    # print("")
    # print("=== Test: test_pipeline_vs_manual ===")
    n_genes, n_tissues, n_grps = 10, 6, 2
    col_names = [
        "muscle_dietM_1", "muscle_dietM_2", "muscle_dietM_3",
        "muscle_dietP_1", "muscle_dietP_2", "muscle_dietP_3",
    ]
    input_matrix = np.empty((n_tissues, n_genes), dtype=np.float64, order="F")
    for i in range(n_genes):
        for j in range(n_tissues):
            input_matrix[j, i] = (i + 1) * 2.0 + (j + 1) * 0.5
    group_s = np.array([1,4], dtype=np.int32)
    group_c = np.array([3,3], dtype=np.int32)
    span = 0.75
    degree = 2

    # Test WITH quantile normalization
    result_pipeline_quantile = tox_normalization_pipeline(
        input_matrix.copy(), group_c, span=span, degree=degree, use_quantile=1
    )

    # Manual with quantile: std_dev → quantile → avg → log2
    temp1 = tox_normalize_by_std_dev(input_matrix.copy(), span=span, degree=degree)
    temp2 = tox_quantile_normalization(temp1)
    temp3 = tox_calculate_tissue_averages(temp2, group_c)
    result_manual_quantile = tox_log2_transformation(temp3)

    # Test WITHOUT quantile normalization
    result_pipeline_no_quantile = tox_normalization_pipeline(
        input_matrix.copy(), group_c, span=span, degree=degree, use_quantile=0
    )

    # Manual without quantile: std_dev → avg → log2
    temp1_no_q = tox_normalize_by_std_dev(input_matrix.copy(), span=span, degree=degree)
    temp2_no_q = tox_calculate_tissue_averages(temp1_no_q, group_c)
    result_manual_no_quantile = tox_log2_transformation(temp2_no_q)

    # Verify shapes
    assert result_pipeline_quantile.shape == (n_grps, n_genes), 1
    assert result_pipeline_no_quantile.shape == (n_grps, n_genes), 2
    assert result_manual_quantile.shape == (n_grps, n_genes), 3
    assert result_manual_no_quantile.shape == (n_grps, n_genes), 4

    # Pipeline should match manual operations
    assert np.allclose(result_pipeline_quantile, result_manual_quantile, rtol=1e-10), \
        "Pipeline with quantile should match manual stepwise operations"
    assert np.allclose(result_pipeline_no_quantile, result_manual_no_quantile, rtol=1e-10), \
        "Pipeline without quantile should match manual stepwise operations"

    # No NaN values
    assert np.all(~np.isnan(result_pipeline_quantile)), "NaN in pipeline output (quantile)"
    assert np.all(~np.isnan(result_pipeline_no_quantile)), "NaN in pipeline output (no quantile)"

    # Results should differ when quantile option changes
    assert np.any(np.abs(result_pipeline_quantile - result_pipeline_no_quantile) > 1e-12), \
        "Quantile and non-quantile results should differ"


if __name__ == '__main__':
    run_all_tests(globals().values())
