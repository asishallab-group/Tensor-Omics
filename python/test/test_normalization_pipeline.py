"""
Test suite for tox_normalization_pipeline (Python)
Mirrors R/Fortran tests: basic, edge case, pipeline vs manual
"""
import numpy as np
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from tensoromics_functions import tox_normalization_pipeline, tox_normalize_by_std_dev, tox_quantile_normalization, tox_calculate_tissue_averages, tox_log2_transformation

def assert_equal(x, y, tol=1e-12, msg=None):
    x = np.asarray(x)
    y = np.asarray(y)
    assert x.shape == y.shape, f"Shape mismatch: {x.shape} vs {y.shape}"
    assert np.all(np.abs(x - y) < tol), f"Values differ: {x} vs {y}"
    if msg:
        print(msg, "✓")

def basic_test():
    n_genes, n_tissues, n_grps = 2, 2, 2
    input_matrix = np.array([[1,3],[2,4]], dtype=np.float64)  # shape (2,2), col-major
    group_s = np.array([1,2], dtype=np.int32)
    group_c = np.array([1,1], dtype=np.int32)
    print(f"[DEBUG] basic_test: n_genes={n_genes} n_tissues={n_tissues} n_grps={n_grps}")
    print(f"[DEBUG] input_matrix=\n{input_matrix}")
    print(f"[DEBUG] group_s={group_s} group_c={group_c}")
    result = tox_normalization_pipeline(input_matrix, group_s, group_c)
    print(f"[DEBUG] buf_log=\n{result}")
    assert np.all(~np.isnan(result)), "NaN in output"
    assert np.all(result >= 0), "Negative values in output"
    assert result.shape == (n_genes, n_grps), f"Output shape mismatch: {result.shape}"
    print("basic_test passed\n")

def edge_case_test():
    n_genes, n_tissues, n_grps = 2, 2, 2
    input_matrix = np.zeros((2,2), dtype=np.float64)
    group_s = np.array([1,2], dtype=np.int32)
    group_c = np.array([1,1], dtype=np.int32)
    print(f"[DEBUG] edge_case_test: n_genes={n_genes} n_tissues={n_tissues} n_grps={n_grps}")
    print(f"[DEBUG] input_matrix=\n{input_matrix}")
    print(f"[DEBUG] group_s={group_s} group_c={group_c}")
    result = tox_normalization_pipeline(input_matrix, group_s, group_c)
    print(f"[DEBUG] buf_log=\n{result}")
    assert np.all(result == 0), "Output not all zeros for zero input"
    print("edge_case_test passed\n")

def pipeline_vs_manual():
    n_genes, n_tissues, n_grps = 2, 2, 2
    col_names = ["muscle_dietM_1", "muscle_dietP_1"]
    input_matrix = np.array([[1,3],[2,4]], dtype=np.float64)
    group_s = np.array([1,2], dtype=np.int32)
    group_c = np.array([1,1], dtype=np.int32)
    print(f"[DEBUG] pipeline_vs_manual: n_genes={n_genes} n_tissues={n_tissues} n_grps={n_grps}")
    print(f"[DEBUG] input_matrix=\n{input_matrix}")
    print(f"[DEBUG] col_names={col_names}")
    print(f"[DEBUG] group_s={group_s} group_c={group_c}")
    buf_stddev = tox_normalize_by_std_dev(input_matrix)
    buf_quant = tox_quantile_normalization(buf_stddev)
    buf_avg = tox_calculate_tissue_averages(buf_quant, group_s, group_c)
    print(f"[DEBUG] buf_avg=\n{buf_avg}")
    manual_out = tox_log2_transformation(buf_avg)
    print(f"[DEBUG] manual_out=\n{manual_out}")
    result = tox_normalization_pipeline(input_matrix, group_s, group_c)
    print(f"[DEBUG] buf_log=\n{result}")
    assert_equal(result, manual_out, msg="pipeline_vs_manual")

if __name__ == "__main__":
    print("=================================================")
    print("    NORMALIZATION PIPELINE FULL PYTHON TESTS")
    print("=================================================\n")
    basic_test()
    edge_case_test()
    pipeline_vs_manual()
    print("=================================================")
    print("             ALL TESTS COMPLETED")
    print("=================================================")
    print("If you see this message, all normalization pipeline Python tests passed! ✓")
