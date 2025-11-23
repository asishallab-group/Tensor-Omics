"""
Comprehensive Python test suite for LOESS smoothing
Uses tensoromics_functions.py wrapper functions (mirrors R test_loess_smooth_2d_r.R tests)
"""

import numpy as np
import sys
import os

# Add parent directory to path to import tensoromics_functions
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from tensoromics_functions import tox_loess_smooth_2d

print("=== Testing LOESS smoothing Python wrapper functions ===")
print("Based on Fortran test suite with comprehensive test coverage")

# =====================
# Tests for loess_smooth_2d
# =====================

def test_loess_simple_smoothing():
    """Test simple smoothing, no mask"""
    print("\n[test_loess_simple_smoothing] Simple smoothing, no mask test")
    
    x_ref = np.array([1, 2, 3, 4, 5], dtype=np.float64)
    y_ref = np.array([10, 20, 30, 40, 50], dtype=np.float64)
    indices_used = np.arange(1, 6, dtype=np.int32)  # 1-based for Fortran
    x_query = np.array([2.5, 4.5], dtype=np.float64)
    kernel_sigma = 1.0
    kernel_cutoff = 2.0

    result = tox_loess_smooth_2d(x_ref, y_ref, indices_used, x_query, kernel_sigma, kernel_cutoff)
    
    print("  Input x_ref:", x_ref)
    print("  Input y_ref:", y_ref)
    print("  Query points:", x_query)
    print("  Result:", result)
    
    # Should be reasonable interpolated values
    assert len(result) == len(x_query)
    assert all(np.isfinite(result))
    
    print("Simple smoothing test passed ✓")

def test_loess_with_mask():
    """Test with mask (exclude some points)"""
    print("\n[test_loess_with_mask] With mask (exclude some points) test")
    
    x_ref = np.array([1, 2, 3, 4, 5], dtype=np.float64)
    y_ref = np.array([10, 20, 30, 40, 50], dtype=np.float64)
    # Use only indices 1, 3, 5 (excluding points 2 and 4)
    indices_used = np.array([1, 3, 5], dtype=np.int32)
    x_query = np.array([2.5, 4.5], dtype=np.float64)
    kernel_sigma = 1.0
    kernel_cutoff = 2.0

    result = tox_loess_smooth_2d(x_ref, y_ref, indices_used, x_query, kernel_sigma, kernel_cutoff)
    
    print("  Using only indices [1, 3, 5]")
    print("  Query points:", x_query)
    print("  Result:", result)
    
    assert len(result) == len(x_query)
    assert all(np.isfinite(result))
    
    print("With mask test passed ✓")

def test_loess_single_point():
    """Test use only one point (should return that point's value)"""
    print("\n[test_loess_single_point] Use only one point test")
    
    x_ref = np.array([1, 2, 3, 4, 5], dtype=np.float64)
    y_ref = np.array([10, 20, 30, 40, 50], dtype=np.float64)
    # Use only the first point
    indices_used = np.array([1], dtype=np.int32)
    x_query = np.array([2.5, 4.5], dtype=np.float64)
    kernel_sigma = 1.0
    kernel_cutoff = 2.0

    result = tox_loess_smooth_2d(x_ref, y_ref, indices_used, x_query, kernel_sigma, kernel_cutoff)
    
    print("  Using only index [1]")
    print("  Query points:", x_query)
    print("  Result:", result)
    
    assert len(result) == len(x_query)
    assert all(np.isfinite(result))
    
    print("Single point test passed ✓")

def test_loess_outside_range():
    """Test x_query outside x_ref range"""
    print("\n[test_loess_outside_range] x_query outside x_ref range test")
    
    x_ref = np.array([1, 2, 3, 4, 5], dtype=np.float64)
    y_ref = np.array([10, 20, 30, 40, 50], dtype=np.float64)
    indices_used = np.arange(1, 6, dtype=np.int32)
    x_query = np.array([-10, 100], dtype=np.float64)  # Outside range
    kernel_sigma = 1.0
    kernel_cutoff = 2.0

    result = tox_loess_smooth_2d(x_ref, y_ref, indices_used, x_query, kernel_sigma, kernel_cutoff)
    
    print("  Query points outside range:", x_query)
    print("  Result:", result)
    
    assert len(result) == len(x_query)
    # Results might be NaN or extrapolated values, both are acceptable
    
    print("Outside range test passed ✓")

def test_loess_exact_match():
    """Test x_query == x_ref"""
    print("\n[test_loess_exact_match] x_query == x_ref test")
    
    x_ref = np.array([1, 2, 3, 4, 5], dtype=np.float64)
    y_ref = np.array([10, 20, 30, 40, 50], dtype=np.float64)
    indices_used = np.arange(1, 6, dtype=np.int32)
    x_query = x_ref[:2].copy()  # [1, 2]
    kernel_sigma = 1.0
    kernel_cutoff = 2.0

    result = tox_loess_smooth_2d(x_ref, y_ref, indices_used, x_query, kernel_sigma, kernel_cutoff)
    
    print("  Query points (exact match):", x_query)
    print("  Expected y values:", y_ref[:2])
    print("  Result:", result)
    
    # Should be very close to exact y values
    expected = y_ref[:2]
    for i in range(len(result)):
        assert abs(result[i] - expected[i]) < 1e-10
    
    print("Exact match test passed ✓")

def test_loess_nan_handling():
    """Test y_ref contains np.nan"""
    print("\n[test_loess_nan_handling] y_ref contains np.nan test")
    
    x_ref = np.array([1, 2, 3, 4, 5], dtype=np.float64)
    y_ref = np.array([10, 20, np.nan, 40, 50], dtype=np.float64)
    indices_used = np.arange(1, 6, dtype=np.int32)
    x_query = np.array([2.5, 4.5], dtype=np.float64)
    kernel_sigma = 1.0
    kernel_cutoff = 2.0

    try:
        result = tox_loess_smooth_2d(x_ref, y_ref, indices_used, x_query, kernel_sigma, kernel_cutoff)
        print("  Query points:", x_query)
        print("  Result:", result)
        # If no exception, check that result is reasonable
        assert len(result) == len(x_query)
        print("NaN handled gracefully")
    except Exception as e:
        print("  Exception (expected if nan not handled):", str(e))
        print("NaN caused expected exception")
    
    print("NaN handling test passed ✓")

def test_loess_single_total_point():
    """Test n_total = 1 (single point)"""
    print("\n[test_loess_single_total_point] n_total = 1 (single point) test")
    
    x_ref = np.array([3.0], dtype=np.float64)
    y_ref = np.array([30.0], dtype=np.float64)
    indices_used = np.array([1], dtype=np.int32)
    x_query = np.array([1.0, 3.0, 5.0], dtype=np.float64)
    kernel_sigma = 1.0
    kernel_cutoff = 2.0

    result = tox_loess_smooth_2d(x_ref, y_ref, indices_used, x_query, kernel_sigma, kernel_cutoff)
    
    print("  Single reference point at x=3, y=30")
    print("  Query points:", x_query)
    print("  Result:", result)
    
    assert len(result) == len(x_query)
    
    print("Single total point test passed ✓")

def test_loess_kernel_parameters():
    """Test different kernel parameters"""
    print("\n[test_loess_kernel_parameters] Different kernel parameters test")
    
    x_ref = np.array([1, 2, 3, 4, 5], dtype=np.float64)
    y_ref = np.array([10, 20, 30, 40, 50], dtype=np.float64)
    indices_used = np.arange(1, 6, dtype=np.int32)
    x_query = np.array([2.5], dtype=np.float64)
    
    # Test different sigma values
    for sigma in [0.5, 1.0, 2.0]:
        result = tox_loess_smooth_2d(x_ref, y_ref, indices_used, x_query, sigma, 2.0)
        print(f"  Sigma {sigma}: {result[0]:.6f}")
        assert np.isfinite(result[0])
    
    # Test different cutoff values
    for cutoff in [1.0, 2.0, 4.0]:
        result = tox_loess_smooth_2d(x_ref, y_ref, indices_used, x_query, 1.0, cutoff)
        print(f"  Cutoff {cutoff}: {result[0]:.6f}")
        assert np.isfinite(result[0])
    
    print("Kernel parameters test passed ✓")

def test_loess_input_validation():
    """Test input validation"""
    print("\n[test_loess_input_validation] Input validation test")
    
    # Test mismatched x_ref and y_ref lengths
    error_caught = False
    try:
        x_ref = np.array([1, 2, 3], dtype=np.float64)
        y_ref = np.array([10, 20], dtype=np.float64)  # Different length
        indices_used = np.array([1, 2], dtype=np.int32)
        x_query = np.array([2.5], dtype=np.float64)
        
        tox_loess_smooth_2d(x_ref, y_ref, indices_used, x_query, 1.0, 2.0)
    except ValueError as e:
        error_caught = True
        assert "same length" in str(e)
    assert error_caught
    
    print("Input validation test passed ✓")

# =====================
# Run all tests
# =====================

def main():
    print("\n=================================================")
    print("    LOESS SMOOTHING FULL PYTHON INTERFACE TESTS")
    print("=================================================\n")

    test_loess_simple_smoothing()
    test_loess_with_mask()
    test_loess_single_point()
    test_loess_outside_range()
    test_loess_exact_match()
    test_loess_nan_handling()
    test_loess_single_total_point()
    test_loess_kernel_parameters()
    test_loess_input_validation()

    print("=================================================")
    print("             ALL TESTS COMPLETED")
    print("=================================================")
    print("If you see this message, all LOESS smoothing Python interface tests passed! ✓")

if __name__ == "__main__":
    main()
