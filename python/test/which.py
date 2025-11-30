"""
'which' utility for Python, like in R/MATLAB.
Includes usage examples.
"""

import numpy as np
import sys
import os

# Add parent directory to path to import tensoromics_functions
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from tensoromics_functions import (
    tox_which
)

def test_tox_which_example_1():
    """Example 1: Simple which test"""
    print("="*50)
    print("TOX_WHICH EXAMPLE 1: Simple mask")
    print("="*50)
    
    mask = np.array([1, 0, 1, 0, 0], dtype=np.int32)
    print(f"Input mask: {mask}")
    
    idx_out, m_out = tox_which(mask)
    print(f"Output indices: {idx_out}")
    print(f"Count: {m_out}")
    
    # Manual verification - should get indices 1 and 3 (1-based)
    expected_count = 2
    expected_indices = [1, 3]
    
    assert m_out == expected_count, f"Count should be {expected_count}, got {m_out}"
    actual_indices = idx_out[:m_out].tolist()
    assert actual_indices == expected_indices, f"Indices should be {expected_indices}, got {actual_indices}"
    
    print("✓ Test passed!")

def test_tox_which_example_2():
    """Example 2: All FALSE case"""
    print("="*50)
    print("TOX_WHICH EXAMPLE 2: All FALSE")
    print("="*50)
    
    mask = np.zeros(5, dtype=np.int32)
    print(f"Input mask: {mask}")
    
    idx_out, m_out = tox_which(mask)
    print(f"Output indices: {idx_out}")
    print(f"Count: {m_out}")
    
    assert m_out == 0, f"Count should be 0, got {m_out}"
    print("✓ Test passed!")

def test_tox_which_example_3():
    """Example 3: All TRUE case"""
    print("="*50)
    print("TOX_WHICH EXAMPLE 3: All TRUE")
    print("="*50)
    
    mask = np.ones(5, dtype=np.int32)
    print(f"Input mask: {mask}")
    
    idx_out, m_out = tox_which(mask)
    print(f"Output indices: {idx_out}")
    print(f"Count: {m_out}")
    
    expected_count = 5
    expected_indices = [1, 2, 3, 4, 5]
    
    assert m_out == expected_count, f"Count should be {expected_count}, got {m_out}"
    actual_indices = idx_out[:m_out].tolist()
    assert actual_indices == expected_indices, f"Indices should be {expected_indices}, got {actual_indices}"
    
    print("✓ Test passed!")

if __name__ == "__main__":
    print("TENSOR-OMICS PYTHON TOX_ FUNCTIONS TEST SUITE")
    print("Testing wrapper functions with tox_ prefix...")
    
    try:

        test_tox_which_example_1()
        test_tox_which_example_2()
        test_tox_which_example_3()
        
        print("\n" + "="*50)
        print("ALL TOX_ FUNCTION TESTS COMPLETED SUCCESSFULLY!")
        print("All utility functions working correctly with tox_ prefix.")
        print("="*50)

    except Exception as e:
        print(f"\nERROR: {e}")
        import traceback
        traceback.print_exc()
        print("Check function implementations and signatures.")
