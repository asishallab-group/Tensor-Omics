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
from test_helpers import run_all_tests


def test_tox_which_example_1():
    """Example 1: Simple which test"""
    mask = np.array([1, 0, 1, 0, 0], dtype=np.int32)
    idx_out = tox_which(mask)
    # Only valid indices are the first nonzero elements
    valid_indices = idx_out[idx_out != 0]
    m_out = len(valid_indices)
    # Manual verification - should get indices 1 and 3 (1-based)
    expected_count = 2
    expected_indices = [1, 3]

    assert m_out == expected_count, f"Count should be {expected_count}, got {m_out}"
    actual_indices = valid_indices.tolist()
    assert actual_indices == expected_indices, f"Indices should be {expected_indices}, got {actual_indices}"


def test_tox_which_example_2():
    """Example 2: All FALSE case"""
    mask = np.zeros(5, dtype=np.int32)
    idx_out = tox_which(mask)
    valid_indices = idx_out[idx_out != 0]
    m_out = len(valid_indices)
    assert m_out == 0, f"Count should be 0, got {m_out}"

def test_tox_which_example_3():
    """Example 3: All TRUE case"""
    mask = np.ones(5, dtype=np.int32)
    idx_out = tox_which(mask)
    valid_indices = idx_out[idx_out != 0]
    m_out = len(valid_indices)
    expected_count = 5
    expected_indices = [1, 2, 3, 4, 5]

    assert m_out == expected_count, f"Count should be {expected_count}, got {m_out}"
    actual_indices = valid_indices.tolist()
    assert actual_indices == expected_indices, f"Indices should be {expected_indices}, got {actual_indices}"


if __name__ == "__main__":
    run_all_tests(globals().values())
