#!/usr/bin/env python3
"""
Test script for the user-facing Python API of the Fortran CSV reader.
"""
import numpy as np
import os
from tensoromics_functions import (
    read_csv_as_strings,
    read_integer_columns,
    read_real_columns,
    read_logical_columns,
    read_character_columns,
    read_complex_columns
)

def test_full_pipeline():
    """
    Runs a comprehensive test of the entire CSV reading pipeline via the API.
    """
    print("=== Testing Full CSV Reading Pipeline via Python API ===")
    
    csv_content = (
        "ID,Status,Score,Name,Coords\n"
        "101,T,95.5,Alpha,(1.1, 2.2)\n"
        "102,F,88.0,Beta,(3.3, -4.4)\n"
        "103,TRUE,72.1,Gamma,(5.5, 6.6)\n"
        "104,FALSE,-5.0,Delta,(7.7, 8.8)\n"
    )
    test_filename = "python_api_test.csv"
    with open(test_filename, "w") as f:
        f.write(csv_content)

    try:
        # Step 1: Read all data as strings using the new API function
        header, data_str = read_csv_as_strings(test_filename)
        
        # Assert read-only flag
        assert not data_str.flags.writeable
        assert not header.flags.writeable
        
        print("Step 1: read_csv_as_strings PASSED")
        np.testing.assert_array_equal(header, ["ID", "Status", "Score", "Name", "Coords"])
        assert data_str[1, 3] == "Beta"

        # Step 2: Convert to specific types
        
        # Integers (Column 1)
        int_data = read_integer_columns(data_str, columns=[1])
        assert not int_data.flags.writeable
        np.testing.assert_array_equal(int_data.flatten(), [101, 102, 103, 104])
        print("Step 2: read_integer_columns PASSED")

        # Reals (Column 3)
        real_data = read_real_columns(data_str, columns=[3])
        assert not real_data.flags.writeable
        np.testing.assert_allclose(real_data.flatten(), [95.5, 88.0, 72.1, -5.0])
        print("Step 3: read_real_columns PASSED")

        # Logicals (Column 2)
        logical_data = read_logical_columns(data_str, columns=[2])
        assert not logical_data.flags.writeable
        np.testing.assert_array_equal(logical_data.flatten(), [True, False, True, False])
        print("Step 4: read_logical_columns PASSED")
        
        # Characters (Column 4)
        char_data = read_character_columns(data_str, columns=[4])
        assert not char_data.flags.writeable
        np.testing.assert_array_equal(char_data.flatten(), ["Alpha", "Beta", "Gamma", "Delta"])
        print("Step 5: read_character_columns PASSED")

        # Complex (Column 5)
        complex_data = read_complex_columns(data_str, columns=[5])
        assert not complex_data.flags.writeable
        np.testing.assert_allclose(complex_data.flatten(), [1.1+2.2j, 3.3-4.4j, 5.5+6.6j, 7.7+8.8j])
        print("Step 6: read_complex_columns PASSED")

    finally:
        if os.path.exists(test_filename):
            os.remove(test_filename)

if __name__ == "__main__":
    test_full_pipeline()
    print("\nAll Python API tests passed! ✓")