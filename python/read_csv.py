#!/usr/bin/env python3
"""
Test script for the Fortran CSV reader library.
This script uses ctypes to call the C-compatible Fortran functions
and verifies their output against expected results.
"""

import ctypes
import numpy as np
import os
import platform
from collections import defaultdict

# --- 1. Library Loading ---

def _load_fortran_library():
    """
    Finds and loads the compiled Fortran shared library (.so, .dll, .dylib)
    created by the build.sh script.
    """
    # The build.sh script creates a symlink in the build directory
    lib_name = "libtensor-omics.so"
        
    project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    # The build.sh script places the library in the root of the build directory
    lib_path = os.path.join(project_root, "build", lib_name)

    if not os.path.exists(lib_path):
        raise OSError(f"Could not find library at '{lib_path}'. "
                      "Please build the Fortran project using the './build.sh' script.")

    # On some systems, dependent libraries like libgomp need to be preloaded.
    try:
        if platform.system() == "Linux":
            ctypes.CDLL("libgomp.so.1", mode=ctypes.RTLD_GLOBAL)
    except OSError:
        print("Warning: Could not preload libgomp.so.1. This may cause issues if OpenMP is used.")

    return ctypes.CDLL(lib_path)

# --- 2. Setup Functions for Fortran Procedures ---

def setup_fortran_functions(lib):
    """
    Sets up the argtypes and restype for all C-wrapped Fortran functions.
    """
    functions = {}
    
    # --- get_csv_dims_c ---
    get_csv_dims = lib.get_csv_dims_c
    get_csv_dims.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.int32), ctypes.c_int,
        ctypes.c_bool, ctypes.c_int,
        ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int)
    ]
    get_csv_dims.restype = None
    functions['get_csv_dims'] = get_csv_dims

    # --- read_csv_to_strings_c ---
    read_csv_to_strings = lib.read_csv_to_strings_c
    read_csv_to_strings.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.int32), ctypes.c_int,
        ctypes.c_bool, ctypes.c_int,
        np.ctypeslib.ndpointer(dtype=np.int32), np.ctypeslib.ndpointer(dtype=np.int32),
        ctypes.POINTER(ctypes.c_int)
    ]
    read_csv_to_strings.restype = None
    functions['read_csv_to_strings'] = read_csv_to_strings

    # --- read_integer_columns_c ---
    read_int = lib.read_integer_columns_c
    read_int.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.int32), ctypes.c_int, ctypes.c_int,
        np.ctypeslib.ndpointer(dtype=np.int32), ctypes.c_int,
        np.ctypeslib.ndpointer(dtype=np.int32), ctypes.POINTER(ctypes.c_int)
    ]
    read_int.restype = None
    functions['read_int'] = read_int
    
    # --- read_real_columns_c ---
    read_real = lib.read_real_columns_c
    read_real.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.int32), ctypes.c_int, ctypes.c_int,
        np.ctypeslib.ndpointer(dtype=np.int32), ctypes.c_int,
        np.ctypeslib.ndpointer(dtype=np.float64), ctypes.POINTER(ctypes.c_int)
    ]
    read_real.restype = None
    functions['read_real'] = read_real

    # --- read_logical_columns_c ---
    read_logical = lib.read_logical_columns_c
    read_logical.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.int32), ctypes.c_int, ctypes.c_int,
        np.ctypeslib.ndpointer(dtype=np.int32), ctypes.c_int,
        np.ctypeslib.ndpointer(dtype=np.bool_), ctypes.POINTER(ctypes.c_int)
    ]
    read_logical.restype = None
    functions['read_logical'] = read_logical
    
    # --- read_character_columns_c ---
    read_char = lib.read_character_columns_c
    read_char.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.int32), ctypes.c_int, ctypes.c_int,
        np.ctypeslib.ndpointer(dtype=np.int32), ctypes.c_int,
        np.ctypeslib.ndpointer(dtype=np.int32), ctypes.POINTER(ctypes.c_int)
    ]
    read_char.restype = None
    functions['read_char'] = read_char

    # --- read_complex_columns_c ---
    read_complex = lib.read_complex_columns_c
    read_complex.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.int32), ctypes.c_int, ctypes.c_int,
        np.ctypeslib.ndpointer(dtype=np.int32), ctypes.c_int,
        np.ctypeslib.ndpointer(dtype=np.complex128), ctypes.POINTER(ctypes.c_int)
    ]
    read_complex.restype = None
    functions['read_complex'] = read_complex
    
    return functions

# --- 3. Helper functions for testing ---

MAX_FIELD_LEN = 512 # This MUST match the Fortran parameter

def _py_string_to_ascii_array(s: str) -> np.ndarray:
    """Converts a Python string to a numpy array of int32 ASCII codes."""
    return np.array([ord(c) for c in s], dtype=np.int32)

def _ascii_to_string_matrix(ascii_flat, rows, cols):
    """Converts a flat ASCII buffer from Fortran into a 2D numpy string array."""
    reshaped = ascii_flat.reshape((rows, cols, MAX_FIELD_LEN), order='C')
    str_matrix = np.empty((rows, cols), dtype=object, order='F')
    for i in range(rows):
        for j in range(cols):
            null_pos = np.where(reshaped[i, j] == 0)[0]
            end_pos = null_pos[0] if len(null_pos) > 0 else MAX_FIELD_LEN
            str_matrix[i, j] = "".join(chr(c) for c in reshaped[i, j, :end_pos])
    return str_matrix.astype(str)

def _string_matrix_to_flat_ascii(str_matrix: np.ndarray):
    """Converts a 2D numpy string array to a flat C-style ASCII integer array."""
    rows, cols = str_matrix.shape
    flat_ascii = np.zeros(rows * cols * MAX_FIELD_LEN, dtype=np.int32)
    for i in range(rows):
        for j in range(cols):
            s = str_matrix[i, j]
            for k, char in enumerate(s):
                if k >= MAX_FIELD_LEN: break
                idx = i * cols * MAX_FIELD_LEN + j * MAX_FIELD_LEN + k
                flat_ascii[idx] = ord(char)
    return flat_ascii

# --- 4. Test Functions ---

def test_full_pipeline(fortran_funcs):
    """
    Runs a comprehensive test of the entire CSV reading pipeline.
    """
    print("=== Testing Full CSV Reading Pipeline ===")
    
    # Arrange: Create a dummy CSV file
    csv_content = (
        "ID,Status,Score,Name,Coords\n"
        "101,T,95.5,Alpha,(1.1, 2.2)\n"
        "102,F,88.0,Beta,(3.3, -4.4)\n"
        "103,TRUE,72.1,Gamma,(5.5, 6.6)\n"
        "104,FALSE,-5.0,Delta,(7.7, 8.8)\n"
    )
    test_filename = "python_pipeline_test.csv"
    with open(test_filename, "w") as f:
        f.write(csv_content)

    try:
        # Act 1: Get dimensions
        fname_ascii = _py_string_to_ascii_array(test_filename)
        num_rows, num_cols, status = ctypes.c_int(), ctypes.c_int(), ctypes.c_int()
        fortran_funcs['get_csv_dims'](fname_ascii, len(test_filename), True, ord(','),
                                      ctypes.byref(num_rows), ctypes.byref(num_cols), ctypes.byref(status))
        
        assert status.value == 0, f"get_csv_dims failed with status {status.value}"
        assert num_rows.value == 4, f"Expected 4 rows, got {num_rows.value}"
        assert num_cols.value == 5, f"Expected 5 cols, got {num_cols.value}"
        print("Step 1: get_csv_dims PASSED")

        # Act 2: Read all data as strings
        nr, nc = num_rows.value, num_cols.value
        # CORRECTED: Allocate a correctly typed and sized buffer for the header.
        header_buffer = np.zeros(nc * MAX_FIELD_LEN, dtype=np.int32)
        data_buffer = np.zeros(nr * nc * MAX_FIELD_LEN, dtype=np.int32)
        fortran_funcs['read_csv_to_strings'](fname_ascii, len(test_filename), True, ord(','),
                                             header_buffer, data_buffer, ctypes.byref(status))
        
        assert status.value == 0, f"read_csv_to_strings failed with status {status.value}"
        data_str_matrix = _ascii_to_string_matrix(data_buffer, nr, nc)
        assert data_str_matrix[1, 3] == "Beta"
        print("Step 2: read_csv_to_strings PASSED")
        
        flat_ascii_input = _string_matrix_to_flat_ascii(data_str_matrix)

        # Act 3 & Assert: Convert to specific types
        
        # Integers (Column 0)
        cols_to_read = np.array([1], dtype=np.int32) # Fortran 1-based index
        int_out = np.zeros(nr * 1, dtype=np.int32)
        fortran_funcs['read_int'](flat_ascii_input, nr, nc, cols_to_read, 1, int_out, ctypes.byref(status))
        assert status.value == 0
        np.testing.assert_array_equal(int_out, [101, 102, 103, 104])
        print("Step 3: read_integer_columns PASSED")

        # Reals (Column 2)
        cols_to_read = np.array([3], dtype=np.int32)
        real_out = np.zeros(nr * 1, dtype=np.float64)
        fortran_funcs['read_real'](flat_ascii_input, nr, nc, cols_to_read, 1, real_out, ctypes.byref(status))
        assert status.value == 0
        np.testing.assert_allclose(real_out, [95.5, 88.0, 72.1, -5.0])
        print("Step 4: read_real_columns PASSED")

        # Logicals (Column 1)
        cols_to_read = np.array([2], dtype=np.int32)
        logical_out = np.zeros(nr * 1, dtype=np.bool_)
        fortran_funcs['read_logical'](flat_ascii_input, nr, nc, cols_to_read, 1, logical_out, ctypes.byref(status))
        assert status.value == 0
        np.testing.assert_array_equal(logical_out, [True, False, True, False])
        print("Step 5: read_logical_columns PASSED")
        
        # Characters (Column 3)
        cols_to_read = np.array([4], dtype=np.int32)
        char_out_ascii = np.zeros(nr * 1 * MAX_FIELD_LEN, dtype=np.int32)
        fortran_funcs['read_char'](flat_ascii_input, nr, nc, cols_to_read, 1, char_out_ascii, ctypes.byref(status))
        assert status.value == 0
        char_out_str = _ascii_to_string_matrix(char_out_ascii, nr, 1).flatten()
        np.testing.assert_array_equal(char_out_str, ["Alpha", "Beta", "Gamma", "Delta"])
        print("Step 6: read_character_columns PASSED")

        # Complex (Column 4)
        cols_to_read = np.array([5], dtype=np.int32)
        complex_out = np.zeros(nr * 1, dtype=np.complex128)
        fortran_funcs['read_complex'](flat_ascii_input, nr, nc, cols_to_read, 1, complex_out, ctypes.byref(status))
        assert status.value == 0
        np.testing.assert_allclose(complex_out, [1.1+2.2j, 3.3-4.4j, 5.5+6.6j, 7.7+8.8j])
        print("Step 7: read_complex_columns PASSED")

    finally:
        # Cleanup
        if os.path.exists(test_filename):
            os.remove(test_filename)

# --- 5. Main Execution Block ---

def main():
    """Load the library and run all tests."""
    print("=================================================")
    print("      PYTHON WRAPPER TESTS FOR FORTRAN CSV READER")
    print("=================================================")
    try:
        fortran_lib = _load_fortran_library()
        fortran_funcs = setup_fortran_functions(fortran_lib)
        
        test_full_pipeline(fortran_funcs)
        
        print("\n=================================================")
        print("           ALL PYTHON TESTS COMPLETED")
        print("=================================================")
        print("If you see this message, all Python interface tests passed! ✓")

    except (ImportError, OSError, AssertionError) as e:
        print("\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        print("                  A TEST FAILED")
        print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        print(f"ERROR: {e}")

if __name__ == "__main__":
    main()
    