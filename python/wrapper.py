import ctypes
import numpy as np
import os
import platform
from collections import defaultdict

# --- 1. Library Loading ---

def _load_fortran_library():
    """
    Finds and loads the compiled Fortran shared library (.so, .dll, .dylib)
    based on the operating system by searching the build directory.
    """
    lib_name = "libfortran_csv_reader"
    
    if platform.system() == "Windows":
        lib_ext = ".dll"
    elif platform.system() == "Darwin": # macOS
        lib_ext = ".dylib"
    else: # Linux
        lib_ext = ".so"
        
    # Assumes this script is in project_root/python/
    project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    build_dir = os.path.join(project_root, "build")
    
    lib_path = None
    # Search for the library file within the build directory
    for root, dirs, files in os.walk(build_dir):
        if f"{lib_name}{lib_ext}" in files:
            lib_path = os.path.join(root, f"{lib_name}{lib_ext}")
            break

    if not lib_path or not os.path.exists(lib_path):
        raise OSError(f"Could not find {lib_name}{lib_ext} in the build directory '{build_dir}'. "
                      "Please build the Fortran project first using 'fpm build --profile release'.")

    return ctypes.CDLL(lib_path)

try:
    fortran_lib = _load_fortran_library()
except (OSError, NameError) as e:
    print(f"Error: {e}")
    fortran_lib = None

# --- 2. Helper Functions for Data Conversion ---

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

# --- 3. Main Python API Function ---

def read_csv(filename: str, dtypes: dict, has_header: bool = True, delimiter: str = ','):
    """
    Reads specified columns from a CSV file into typed numpy arrays.
    """
    if not fortran_lib:
        raise ImportError("Fortran library is not loaded.")

    fname_ascii = _py_string_to_ascii_array(filename)
    num_rows = ctypes.c_int()
    num_cols = ctypes.c_int()
    status = ctypes.c_int()
    
    fortran_lib.get_csv_dims_c.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.int32), ctypes.c_int,
        ctypes.c_bool, ctypes.c_int,
        ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int)
    ]
    fortran_lib.get_csv_dims_c(fname_ascii, len(filename), has_header, ord(delimiter),
                               ctypes.byref(num_rows), ctypes.byref(num_cols), ctypes.byref(status))
    if status.value != 0:
        raise IOError(f"Fortran failed to get dimensions for '{filename}' (status: {status.value})")
    
    nr, nc = num_rows.value, num_cols.value
    if nr == 0:
        return {}

    header_buffer = np.zeros(nc * MAX_FIELD_LEN, dtype=np.int32)
    data_buffer = np.zeros(nr * nc * MAX_FIELD_LEN, dtype=np.int32)

    fortran_lib.read_csv_to_strings_c.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.int32), ctypes.c_int,
        ctypes.c_bool, ctypes.c_int,
        np.ctypeslib.ndpointer(dtype=np.int32), np.ctypeslib.ndpointer(dtype=np.int32),
        ctypes.POINTER(ctypes.c_int)
    ]
    fortran_lib.read_csv_to_strings_c(fname_ascii, len(filename), has_header, ord(delimiter),
                                      header_buffer, data_buffer, ctypes.byref(status))
    if status.value != 0:
        raise IOError(f"Fortran failed to read file '{filename}' (status: {status.value})")

    data_str_matrix = _ascii_to_string_matrix(data_buffer, nr, nc)
    
    cols_by_type = defaultdict(list)
    for col_idx, dtype in dtypes.items():
        if 0 <= col_idx < nc:
            cols_by_type[dtype.lower()].append(col_idx)

    results = {}
    
    type_map = {
        'int': ('read_integer_columns_c', np.int32),
        'real': ('read_real_columns_c', np.float64),
        'logical': ('read_logical_columns_c', np.bool_),
        'char': ('read_character_columns_c', None),
        'complex': ('read_complex_columns_c', np.complex128)
    }

    flat_ascii_input = _string_matrix_to_flat_ascii(data_str_matrix)

    for dtype, cols in cols_by_type.items():
        if dtype not in type_map:
            print(f"Warning: Unknown dtype '{dtype}' requested, skipping.")
            continue
            
        func_name, np_type = type_map[dtype]
        cols_to_read = np.array(cols, dtype=np.int32) + 1
        num_cols_to_read = len(cols)

        if dtype == 'char':
            result_matrix = data_str_matrix[:, cols]
        else:
            output_buffer = np.zeros(nr * num_cols_to_read, dtype=np_type)
            
            func = getattr(fortran_lib, func_name)
            func.argtypes = [
                np.ctypeslib.ndpointer(dtype=np.int32), ctypes.c_int, ctypes.c_int,
                np.ctypeslib.ndpointer(dtype=np.int32), ctypes.c_int,
                np.ctypeslib.ndpointer(dtype=np_type), ctypes.POINTER(ctypes.c_int)
            ]
            
            func(flat_ascii_input, nr, nc,
                 cols_to_read, num_cols_to_read,
                 output_buffer, ctypes.byref(status))
            
            if status.value != 0:
                print(f"Warning: Fortran conversion for dtype '{dtype}' returned status {status.value}")

            result_matrix = output_buffer.reshape((nr, num_cols_to_read), order='F')

        for i, col_idx in enumerate(cols):
            results[col_idx] = result_matrix[:, i]

    return results

# --- 4. Python Test Script ---
if __name__ == '__main__':
    csv_content = (
        "ID,Status,Score,Name,Coords\n"
        "101,T,95.5,Alpha,(1.1, 2.2)\n"
        "102,F,88.0,Beta,(3.3, -4.4)\n"
        "103,TRUE,72.1,Gamma,(5.5, 6.6)\n"
        "104,FALSE,-5.0,Delta,(7.7, 8.8)\n"
    )
    test_filename = "python_test.csv"
    with open(test_filename, "w") as f:
        f.write(csv_content)

    print(f"--- Testing Python Wrapper for '{test_filename}' ---")
    
    try:
        desired_types = {
            0: 'int', 1: 'logical', 2: 'real', 3: 'char', 4: 'complex'
        }
        
        data = read_csv(test_filename, dtypes=desired_types)
        
        print("\nSuccessfully read and converted data:")
        for col, arr in data.items():
            print(f"  Column {col} (type: {desired_types[col]}): {arr}")

        assert np.array_equal(data[0], np.array([101, 102, 103, 104], dtype=np.int32))
        assert np.array_equal(data[1], np.array([True, False, True, False], dtype=np.bool_))
        assert np.allclose(data[2], np.array([95.5, 88.0, 72.1, -5.0], dtype=np.float64))
        assert np.array_equal(data[3], np.array(["Alpha", "Beta", "Gamma", "Delta"]))
        assert np.allclose(data[4], np.array([1.1+2.2j, 3.3-4.4j, 5.5+6.6j, 7.7+8.8j], dtype=np.complex128))

        print("\nAll Python wrapper tests PASSED!")

    except Exception as e:
        print(f"\nAn error occurred during the test: {e}")
        
    finally:
        if os.path.exists(test_filename):
            os.remove(test_filename)
            print(f"Cleaned up {test_filename}")
