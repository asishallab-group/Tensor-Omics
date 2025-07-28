import ctypes
import numpy as np
import os
import platform

# Determine library extension based on OS
if platform.system() == "Windows":
    lib_ext = ".dll"
elif platform.system() == "Darwin": # macOS
    lib_ext = ".dylib"
else: # Linux
    lib_ext = ".so"

# Find and load the Fortran shared library
# Assumes this script is in project_root/python/
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# The library name is based on the fpm.toml 'name' field
lib_path = os.path.join(project_root, "build", "fpm_build", f"libfortran_csv_reader{lib_ext}")

try:
    fortran_lib = ctypes.CDLL(lib_path)
except OSError as e:
    print(f"Error loading library: {e}")
    print(f"Searched for library at: {lib_path}")
    print("Please ensure the Fortran code has been built with 'fpm build --profile release'")
    exit(1)

# --- Helper Functions to convert between Python and Fortran C types ---

def _py_string_to_ascii_array(s: str) -> np.ndarray:
    """Converts a Python string to a null-terminated numpy array of int32 ASCII codes."""
    return np.array([ord(c) for c in s] + [0], dtype=np.int32)

def _ascii_array_to_py_string(arr: np.ndarray) -> str:
    """Converts a null-terminated numpy array of int32 ASCII codes to a Python string."""
    return "".join(chr(c) for c in arr if c > 0)

# --- Main Python API ---

def read_csv_to_strings(filename: str, has_header: bool = True, delimiter: str = ',') -> tuple:
    """
    Reads a CSV file into a numpy array of strings for the header and a 2D array for the data.

    Returns:
        A tuple containing (header_array, data_array).
        header_array is 1D, data_array is 2D. Both are Fortran-ordered numpy arrays.
    """
    # 1. Get dimensions from Fortran
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
        raise IOError(f"Fortran failed to get dimensions for '{filename}' with status {status.value}")

    # 2. Allocate buffers in Python based on the dimensions
    MAX_FIELD_LEN = 512 # This must match the Fortran parameter
    header_buffer = np.zeros(num_cols.value * MAX_FIELD_LEN, dtype=np.int32)
    data_buffer = np.zeros(num_rows.value * num_cols.value * MAX_FIELD_LEN, dtype=np.int32)

    # 3. Call Fortran again to fill the buffers
    fortran_lib.read_csv_to_strings_c.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.int32), ctypes.c_int,
        ctypes.c_bool, ctypes.c_int,
        np.ctypeslib.ndpointer(dtype=np.int32), np.ctypeslib.ndpointer(dtype=np.int32),
        ctypes.POINTER(ctypes.c_int)
    ]
    fortran_lib.read_csv_to_strings_c(fname_ascii, len(filename), has_header, ord(delimiter),
                                      header_buffer, data_buffer, ctypes.byref(status))
    if status.value != 0:
        raise IOError(f"Fortran failed to read file '{filename}' with status {status.value}")

    # 4. Convert the flat ASCII buffers back into numpy string arrays
    header = np.array([_ascii_array_to_py_string(s) for s in np.split(header_buffer, num_cols.value)])
    
    data_flat = [_ascii_array_to_py_string(s) for s in np.split(data_buffer, num_rows.value * num_cols.value)]
    data = np.array(data_flat).reshape((num_rows.value, num_cols.value), order='F')
    
    return header, data

# --- Python Test Script ---
if __name__ == '__main__':
    # Create a dummy CSV file to test with
    csv_content = "ID,Name,Value\n1,Alpha,10.5\n2,Beta,20.5"
    test_filename = "test_data.csv"
    with open(test_filename, "w") as f:
        f.write(csv_content)

    print(f"--- Testing Python Wrapper for '{test_filename}' ---")
    
    try:
        # 1. Read the CSV into string arrays
        header, data = read_csv_to_strings(test_filename, has_header=True)
        
        print("Header read from Fortran:")
        print(header)
        assert np.array_equal(header, np.array(["ID", "Name", "Value"]))
        
        print("\nData read from Fortran:")
        print(data)
        assert data[1, 1] == "Beta"

        # 2. Convert integer column
        # In a real application, you would create functions like get_integer_columns etc.
        # For this test, we demonstrate the conversion directly
        id_col = data[:, 0].astype(np.int32)
        print("\nID column converted to integers:")
        print(id_col)
        assert np.array_equal(id_col, np.array([1, 2]))

        print("\nAll Python wrapper tests PASSED!")

    finally:
        # Clean up the dummy file
        if os.path.exists(test_filename):
            os.remove(test_filename)
            print(f"\nCleaned up {test_filename}")

