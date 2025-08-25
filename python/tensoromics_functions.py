# tensoromics_functions.py

import ctypes
import numpy as np
import os
from pathlib import Path

# --- Library Loading ---
# Assumes the library is in a 'build' directory parallel to the 'python' directory
try:
    lib_path = Path(__file__).parent.parent / "build" / "libtensor-omics.so"
    fortran_lib = ctypes.CDLL(str(lib_path))
except OSError as e:
    raise ImportError(f"Could not load the Fortran shared library at '{lib_path}'. "
                      "Please ensure the project is compiled. Original error: {e}")

# This constant must match the Fortran parameter
MAX_FIELD_LEN = 512

# --- Error Handling ---
def _handle_fortran_errors(status_code: int, context: str):
    """Translates Fortran status codes into Python exceptions."""
    errors = {
        1: "No data to process.",
        2: "Column index out of bounds.",
        3: "At least one data conversion error occurred.",
        4: "The specified file is empty.",
        5: "File has a header but no data rows.",
        10: "File not found or cannot be opened.",
    }
    if status_code != 0:
        message = errors.get(status_code, f"An unknown Fortran error occurred (code: {status_code}).")
        raise ValueError(f"Error in {context}: {message}")

# --- Utilities ---
def _string_matrix_to_flat_ascii(str_matrix: np.ndarray) -> np.ndarray:
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

def _ascii_to_string_matrix(ascii_flat: np.ndarray, rows: int, cols: int) -> np.ndarray:
    """Converts a flat ASCII buffer from Fortran into a 2D numpy string array."""
    reshaped = ascii_flat.reshape((rows, cols, MAX_FIELD_LEN))
    str_matrix = np.empty((rows, cols), dtype=object)
    for i in range(rows):
        for j in range(cols):
            null_pos = np.where(reshaped[i, j] == 0)[0]
            end_pos = null_pos[0] if len(null_pos) > 0 else MAX_FIELD_LEN
            str_matrix[i, j] = "".join(chr(c) for c in reshaped[i, j, :end_pos])
    
    # Set the final array to read-only as required
    final_matrix = np.array(str_matrix, dtype=str)
    final_matrix.flags.writeable = False
    # NOTE: Returned NumPy arrays are read-only for safety.
    # If you need to modify them (e.g., for plotting), use .copy().
    return final_matrix

# --- User-Facing API Functions ---

def read_csv_dimensions(filepath: str, has_header: bool = True, delimiter: str = ',') -> tuple[int, int]:
    """
    Gets the dimensions (rows, columns) of a CSV file without reading its data.

    Parameters
    ----------
    filepath : str
        Path to the CSV file.
    has_header : bool, optional
        Whether the file has a header row, by default True.
    delimiter : str, optional
        The column delimiter, by default ','.

    Returns
    -------
    tuple[int, int]
        A tuple containing (number_of_data_rows, number_of_columns).
    """
    if not Path(filepath).exists():
        raise FileNotFoundError(f"File not found at '{filepath}'")
    if len(delimiter) != 1:
        raise ValueError("Delimiter must be a single character.")

    fname_ascii = np.array([ord(c) for c in filepath], dtype=np.int32)
    num_rows, num_cols, status = ctypes.c_int(0), ctypes.c_int(0), ctypes.c_int(0)

    f = fortran_lib.get_csv_dims_c
    f.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
        ctypes.c_int, ctypes.c_bool, ctypes.c_int,
        ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int)
    ]
    f.restype = None
    
    f(fname_ascii, len(filepath), has_header, ord(delimiter),
      ctypes.byref(num_rows), ctypes.byref(num_cols), ctypes.byref(status))
    
    _handle_fortran_errors(status.value, "read_csv_dimensions")
    return num_rows.value, num_cols.value

def read_csv_as_strings(filepath: str, has_header: bool = True, delimiter: str = ',') -> tuple[np.ndarray | None, np.ndarray]:
    """
    Reads a CSV file into a NumPy array of strings.

    Parameters
    ----------
    filepath : str
        Path to the CSV file.
    has_header : bool, optional
        Whether the file has a header row, by default True.
    delimiter : str, optional
        The column delimiter, by default ','.

    Returns
    -------
    tuple[np.ndarray | None, np.ndarray]
        A tuple containing (header, data). Header is None if has_header is False.
        Data is a 2D NumPy array of strings. Returned arrays are read-only.
    """
    num_rows, num_cols = read_csv_dimensions(filepath, has_header, delimiter)
    fname_ascii = np.array([ord(c) for c in filepath], dtype=np.int32)
    status = ctypes.c_int(0)

    header_buffer = np.zeros(num_cols * MAX_FIELD_LEN, dtype=np.int32) if has_header else np.array([], dtype=np.int32)
    data_buffer = np.zeros(num_rows * num_cols * MAX_FIELD_LEN, dtype=np.int32)
    
    f = fortran_lib.read_csv_to_strings_c
    f.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.int32), ctypes.c_int,
        ctypes.c_bool, ctypes.c_int,
        np.ctypeslib.ndpointer(dtype=np.int32), np.ctypeslib.ndpointer(dtype=np.int32),
        ctypes.POINTER(ctypes.c_int)
    ]
    f.restype = None
    
    f(fname_ascii, len(filepath), has_header, ord(delimiter),
      header_buffer, data_buffer, ctypes.byref(status))
    
    _handle_fortran_errors(status.value, "read_csv_as_strings")
    
    data_out = _ascii_to_string_matrix(data_buffer, num_rows, num_cols)
    header_out = _ascii_to_string_matrix(header_buffer, 1, num_cols).flatten() if has_header else None
    
    if header_out is not None:
        header_out.flags.writeable = False
        # NOTE: Returned NumPy arrays are read-only for safety.
        # If you need to modify them (e.g., for plotting), use .copy().
    
    return header_out, data_out

def _read_typed_columns(string_data: np.ndarray, columns: list[int], dtype: np.dtype, fortran_func_name: str) -> np.ndarray:
    """Generic internal function to convert string data to a numeric type."""
    if not isinstance(string_data, np.ndarray) or string_data.ndim != 2:
        raise TypeError("string_data must be a 2D NumPy array.")
    if not all(isinstance(c, int) and c > 0 for c in columns):
        raise ValueError("columns must be a list of 1-based positive integers.")

    num_rows, num_cols_in = string_data.shape
    num_cols_to_read = len(columns)
    
    flat_ascii_input = _string_matrix_to_flat_ascii(string_data)
    cols_to_read_f = np.array(columns, dtype=np.int32)
    output_buffer = np.zeros(num_rows * num_cols_to_read, dtype=dtype)
    status = ctypes.c_int(0)
    
    f = getattr(fortran_lib, fortran_func_name)
    f.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.int32), ctypes.c_int, ctypes.c_int,
        np.ctypeslib.ndpointer(dtype=np.int32), ctypes.c_int,
        np.ctypeslib.ndpointer(dtype=dtype), ctypes.POINTER(ctypes.c_int)
    ]
    f.restype = None

    f(flat_ascii_input, num_rows, num_cols_in, cols_to_read_f, num_cols_to_read, output_buffer, ctypes.byref(status))
    _handle_fortran_errors(status.value, fortran_func_name)
    
    # Reshape to Fortran order, then transpose to get C order semantics
    reshaped_output = output_buffer.reshape((num_cols_to_read, num_rows)).T
    
    reshaped_output.flags.writeable = False
    # NOTE: Returned NumPy arrays are read-only for safety.
    # If you need to modify them (e.g., for plotting), use .copy().
    return reshaped_output

def read_integer_columns(string_data: np.ndarray, columns: list[int]) -> np.ndarray:
    """Extracts and converts specified columns to integers."""
    return _read_typed_columns(string_data, columns, np.int32, 'read_integer_columns_c')

def read_real_columns(string_data: np.ndarray, columns: list[int]) -> np.ndarray:
    """Extracts and converts specified columns to double precision reals."""
    return _read_typed_columns(string_data, columns, np.float64, 'read_real_columns_c')

def read_logical_columns(string_data: np.ndarray, columns: list[int]) -> np.ndarray:
    """Extracts and converts specified columns to booleans."""
    return _read_typed_columns(string_data, columns, np.bool_, 'read_logical_columns_c')
    
def read_complex_columns(string_data: np.ndarray, columns: list[int]) -> np.ndarray:
    """Extracts and converts specified columns to double precision complex numbers."""
    return _read_typed_columns(string_data, columns, np.complex128, 'read_complex_columns_c')

def read_character_columns(string_data: np.ndarray, columns: list[int]) -> np.ndarray:
    """Extracts specified string columns without type conversion."""
    if not isinstance(string_data, np.ndarray) or string_data.ndim != 2:
        raise TypeError("string_data must be a 2D NumPy array.")
    if not all(isinstance(c, int) and c > 0 for c in columns):
        raise ValueError("columns must be a list of 1-based positive integers.")

    # Fortran columns are 1-based, Python is 0-based.
    py_indices = [c - 1 for c in columns]
    output_array = string_data[:, py_indices]
    
    output_array.flags.writeable = False
    # NOTE: Returned NumPy arrays are read-only for safety.
    # If you need to modify them (e.g., for plotting), use .copy().
    return output_array