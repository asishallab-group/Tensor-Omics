import numpy as np
import sys
import os

# Add parent directory to path to import tensoromics_functions
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from tensoromics_functions import (tox_serialize_char_nd, tox_serialize_int_nd, tox_serialize_real_nd,
                                   tox_deserialize_char_nd, tox_deserialize_int_nd, tox_deserialize_real_nd,
                                   tox_serialize_logical_nd, tox_serialize_complex_nd,
                                   tox_deserialize_logical_nd, tox_deserialize_complex_nd)

# Tests for integer
def int_test():
    array = np.array([1, 2, 3, 4, 5], dtype=np.int32, order='F')
    filename = "test_int_1d.bin"
    tox_serialize_int_nd(array, filename)
    print(f"Serialized array to {filename}")
    res = tox_deserialize_int_nd(filename)
    assert np.array_equal(res, array), "Deserialized array does not match original"
    print(res)

    array_2d = np.array([[1, 2, 3], [4, 5, 6]], dtype=np.int32, order='F')
    filename_2d = "test_int_2d.bin"
    tox_serialize_int_nd(array_2d, filename_2d)
    print(f"Serialized 2D array to {filename_2d}")
    res_2d = tox_deserialize_int_nd(filename_2d)
    print(res_2d)
    assert np.array_equal(res_2d, array_2d), "Deserialized 2D array does not match original"

    array_3d = np.array([[[1, 2], [3, 4]], [[5, 6], [7, 8]]], dtype=np.int32, order='F')
    filename_3d = "test_int_3d.bin"
    tox_serialize_int_nd(array_3d, filename_3d)
    print(f"Serialized 3D array to {filename_3d}")
    res_3d = tox_deserialize_int_nd(filename_3d)
    print(res_3d)
    assert np.array_equal(res_3d, array_3d), "Deserialized 3D array does not match original"
    
# Tests for real
def real_test():
    array = np.array([1.5, 2.3, 3.2, 4.0, 5.0], dtype=np.float64, order='F')
    filename = "test_real_1d.bin"
    tox_serialize_real_nd(array, filename)
    print(f"Serialized array to {filename}")
    res = tox_deserialize_real_nd(filename)
    assert np.array_equal(res, array), "Deserialized array does not match original"
    print(res)

    array_2d = np.asfortranarray([[1.0, 2.0, 3.0], [4.0, 5.7, 6.0]], dtype=np.float64)
    filename_2d = "test_real_2d.bin"
    tox_serialize_real_nd(array_2d, filename_2d)
    print(f"Serialized 2D array to {filename_2d}")
    res_2d = tox_deserialize_real_nd(filename_2d)
    print(res_2d)
    assert np.array_equal(res_2d, array_2d), "Deserialized 2D array does not match original"

    array_3d = np.array([[[1.0, 2.0], [3.3, 4.0]], [[5.0, 6.8], [7.0, 8.0]]], dtype=np.float64, order='F')
    filename_3d = "test_real_3d.bin"
    tox_serialize_real_nd(array_3d, filename_3d)
    print(f"Serialized 3D array to {filename_3d}")
    res_3d = tox_deserialize_real_nd(filename_3d)
    print(res_3d)
    assert np.array_equal(res_3d, array_3d), "Deserialized 3D array does not match original"

    empty_array = np.array([], dtype=np.float64, order='F')
    empty_filename = "test_real_empty.bin"
    tox_serialize_real_nd(empty_array, empty_filename)
    print(f"Serialized empty array to {empty_filename}")
    res_empty = tox_deserialize_real_nd(empty_filename)
    assert res_empty.size == 0, "Deserialized empty array should be empty"

# Tests for chars
def char_test():
    array = np.asfortranarray(["hello", "world"], dtype='U5')
    filename = "test_char_1d.bin"
    tox_serialize_char_nd(array, filename)
    print(f"Serialized array to {filename}")
    res = tox_deserialize_char_nd(filename)
    print(res)
    assert np.array_equal(res, array), "Deserialized array does not match original"
    

    array_2d = np.asfortranarray([["foo", "bar"], ["baz", "qux"]], dtype='U5')
    filename_2d = "test_char_2d.bin"
    tox_serialize_char_nd(array_2d, filename_2d)
    print(f"Serialized 2D array to {filename_2d}")
    res_2d = tox_deserialize_char_nd(filename_2d)
    print(res_2d)
    assert np.array_equal(res_2d, array_2d), "Deserialized 2D array does not match original"

    array_3d = np.array([[["abb", "bbbbbbb"], ["cfs", "d"]], [["e", ""], ["g", "h"]]], dtype='U5', order='F')
    filename_3d = "test_char_3d.bin"
    tox_serialize_char_nd(array_3d, filename_3d)
    print(f"Serialized 3D array to {filename_3d}")
    res_3d = tox_deserialize_char_nd(filename_3d)
    print(res_3d)
    assert np.array_equal(res_3d, array_3d), "Deserialized 3D array does not match original"

# Tests for logical
def logical_test():
    # 1D logical array
    array = np.array([True, False, True, False, True], dtype=np.bool_, order='F')
    filename = "test_logical_1d.bin"
    tox_serialize_logical_nd(array, filename)
    print(f"Serialized logical array to {filename}")
    res = tox_deserialize_logical_nd(filename)
    print("Original:", array)
    print("Deserialized:", res)
    assert np.array_equal(res, array), "Deserialized logical array does not match original"
    
    # 2D logical array
    array_2d = np.array([[True, False, True], [False, True, False]], dtype=np.bool_, order='F')
    filename_2d = "test_logical_2d.bin"
    tox_serialize_logical_nd(array_2d, filename_2d)
    print(f"Serialized 2D logical array to {filename_2d}")
    res_2d = tox_deserialize_logical_nd(filename_2d)
    print("Original 2D:")
    print(array_2d)
    print("Deserialized 2D:")
    print(res_2d)
    assert np.array_equal(res_2d, array_2d), "Deserialized 2D logical array does not match original"

    # 3D logical array
    array_3d = np.array([[[True, False], [False, True]], [[False, True], [True, False]]], dtype=np.bool_, order='F')
    filename_3d = "test_logical_3d.bin"
    tox_serialize_logical_nd(array_3d, filename_3d)
    print(f"Serialized 3D logical array to {filename_3d}")
    res_3d = tox_deserialize_logical_nd(filename_3d)
    print("Original 3D shape:", array_3d.shape)
    print("Deserialized 3D shape:", res_3d.shape)
    assert np.array_equal(res_3d, array_3d), "Deserialized 3D logical array does not match original"
    
    # Edge cases: all True and all False
    all_true = np.array([True, True, True, True], dtype=np.bool_, order='F')
    all_false = np.array([False, False, False, False], dtype=np.bool_, order='F')
    
    tox_serialize_logical_nd(all_true, "test_logical_all_true.bin")
    tox_serialize_logical_nd(all_false, "test_logical_all_false.bin")
    
    res_all_true = tox_deserialize_logical_nd("test_logical_all_true.bin")
    res_all_false = tox_deserialize_logical_nd("test_logical_all_false.bin")
    
    assert np.array_equal(res_all_true, all_true), "All True array mismatch"
    assert np.array_equal(res_all_false, all_false), "All False array mismatch"
    print("All True/False edge cases passed")

# Tests for complex
def complex_test():
    # 1D complex array
    array = np.array([1+2j, 3+4j, 5+6j, 7+8j], dtype=np.complex128, order='F')
    filename = "test_complex_1d.bin"
    tox_serialize_complex_nd(array, filename)
    print(f"Serialized complex array to {filename}")
    res = tox_deserialize_complex_nd(filename)
    print("Original:", array)
    print("Deserialized:", res)
    assert np.allclose(res, array), "Deserialized complex array does not match original"
    
    # 2D complex array
    array_2d = np.array([[1+1j, 2+2j], [3+3j, 4+4j]], dtype=np.complex128, order='F')
    filename_2d = "test_complex_2d.bin"
    tox_serialize_complex_nd(array_2d, filename_2d)
    print(f"Serialized 2D complex array to {filename_2d}")
    res_2d = tox_deserialize_complex_nd(filename_2d)
    print("Original 2D:")
    print(array_2d)
    print("Deserialized 2D:")
    print(res_2d)
    assert np.allclose(res_2d, array_2d), "Deserialized 2D complex array does not match original"

    # 3D complex array with various values
    array_3d = np.array([[[1+0j, 0+1j], [1+1j, 0+0j]], 
                         [[2+0j, 0+2j], [2+2j, 0+0j]]], dtype=np.complex128, order='F')
    filename_3d = "test_complex_3d.bin"
    tox_serialize_complex_nd(array_3d, filename_3d)
    print(f"Serialized 3D complex array to {filename_3d}")
    res_3d = tox_deserialize_complex_nd(filename_3d)
    print("Original 3D shape:", array_3d.shape)
    print("Deserialized 3D shape:", res_3d.shape)
    assert np.allclose(res_3d, array_3d), "Deserialized 3D complex array does not match original"
    
    # Edge cases: pure real, pure imaginary, mixed
    pure_real = np.array([1+0j, 2+0j, 3+0j], dtype=np.complex128, order='F')
    pure_imag = np.array([0+1j, 0+2j, 0+3j], dtype=np.complex128, order='F')
    mixed = np.array([1+2j, -1-2j, 3-4j, -3+4j], dtype=np.complex128, order='F')
    
    tox_serialize_complex_nd(pure_real, "test_complex_pure_real.bin")
    tox_serialize_complex_nd(pure_imag, "test_complex_pure_imag.bin")
    tox_serialize_complex_nd(mixed, "test_complex_mixed.bin")
    
    res_pure_real = tox_deserialize_complex_nd("test_complex_pure_real.bin")
    res_pure_imag = tox_deserialize_complex_nd("test_complex_pure_imag.bin")
    res_mixed = tox_deserialize_complex_nd("test_complex_mixed.bin")
    
    assert np.allclose(res_pure_real, pure_real), "Pure real complex array mismatch"
    assert np.allclose(res_pure_imag, pure_imag), "Pure imaginary complex array mismatch"
    assert np.allclose(res_mixed, mixed), "Mixed complex array mismatch"
    print("Complex edge cases passed")

print("============ INT TEST ============")
int_test()
print("============ REAL TEST ============")
real_test()
print("============ CHAR TEST ============")
char_test()
print("============ LOGICAL TEST ============")
logical_test()
print("============ COMPLEX TEST ============")
complex_test()