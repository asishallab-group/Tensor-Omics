import numpy as np
import sys
import os

# Add parent directory to path to import tensoromics_functions
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from tensoromics_functions import tox_serialize_char_nd, tox_serialize_int_nd, tox_serialize_real_nd 
from tensoromics_functions import tox_deserialize_char_nd, tox_deserialize_int_nd, tox_deserialize_real_nd

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

print("============ INT TEST ============")
int_test()
print("============ REAL TEST ============")
real_test()
print("============ CHAR TEST ============")
char_test()