import sys
import os
import numpy as np
import ctypes

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
# Load library
dll_path = os.path.abspath("build/libtensor-omics.so")
ctypes.CDLL("libgomp.so.1", mode=ctypes.RTLD_GLOBAL)
ctypes.CDLL("libzip.so", mode=ctypes.RTLD_GLOBAL)
ctypes.CDLL("libxxhash.so", mode=ctypes.RTLD_GLOBAL)
lib = ctypes.CDLL(dll_path)


from tensoromics_functions_tox_data import (
    read_expression_vectors_tsv,
    read_orthofinder_file,
    filter_unassigned_genes,
    read_gene_ids_from_tsv_file,
    validate_all_data,
    validate_data_structure,
    validate_expression_data,
    validate_family_centroids,
    validate_string_array_uniqueness,
    validate_gene_to_family_mapping,
    validate_shift_vectors,
    save_tox_data,
    read_tox_data,
    create_zip_archive,
    extract_zip_archive
)
from tensoromics_functions import (
    tox_group_centroid,
    tox_compute_shift_vector_field,
    tox_serialize_int_nd,
    tox_deserialize_int_nd,
    tox_serialize_real_nd,
    tox_deserialize_real_nd,
    tox_serialize_char_nd,
    tox_deserialize_char_nd
)

# ---- Example: replicate Fortran test logic in Python ----

# Define your file lists (replace with your actual file paths)
file = ["material/kallisto_sex_data_no_na.tsv"]
# Parameters
n_genes = 88327
n_families = 15512
gene_len = 32
family_len = 32
n_samples = 67
value_cols = np.arange(2, 68)

# Allocate result matrices
kallisto_expr = np.zeros((n_samples, n_genes), dtype=np.float64, order='F')

# Read gene IDs from first file
print("Reading gene IDs...")
gene_ids = read_gene_ids_from_tsv_file(file, n_genes, gene_len, n_header_rows=1, gene_col=1)
print("gene_ids sample:", gene_ids[:5])
print("gene_ids type:", type(gene_ids))

# Read 6-replicate files
print("Reading 6-replicate files...")

kallisto_expr = read_expression_vectors_tsv(
    file_list=file,
    gene_ids=gene_ids,  # Now accepts numpy array
    n_samples=67,
    n_header_rows=1,
    gene_col=1,
    value_cols=value_cols,
    delimiter="\t"
)
print("expression shape:", kallisto_expr.shape)

# Read family mapping
print("Reading family file...")
family_result = read_orthofinder_file("material/Orthogroups.tsv", gene_ids, family_len, n_families)
family_ids = family_result['family_ids']
gene_to_fam = family_result['gene_to_fam']
print("family_ids sample:", family_ids[:5])
print("family_ids type:", type(family_ids))
print("gene_to_fam sample:", gene_to_fam[:10])

# Filter out genes without family assignments
print("Filtering unassigned genes...")
filter_result = filter_unassigned_genes(gene_to_fam)
mask = np.array(filter_result['mask'], dtype=bool)
n_genes_kept = filter_result['n_genes_kept']

print("n_genes_kept:", n_genes_kept)
print("mask sample:", mask[:10])

print("First gene expression: ", kallisto_expr[:, 0])

# Filter arrays using mask
filtered_gene_ids = np.array(gene_ids)[mask]
filtered_kallisto_expr = kallisto_expr[:, mask]
filtered_gene_to_fam = np.array(gene_to_fam)[mask]


print("filtered_gene_ids sample:", filtered_gene_ids[:5])
print("filtered_kallisto_expr shape:", filtered_kallisto_expr.shape)
print("filtered_gene_to_fam sample:", filtered_gene_to_fam[:10])
validate_string_array_uniqueness(gene_ids)
validate_string_array_uniqueness(family_ids)
validate_expression_data(kallisto_expr, True)

ortholog_set = np.array([True for i in range(n_genes_kept)])

centroids = tox_group_centroid(filtered_kallisto_expr, filtered_gene_to_fam, n_families, "all", ortholog_set)

print("Validating centroids...")
validate_family_centroids(centroids)
print("Centroids validated")

shift_vectors_result = tox_compute_shift_vector_field(filtered_kallisto_expr, centroids, filtered_gene_to_fam)
shift_vectors = shift_vectors_result

print("Validating shift vectors...")
validate_shift_vectors(shift_vectors, filtered_kallisto_expr, centroids, filtered_gene_to_fam, 
                      n_genes_kept, n_samples, n_families)
print("Shift vectors validated")

validate_gene_to_family_mapping(filtered_gene_to_fam, n_families)

print("Number of genes kept after filtering:", n_genes_kept)
print("gene_family_size:", len(filtered_gene_to_fam))
print("filtered gene ids size:", len(filtered_gene_ids))

print("Validating data structure...")
validate_data_structure(n_genes_kept, n_families, n_samples, filtered_gene_ids, 
                        family_ids, filtered_gene_to_fam, filtered_kallisto_expr, centroids, shift_vectors)
print("Data structure validated")

print("Validating all data...")
validate_all_data(n_genes_kept, n_families, n_samples, filtered_gene_ids, family_ids, 
                  filtered_gene_to_fam, filtered_kallisto_expr, centroids, shift_vectors)
print("All data validated")

print("Testing save function...")
# Clean up any existing test archives
for archive_name in ["test_archive_1_py.zip", "test_archive_2_py.zip", "test_archive_3_py.zip", "test_archive_4_py.zip"]:
    if os.path.exists(archive_name):
        os.remove(archive_name)
        print(f"Removed existing {archive_name}")

save_tox_data("test_archive_1_py.zip", gene_ids=filtered_gene_ids, gene_ids_name="gene_ids_v1.bin",
                      expression_vectors=filtered_kallisto_expr, expression_vectors_name="kallisto_data_v1.bin")

save_tox_data("test_archive_2_py.zip", family_centroids=centroids, family_centroids_name="centroids.bin")
try:
    save_tox_data("test_archive_3_py.zip", family_centroids=centroids, gene_ids=kallisto_expr)
except:
    print("Successfully threw exception")

save_tox_data("test_archive_4_py.zip", gene_ids=filtered_gene_ids, gene_ids_name="gene_ids_v1.bin",
                      expression_vectors=filtered_kallisto_expr, expression_vectors_name="kallisto_data_v1.bin",
                      gene_to_fam=filtered_gene_to_fam, gene_to_fam_name="gene_to_fam_v1.bin",
                      family_ids=family_ids, family_ids_name="family_ids_v1.bin",
                      family_centroids=centroids, family_centroids_name="family_centroids_v1.bin",
                      shift_vectors=shift_vectors, shift_vectors_name="shift_vectors_v1.bin")

result_1 = read_tox_data("test_archive_4_py.zip", load_gene_ids=True, load_expression_vectors=True, 
                                 load_gene_to_fam=True, load_family_ids=True, load_family_centroids=True, load_shift_vectors=True)
result_2 = read_tox_data("test_archive_4_py.zip", load_gene_ids=True, load_gene_to_fam=True)
try:
    result_3 = read_tox_data("test_archive_3_py.zip", load_gene_ids=True)
except:
    result_3 = None
    print("Successfully threw exception")
try:
    result_f = read_tox_data("test_archive_1_f.zip", True, True, True, True, True, True)
except:
    result_f = None 
    print("Could not read fortran archive. Does the file exist?")

try: 
    result_r = read_tox_data(zip_filename="test_archive_1_R.zip", gene_ids=True, expression_vectors=True, 
                                 gene_to_fam=True, family_ids=True, family_centroids=True, shift_vectors=True)
except:
    print("Could not read R archive. Does the file exist?")

# ---- NEW TESTS: Non-standard arrays and direct create_zip_archive calls ----

print("\n" + "="*50)
print("TESTING NON-STANDARD ARRAYS AND DIRECT create_zip_archive")
print("="*50)

def test_non_standard_arrays():
    """Test creating and saving non-standard arrays using direct create_zip_archive calls"""
    
    # Clean up any existing test archives from previous runs
    for archive_name in ["test_non_standard_1.zip", "test_non_standard_2.zip", "test_non_standard_3.zip"]:
        if os.path.exists(archive_name):
            os.remove(archive_name)
            print(f"Removed existing {archive_name}")
    
    # Create various non-standard arrays
    print("Creating non-standard arrays...")
    
    # 1. 3D array of integers
    array_3d_int = np.random.randint(0, 100, size=(5, 10, 3), dtype=np.int32)
    print(f"3D int array shape: {array_3d_int.shape}, dtype: {array_3d_int.dtype}")
    
    # 2. 1D array of floats with unusual values
    array_1d_float = np.array([np.nan, np.inf, -np.inf, 0.0, 1.5, -2.3], dtype=np.float64)
    print(f"1D float array: {array_1d_float}, dtype: {array_1d_float.dtype}")
    
    # 3. 2D array of strings with different lengths
    array_2d_char = np.array([
        ["short", "medium_length", "very_long_string_here"],
        ["a", "bb", "ccc"],
        ["test1", "test2", "test3"]
    ], dtype='U')
    print(f"2D char array shape: {array_2d_char.shape}, dtype: {array_2d_char.dtype}")
    
    # 4. Large 1D boolean array
    array_1d_bool = np.random.choice([True, False], size=1000)
    print(f"1D bool array shape: {array_1d_bool.shape}, dtype: {array_1d_bool.dtype}")
    
    # 5. Complex number array (converted to two real arrays for serialization)
    array_complex_real = np.random.rand(4, 4).astype(np.float64)
    array_complex_imag = np.random.rand(4, 4).astype(np.float64)
    print(f"Complex real part shape: {array_complex_real.shape}")
    print(f"Complex imag part shape: {array_complex_imag.shape}")
    
    # Serialize all arrays to temporary files
    print("\nSerializing arrays to temporary files...")
    
    tox_serialize_int_nd(array_3d_int, "temp_3d_int.bin")
    tox_serialize_real_nd(array_1d_float, "temp_1d_float.bin")
    tox_serialize_char_nd(array_2d_char, "temp_2d_char.bin")
    
    # For boolean array, convert to int32 for serialization
    array_1d_bool_int = array_1d_bool.astype(np.int32)
    tox_serialize_int_nd(array_1d_bool_int, "temp_1d_bool.bin")
    
    tox_serialize_real_nd(array_complex_real, "temp_complex_real.bin")
    tox_serialize_real_nd(array_complex_imag, "temp_complex_imag.bin")
    
    # Test 1: Direct call to create_zip_archive with non-standard arrays
    print("\nTest 1: Direct create_zip_archive call with non-standard arrays")
    keys = [
        "custom_3d_int_data",
        "special_float_array", 
        "string_matrix",
        "boolean_mask",
        "complex_real_part",
        "complex_imag_part"
    ]
    
    filenames = [
        "temp_3d_int.bin",
        "temp_1d_float.bin", 
        "temp_2d_char.bin",
        "temp_1d_bool.bin",
        "temp_complex_real.bin",
        "temp_complex_imag.bin"
    ]
    
    create_zip_archive("test_non_standard_1.zip", keys, filenames)
    print("Successfully created archive with non-standard arrays")
    
    # Test 2: Mixed standard and non-standard arrays
    print("\nTest 2: Mixed standard and non-standard arrays")
    
    # Add some standard arrays to the mix
    if 'filtered_gene_ids' in locals():
        tox_serialize_char_nd(filtered_gene_ids, "temp_standard_gene_ids.bin")
        tox_serialize_real_nd(filtered_kallisto_expr, "temp_standard_expr.bin")
        
        keys_mixed = keys + ["standard_gene_ids", "standard_expression"]
        filenames_mixed = filenames + ["temp_standard_gene_ids.bin", "temp_standard_expr.bin"]
        
        create_zip_archive("test_mixed_arrays.zip", keys_mixed, filenames_mixed)
        print("Successfully created archive with mixed standard and non-standard arrays")
    
    # Test 3: Read back and verify non-standard arrays
    print("\nTest 3: Reading back non-standard arrays")
    
    
    zip_filename = "test_non_standard_1.zip"

    file_mapping = extract_zip_archive(zip_filename)
    # Deserialize and verify some arrays
    if "custom_3d_int_data" in file_mapping:
        loaded_3d_int = tox_deserialize_int_nd(file_mapping["custom_3d_int_data"])
        print(f"Loaded 3D int array shape: {loaded_3d_int.shape}")
        print(f"Original and loaded arrays equal: {np.array_equal(array_3d_int, loaded_3d_int)}")
    
    if "special_float_array" in file_mapping:
        loaded_1d_float = tox_deserialize_real_nd(file_mapping["special_float_array"])
        print(f"Loaded 1D float array: {loaded_1d_float}")
        # For NaN/inf comparison, we need special handling
        print(f"Float arrays equal (considering NaN): {np.array_equal(array_1d_float, loaded_1d_float, equal_nan=True)}")
    
    if "string_matrix" in file_mapping:
        loaded_2d_char = tox_deserialize_char_nd(file_mapping["string_matrix"])
        print(f"Loaded 2D char array shape: {loaded_2d_char.shape}")
        print(f"String arrays equal: {np.array_equal(array_2d_char, loaded_2d_char)}")
    
    # Cleanup extracted files
    for filename in file_mapping.values():
        if os.path.exists(filename):
            os.remove(filename)
    if os.path.exists("manifest.txt"):
        os.remove("manifest.txt")
    
    # Test 4: Error handling - mismatched keys and filenames
    print("\nTest 4: Error handling with mismatched arrays")
    try:
        create_zip_archive("test_error.zip", ["key1", "key2"], ["file1.bin"])  # Mismatched lengths
        print("Should have thrown an error!")
    except ValueError as e:
        print(f"Correctly caught error: {e}")
    
    # Cleanup temporary files
    print("\nCleaning up temporary files...")
    temp_files = [
        "temp_3d_int.bin", "temp_1d_float.bin", "temp_2d_char.bin",
        "temp_1d_bool.bin", "temp_complex_real.bin", "temp_complex_imag.bin"
    ]
    
    if 'filtered_gene_ids' in locals():
        temp_files.extend(["temp_standard_gene_ids.bin", "temp_standard_expr.bin"])
    
    for temp_file in temp_files:
        if os.path.exists(temp_file):
            os.remove(temp_file)
            print(f"Removed: {temp_file}")
    
    print("All non-standard array tests completed!")

# Run the non-standard array tests
test_non_standard_arrays()

print("\n" + "="*50)
print("ALL TESTS COMPLETED SUCCESSFULLY!")
print("="*50)