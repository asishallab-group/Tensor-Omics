import ctypes
import os

# 1. Define constants and library path
LIB_NAME = "libF42CsvReader.so"
if os.name == "nt":
    LIB_NAME = "libF42CsvReader.dll"
elif os.uname().sysname == "Darwin":
    LIB_NAME = "libF42CsvReader.dylib"

SHARED_LIB_PATH = os.path.join("./build", LIB_NAME)
CSV_FILE = b"data.csv" # Pass strings as bytes to C
MAX_STR_LEN = 512

# 2. Load the shared library
try:
    lib = ctypes.CDLL(SHARED_LIB_PATH)
    print(f"Successfully loaded library: {SHARED_LIB_PATH}")
except OSError as e:
    print(f"Error loading library: {e}")
    exit()

# 3. Define function prototypes for type safety
# This tells ctypes what argument and return types the C functions expect.

# void read_csv_file_C(char*, int, bool, char, int*);
read_csv = lib.read_csv_file_C
read_csv.argtypes = [ctypes.c_char_p, ctypes.c_int, ctypes.c_bool, ctypes.c_char, ctypes.POINTER(ctypes.c_int)]

# int get_num_rows_C();
get_rows = lib.get_num_rows_C
get_rows.restype = ctypes.c_int

# int get_num_cols_C();
get_cols = lib.get_num_cols_C
get_cols.restype = ctypes.c_int

# void get_header_C(int, char*, int);
get_header = lib.get_header_C
get_header.argtypes = [ctypes.c_int, ctypes.c_char_p, ctypes.c_int]

# void get_cell_C(int, int, int*, int64_t*, double*, char*, int);
get_cell = lib.get_cell_C
get_cell.argtypes = [
    ctypes.c_int, ctypes.c_int,
    ctypes.POINTER(ctypes.c_int),     # data_type
    ctypes.POINTER(ctypes.c_int64),   # i_val
    ctypes.POINTER(ctypes.c_double),  # r_val
    ctypes.c_char_p,                  # c_val_c
    ctypes.c_int
]

# void cleanup_csv_data_C();
cleanup = lib.cleanup_csv_data_C

# --- Now, execute the test workflow ---

# 4. Call the CSV reader
io_status = ctypes.c_int()
read_csv(CSV_FILE, len(CSV_FILE), True, b',', ctypes.byref(io_status))

if io_status.value != 0:
    raise RuntimeError(f"Error reading CSV file. IO Status: {io_status.value}")
print("CSV file read successfully.")

# 5. Get dimensions
num_rows = get_rows()
num_cols = get_cols()
print(f"Dimensions: {num_rows} rows x {num_cols} cols")

# 6. Get and print header
header_parts = []
for j in range(1, num_cols + 1):
    buffer = ctypes.create_string_buffer(MAX_STR_LEN)
    get_header(j, buffer, MAX_STR_LEN)
    header_parts.append(buffer.value.decode('utf-8'))
print("Header:")
print("\t".join(header_parts))
print("-" * 50)

# 7. Get and print cell data
print("Data:")
for i in range(1, num_rows + 1):
    row_parts = []
    for j in range(1, num_cols + 1):
        # Create pointers for the output values
        data_type = ctypes.c_int()
        i_val = ctypes.c_int64()
        r_val = ctypes.c_double()
        c_val = ctypes.create_string_buffer(MAX_STR_LEN)

        get_cell(i, j, ctypes.byref(data_type), ctypes.byref(i_val), ctypes.byref(r_val), c_val, MAX_STR_LEN)
        
        cell_str = ""
        if data_type.value == 1: # Integer
            cell_str = str(i_val.value)
        elif data_type.value == 2: # Real
            cell_str = str(r_val.value)
        elif data_type.value == 3: # Character
            cell_str = c_val.value.decode('utf-8')
        row_parts.append(cell_str)
    print("\t".join(row_parts))

# 8. IMPORTANT: Call the cleanup routine
cleanup()
print("Cleanup function called. Memory released.")