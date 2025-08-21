#!/usr/bin/env python3

"""
This script performs an optimized SVD on a large, sparse protein similarity matrix.

The key improvement in this version is the use of a `scipy.sparse.linalg.LinearOperator`.
This allows us to implicitly "center" the sparse matrix without converting it to a dense,
memory-intensive format. The centering operation is performed as a part of the
matrix-vector product, which is the core operation of the SVD solver.

This approach addresses two critical points:
1. It correctly centers the data, which is essential for a true Principal Component Analysis.
2. It prevents memory overflow by avoiding the creation of a dense matrix.
"""

import os
import pickle
import numpy as np
import pandas as pd
from Bio import SeqIO
from scipy.sparse import coo_matrix, csr_matrix
from scipy.sparse.linalg import svds, LinearOperator

# === CONFIGURATION ===
CONFIG = {
    "FASTA_FILE": os.path.join("material", "uniref50_morethan1.fasta"),
    "PKL_FILE": os.path.join("results", "geometric_mean_pairs.pkl"),
    "REDUCED_CSV": os.path.join("results", "centered_svd_reduced_matrix.csv"),
    "EXPLAINED_VARIANCE_CSV": os.path.join("results", "centered_svd_explained_variance.csv"),
    "MEAN_VECTOR_PKL": os.path.join("results", "mean_vector.pkl"),
    "MIN_EXPLAINED_VARIANCE": 0.80,
    "MAX_COMPONENTS": 100 # Maximum number of components for the truncated SVD
}

print("Starting script with the following configurations:")
for key, value in CONFIG.items():
    print(f"- {key}: {value}")

# --- Step 1: Load Protein IDs and Similarity Data ---
print("\n[1/5] Loading protein IDs and similarity data...")
try:
    with open(CONFIG["FASTA_FILE"]) as handle:
        protein_ids = [record.id for record in SeqIO.parse(handle, "fasta")]
    
    with open(CONFIG["PKL_FILE"], 'rb') as f:
        sim_dict = pickle.load(f)
except FileNotFoundError as e:
    print(f"Error: Required file not found. Please check your config. ({e})")
    exit()

N = len(protein_ids)
print(f"Total proteins found: {N}")
protein_index = {p: i for i, p in enumerate(protein_ids)}

# --- Step 2: Build the Sparse Matrix ---
print("\n[2/5] Building sparse similarity matrix...")
rows = []
cols = []
data = []
for (p1, p2), val in sim_dict.items():
    if p1 in protein_index and p2 in protein_index:
        rows.append(protein_index[p1])
        cols.append(protein_index[p2])
        data.append(val)
        # Add symmetric entry as the matrix is symmetric
        if p1 != p2:
            rows.append(protein_index[p2])
            cols.append(protein_index[p1])
            data.append(val)
# Add diagonal entries for self-similarity (score of 1.0)
for i in range(N):
    rows.append(i)
    cols.append(i)
    data.append(1.0)

# Create the sparse matrix in CSR format for efficient row-wise operations
M_sparse = coo_matrix((data, (rows, cols)), shape=(N, N), dtype=np.float32).tocsr()
print(f"Sparse matrix built with shape {M_sparse.shape} and {M_sparse.nnz} non-zero entries.")

# --- Step 3: Calculate and Save the Mean Vector for Centering ---
print("\n[3/5] Calculating and saving the mean vector...")
# To center the matrix, we need the mean of each row. The true mean
# of a row is its sum divided by the total number of columns (N).
mean_vector = M_sparse.sum(axis=0) / N
mean_vector = np.array(mean_vector).flatten() # Convert to a 1D numpy array

# Save the mean vector to a file for future use in inference
with open(CONFIG["MEAN_VECTOR_PKL"], 'wb') as f:
    pickle.dump(mean_vector, f)
print(f"Mean vector saved to {CONFIG['MEAN_VECTOR_PKL']}")

# --- Step 4: Define the Linear Operator and Run SVD ---
print(f"\n[4/5] Running SVD using a Linear Operator...")

# This function implicitly centers the matrix during a matrix-vector product
def matvec_centered_op(x):
    # The operation is: (M_sparse @ x) - (mean_vector @ x)
    # This is equivalent to (M_sparse - M_mean_matrix) @ x
    result_sparse = M_sparse.dot(x)
    result_mean = mean_vector.dot(x)
    return result_sparse - result_mean

# Create the LinearOperator object
A_centered_op = LinearOperator(
    shape=(N, N),
    matvec=matvec_centered_op,
    rmatvec=matvec_centered_op # SVD for a symmetric matrix uses the same operator for both
)

# Run the truncated SVD on the LinearOperator
try:
    U, S, VT = svds(A_centered_op, k=min(CONFIG["MAX_COMPONENTS"], N - 1))
    # svds returns components in ascending order of singular values; reverse them
    idx = np.argsort(S)[::-1]
    S = S[idx]
    U = U[:, idx]
    VT = VT[idx, :]
except Exception as e:
    print(f"Error during SVD: {e}")
    print("Please check the matrix shape and the value of MAX_COMPONENTS.")
    exit()

# --- Step 5: Calculate Explained Variance and Save Results ---
print("\n[5/5] Calculating explained variance and saving results...")
total_variance = np.sum(S ** 2)
explained_var = (S ** 2) / total_variance
cumulative_var = np.cumsum(explained_var)
n_components = np.searchsorted(cumulative_var, CONFIG["MIN_EXPLAINED_VARIANCE"]) + 1

print(f"Components to reach {CONFIG['MIN_EXPLAINED_VARIANCE']*100}% explained variance: {n_components}")

# Save the explained variance data
pd.DataFrame({
    "singular_value": S,
    "explained_variance_ratio": explained_var,
    "cumulative_variance": cumulative_var
}).to_csv(CONFIG["EXPLAINED_VARIANCE_CSV"], index=False)
print(f"Explained variance data saved to {CONFIG['EXPLAINED_VARIANCE_CSV']}")

# Save the reduced matrix with the selected number of components
# M_reduced = U[:, :n_components] @ np.diag(S[:n_components])
# Note: In a centered SVD (PCA), the reduced matrix is typically just U @ sqrt(S)
# or just U. We'll save the projection U[:, :n_components]
M_reduced = U[:, :n_components]
reduced_df = pd.DataFrame(M_reduced, index=protein_ids)
reduced_df.to_csv(CONFIG["REDUCED_CSV"])
print(f"Reduced matrix with {n_components} components saved to {CONFIG['REDUCED_CSV']}")

print("\nDone! ✅")
