#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import math
import json
import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix, csr_matrix, load_npz, save_npz
from scipy.sparse.linalg import LinearOperator, eigsh

# --------------------------------------------------------------------
# 1) Build sparse similarity matrix from Diamond output
# --------------------------------------------------------------------
def build_similarity_matrix_from_diamond(tsv_path, out_npz_path, ids_json_path,
                                         geom_formula='sqrt', min_weight=1e-12):
    """
    Parse Diamond TSV (qseqid, sseqid, pident, overlap) and create
    symmetric sparse similarity matrix (CSR).
    Also save mapping {protein_id: index}.
    """
    ids = {}
    idx_counter = 0
    rows, cols, data = [], [], []

    for chunk in pd.read_csv(tsv_path, sep='\t', header=0, chunksize=2_000_000,
                             usecols=['qseqid','sseqid','pident','overlap']):
        for _, r in chunk.iterrows():
            q = r['qseqid']; s = r['sseqid']
            if q not in ids:
                ids[q] = idx_counter; idx_counter += 1
            if s not in ids:
                ids[s] = idx_counter; idx_counter += 1

            i = ids[q]; j = ids[s]
            pident = float(r['pident']) / 100.0
            overlap = float(r['overlap']) / 100.0
            if geom_formula == 'sqrt':
                w = math.sqrt(max(pident, 0.0) * max(overlap, 0.0))
            else:
                w = math.exp((math.log(max(pident, 1e-300)) + math.log(max(overlap,1e-300))) / 2.0)

            if w < min_weight:
                continue

            rows.append(i); cols.append(j); data.append(w)
            if i != j:
                rows.append(j); cols.append(i); data.append(w)

    n = idx_counter
    print(f"Unique proteins: {n}, non-zeros: {len(data)}")

    M = coo_matrix((data, (rows, cols)), shape=(n, n))
    M.sum_duplicates()
    M = M.tocsr()

    save_npz(out_npz_path, M)
    with open(ids_json_path, 'w') as fh:
        json.dump(ids, fh)

    return out_npz_path, ids_json_path

# --------------------------------------------------------------------
# 2) Compute and save centering means
# --------------------------------------------------------------------
def compute_and_save_means(K_csr, means_npz_path):
    n = K_csr.shape[0]
    col_mean = np.asarray(K_csr.mean(axis=0)).ravel()
    overall_mean = float(K_csr.sum()) / (n * n)
    np.savez_compressed(means_npz_path, col_mean=col_mean, overall_mean=overall_mean)
    print(f"Means saved to {means_npz_path}")
    return col_mean, overall_mean

# --------------------------------------------------------------------
# 3) Centered LinearOperator and trace
# --------------------------------------------------------------------
def make_centered_kernel_operator(K_csr):
    def matvec(v):
        mean_v = v.mean()
        r = v - mean_v
        t = K_csr.dot(r)
        return t - t.mean()
    return LinearOperator(shape=K_csr.shape, matvec=matvec, rmatvec=matvec, dtype=np.float64)

def compute_trace_Kc(K_csr):
    n = K_csr.shape[0]
    trace_K = K_csr.diagonal().sum()
    total_sum = K_csr.sum()
    return float(trace_K - (1.0 / n) * total_sum)

# --------------------------------------------------------------------
# 4) Eigen decomposition until target variance
# --------------------------------------------------------------------
def compute_top_eigenpairs_until_variance(K_csr, target_variance=0.80,
                                          k0=50, step=50, k_max=2000):
    Kc_op = make_centered_kernel_operator(K_csr)
    trace_Kc = compute_trace_Kc(K_csr)
    if trace_Kc <= 0:
        raise ValueError("trace(Kc) <= 0, check data.")

    k = k0
    while True:
        print(f"Computing top {k} eigenpairs...")
        eigvals, eigvecs = eigsh(Kc_op, k=k, which='LA')
        order = np.argsort(eigvals)[::-1]
        eigvals = eigvals[order]
        eigvecs = eigvecs[:, order]
        cumvar = np.cumsum(eigvals)
        frac = cumvar / trace_Kc
        e = np.searchsorted(frac, target_variance) + 1
        print(f"First {e} comps explain {frac[e-1]:.4f} variance")

        if frac[e-1] >= target_variance or k >= k_max:
            return eigvals[:e], eigvecs[:, :e]
        k = min(k + step, k_max)

# --------------------------------------------------------------------
# 5) Main script
# --------------------------------------------------------------------
if __name__ == "__main__":
    # Files
    diamond_tsv = "diamond_results.tsv"
    matrix_npz = "similarity_matrix.npz"
    ids_json = "ids.json"
    means_npz = "centering_means.npz"
    reduced_csv = "kpca_reduced.csv"

    # Step 1: Build matrix if not yet done
    if not os.path.exists(matrix_npz):
        build_similarity_matrix_from_diamond(diamond_tsv, matrix_npz, ids_json)

    # Step 2: Load matrix and compute/save means
    K = load_npz(matrix_npz).tocsr()
    col_mean, overall_mean = compute_and_save_means(K, means_npz)

    # Step 3: Eigen decomposition
    eigvals, eigvecs = compute_top_eigenpairs_until_variance(K, target_variance=0.80)

    # Step 4: Koordinaten im reduzierten Raum berechnen
    sqrt_lambda = np.sqrt(np.clip(eigvals, a_min=1e-300, a_max=None))
    coords = eigvecs * sqrt_lambda[np.newaxis, :]  # shape (n, e)

    # Step 5: IDs laden und als CSV speichern
    with open(ids_json, 'r') as fh:
        id_map = json.load(fh)
    index_to_id = {v: k for k, v in id_map.items()}
    ids_ordered = [index_to_id[i] for i in range(len(index_to_id))]

    df = pd.DataFrame(coords, index=ids_ordered,
                      columns=[f"PC{i+1}" for i in range(coords.shape[1])])
    df.to_csv(reduced_csv)
    print(f"Reduced matrix saved to {reduced_csv}")
