"""
Synthetic benchmark generator for the STC/simulated-annealing experiments described in
`misc/mod_STC.md`'s "Simulated annealing CLI in C" section.

Adapted from `generate_datasets.py` -- same SEED, same N_POINTS, same low/medium/high noise
levels, and several of the same shapes (a line, a circular arc, a sine curve, a two-branch "Y"
bifurcation), plus one new shape (a Gaussian hill). The output is completely different, though:
one JSON file per dataset in `mod_STC.md`'s benchmark schema, not a CSV. Alongside its noisy
data points, every dataset carries a dense noise-free reference point cloud and a true
orthonormal tangent basis at every reference point, for every true manifold -- exactly what the
objective function in `mod_STC.md` needs to score a candidate STC run against ground truth.
Every shape here is defined by an exact, closed-form embedding of a low-dimensional domain into
`dim`-dimensional ambient space, so the reference cloud and tangent bases are computed once from
that embedding (central-difference Jacobian, orthonormalized via QR) rather than estimated from
noisy samples.

For now, every dataset uses one fixed, isotropic Gaussian noise standard deviation applied to
every ambient coordinate (heteroscedastic/per-point noise is future work, not implemented here).
`mixed_line_circle_noise` additionally includes a deliberately unstructured third segment
(`true_manifold_id = -1`, no corresponding entry in `manifolds[]`) -- points with no true
manifold at all, which a correct STC run should leave uncovered rather than absorb into some
super-ensemble.

`REF_POINTS_1D`/`REF_GRID_RES_2D` trade reference-cloud density (and so file size) against how
finely the ground truth resolves relative to the noise scale; both are kept well above the
lowest noise level (`NOISE_LEVELS["low"]`) so the nearest-reference-point approximation stays
accurate at every noise level generated here.

Run: python3 generate_sa_benchmark_datasets.py
"""

import json
import os

import numpy as np

SEED = 42
N_POINTS = 500
NOISE_LEVELS = {
    "low": 0.02,
    "medium": 0.08,
    "high": 0.2,
}
REF_POINTS_1D = 1000
REF_GRID_RES_2D = 150
TANGENT_FD_STEP = 1e-6
GAUSSIAN_HILL_AMPLITUDE = 1.0
GAUSSIAN_HILL_WIDTH = 0.3

HERE = os.path.dirname(os.path.abspath(__file__))


def dense_1d_grid(lo, hi, n=REF_POINTS_1D):
    return np.linspace(lo, hi, n).reshape(-1, 1)


def dense_2d_grid(lo, hi, res=REF_GRID_RES_2D):
    axis = np.linspace(lo, hi, res)
    uu, vv = np.meshgrid(axis, axis)
    return np.column_stack([uu.ravel(), vv.ravel()])


def tangent_basis(embed_fn, t0, domain_dim, ambient_dim, h=TANGENT_FD_STEP):
    """Central-difference Jacobian of embed_fn at t0, orthonormalized via QR."""
    jac = np.zeros((domain_dim, ambient_dim))
    for i in range(domain_dim):
        t_plus, t_minus = t0.copy(), t0.copy()
        t_plus[i] += h
        t_minus[i] -= h
        jac[i] = (embed_fn(t_plus) - embed_fn(t_minus)) / (2 * h)
    q, _ = np.linalg.qr(jac.T)
    return q.T


def manifold_geometry(manifold_id, embed_fn, domain_dim, ambient_dim, ref_domain):
    """Reference cloud + tangent bases for one manifold -- independent of noise level, so this
    is computed once per (shape, dim) and reused across that dataset's low/medium/high variants."""
    ref_points = np.array([embed_fn(t) for t in ref_domain])
    ref_bases = np.array([tangent_basis(embed_fn, t, domain_dim, ambient_dim) for t in ref_domain])
    return {
        "manifold_id": manifold_id,
        "intrinsic_dim": domain_dim,
        "reference_points": ref_points,
        "reference_tangent_bases": ref_bases,
    }


def manifold_json(geometry, noise_sd):
    return {
        "manifold_id": geometry["manifold_id"],
        "intrinsic_dim": geometry["intrinsic_dim"],
        "noise_sd": noise_sd,
        "reference_points": geometry["reference_points"].tolist(),
        "reference_tangent_bases": geometry["reference_tangent_bases"].tolist(),
    }


def save_json_dataset(dataset_id, ambient_dim, manifolds, coords, true_ids, filename, output_dir):
    payload = {
        "dataset_id": dataset_id,
        "ambient_dim": ambient_dim,
        "generation_seed": SEED,
        "manifolds": manifolds,
        "points": {
            "coordinates": coords.tolist(),
            "true_manifold_id": [int(v) for v in true_ids],
        },
    }
    with open(os.path.join(output_dir, filename), "w") as fh:
        json.dump(payload, fh)
    print(f"Saved: {filename} ({coords.shape[0]} points, {len(manifolds)} manifold(s))")


# =========================
# Shape embeddings: domain -> ambient point. Each takes `dim` (2 or 3) and returns a
# callable(t: np.ndarray) -> np.ndarray of length `dim`.
# =========================

def line_embed(dim):
    if dim == 2:
        return lambda t: np.array([t[0], 2 * t[0] + 1])
    return lambda t: np.array([t[0], 2 * t[0] + 1, 0.3 * t[0]])


def circular_arc_embed(dim):
    if dim == 2:
        return lambda t: np.array([np.cos(t[0]), np.sin(t[0])])
    return lambda t: np.array([np.cos(t[0]), np.sin(t[0]), 0.3 * t[0]])


def s_curve_embed(dim):
    if dim == 2:
        return lambda t: np.array([t[0], np.sin(2 * np.pi * t[0])])
    return lambda t: np.array([t[0], np.sin(2 * np.pi * t[0]), 0.3 * t[0]])


def gaussian_hill_embed(dim):
    if dim == 2:
        return lambda t: np.array([
            t[0],
            GAUSSIAN_HILL_AMPLITUDE * np.exp(-(t[0] ** 2) / (2 * GAUSSIAN_HILL_WIDTH ** 2)),
        ])
    return lambda t: np.array([
        t[0], t[1],
        GAUSSIAN_HILL_AMPLITUDE * np.exp(-(t[0] ** 2 + t[1] ** 2) / (2 * GAUSSIAN_HILL_WIDTH ** 2)),
    ])


def y_branch_embed(dim, direction):
    """`direction` = (dx, dy) the branch travels away from the shared vertex at the origin,
    t in [0, 1]."""
    dx, dy = direction
    if dim == 2:
        return lambda t: np.array([dx * t[0], dy * t[0]])
    return lambda t: np.array([dx * t[0], dy * t[0], 0.3 * t[0]])


# =========================
# Per-shape dataset builders. Each returns (geometries, pieces): `geometries` is the list of
# `manifolds[]` entries (independent of noise level); `pieces` is a list of
# (manifold_id_or_None, embed_fn, domain_dim, domain_lo, domain_hi, n_points), one per segment
# to realize at each noise level -- `manifold_id=None` marks a deliberately unstructured piece
# (no true manifold, sampled directly over the ambient box, `true_manifold_id = -1`).
# =========================

def single_manifold_pieces(dim, embed_fn, domain_dim, domain_lo, domain_hi):
    ref_domain = dense_1d_grid(domain_lo, domain_hi) if domain_dim == 1 else dense_2d_grid(domain_lo, domain_hi)
    geometries = [manifold_geometry(0, embed_fn, domain_dim, dim, ref_domain)]
    pieces = [(0, embed_fn, domain_dim, domain_lo, domain_hi, N_POINTS)]
    return geometries, pieces


def y_bifurcation_pieces(dim):
    n1 = N_POINTS // 3
    n2 = N_POINTS // 3
    n3 = N_POINTS - n1 - n2
    branches = [
        (0, y_branch_embed(dim, (0.0, -0.8)), n1),
        (1, y_branch_embed(dim, (-0.6, 0.8)), n2),
        (2, y_branch_embed(dim, (0.6, 0.8)), n3),
    ]
    geometries, pieces = [], []
    for manifold_id, embed_fn, n in branches:
        geometries.append(manifold_geometry(manifold_id, embed_fn, 1, dim, dense_1d_grid(0.0, 1.0)))
        pieces.append((manifold_id, embed_fn, 1, 0.0, 1.0, n))
    return geometries, pieces


def mixed_line_circle_noise_pieces(dim):
    n1 = N_POINTS // 3
    n2 = N_POINTS // 3
    n3 = N_POINTS - n1 - n2
    line_fn = line_embed(dim)
    circle_fn = circular_arc_embed(dim)
    geometries = [
        manifold_geometry(0, line_fn, 1, dim, dense_1d_grid(-1.0, 1.0)),
        manifold_geometry(1, circle_fn, 1, dim, dense_1d_grid(-np.pi / 4, np.pi / 4)),
    ]
    pieces = [
        (0, line_fn, 1, -1.0, 1.0, n1),
        (1, circle_fn, 1, -np.pi / 4, np.pi / 4, n2),
        (None, None, None, -1.0, 1.0, n3),  # unstructured background, true_manifold_id = -1
    ]
    return geometries, pieces


SHAPES = {
    "line": lambda dim: single_manifold_pieces(dim, line_embed(dim), 1, -1.0, 1.0),
    "circular_arc": lambda dim: single_manifold_pieces(dim, circular_arc_embed(dim), 1, -np.pi / 4, np.pi / 4),
    "s_curve": lambda dim: single_manifold_pieces(dim, s_curve_embed(dim), 1, -1.0, 1.0),
    "gaussian_hill": lambda dim: single_manifold_pieces(dim, gaussian_hill_embed(dim), dim - 1, -1.0, 1.0),
    "y_bifurcation": y_bifurcation_pieces,
    "mixed_line_circle_noise": mixed_line_circle_noise_pieces,
}


def realize_pieces(pieces, dim, sigma):
    """Draw one noisy instance of every piece for a given noise level."""
    coords, true_ids = [], []
    for manifold_id, embed_fn, domain_dim, lo, hi, n in pieces:
        if manifold_id is None:
            coords.append(np.random.uniform(lo, hi, size=(n, dim)))
            true_ids.append(np.full(n, -1, dtype=int))
            continue
        t = np.random.uniform(lo, hi, size=(n, domain_dim))
        clean = np.array([embed_fn(row) for row in t])
        coords.append(clean + np.random.normal(0.0, sigma, size=clean.shape))
        true_ids.append(np.full(n, manifold_id, dtype=int))
    return np.concatenate(coords), np.concatenate(true_ids)


def generate(dim):
    """Generate every shape's low/medium/high-noise JSON datasets at ambient dimension `dim`."""
    output_dir = os.path.join(HERE, "data_json", f"{dim}d")
    os.makedirs(output_dir, exist_ok=True)
    np.random.seed(SEED)

    for name, builder in SHAPES.items():
        geometries, pieces = builder(dim)
        for label, sigma in NOISE_LEVELS.items():
            coords, true_ids = realize_pieces(pieces, dim, sigma)
            manifolds = [manifold_json(g, sigma) for g in geometries]
            dataset_id = f"{name}_{dim}d_noise_{label}"
            save_json_dataset(dataset_id, dim, manifolds, coords, true_ids, f"{dataset_id}.json", output_dir)


if __name__ == "__main__":
    for dim in (2, 3):
        generate(dim)
