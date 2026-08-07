"""
Synthetic point-cloud generator for the STC experiments in this directory.

Adapted from `python/generate_smoothing_datasets.py` on branch `smoothing` (the LoManLe
manifold-learning experiments) -- same generators, same seed, same noise levels, same output
CSV schema (`x1, x2, ..., y_original`), so the exact same datasets that exercised LoManLe's
skeleton-fitting also exercise STC's clustering. The only real change: this writes both the
2D and 3D variant of every dataset in one run (DIM used to be a global you edited and reran
for), and writes into `misc/STC-experiments/data/` instead of `results/data/`.

Run: python3 generate_datasets.py
"""

import numpy as np
import pandas as pd
import os

SEED = 42
N_POINTS = 500
NOISE_LEVELS = {
    "low": 0.02,
    "medium": 0.08,
    "high": 0.2,
}

HERE = os.path.dirname(os.path.abspath(__file__))


def save_csv_nd(X, y, filename, output_dir):
    """Saves an N-dimensional dataset to CSV: X is (N_POINTS, n_features), y is (N_POINTS,)."""
    if X.ndim == 1:
        X = X.reshape(-1, 1)
    data = {f"x{i + 1}": X[:, i] for i in range(X.shape[1])}
    data["y_original"] = y
    df = pd.DataFrame(data)
    df.to_csv(os.path.join(output_dir, filename), index=False)
    print(f"Saved: {filename} (Shape: {df.shape})")


def generate(dim):
    """Generate every dataset variant at ambient dimension `dim` (2 or 3)."""
    output_dir = os.path.join(HERE, "data", f"{dim}d")
    os.makedirs(output_dir, exist_ok=True)
    rng_state = np.random.RandomState(SEED)
    # np.random.* (module-level, not rng_state) is what the original generator used
    # throughout -- reseeding here keeps every dim's run independently reproducible from
    # the same SEED, matching what re-running the original script with a different DIM
    # and the same SEED did.
    np.random.seed(SEED)
    del rng_state

    n_features = dim - 1

    # =========================
    # Functional datasets
    # =========================
    X_func = np.random.uniform(-1, 1, size=(N_POINTS, n_features))

    def generate_functional_nd(name, f):
        for label, sigma in NOISE_LEVELS.items():
            noise = np.random.normal(0, sigma, N_POINTS)
            y = f(X_func) + noise
            save_csv_nd(X_func, y, f"{name}_{dim}d_noise_{label}.csv", output_dir)

    generate_functional_nd("linear", lambda X: np.sum(2 * X, axis=1) + 1)
    generate_functional_nd(
        "cubic",
        lambda X: 4 * X[:, 0] ** 3 - (2 * np.sum(X[:, 1:] ** 2, axis=1) if X.shape[1] > 1 else 0) + X[:, 0],
    )
    generate_functional_nd("exponential", lambda X: np.exp(np.sum(X, axis=1)))

    # =========================
    # Non-functional datasets
    # =========================
    t_s = np.random.uniform(-1, 1, N_POINTS)
    X_s = np.random.uniform(-1, 1, size=(N_POINTS, n_features))
    X_s[:, 0] = t_s
    y_s_clean = np.sin(2 * np.pi * t_s)
    for label, sigma in NOISE_LEVELS.items():
        y_s = y_s_clean + np.random.normal(0, sigma, N_POINTS)
        save_csv_nd(X_s, y_s, f"s_curve_{dim}d_noise_{label}.csv", output_dir)

    X_arc = np.random.uniform(-1, 1, size=(N_POINTS, n_features))
    y_arc_clean = np.sin(np.cos(X_arc[:, 0]))
    for label, sigma in NOISE_LEVELS.items():
        y_arc = y_arc_clean + np.random.normal(0, sigma, N_POINTS)
        save_csv_nd(X_arc, y_arc, f"circular_arc_{dim}d_noise_{label}.csv", output_dir)

    # =========================
    # Bad cases
    # =========================
    X_kink = np.random.uniform(-1, 1, size=(N_POINTS, n_features))
    mask1 = X_kink[:, 0] <= 0.5
    mask2 = ~mask1
    y_kink_clean = np.zeros(N_POINTS)
    y_kink_clean[mask1] = 2 * X_kink[mask1, 0]
    y_kink_clean[mask2] = 1.0
    for label, sigma in NOISE_LEVELS.items():
        y_kink = y_kink_clean + np.random.normal(0, sigma, N_POINTS)
        save_csv_nd(X_kink, y_kink, f"kinked_curve_{dim}d_noise_{label}.csv", output_dir)

    n_trunk = N_POINTS // 3
    n_branch = N_POINTS // 3
    n_right = N_POINTS - n_trunk - n_branch
    X_bif2 = np.random.uniform(-1, 1, size=(N_POINTS, n_features))
    X_bif2[:n_trunk, 0] = 0
    y_trunk = np.linspace(-0.8, 0, n_trunk)
    t_left = np.linspace(0, 1, n_branch)
    X_bif2[n_trunk:n_trunk + n_branch, 0] = -0.6 * t_left
    y_left = 0.8 * t_left
    t_right = np.linspace(0, 1, n_right)
    X_bif2[n_trunk + n_branch:, 0] = 0.6 * t_right
    y_right = 0.8 * t_right
    y_bifurc2_clean = np.concatenate([y_trunk, y_left, y_right])
    for label, sigma in NOISE_LEVELS.items():
        X_bif2_noisy = X_bif2 + np.random.normal(0, sigma, size=X_bif2.shape)
        y_bifurc2 = y_bifurc2_clean + np.random.normal(0, sigma, N_POINTS)
        save_csv_nd(X_bif2_noisy, y_bifurc2, f"bifurcation_2way_{dim}d_noise_{label}.csv", output_dir)

    n_trunk3 = N_POINTS // 4
    n_branch3 = N_POINTS // 4
    n_right3 = N_POINTS - n_trunk3 - 2 * n_branch3
    X_bif3 = np.random.uniform(-1, 1, size=(N_POINTS, n_features))
    X_bif3[:n_trunk3, 0] = 0
    y_trunk3 = np.linspace(-0.8, 0, n_trunk3)
    t_left3 = np.linspace(0, 1, n_branch3)
    X_bif3[n_trunk3:n_trunk3 + n_branch3, 0] = -0.7 * t_left3
    y_left3 = 0.7 * t_left3
    t_center3 = np.linspace(0, 1, n_branch3)
    X_bif3[n_trunk3 + n_branch3:n_trunk3 + 2 * n_branch3, 0] = 0
    y_center3 = 0.9 * t_center3
    t_right3 = np.linspace(0, 1, n_right3)
    X_bif3[n_trunk3 + 2 * n_branch3:, 0] = 0.7 * t_right3
    y_right3 = 0.7 * t_right3
    y_bifurc3_clean = np.concatenate([y_trunk3, y_left3, y_center3, y_right3])
    for label, sigma in NOISE_LEVELS.items():
        X_bif3_noisy = X_bif3 + np.random.normal(0, sigma, size=X_bif3.shape)
        y_bifurc3 = y_bifurc3_clean + np.random.normal(0, sigma, N_POINTS)
        save_csv_nd(X_bif3_noisy, y_bifurc3, f"bifurcation_3way_{dim}d_noise_{label}.csv", output_dir)

    X_hetero = np.random.uniform(0, 1, size=(N_POINTS, n_features))
    x_mean = np.mean(X_hetero, axis=1)
    for label, scale in NOISE_LEVELS.items():
        sigma_x = 0.05 + 2 * scale * x_mean
        noise_hetero = np.random.normal(0, sigma_x, size=N_POINTS)
        y_hetero = np.sum(2 * X_hetero, axis=1) + 1 + noise_hetero
        save_csv_nd(X_hetero, y_hetero, f"heteroscedastic_{dim}d_noise_{label}.csv", output_dir)

    # =========================
    # Mixed generators
    # =========================
    n1 = N_POINTS // 3
    n2 = N_POINTS // 3
    n3 = N_POINTS - n1 - n2
    X_mix_base = np.random.uniform(-1, 1, size=(N_POINTS, n_features))
    x_a = np.linspace(-1, 1, n1)
    y_a_clean = 2 * x_a + 1
    X_mix_base[:n1, 0] = x_a
    t_mix = np.linspace(-np.pi / 4, np.pi / 4, n2)
    x_b = np.cos(t_mix)
    y_b_clean = np.sin(t_mix)
    X_mix_base[n1:n1 + n2, 0] = x_b
    y_c = np.random.uniform(-1, 1, n3)
    for label, sigma in NOISE_LEVELS.items():
        y_a = y_a_clean + np.random.normal(0, sigma, n1)
        y_b = y_b_clean + np.random.normal(0, sigma, n2)
        y_mix = np.concatenate([y_a, y_b, y_c])
        save_csv_nd(X_mix_base, y_mix, f"mixed_generators_{dim}d_noise_{label}.csv", output_dir)


if __name__ == "__main__":
    for dim in (2, 3):
        generate(dim)
