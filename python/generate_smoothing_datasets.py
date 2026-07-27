import numpy as np
import pandas as pd
import os

# =========================
# GLOBAL CONFIGURATION
# =========================
SEED = 42
N_POINTS = 500
DIM = 3  # <--- Set any total dimension here (e.g., 3 means X1, X2, and target Y)
OUTPUT_DIR = f"results/data/{DIM}d"

NOISE_LEVELS = {
    "low": 0.02,
    "medium": 0.08,
    "high": 0.2
}

np.random.seed(SEED)
os.makedirs(OUTPUT_DIR, exist_ok=True)

# =========================
# UTILITY TO SAVE CSV (N-DIMENSIONAL)
# =========================
def save_csv_nd(X, y, filename):
    """
    Saves an N-dimensional dataset to CSV.
    X: feature matrix of shape (N_POINTS, DIM - 1)
    y: target vector of shape (N_POINTS,)
    """
    if X.ndim == 1:
        X = X.reshape(-1, 1)
        
    num_features = X.shape[1]
    
    # Dynamically build feature names: {'x1': ..., 'x2': ...}
    data = {f"x{i+1}": X[:, i] for i in range(num_features)}
    data["y_original"] = y
    
    df = pd.DataFrame(data)
    df.to_csv(os.path.join(OUTPUT_DIR, filename), index=False)
    print(f"Saved: {filename} (Shape: {df.shape})")

# =========================
# 2.2 FUNCTIONAL DATASETS
# =========================
n_features = DIM - 1
X_func = np.random.uniform(-1, 1, size=(N_POINTS, n_features))

def generate_functional_nd(name, f):
    for label, sigma in NOISE_LEVELS.items():
        noise = np.random.normal(0, sigma, N_POINTS)
        y = f(X_func) + noise
        save_csv_nd(X_func, y, f"{name}_{DIM}d_noise_{label}.csv")

# 1) Linear: combination of all features + 1
generate_functional_nd("linear", lambda X: np.sum(2 * X, axis=1) + 1)

# 2) Cubic: 4*x1^3 - 2*x2^2 + x3... (adapted gracefully for any N-D)
generate_functional_nd("cubic", lambda X: 4 * X[:, 0]**3 - (2 * np.sum(X[:, 1:]**2, axis=1) if X.shape[1] > 1 else 0) + X[:, 0])

# 3) Exponential: exp of the sum of features
generate_functional_nd("exponential", lambda X: np.exp(np.sum(X, axis=1)))


# =========================
# 2.3 NON-FUNCTIONAL DATASETS
# =========================

# 1) S-Curve (Manifold-like behavior mapped across all dimensions)
t_s = np.random.uniform(-1, 1, N_POINTS)
X_s = np.random.uniform(-1, 1, size=(N_POINTS, n_features))
X_s[:, 0] = t_s  # Primary driving component is structural
y_s_clean = np.sin(2 * np.pi * t_s)

for label, sigma in NOISE_LEVELS.items():
    y_s = y_s_clean + np.random.normal(0, sigma, N_POINTS)
    save_csv_nd(X_s, y_s, f"s_curve_{DIM}d_noise_{label}.csv")

# 2) Circular arc
X_arc = np.random.uniform(-1, 1, size=(N_POINTS, n_features))
y_arc_clean = np.sin(np.cos(X_arc[:, 0]))

for label, sigma in NOISE_LEVELS.items():
    y_arc = y_arc_clean + np.random.normal(0, sigma, N_POINTS)
    save_csv_nd(X_arc, y_arc, f"circular_arc_{DIM}d_noise_{label}.csv")


# =========================
# 2.4 BAD CASES
# =========================

# 1) Curve with a kink
X_kink = np.random.uniform(-1, 1, size=(N_POINTS, n_features))
# Process structural part based on the first feature split
mask1 = X_kink[:, 0] <= 0.5
mask2 = ~mask1

y_kink_clean = np.zeros(N_POINTS)
y_kink_clean[mask1] = 2 * X_kink[mask1, 0]
y_kink_clean[mask2] = 1.0

for label, sigma in NOISE_LEVELS.items():
    y_kink = y_kink_clean + np.random.normal(0, sigma, N_POINTS)
    save_csv_nd(X_kink, y_kink, f"kinked_curve_{DIM}d_noise_{label}.csv")

# 2) Bifurcation dataset (2-way split, Y-shaped structure)
n_trunk = N_POINTS // 3
n_branch = N_POINTS // 3
n_right = N_POINTS - n_trunk - n_branch

# Initialize a base noise matrix for background features
X_bif2 = np.random.uniform(-1, 1, size=(N_POINTS, n_features))

# Trunk features
X_bif2[:n_trunk, 0] = 0
y_trunk = np.linspace(-0.8, 0, n_trunk)

# Left branch features
t_left = np.linspace(0, 1, n_branch)
X_bif2[n_trunk:n_trunk+n_branch, 0] = -0.6 * t_left
y_left = 0.8 * t_left

# Right branch features
t_right = np.linspace(0, 1, n_right)
X_bif2[n_trunk+n_branch:, 0] = 0.6 * t_right
y_right = 0.8 * t_right

y_bifurc2_clean = np.concatenate([y_trunk, y_left, y_right])

for label, sigma in NOISE_LEVELS.items():
    # Inject spatial noise safely across variables
    X_bif2_noisy = X_bif2 + np.random.normal(0, sigma, size=X_bif2.shape)
    y_bifurc2 = y_bifurc2_clean + np.random.normal(0, sigma, N_POINTS)
    save_csv_nd(X_bif2_noisy, y_bifurc2, f"bifurcation_2way_{DIM}d_noise_{label}.csv")

# 3) Trifurcation dataset (3-way split)
n_trunk3 = N_POINTS // 4
n_branch3 = N_POINTS // 4
n_right3 = N_POINTS - n_trunk3 - 2 * n_branch3

X_bif3 = np.random.uniform(-1, 1, size=(N_POINTS, n_features))

# Trunk
X_bif3[:n_trunk3, 0] = 0
y_trunk3 = np.linspace(-0.8, 0, n_trunk3)

# Left branch
t_left3 = np.linspace(0, 1, n_branch3)
X_bif3[n_trunk3:n_trunk3+n_branch3, 0] = -0.7 * t_left3
y_left3 = 0.7 * t_left3

# Center branch
t_center3 = np.linspace(0, 1, n_branch3)
X_bif3[n_trunk3+n_branch3:n_trunk3+2*n_branch3, 0] = 0
y_center3 = 0.9 * t_center3

# Right branch
t_right3 = np.linspace(0, 1, n_right3)
X_bif3[n_trunk3+2*n_branch3:, 0] = 0.7 * t_right3
y_right3 = 0.7 * t_right3

y_bifurc3_clean = np.concatenate([y_trunk3, y_left3, y_center3, y_right3])

for label, sigma in NOISE_LEVELS.items():
    X_bif3_noisy = X_bif3 + np.random.normal(0, sigma, size=X_bif3.shape)
    y_bifurc3 = y_bifurc3_clean + np.random.normal(0, sigma, N_POINTS)
    save_csv_nd(X_bif3_noisy, y_bifurc3, f"bifurcation_3way_{DIM}d_noise_{label}.csv")

# 4) Heteroscedastic noise (scaled by level)
X_hetero = np.random.uniform(0, 1, size=(N_POINTS, n_features))
x_mean = np.mean(X_hetero, axis=1)

for label, scale in NOISE_LEVELS.items():
    sigma_x = 0.05 + 2 * scale * x_mean
    noise_hetero = np.random.normal(0, sigma_x, size=N_POINTS)
    y_hetero = np.sum(2 * X_hetero, axis=1) + 1 + noise_hetero
    save_csv_nd(X_hetero, y_hetero, f"heteroscedastic_{DIM}d_noise_{label}.csv")


# =========================
# 3) Mixed Generators
# =========================
n1 = N_POINTS // 3
n2 = N_POINTS // 3
n3 = N_POINTS - n1 - n2

# Initialize background features
X_mix_base = np.random.uniform(-1, 1, size=(N_POINTS, n_features))

# Line structural domain
x_a = np.linspace(-1, 1, n1)
y_a_clean = 2 * x_a + 1
X_mix_base[:n1, 0] = x_a

# Arc structural domain
t_mix = np.linspace(-np.pi/4, np.pi/4, n2)
x_b = np.cos(t_mix)
y_b_clean = np.sin(t_mix)
X_mix_base[n1:n1+n2, 0] = x_b

# Random domain (Uniform noise natively across all components)
y_c = np.random.uniform(-1, 1, n3)

for label, sigma in NOISE_LEVELS.items():
    y_a = y_a_clean + np.random.normal(0, sigma, n1)
    y_b = y_b_clean + np.random.normal(0, sigma, n2)
    
    y_mix = np.concatenate([y_a, y_b, y_c])
    save_csv_nd(X_mix_base, y_mix, f"mixed_generators_{DIM}d_noise_{label}.csv")


print(f"\nALL {DIM}D datasets (low / medium / high) were generated successfully.")