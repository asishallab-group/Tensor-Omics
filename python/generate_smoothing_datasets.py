import numpy as np
import pandas as pd
import os

# =========================
# GLOBAL CONFIGURATION
# =========================
SEED = 42
N_POINTS = 500
OUTPUT_DIR = "results/data"

NOISE_LEVELS = {
    "low": 0.02,
    "medium": 0.08,
    "high": 0.2
}

np.random.seed(SEED)
os.makedirs(OUTPUT_DIR, exist_ok=True)

# =========================
# UTILITY TO SAVE CSV
# =========================
def save_csv(x, y, filename):
    df = pd.DataFrame({
        "x": x,
        "y_original": y
    })
    df.to_csv(os.path.join(OUTPUT_DIR, filename), index=False)
    print(f"Saved: {filename}")

# =========================
# 2.2 FUNCTIONAL DATASETS
# =========================

x = np.linspace(-1, 1, N_POINTS)  

def generate_functional(name, f):
    for label, sigma in NOISE_LEVELS.items():
        noise = np.random.normal(0, sigma, N_POINTS)
        y = f(x) + noise
        save_csv(x, y, f"{name}_noise_{label}.csv")

generate_functional("linear", lambda x: 2*x + 1)
generate_functional("cubic", lambda x: 4*x**3 - 2*x**2 + x)
generate_functional("exponential", lambda x: np.exp(2*x))

# =========================
# 2.3 NON-FUNCTIONAL DATASETS
# =========================

# S-Curve with noise
t = np.linspace(-1, 1, N_POINTS)  
x_s = t
y_s_clean = np.sin(2 * np.pi * t)

for label, sigma in NOISE_LEVELS.items():
    y_s = y_s_clean + np.random.normal(0, sigma, N_POINTS)
    save_csv(x_s, y_s, f"s_curve_noise_{label}.csv")

# Circular arc with noise
t = np.linspace(-np.pi, np.pi, N_POINTS)  # Adjusted to cover the full range of -1 to 1 in x_arc
x_arc = np.linspace(-1, 1, N_POINTS)  # Now x_arc goes from -1 to 1
y_arc_clean = np.sin(np.cos(x_arc))  # Adjusted so t is cos(x_arc)

for label, sigma in NOISE_LEVELS.items():
    y_arc = y_arc_clean + np.random.normal(0, sigma, N_POINTS)
    save_csv(x_arc, y_arc, f"circular_arc_noise_{label}.csv")

# =========================
# 2.4 BAD CASES
# =========================

# 1) Curve with a kink
x1 = np.linspace(0, 0.5, N_POINTS//2)  # First segment from 0,0 to 0.5,1
y1 = 2 * x1

x2 = np.linspace(0.5, 1, N_POINTS//2)  # Second segment from 0.5,1 to 1,1
y2 = np.ones_like(x2)

x_kink = np.concatenate([x1, x2])
y_kink_clean = np.concatenate([y1, y2])

for label, sigma in NOISE_LEVELS.items():
    y_kink = y_kink_clean + np.random.normal(0, sigma, N_POINTS)
    save_csv(x_kink, y_kink, f"kinked_curve_noise_{label}.csv")

x_aux = np.linspace(0, 1, N_POINTS)  

# 2) Heteroscedastic noise (scaled by level)
for label, scale in NOISE_LEVELS.items():
    sigma_x = 0.05 + 2*scale * x_aux  # Define σ(x) as 0.05 + scale * x
    noise_hetero = np.random.normal(0, sigma_x, size=len(x_aux))  # Generate noise with variance σ(x)^2
    y_hetero = 2 * x_aux + 1 + noise_hetero  # Calculate y as f(x) + noise
    save_csv(x_aux, y_hetero, f"heteroscedastic_noise_{label}.csv")

# =========================
# 3) Mixed Generators (CORRECTED)
# =========================

n1 = N_POINTS // 3
n2 = N_POINTS // 3
n3 = N_POINTS - n1 - n2

# Line
x_a = np.linspace(-1, 1, n1)  
y_a_clean = 2*x_a + 1

# Arc
t = np.linspace(-np.pi/4, np.pi/4, n2)  
x_b = np.cos(t)
y_b_clean = np.sin(t)

# Pure scatter (NO extra noise added)
x_c = np.random.uniform(-1, 1, n3)  
y_c = np.random.uniform(-1, 1, n3)  

for label, sigma in NOISE_LEVELS.items():
    # Noise ONLY in the structured parts
    y_a = y_a_clean + np.random.normal(0, sigma, n1)
    y_b = y_b_clean + np.random.normal(0, sigma, n2)

    x_mix = np.concatenate([x_a, x_b, x_c])
    y_mix = np.concatenate([y_a, y_b, y_c])

    save_csv(x_mix, y_mix, f"mixed_generators_noise_{label}.csv")


print("\n✅ ALL datasets (low / medium / high) were generated successfully.")
