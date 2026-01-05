import numpy as np
import pandas as pd
import os

# =========================
# CONFIGURACIÓN GLOBAL
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
# UTILIDAD PARA GUARDAR CSV
# =========================
def save_csv(x, y, filename):
    df = pd.DataFrame({
        "x": x,
        "y_original": y
    })
    df.to_csv(os.path.join(OUTPUT_DIR, filename), index=False)
    print(f"Saved: {filename}")

# =========================
# 2.2 DATASETS FUNCIONALES
# =========================

x = np.linspace(-1, 1, N_POINTS)  # Cambiado para incluir valores negativos

def generate_functional(name, f):
    for label, sigma in NOISE_LEVELS.items():
        noise = np.random.normal(0, sigma, N_POINTS)
        y = f(x) + noise
        save_csv(x, y, f"{name}_noise_{label}.csv")

generate_functional("linear", lambda x: 2*x + 1)
generate_functional("cubic", lambda x: 4*x**3 - 2*x**2 + x)
generate_functional("exponential", lambda x: np.exp(2*x))

# =========================
# 2.3 NO FUNCIONALES
# =========================

# S-Curve con ruido
t = np.linspace(-1, 1, N_POINTS)  # Cambiado para incluir valores negativos
x_s = t
y_s_clean = np.sin(2 * np.pi * t)

for label, sigma in NOISE_LEVELS.items():
    y_s = y_s_clean + np.random.normal(0, sigma, N_POINTS)
    save_csv(x_s, y_s, f"s_curve_noise_{label}.csv")

# Arco circular con ruido
t = np.linspace(-np.pi, np.pi, N_POINTS)  # Ajustado para cubrir todo el rango de -1 a 1 en x_arc
x_arc = np.linspace(-1, 1, N_POINTS)  # Ahora x_arc va de -1 a 1
y_arc_clean = np.sin(np.cos(x_arc))  # Ajustado para que t sea cos(x_arc)

for label, sigma in NOISE_LEVELS.items():
    y_arc = y_arc_clean + np.random.normal(0, sigma, N_POINTS)
    save_csv(x_arc, y_arc, f"circular_arc_noise_{label}.csv")

# =========================
# 2.4 BAD CASES
# =========================

# 1) Curva con quiebre (kink)
x1 = np.linspace(0, 0.5, N_POINTS//2)  # Primer segmento de 0,0 a 0.5,1
y1 = 2 * x1

x2 = np.linspace(0.5, 1, N_POINTS//2)  # Segundo segmento de 0.5,1 a 1,1
y2 = np.ones_like(x2)

x_kink = np.concatenate([x1, x2])
y_kink_clean = np.concatenate([y1, y2])

for label, sigma in NOISE_LEVELS.items():
    y_kink = y_kink_clean + np.random.normal(0, sigma, N_POINTS)
    save_csv(x_kink, y_kink, f"kinked_curve_noise_{label}.csv")

x_aux = np.linspace(0, 1, N_POINTS)  # Cambiado para incluir valores negativos

# 2) Ruido heteroscedástico (escalado por nivel)
for label, scale in NOISE_LEVELS.items():
    sigma_x = 0.05 + 2*scale * x_aux  # Definir σ(x) como 0.05 + scale * x
    noise_hetero = np.random.normal(0, sigma_x, size=len(x_aux))  # Generar ruido con varianza σ(x)^2
    y_hetero = 2 * x_aux + 1 + noise_hetero  # Calcular y como f(x) + ruido
    save_csv(x_aux, y_hetero, f"heteroscedastic_noise_{label}.csv")

# =========================
# 3) Generadores mixtos (CORREGIDO)
# =========================

n1 = N_POINTS // 3
n2 = N_POINTS // 3
n3 = N_POINTS - n1 - n2

# Línea
x_a = np.linspace(-1, 1, n1)  # Cambiado para incluir valores negativos
y_a_clean = 2*x_a + 1

# Arco
t = np.linspace(-np.pi/4, np.pi/4, n2)  # Cambiado para incluir valores negativos
x_b = np.cos(t)
y_b_clean = np.sin(t)

# Scatter puro (NO se le añade ruido extra)
x_c = np.random.uniform(-1, 1, n3)  # Cambiado para incluir valores negativos
y_c = np.random.uniform(-1, 1, n3)  # Cambiado para incluir valores negativos

for label, sigma in NOISE_LEVELS.items():
    # Ruido SOLO en las partes con estructura
    y_a = y_a_clean + np.random.normal(0, sigma, n1)
    y_b = y_b_clean + np.random.normal(0, sigma, n2)

    x_mix = np.concatenate([x_a, x_b, x_c])
    y_mix = np.concatenate([y_a, y_b, y_c])

    save_csv(x_mix, y_mix, f"mixed_generators_noise_{label}.csv")


print("\n✅ TODOS los datasets (low / medium / high) fueron generados correctamente.")
