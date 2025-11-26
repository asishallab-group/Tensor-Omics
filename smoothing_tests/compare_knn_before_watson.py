import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Datos proporcionados en la salida
df = pd.read_csv('benchmark_results_full.tsv', sep='\t')


# Separar los regímenes para graficar
df_fixed = df[df['k_regime'] == 'FIXED'].sort_values(by='N')
df_relative = df[df['k_regime'] == 'RELATIVE'].sort_values(by='N')
df_root = df[df['k_regime'] == 'ROOT'].sort_values(by='N')

plt.figure(figsize=(10, 6))

# Trazar los puntos en escala log-log
# ETIQUETA CORREGIDA para reflejar el costo de ordenación/indexing en 1D
plt.loglog(df_fixed['N'], df_fixed['seconds'], 'o-', label=r'FIXED ($k=30$): $O(N \log N)$ (Indexing 1D)', color='skyblue')
plt.loglog(df_root['N'], df_root['seconds'], 's-', label=r'ROOT ($k \approx \sqrt{N}$): $O(N^{1.5})$', color='orange')
plt.loglog(df_relative['N'], df_relative['seconds'], '^-', label=r'RELATIVE ($k=0.1\%\cdot N$): $O(N^2)$', color='darkred')

# Añadir líneas de referencia de complejidad
N_ref = np.array([10000, 1000000])

# O(N log N) - Usando el punto de 1M para calibrar
ref_nlogn = df_fixed[df_fixed['N'] == 1000000]['seconds'].iloc[0] / (1000000 * np.log(1000000))
T_nlogn = ref_nlogn * N_ref * np.log(N_ref)
plt.loglog(N_ref, T_nlogn, '--', color='skyblue', alpha=0.6, linewidth=1, label=r'Ref: $O(N \log N)$')

# O(N^2) - Usando el punto de 1M para calibrar
ref_n2 = df_relative[df_relative['N'] == 1000000]['seconds'].iloc[0] / (1000000**2)
T_n2 = ref_n2 * N_ref**2
plt.loglog(N_ref, T_n2, ':', color='darkred', alpha=0.6, linewidth=1, label=r'Ref: $O(N^2)$')

# O(N^1.5) - Usando el punto de 1M para calibrar
ref_n15 = df_root[df_root['N'] == 1000000]['seconds'].iloc[0] / (1000000**1.5)
T_n15 = ref_n15 * N_ref**1.5
plt.loglog(N_ref, T_n15, '-.', color='orange', alpha=0.6, linewidth=1, label=r'Ref: $O(N^{1.5})$')


plt.title('KNN Smoothing Performance Scaling (1D)', fontsize=14)
plt.xlabel(r'Number of Points ($N$) [Log Scale]', fontsize=12)
plt.ylabel('Execution Time (seconds) [Log Scale]', fontsize=12)
plt.legend(loc='upper left', frameon=True, shadow=True)
plt.grid(True, which="both", ls="--", linewidth=0.5)
plt.tight_layout()

plt.savefig("benchmark_scaling.png")
print("benchmark_scaling.png")