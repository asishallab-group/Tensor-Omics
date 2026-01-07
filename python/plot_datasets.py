import pandas as pd
import matplotlib.pyplot as plt
import glob
import os

DATA_DIR = "results/data"
PLOT_DIR = "results/plots"

os.makedirs(PLOT_DIR, exist_ok=True)

csv_files = glob.glob(os.path.join(DATA_DIR, "*.csv"))

for file in csv_files:
    name = os.path.splitext(os.path.basename(file))[0]

    df = pd.read_csv(file)

    plt.figure(figsize=(6, 5))
    plt.scatter(df["x"], df["y_original"], s=10)
    plt.xlabel("x")
    plt.ylabel("y_original")
    plt.title(name.replace("_", " "))
    plt.grid(True)

    output_path = os.path.join(PLOT_DIR, f"{name}.png")
    plt.savefig(output_path, dpi=150)
    plt.close()

    print(f"Plot saved: {output_path}")

print("\nAll plots were generated.")
