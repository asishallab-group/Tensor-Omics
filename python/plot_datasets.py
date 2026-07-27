import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # Necessary for 3D projection
import glob
import os
from sklearn.decomposition import PCA

DIM = 2  # <--- Set the dimension you want to plot (e.g., 3, 5, etc.)

# Dynamic paths based on the dimension
DATA_DIR = f"results/data/{DIM}d"
PLOT_DIR = f"results/plots/examples/{DIM}d"

os.makedirs(PLOT_DIR, exist_ok=True)

csv_files = glob.glob(os.path.join(DATA_DIR, "*.csv"))

for file in csv_files:
    name = os.path.splitext(os.path.basename(file))[0]
    
    if "anwil" in name:
        continue
        
    df = pd.read_csv(file)
    
    # Check if 'y_original' exists in the dataset
    if "y_original" not in df.columns:
        print(f"Skipping {name}: 'y_original' column not found.")
        continue

    # Identify feature columns (everything except the target y)
    features = [col for col in df.columns if col != "y_original"]
    num_features = len(features)

    # --- CASE 1: STANDARD 2D DATASET (1 Feature + 1 Target) ---
    if num_features == 1:
        plt.figure(figsize=(6, 5))
        plt.scatter(df[features[0]], df["y_original"], s=10, alpha=0.7)
        plt.xlabel(features[0])
        plt.ylabel("y_original")
        plt.title(name.replace("_", " "))
        plt.grid(True)

    # --- NEW CASE 2: 3D DATASET (2 Features + 1 Target) ---
    # We can visualize this directly using Matplotlib's 3D toolkit
    elif num_features == 2:
        # Initialize a figure with 3D projection
        fig = plt.figure(figsize=(8, 7))
        ax = fig.add_subplot(111, projection='3d')
        
        # Extract features and target
        x1 = df[features[0]]
        x2 = df[features[1]]
        y = df["y_original"]
        
        # Plot points. We use a color map based on Y for better depth perception
        scatter = ax.scatter(x1, x2, y, c=y, cmap='plasma', s=10, alpha=0.6)
        
        # Set proper labels
        ax.set_xlabel(features[0])
        ax.set_ylabel(features[1])
        ax.set_zlabel("y_original")
        ax.set_title(name.replace("_", " ") + " (3D View)")
        
        # Add a color bar to show the mapping of Y values
        fig.colorbar(scatter, ax=ax, label='Target Value')
        
        # Adjust the starting viewing angle (elevation, azimuth)
        # 30, -60 is a good default to show both features and the target
        ax.view_init(elev=30, azim=-60)
        
        plt.tight_layout()

    # --- CASE 3: HIGH-DIMENSIONAL DATASET (D > 3, or Features > 2) ---
    # Continue using paired plots and PCA for D > 3
    else:
        # Create a grid of subplots: One row for Each Feature vs Y + One for PCA
        fig, axes = plt.subplots(1, num_features + 1, figsize=(4 * (num_features + 1), 4))
        fig.suptitle(name.replace("_", " "), fontsize=14, fontweight='bold')
        
        # 1. Plot each feature individual relation against Y
        for i, col in enumerate(features):
            axes[i].scatter(df[col], df["y_original"], s=8, alpha=0.6, color='blue')
            axes[i].set_xlabel(col)
            axes[i].set_ylabel("y_original")
            axes[i].set_title(f"{col} vs Target")
            axes[i].grid(True)
            
        # 2. PCA Trick: Compress all features into 1D, and plot against Y
        pca = PCA(n_components=1)
        X_compressed = pca.fit_transform(df[features])
        
        axes[-1].scatter(X_compressed, df["y_original"], s=8, alpha=0.6, color='purple')
        axes[-1].set_xlabel("First Principal Component (PCA 1)")
        axes[-1].set_ylabel("y_original")
        axes[-1].set_title("Compressed Space vs Target")
        axes[-1].grid(True)
        
        plt.tight_layout()

    # Save the final image
    output_path = os.path.join(PLOT_DIR, f"{name}.png")
    plt.savefig(output_path, dpi=150)
    plt.close()

    print(f"Plot saved: {output_path} (Detected structure: {num_features}D features + 1D target)")

print("\nAll plots were generated successfully.")