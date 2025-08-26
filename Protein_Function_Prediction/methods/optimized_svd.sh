#!/bin/bash
#SBATCH --job-name=optimized_svd
#SBATCH --output=optimized_svd_%j.out
#SBATCH --error=optimized_svd_%j.err
#SBATCH --cpus-per-task=64
#SBATCH --mem=500G

echo "Starting optimized SVD script"
python3 methods/optimized_svd.py

echo "Finished!"
