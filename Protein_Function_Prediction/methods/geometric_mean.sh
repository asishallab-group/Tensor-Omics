#!/bin/bash
#SBATCH --job-name=geomean
#SBATCH --output=geomean_%j.out
#SBATCH --error=geomean_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=100G

echo "Starting script"
python3 methods/geometric_mean.py

echo "Finished!"
