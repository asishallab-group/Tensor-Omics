#!/bin/bash
#SBATCH --job-name=svd
#SBATCH --output=svd_%j.out
#SBATCH --error=svd_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=100G

echo "Starting script"
python3 methods/svd.py

echo "Finished!"
