#!/bin/bash
#SBATCH --job-name=calculate_kidera
#SBATCH --output=calculate_kidera_%j.out
#SBATCH --error=calculate_kidera_%j.err
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=10
#SBATCH --mem=200G

Rscript /media/BioNAS2/Tensor_Omics/Protein_Function_Prediction/methods/calculate_kidera.R
