#!/bin/bash
#SBATCH --job-name=create_lexical_axes
#SBATCH --output=create_lexical_axes_%j.out
#SBATCH --error=create_lexical_axes_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=50G

Rscript /media/BioNAS2/Tensor_Omics/Protein_Function_Prediction/methods/create_lexical_axes.R
