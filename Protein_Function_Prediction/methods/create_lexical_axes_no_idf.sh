#!/bin/bash
#SBATCH --job-name=create_lexical_axes_no_idf
#SBATCH --output=create_lexical_axes_no_idf_%j.out
#SBATCH --error=create_lexical_axes_no_idf_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=100G

Rscript /media/BioNAS2/Tensor_Omics/Protein_Function_Prediction/methods/create_lexical_axes_no_idf.R
