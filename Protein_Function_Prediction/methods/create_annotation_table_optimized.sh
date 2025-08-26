#!/bin/bash
#SBATCH --job-name=create_tables
#SBATCH --output=create_tables_%j.out
#SBATCH --error=create_tables_%j.err
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=10
#SBATCH --mem=400G

Rscript /media/BioNAS2/Tensor_Omics/Protein_Function_Prediction/methods/create_annotation_table_optimized.R
Rscript /media/BioNAS2/Tensor_Omics/Protein_Function_Prediction/methods/create_lexical_axes.R
