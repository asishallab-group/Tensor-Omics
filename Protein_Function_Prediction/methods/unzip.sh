#!/bin/bash
#SBATCH --job-name=unzip
#SBATCH --output=unzip_%j.out
#SBATCH --error=unzip_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5

gunzip -c /media/BioNAS/ag_hallab/Tensor_Omics/Protein_Function_Prediction/material/20250604/protein2ipr.dat.gz \
	> /media/BioNAS/ag_hallab/Tensor_Omics/Protein_Function_Prediction/material/protein2ipr.dat

# gzip -dc /media/BioNAS/ag_hallab/Tensor_Omics/Protein_Function_Prediction/material/uniref90.fasta.gz \
#         > /media/BioNAS/ag_hallab/Tensor_Omics/Protein_Function_Prediction/material/uniref90.fasta

