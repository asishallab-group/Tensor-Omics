#!/bin/bash
#SBATCH --job-name=create_chemical_space
#SBATCH --output=create_chemical_space_%j.out
#SBATCH --error=create_chemical_space_%j.err
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=8
#SBATCH --mem=150G

Rscript create_chemical_space.R \
	--svd ../results/svd_reduced_matrix.csv \
	--kidera ../material/kidera_factors_uniref50_morethan1_ref_prots.rds \
	--out ../results/chemical_space.rds
