#!/bin/bash
#SBATCH --job-name=run_pipeline
#SBATCH --output=run_pipeline_%j.out
#SBATCH --error=run_pipeline_%j.err
#SBATCH --cpus-per-task=10
#SBATCH --mem=250G

Rscript methods/create_lexical_axes.R

Rscript methods/run_query.R \
        --blast=results/p_coccineus_test_data.txt \
        --chemical_space=results/chemical_space.rds \
        --query_kidera=results/query_kidera.rds \
	--mean_vector=results/protein_mean_vector.csv \
        --go=material/lexical_axes_GO.rds \
        --ipr=material/lexical_axes_IPR.rds \
        --hrd=material/lexical_axes_HRD.rds \
        --out=results/ \
        --kcap=50 \
        --threshold=0.7 \
        --combine=sum
