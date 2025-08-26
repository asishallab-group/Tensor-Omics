#!/bin/bash
#SBATCH --job-name=run_pipeline_no_idf
#SBATCH --output=run_pipeline_no_idf%j.out
#SBATCH --error=run_pipeline_no_idf%j.err
#SBATCH --cpus-per-task=10
#SBATCH --mem=250G

Rscript methods/create_lexical_axes_no_idf.R

Rscript methods/run_query.R \
        --blast=results/p_coccineus_test_data.txt \
        --chemical_space=results/chemical_space.rds \
        --query_kidera=results/query_kidera.rds \
	--mean_vector=results/protein_mean_vector.csv \
        --go=material/lexical_axes_GO_no_idf.rds \
        --ipr=material/lexical_axes_IPR_no_idf.rds \
        --hrd=material/lexical_axes_HRD_no_idf.rds \
        --out=results_no_idf/ \
        --kcap=50 \
        --threshold=0.7 \
        --combine=sum
