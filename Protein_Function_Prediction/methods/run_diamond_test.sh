#!/bin/bash
#SBATCH --job-name=diamond
#SBATCH --output=diamond%j.out
#SBATCH --cpus-per-task=20
#SBATCH --mem=120G

echo "Starting diamond script"

# echo "creating protein database"
# diamond makedb --in material/uniref50_morethan1.fasta -d results/db_uniref50_morethan1

echo "executing diamond analysis"
diamond blastp -q material/p_coccineus_different_root.faa -d results/db_uniref50_morethan1.dmnd -o results/p_coccineus_test_data.txt -f6 qseqid sseqid qstart qend qlen sstart send slen pident evalue bitscore --more-sensitive -p 20 --quiet --evalue 1e-5

echo  "diamond finished!"
