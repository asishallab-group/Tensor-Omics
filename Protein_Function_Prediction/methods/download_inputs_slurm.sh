#!/bin/bash
#SBATCH --job-name=uniref90_download
#SBATCH --output=uniref90_download_%j.log
#SBATCH --error=uniref90_download_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G

echo "Download started: $(date)"

wget -P /media/BioNAS/ag_hallab/Tensor_Omics/Protein_Function_Prediction/material https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz
wget -P /media/BioNAS/ag_hallab/Tensor_Omics/Protein_Function_Prediction/material https://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_all.gaf.gz
echo "Download finished: $(date)"

# unzip files
gzip -dc /media/BioNAS/ag_hallab/Tensor_Omics/Protein_Function_Prediction/material/uniref90.fasta.gz \
	> /media/BioNAS/ag_hallab/Tensor_Omics/Protein_Function_Prediction/material/uniref90.fasta

gzip -dc /media/BioNAS/ag_hallab/Tensor_Omics/Protein_Function_Prediction/material/goa_uniprot_all.gaf.gz \
	> /media/BioNAS/ag_hallab/Tensor_Omics/Protein_Function_Prediction/material/goa_uniprot_all.gaf

