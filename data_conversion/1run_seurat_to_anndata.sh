#!/bin/bash
#SBATCH --job-name=seurat_to_anndata
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32gb
#SBATCH --time=72:00:00
#SBATCH --array=0
#SBATCH --output=logs/seurat_to_anndata_%j.log

# Load modules 
module load singularity 
module load R-bundle-Bioconductor

# Singularity komutu
sing_exec="singularity exec -B /mnt/beegfs/ /mnt/beegfs/amitjavila/singularity/260323_xenium_cellsp.sif"

# CALISMA KLASORUN
WORKDIR=~/my_cellsp_analysis
cd $WORKDIR

# ORNEK ADI
SAMPLE="TEST4" 

# 1. Adim: Seurat objesini CSV'lere cevir
Rscript $WORKDIR/seurat_to_csv.R $SAMPLE

# 2. Adim: CSV'leri h5ad anndata objesine cevir
$sing_exec python $WORKDIR/csv_to_anndata.py --sample $SAMPLE
