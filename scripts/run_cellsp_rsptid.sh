#!/bin/bash
#SBATCH --job-name=cellsp_RSptid
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128gb
#SBATCH --time=72:00:00
#SBATCH --output=logs/cellsp_RSptid_%j.log

set -e
module load singularity
sing_exec="singularity exec -B /mnt/beegfs/amitjavila/,/mnt/beegfs/ksaglam/ /mnt/beegfs/amitjavila/singularity/260323_xenium_cellsp.sif"

cd ~/my_cellsp_analysis

SAMPLE="TEST4"
CELLTYPE="RSptid"                              
INSTANT_DIST=2

mkdir -p cellsp_anndata/$SAMPLE

echo "--- Starting CellSP Analysis for $SAMPLE and $CELLTYPE ---"

echo ">>> STEP 1: Importing Data..."
$sing_exec python cellsp_1import_data.py --sample $SAMPLE --cellType $CELLTYPE

echo ">>> STEP 3: Running Instant..."
$sing_exec python cellsp_3run_instant.py --sample $SAMPLE --cellType $CELLTYPE --instant_dist $INSTANT_DIST --threads 16

echo ">>> STEP 4: Running Sprawl..."
$sing_exec python cellsp_4run_sprawl.py --sample $SAMPLE --cellType $CELLTYPE --threads 16

echo ">>> STEP 5: Running Bicluster Instant (UMAP BURADA OLABILIR)..."
$sing_exec python cellsp_5bicluster_instant.py --sample $SAMPLE --cellType $CELLTYPE --instant_dist $INSTANT_DIST --threads 16

echo ">>> STEP 6: Running Bicluster Sprawl (VEYA BURADA)..."
$sing_exec python cellsp_6bicluster_sprawl.py --sample $SAMPLE --cellType $CELLTYPE --threads 16

echo ">>> STEP 8: Report Creation..."
$sing_exec python cellsp_8report_creation.py --sample $SAMPLE --cellType $CELLTYPE

echo "--- Analysis Completed Successfully ---"
