#!/bin/bash
#SBATCH --job-name=run_cellsp
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128gb
#SBATCH --time=72:00:00
#SBATCH --output=logs/run_cellsp_cell_type_esptid_%j.log 

# Load modules 
module load singularity

# Define singularity commands
sing_exec="singularity exec -B /mnt/beegfs/amitjavila/,/mnt/beegfs/ksaglam/ /mnt/beegfs/amitjavila/singularity/260323_xenium_cellsp.sif"

# Define SAMPLE and CELLTYPE
SAMPLE="TEST4"
CELLTYPE="ESptid"
INSTANT_DIST=2
SPRAWL_METHODS=('Peripheral' 'Radial' 'Punctate' 'Central')
THREADS=16

# Define WORKDIR
WORKDIR=/mnt/beegfs/ksaglam/my_cellsp_analysis
cd $WORKDIR

# Create OUTDIR
OUTDIR=cellsp_anndata
mkdir -p cellsp_anndata/$SAMPLE

# TUM ADIMLAR AKTIF EDILDI (Python dosyalari ana klasorde)
echo "Adim 1: Import Data basliyor..."
$sing_exec python $WORKDIR/cellsp_1import_data.py --sample $SAMPLE --cellType $CELLTYPE

echo "Adim 3: InSTAnT basliyor..."
$sing_exec python $WORKDIR/cellsp_3run_instant.py --sample $SAMPLE --cellType $CELLTYPE --instant_dist $INSTANT_DIST --threads $THREADS

echo "Adim 4: SPRAWL basliyor..."
$sing_exec python $WORKDIR/cellsp_4run_sprawl.py --sample $SAMPLE --cellType $CELLTYPE --methods ${SPRAWL_METHODS[@]} --threads $THREADS

echo "Adim 5: Bicluster InSTAnT basliyor..."
$sing_exec python $WORKDIR/cellsp_5bicluster_instant.py --sample $SAMPLE --cellType $CELLTYPE --instant_dist $INSTANT_DIST --threads $THREADS

echo "Adim 6: Bicluster SPRAWL basliyor..."
$sing_exec python $WORKDIR/cellsp_6bicluster_sprawl.py --sample $SAMPLE --cellType $CELLTYPE --methods ${SPRAWL_METHODS[@]} --threads $THREADS

echo "Adim 8: Rapor Olusturma basliyor..."
$sing_exec python $WORKDIR/cellsp_8report_creation.py --sample $SAMPLE --cellType $CELLTYPE

echo "--- TUM CELLSP ANALIZI TAMAMLANDI ---"
