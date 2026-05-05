#!/bin/bash
#SBATCH --job-name=run_ella
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128gb
#SBATCH --time=72:00:00
#SBATCH --output=logs/run_ella_%j.log

module load singularity

WORKDIR=/mnt/beegfs/ksaglam/my_cellsp_analysis
cd $WORKDIR

sing_exec="singularity exec -B /mnt/beegfs/amitjavila/,/mnt/beegfs/ksaglam/ /mnt/beegfs/amitjavila/singularity/260420_xenium_ella.sif"

echo "========================================="
echo "Adim 1: Veri ELLA formatina cevriliyor (Prepare Data)..."
echo "========================================="
$sing_exec python -m ella.data.prepare_data -i ella_inputs/data.pkl -o prepared_data

echo "========================================="
echo "Adim 2: Yapay Zeka Modeli Egitiliyor (Train)..."
echo "========================================="
$sing_exec ella-train --config-name mini_demo

echo "========================================="
echo "Adim 3: Sonuclar Cekiliyor (Estimate)..."
echo "========================================="
# Sonuclarin kaydedilecegi klasoru olusturalim
mkdir -p ella_results

# Kılavuzdaki tahmın komutunu calistiriyoruz
$sing_exec ella-estimate -d lightning_logs/debug -p "gene_0-kernel_.*" -b 10 -o ella_results/out.json

echo "--- ELLA ANALIZI TAMAMLANDI! (Tebrikler!) ---"
