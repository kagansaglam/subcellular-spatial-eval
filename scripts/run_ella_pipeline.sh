#!/bin/bash
#SBATCH --job-name=ella_TEST4
#SBATCH --partition=test-gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64gb
#SBATCH --time=24:00:00
#SBATCH --gres=gpu:1
#SBATCH --output=logs/ella_TEST4_%j.log

mkdir -p logs

module load singularity

SIF_IMAGE="/mnt/beegfs/amitjavila/singularity/260427_xenium_ella.sif"
SING_BIND="-B /mnt/beegfs/amitjavila/,/mnt/beegfs/ksaglam/"

# COK KRITIK: Singularity'nin GPU'yu kullanabilmesi icin --nv bayragini ekledik!
EXEC="singularity exec --nv $SING_BIND $SIF_IMAGE"

BASE_DIR=$(pwd)

echo "--- Starting ELLA Pipeline (GPU Accelerated) ---"

# STEP 1: Data Preparation
echo ">>> Step 1: Preparing data..."
$EXEC python -m ella.data.prepare_data -i $BASE_DIR/ella_inputs/fixed_data.pkl -o prepared_data

# STEP 2: Training (GPU + 3 Workers)
echo ">>> Step 2: Running training model..."
$EXEC ella-train --config-name mini_demo \
    data.data_path="$BASE_DIR/prepared_data/training_data.jsonl" \
    log.save_dir="$BASE_DIR/lightning_logs/debug" \
    estimation.search_dir="$BASE_DIR/lightning_logs/debug" \
    estimation.output_path="$BASE_DIR/lightning_logs/debug/gene_0_estimation_result.json" \
    data.num_workers=3

# STEP 3: Estimation
echo ">>> Step 3: Running estimation..."
$EXEC ella-estimate -d lightning_logs/debug -p "gene_0-kernel_.*" -b 10 -o lightning_logs/debug/out.json

# STEP 4: Postprocessing
echo ">>> Step 4: Running postprocessing..."
$EXEC python ella_postprocess.py

echo "--- Pipeline Execution Completed ---"
