# Comparative Evaluation of Image-Based Spatial Transcriptomics Methods

This repository contains the computational pipelines, data conversion scripts, and customized Python modules used for the thesis: **"Comparative evaluation of image-based spatial transcriptomics methods for subcellular RNA localization patterns."**

## Project Overview
The study focuses on evaluating and mapping the subcellular spatial transcriptomic landscape of mouse spermatogenesis using high-resolution in situ sequencing (10x Genomics Xenium). 
Target cell types include Round Spermatids (`RSptid`), Elongating Spermatids (`ESptid`), and Spermatogonia (`Spgonia`).

## Computational Workflows

### 1. Data Conversion (R to Python)
* Scripts: `data_conversion/seurat_to_csv.R` and `1run_seurat_to_anndata.sh`
* Bridging the ecosystem gap by converting pre-processed spatial Seurat objects into Scanpy-compatible AnnData (`.h5ad`) formats for downstream deep learning and spatial mapping.

### 2. ELLA (Deep Learning Spatial Model)
* Script: `scripts/run_ella_pipeline.sh`
* Accelerated using Nvidia A40 GPUs via Slurm workload manager to estimate high-dimensional spatial gene expression patterns.

### 3. cellSP (Subcellular Spatial Patterning)
* Script: `scripts/cellsp_8report_creation.py`
* **Methodological Fix:** Addressed strict `64-dimension` limit errors (`ValueError: setting an array element with a sequence`) inherent to sparse Xenium outputs during UMAP generation. A custom "nuking" function (`nuke_array`) was developed to forcefully flatten nested structures into continuous 2D `float32` arrays, successfully bypassing algorithmic limitations.
