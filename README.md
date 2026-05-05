# Comparative Evaluation of Image-Based Spatial Transcriptomics Methods

This repository contains the computational pipelines, data conversion scripts, and customized Python modules used for the thesis: **"Comparative evaluation of image-based spatial transcriptomics methods for subcellular RNA localization patterns."**

## Project Overview
The study focuses on evaluating and mapping the subcellular spatial transcriptomic landscape of mouse spermatogenesis using high-resolution in situ sequencing (10x Genomics Xenium). 
Target cell types include Round Spermatids (`RSptid`), Elongating Spermatids (`ESptid`), and Spermatogonia (`Spgonia`).

## Computational Workflow Diagram

```mermaid
graph TD
    %% Phase 1: Data Integration
    subgraph Data Harmonization
        A[Raw Xenium Data & Seurat Object] -->|1run_seurat_to_anndata.sh| B(Extract to CSV via R)
        B -->|csv_to_anndata.py| C{Unified AnnData: .h5ad}
    end

    %% Phase 2: ELLA Pipeline
    subgraph ELLA Pipeline: Deep Learning Estimation
        C -->|run_ella_pipeline.sh| K[ELLA: Neural Network Training]
        K --> L((ELLA Spatial Predictions))
    end

    %% Phase 3: cellSP Pipeline
    subgraph cellSP Pipeline: Subcellular Patterning
        C -->|cellsp_1import_data.py| D[Filter Target Cells e.g., RSptid]
        
        D -->|cellsp_3run_instant.py| E[InSTAnT: Transcript Localization]
        D -->|cellsp_4run_sprawl.py| F[SPRAWL: Spatial Clustering]
        
        E -->|cellsp_5bicluster_instant.py| G[Bicluster InSTAnT Modules]
        F -->|cellsp_6bicluster_sprawl.py| H[Bicluster SPRAWL Modules]
        
        G --> I[cellsp_8report_creation.py]
        H --> I
        
        I --> J((Subcellular Heatmaps & Colocalization Plots))
    end

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
