# python

#############################################
## run python from singularity with image with cellSP 
# module load singularity
# singularity exec -B "/mnt/beegfs/amitjavila/" /mnt/beegfs/amitjavila/singularity/260323_xenium_cellsp.sif python
#############################################

## this script runs cellSP using custom anndata created from Seurat as input
## see cellSP tutorial https://github.com/bhavaygg/CellSP/blob/main/figures/tutorial.ipynb
## refer to the tutorial for more info about the different functions

# SET UP -----------------------

## import packages
import os
import cellSP
import argparse

## add command line arguments

### sample
### cellType to filter the data... later I will add an option to filter by cluster
### distance threshold in InSTAnT

### create parser
parser = argparse.ArgumentParser(description="Run SPRAWL biclustering for cellSP")

### add arguments
parser.add_argument("--sample", type=str, help="Sample to analyse, only used for printing")
parser.add_argument("--cellType", type=str, help="Cell Type or Seurat cluster to filter the data")
parser.add_argument("--methods", type=list, nargs="+", default=['Peripheral', 'Radial', 'Punctate', 'Central'], 
                    help="Patters SPRAWL looks for")
parser.add_argument("--threads", type=int, help="Number of threads used in biclustering")

### parse arguments
args = parser.parse_args()
sample = args.sample
cellType = args.cellType
sprawl_methods = args.methods
nthreads = args.threads

print(f"====== Biclustering SPRAWL for sample {sample} and cell type {cellType}... ======")

# IMPORT ANNDATA -------------------

print(f"Importing anndata for sample {sample} and cell type {cellType}...")

## get anndata object filename
filename=f"cellsp_anndata/{sample}/adata_st_{sample}_{cellType}.h5ad"

## import spatial tx anndata object
adata_st = cellSP.ds.load_data(st_adata = filename)

# MODULE DISCOVERY ------------------------

print(f"Running SPRAWL biclustering for sample {sample} and cell type {cellType}...") 

## Perform biclustering of SPRAWL output gene subcellular pattern matrices
adata_st = cellSP.ch.bicluster_sprawl(adata_st, methods=sprawl_methods, threads=nthreads, num_biclusters = 10, randomized_searches=5000)

# ## We can inspect the results of the biclustering looking at adata_st.uns
# adata_st.uns['sprawl_biclustering']

# Save updated anndata object
print(f"Writing anndata with SPRAWL biclustering results for sample {sample} and cell type {cellType}")
cellSP.io.write_h5ad(adata_st, filename)

