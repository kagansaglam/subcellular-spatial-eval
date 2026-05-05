# python

#############################################
## run python from singularity with image with cellSP 
# module load singularity
# singularity exec -B "/mnt/beegfs/amitjavila/" /mnt/beegfs/amitjavila/singularity/260323_xenium_cellsp.sif python
#############################################

## this script runs cellSP using custom anndata created from Seurat as input
## see cellSP tutorial https://github.com/bhavaygg/CellSP/blob/main/figures/tutorial.ipynb
## refer to the tutorial for more info about the different functions

# SET UP -------------------------

## import packages
import os
import cellSP
import argparse

## add command line arguments

### sample
### cellType to filter the data... later I will add an option to filter by cluster
### distance threshold in InSTAnT

### create parser
parser = argparse.ArgumentParser(description="Run InSTAnT and SPRAWL for cellSP")

### add arguments
parser.add_argument("--sample", type=str, help="Sample to analyse, only used for printing")
parser.add_argument("--cellType", type=str, help="Cell Type or Seurat cluster to filter the data")
parser.add_argument("--instant_dist", type=int, help="Distance threshold for pair-pair collocalization in InSTAnT")
parser.add_argument("--threads", type=int, help="Threads to use in run_instant() and run_sprawl()")

### parse arguments
args = parser.parse_args()
sample = args.sample
cellType = args.cellType
instant_dist = args.instant_dist
nthreads = args.threads

print(f"====== Running InSTAnT for sample {sample} and cell type {cellType}... ======")


# IMPORT ANNDATA -------------------

print(f"Importing anndata for sample {sample} and cell type {cellType}...")

## get anndata object filename
filename=f"cellsp_anndata/{sample}/adata_st_{sample}_{cellType}.h5ad"

## import spatial tx anndata object
adata_st = cellSP.ds.load_data(st_adata = filename)

# SUBCELLULAR PATTERN DISCOVERY ----------------

## Add "absZ_raw" in adata.uns["transcripts"], 
##  because it internally converts absZ_raw to absZ when is_sliced=False
##  although I already have absZ, but if not, it gives error
adata_st.uns["transcripts"]["absZ_raw"] = adata_st.uns["transcripts"]["absZ"]

## Run InSTAnT ==> OK
## - set is_sliced=False, because if True, it expects sliced Z axis, and I have a continuous Z axis
print("Running InSTAnT...")
adata_st = cellSP.ch.run_instant(adata_st = adata_st, distance_threshold=instant_dist, threads=nthreads, is_sliced=False)

## save updated undata
print(f"Writing anndata with InSTAnT results for sample {sample} and cell type {cellType}")
cellSP.io.write_h5ad(adata_st, filename)

