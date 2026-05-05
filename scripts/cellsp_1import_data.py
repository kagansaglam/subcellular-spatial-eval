# python

#############################################
## run python from singularity with image with cellSP
# module load singularity
# singularity exec -B "/mnt/beegfs/amitjavila/" /mnt/beegfs/amitjavila/singularity/260323_xenium_cellsp.sif python
#############################################

## this script runs cellSP using custom anndata created from Seurat as input
## see cellSP tutorial https://github.com/bhavaygg/CellSP/blob/main/figures/tutorial.ipynb
## refer to the tutorial for more info about the different functions

## import packages
import os
import cellSP
import argparse

## get current directory
os.getcwd()

## add command line arguments

### sample
### cellType to filter the data... later I will add an option to filter by cluster

### create parser
parser = argparse.ArgumentParser(description="import and filter anndata for cellSP")

### add arguments
parser.add_argument("--sample", type=str, help="Sample to analyse, only used for printing")
parser.add_argument("--cellType", type=str, help="Cell Type or Seurat cluster to filter the data")

### parse arguments
args = parser.parse_args()
sample = args.sample
cellType = args.cellType

# IMPORT ANNDATA --------------

print(f"Importing anndata for sample {sample} and filtering by cellType/cluster {cellType}...")

## define input filename
infilename = f"seurat_to_anndata/{sample}/{sample}_anndata.h5ad"

## import spatial tx anndata object
adata_st = cellSP.ds.load_data(st_adata = infilename)

# SUBSET ANNDATA BY CELL TYPE -------------------

## subset by cell type
subsets = {
    ct: adata_st[adata_st.obs["cell_type"] == ct].copy()
    for ct in adata_st.obs["cell_type"].unique()
}

## take desired subset
subset = subsets[cellType]

## get the cells in the subset
cells = subset.obs["uID"]

## filter the distinct adata layers by cells present in subset
subset.uns["transcripts"] = adata_st.uns["transcripts"].loc[adata_st.uns["transcripts"]['uID'].isin(cells)]
subset.uns["cell_boundary"] = adata_st.uns["cell_boundary"].loc[adata_st.uns["cell_boundary"].index.isin(cells)]

# WRITE NEW ADATA TO FILE -----------------------------

print(f"Writting filtered anndata for sample {sample} and cell type {cellType}")

# ## create directory
# os.mkdir(f"cellsp_anndata/{sample}")

## write adata
filename=f"cellsp_anndata/{sample}/adata_st_{sample}_{cellType}.h5ad"
cellSP.io.write_h5ad(subset, filename)
