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

# IMPORT ANNDATA
print("Importing anndata...")

## import spatial tx anndata object
adata_st = cellSP.ds.load_data(st_adata = f"adata_st.h5ad")

#############################################
# # SINGLE CELL IMPUTATION
# ### we cannot import single cell data, as we don't have it... it would be nice to have though 
# ### we can download public data and analyse it with seurat and then transform it to anndata as it has been done with the spatial tx
# ### until then, we cannot run the denoising of spatialdata and imputation of genes not in st data
# ## import single-cell data
# adata_sc = cellSP.ds.load_data(sc_adata= 'adata_sc.h5ad')
# ## denoise single-cell data
# adata_sc = cellSP.pp.impute(adata_sc, t="auto")
# ## impute genes not in st data
# adata_st = cellSP.pp.run_tangram(adata_sc, adata_st, device='cuda')
# ## save updated anndata object
# cellSP.io.write_h5ad(adata_st, f"adata_st.h5ad")
#############################################
