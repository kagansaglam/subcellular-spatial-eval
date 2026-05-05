# -*- coding: utf-8 -*-
import os
import cellSP
import argparse
import numpy as np
import pandas as pd
import scipy.sparse as sp

parser = argparse.ArgumentParser()
parser.add_argument("--sample", type=str)
parser.add_argument("--cellType", type=str)
args = parser.parse_args()
sample = args.sample
cellType = args.cellType

print(f"Importing anndata for {sample} - {cellType}...")
filename = f"cellsp_anndata/{sample}/adata_st_{sample}_{cellType}.h5ad"
adata_st = cellSP.ds.load_data(st_adata=filename)

print(">>> NUKING arrays for UMAP (Absolute Flattening) <<<")

# NUKLEER FONKSIYON: Ne gelirse gelsin 2D Float32 matrisine zorla cevirir
def nuke_array(arr):
    if sp.issparse(arr): 
        arr = arr.toarray()
    clean_rows = []
    for row in arr:
        # Eger satirin icinde baska liste varsa onu da dumduz et (np.ravel)
        clean_rows.append(np.ravel(row))
    return np.vstack(clean_rows).astype(np.float32)

# 1. Asil Matrisi (X) Ez
adata_st.X = nuke_array(adata_st.X)

# 2. Diger Matrisleri (obsm) Ez - Hata vereni atla
for k in list(adata_st.obsm.keys()):
    try:
        adata_st.obsm[k] = nuke_array(adata_st.obsm[k])
    except Exception:
        pass

# 3. Metadata'daki nesneleri string'e cevir
for col in adata_st.obs.columns:
    if adata_st.obs[col].dtype == 'object':
        adata_st.obs[col] = adata_st.obs[col].astype(str)

print(">>> Running module visualization... <<<")
cellSP.vs.create_report(adata_st)

print("--- Report Done! ---")
