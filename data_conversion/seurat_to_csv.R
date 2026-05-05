# R

# This script takes Xenium data from a Seurat object and raw files (transcripts, cells and cell boundaries)
# and writes expression, cell metadata, gene metadata, transcripts coords, cell boundaries and cell centroids to CSV
# R versiyonuna uygun olan tek bir Adria kutuphanesini seciyoruz
current_r_version <- paste0(R.Version()$major, ".", strsplit(R.Version()$minor, "\\.")[[1]][1])
adria_path <- paste0("/mnt/beegfs/amitjavila/R/x86_64-pc-linux-gnu-library/", current_r_version)
.libPaths(c(adria_path, .libPaths()))

library(dplyr)
library(purrr)
library(stringr)
library(Seurat)
library(SeuratObject)
library(arrow)
library(data.table)

args = commandArgs(T)
name <- args[1] #name <- "TEST4"
# seurat_object <- args[2] # seurat_object <- "~/projects/pic_variation/5nxf1_iap_B6CAST/3xenium_B6CAST/output/samples/seurat_objects/TEST4.seuratObject.RDS"
# transcripts_parquet <- args[3] # transcripts_parquet <- "~/projects/pic_variation/0data/private/Xenium_B6CAST/output-XETG00289__0052444__TEST4__20250213__120501//transcripts.parquet"
# cell_boundaries_parquet <- args[4] # cell_boundaries_parquet <- "~/projects/pic_variation/0data/private/Xenium_B6CAST/output-XETG00289__0052444__TEST4__20250213__120501//cell_boundaries.parquet"
# cells_parquet <- args[5] #cells_parquet <- "~/projects/pic_variation/0data/private/Xenium_B6CAST/output-XETG00289__0052444__TEST4__20250213__120501/cells.parquet"
# outdir <- args[6] #outdir <- "~/projects/pic_variation/5nxf1_iap_B6CAST/3xenium_B6CAST/test/cellsp/seurat_to_anndata"

# define input files based on name 
seurat_object <- paste0("/mnt/beegfs/amitjavila/projects/pic_variation/5nxf1_iap_B6CAST/3xenium_B6CAST/output/samples/seurat_objects/",name,".seuratObject.RDS")
transcripts_parquet <- list.files("/mnt/beegfs/amitjavila/projects/pic_variation/0data/private/Xenium_B6CAST/", "transcripts.parquet", recursive=T, full.names=T) %>% purrr::keep(stringr::str_detect(., paste0(name,"__")))
cell_boundaries_parquet <- list.files("/mnt/beegfs/amitjavila/projects/pic_variation/0data/private/Xenium_B6CAST/", "cell_boundaries.parquet", recursive=T, full.names=T) %>% purrr::keep(stringr::str_detect(., paste0(name,"__")))
cells_parquet <- list.files("/mnt/beegfs/amitjavila/projects/pic_variation/0data/private/Xenium_B6CAST/", "cells.parquet", recursive=T, full.names=T) %>% purrr::keep(stringr::str_detect(., paste0(name,"__")))

outdir <- paste0("~/my_cellsp_analysis/seurat_to_anndata/",name)
print(paste0("Starting Seurat to CSV for sample :", name))
print(paste0("Seurat object: ", seurat_object))
print(paste0("Raw transcripts: ", transcripts_parquet))
print(paste0("Cell boundaries: ", cell_boundaries_parquet))
print(paste0("Cell coordinates: ", cells_parquet))
print(paste0("Output directory: ", outdir))


# IMPORT SEURAT ======================================

print("Importing Seurat object")

## Import Seurat Object
## Updating SeuratObject because it was serialized with an older version of SeuratObject
## Set active identity (cell type)
x_seu <- readRDS(seurat_object)
x_seu <- SeuratObject::UpdateSeuratObject(x_seu)
Seurat::Idents(x_seu) <- "cell.type"
head(x_seu@active.ident)

# # CROP SEURAT OBJECT =================================
# ### ONLY FOR TESTING PURPOSES, BECAUSE I WILL TAKE THE FULL OBJECT RUN CELLSP IN EACH CELL TYPE OR CLUSTER SEPARATELY
# print("Cropping Seurat object to reduce size of object")
#
# ## increase max gb for parallelization
# options(future.globals.maxSize = 5 * 1024^3)  # 5 GB
#
# ## crop using plot coordinates
# # cropped.coords <- SeuratObject::Crop(x_seu[["fov"]], x = c(1000, 2000), y = c(2000, 3000), coords = "plot")
# cropped.coords <- SeuratObject::Crop(x_seu[["fov"]], x = c(0, 4000), y = c(0, 4500), coords = "plot")
#
# ## restart max gb for parallelization
# options(future.globals.maxSize = 1 * 1024^3)  # 1 GB
#
# ## add a new fov
# x_seu[["zoom"]] <- cropped.coords
#
# # ## plot imgdimplot from the new zoom fov
# # Seurat::ImageDimPlot(x_seu, fov = "zoom", group.by = "ident", axes = T, border.color = NA, border.size = .01)
#
# # Get barcodes of cells in cropped plot to filter the other 
# crop_cell_ids <- Seurat::Cells(cropped.coords)
# crop_barcodes <- sub("_.*","",crop_cell_ids)

# GET EXPRESSION DATA ========================

print("Getting expression data")

## export expression data from seurat object
## format cell names to include only barcodes withouth the sample
# ## get the expression from only the cropped cells
## transpose to get cells in rows and genes in cols
## change dashes and dots for underscores in gene names because it can cause issues in instant 
# expr <- Seurat::GetAssayData(x_seu, layer = "counts")  # or "data"
expr <- Seurat::GetAssayData(x_seu, layer = "data")  # or "data"
colnames(expr) <- sub("_.*","",colnames(expr))
# expr <- expr[, crop_barcodes]
expr <- t(expr)
colnames(expr) <- gsub("-","_",gsub("\\.","_",gsub("\\/\\/","__",colnames(expr))))
head(expr)

# GET GENE AND CELL METADATA =================

print("Getting gene and cell metadata")

## cell metadata
meta_cell <- x_seu@meta.data %>% 
  dplyr::select(uID=barcodes, umap1, umap2, cell_type=cell.type, sample_cluster_res_1=seurat_clusters) %>% 
  dplyr::filter(uID %in% rownames(expr))
rownames(meta_cell) <- meta_cell$uID

## gene metadata
## change dashes and dots for underscores in gene names, since they can cause issues in InSTAnT
meta_gene <- data.frame(gene=rownames(x_seu)) %>% dplyr::arrange(gene)
meta_gene$gene <- gsub("-","_",gsub("\\.","_",gsub("\\/\\/","__",meta_gene$gene)))
rownames(meta_gene) <- meta_gene$gene

# GET RAW TRANSCRIPTS ==========================

print("Getting raw transcripts")

## Import raw transcript
## Filter trans pqt by quality and cell_id present in expression data (if coords are cropped, only transcripts within cropped coords will be present)
# ## Add transcript id to rownames
## Format it as explained in cellsp tutorial
## Add X at the start of the gene name to match the expression data and gene metadata
## Filter out "*Control*" and "*Unassigned*" transcripts
## Change dashes and dots in gene names, since they can cause issues in InSTAnT
trans_pqt <- arrow::read_parquet(transcripts_parquet)
trans_pqt <- trans_pqt %>% dplyr::filter(qv >= 20, cell_id %in% rownames(expr))
# rownames(trans_pqt) <- trans_pqt$transcript_id
trans_pqt <- trans_pqt %>% dplyr::transmute(gene=paste0("X",feature_name), uID=cell_id, absX=x_location, absY=y_location, absZ=z_location)
trans_pqt <- trans_pqt %>% dplyr::filter(!stringr::str_detect(gene, "Control|Unassigned"))
trans_pqt$gene <- gsub("-","_",gsub("\\.","_",gsub("\\/\\/","__",trans_pqt$gene)))


# GET CELL BOUNDARIES ==========================

print("Getting cell boundaries")

## Imort cell boundaries
## Filter cells for those present in the seurat object
## Format it as explained in cellsp tutorial
cellb_pqt <- arrow::read_parquet(cell_boundaries_parquet)
cellb_pqt <- cellb_pqt %>% dplyr::filter(cell_id %in% rownames(expr)) 
cellb_pqt <- cellb_pqt %>% dplyr::transmute(uID=cell_id, vertex_x, vertex_y)
print(head(cellb_pqt))


# GET CELL COORDINATES ========================

print("Getting cell coordinates")

## Import cell coordinates
## Filter cells for those present in the seurat object
## Format it as explained in cellsp tutorial
## Add cell id as rownames
cell_coord_pqt <- arrow::read_parquet(cells_parquet)
cell_coord_pqt <- cell_coord_pqt %>% dplyr::filter(cell_id %in% rownames(expr)) 
cell_coord_pqt <- cell_coord_pqt %>% dplyr::transmute(uID=cell_id, x=x_centroid, y=y_centroid) #%>% tibble::column_to_rownames("uID")
rownames(cell_coord_pqt) <- cell_coord_pqt$uID

# WRITE DATA TO FILE ==========================

print("Writting data to file")

## create outdir
dir.create(outdir,F,T) # the outdir already contains the sample name here

## expression
out_expr <- paste0(outdir, "/", name, "_expr.csv"); print(out_expr)
data.table::fwrite(as.matrix(expr), out_expr, quote = F, sep = ",", row.names = T, col.names = T)

## gene metadata
out_meta_gene <- paste0(outdir,"/", name, "_gene_metadata.csv"); print(out_meta_gene)
data.table::fwrite(meta_gene, out_meta_gene, quote = F, sep = ",", row.names = T, col.names = T)

## cell metadata
out_meta_cell <- paste0(outdir,"/", name, "_cell_metadata.csv")
print(out_meta_cell)
data.table::fwrite(meta_cell, out_meta_cell, quote = F, sep = ",", row.names = T, col.names = T)

## raw transcripts
out_trans <- paste0(outdir, "/", name, "_transcripts.csv")
print(out_trans)
data.table::fwrite(trans_pqt, out_trans, quote = F, sep = ",", row.names = T, col.names = T)

## cell boundaries
## - set row.names=F because we can't set cell id as rownames because it would be duplicated
## - when imported in anndata, set 
out_cellb <- paste0(outdir, "/", name, "_cell_boundaries.csv")
print(out_cellb)
data.table::fwrite(cellb_pqt, out_cellb, quote = F, sep = ",", row.names = F, col.names = T)

## cell coordinates
out_cell_coords <- paste0(outdir, "/", name, "_cell_coords.csv")
print(out_cell_coords)
data.table::fwrite(cell_coord_pqt, out_cell_coords, quote = F, sep = ",", row.names = T, col.names = T)

