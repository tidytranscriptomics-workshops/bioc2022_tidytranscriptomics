
library(tidyverse)
library(glue)
library(Seurat)
library(tidyseurat)
library(tidySingleCellExperiment)

seurat_obj <- readRDS("/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/mangiola.s/PostDoc/oligo_breast/expanded_analyses_with_control/cancer_only_analyses/lymphoid/cancer_lymphoid_cell_type_curated.rds")

set.seed(123)
seurat_obj = seurat_obj |> RunPCA(npcs = 20) |> select(-contains("UMAP")) |> RunUMAP(dims=1:20)
seurat_obj = seurat_obj |> select(.cell, file, 3, 8, 9, S.Score, G2M.Score , Phase , curated_cell_type , contains("UMAP"))

seurat_obj = seurat_obj %>% filter(.cell %in% (seurat_obj %>% sample_n(3000) %>%  pull(.cell) %>% c(seurat_obj %>% filter(grepl("Delta", curated_cell_type)) %>% pull(.cell)) %>% unique))

seurat_obj = seurat_obj %>% FindVariableFeatures(assay="RNA", nfeatures = 500)
seurat_obj = seurat_obj[VariableFeatures(seurat_obj, assay="RNA") %>% c("CD3D", "TRDC", "TRGC1", "TRGC2", "CD8A", "CD8B"),]
seurat_obj = seurat_obj %>% mutate(
	file = factor(file), Barcode = factor(Barcode), batch= factor(batch), BCB= factor(BCB), Phase= factor(Phase), curated_cell_type= factor(curated_cell_type),	
	nCount_RNA = as.integer(nCount_RNA), nFeature_RNA= as.integer(nFeature_RNA), nCount_SCT= as.integer(nCount_SCT), nFeature_SCT= as.integer(nFeature_SCT)
) 
#seurat_obj[["SCT"]]@scale.data = seurat_obj[["SCT"]]@scale.data[c("CD3D", "TRDC", "TRGC1", "TRGC2", "CD8A", "CD8B"),]
#seurat_obj[["SCT"]]@data = seurat_obj[["SCT"]]@data[c("CD3D", "TRDC", "TRGC1", "TRGC2", "CD8A", "CD8B"),]
DefaultAssay(seurat_obj) = "SCT"
seurat_obj[["integrated"]] = NULL

sce_obj = seurat_obj %>%
	as.SingleCellExperiment() |>
	
	# Add factor of interest
	nest(data = -file) |> 
	mutate(condition = sample(c("treated", "untreated"), n(), replace = TRUE)) |>
	unnest(data)

# Parse
sce_obj = 
	sce_obj |> 
	select(-condition) |>  
	left_join(readRDS("~/metadata_oligo.rds")) |> 
	rename(treatment = type) |> 
	mutate(treatment = if_else(treatment=="OMBC", "treated", "untreated")) |> 
	#select(-file) |> 
	mutate(sample = glue("S{as.integer(as.factor(sample))}")) |> 
	rename(cell_type = curated_cell_type) |> 

	# filtering because of to few samples per cell types
	filter(cell_type !="CD8+_Tem") |> 
	
	# Replace file path
	mutate(file = file |> str_replace("bhupinder_10X_260819", "single_cell"))

# job::job({
	save(sce_obj , file="data/sce_obj.rda", compress = "xz")
# })
