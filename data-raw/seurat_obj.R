seurat_obj <- readRDS("/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/mangiola.s/PostDoc/oligo_breast/expanded_analyses_with_control/cancer_only_analyses/lymphoid/cancer_lymphoid_cell_type_curated.rds")

seurat_obj = seurat_obj |> RunPCA() |> select(-contains("UMAP")) |> RunUMAP(dims=1:20)
seurat_obj = seurat_obj |> select(cell, file, 3, 8, 9,S.Score, G2M.Score , Phase , curated_cell_type , contains("UMAP"))
seurat_obj = seurat_obj %>% filter(.cell %in% (seurat_obj %>% sample_frac(0.5) %>%  pull(.cell) %>% c(seurat_obj %>% filter(grepl("Delta", curated_cell_type)) %>% pull(.cell)) %>% unique))
seurat_obj = seurat_obj[VariableFeatures(seurat_obj),]
seurat_obj[["SCT"]]@scale.data = seurat_obj[["SCT"]]@scale.data[c("CD3D", "TRDC", "TRGC1", "TRGC2", "CD8A", "CD8B"),]
seurat_obj[["SCT"]]@data = seurat_obj[["SCT"]]@data[c("CD3D", "TRDC", "TRGC1", "TRGC2", "CD8A", "CD8B"),]
DefaultAssay(seurat_obj) = "SCT"
seurat_obj[["integrated"]] = NULL
save(seurat_obj , file="data/seurat_obj.rda", compress = "xz")