seurat_obj_UMAP3 = 
  seurat_obj_for_BioCAsia2021 %>% 
  RunUMAP(dims = 1:30, n.components = 3L, spread    = 0.5,min.dist  = 0.01, n.neighbors = 10L)

seurat_obj_UMAP3[["RNA"]] = NULL
seurat_obj_UMAP3[["SCT"]] = NULL
seurat_obj_UMAP3 = seurat_obj_UMAP3[1,] 

seurat_obj_UMAP3 %>% saveRDS("~/PostDoc/workshops/bioc2022_tidytranscriptomics/dev/seurat_obj_UMAP3.rds", compress = "xz")
save(seurat_obj_UMAP3, file="data/seurat_obj_UMAP3.rda", compress = "xz")


seurat_obj_UMAP3 |>
  
  plot_ly(
    x = ~`UMAP_1`,
    y = ~`UMAP_2`,
    z = ~`UMAP_3`,
    color = ~curated_cell_type,
    size=0.05
  )
