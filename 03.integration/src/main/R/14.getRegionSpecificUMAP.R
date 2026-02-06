library(Seurat)
library(tidyverse)

# * meta
projd <- here::here()
workd <- file.path(projd, "03.integration")
outd <- file.path(workd, "out", "tfneu_vf_region_cca_k5")

seud <- file.path("/tscc/projects/ps-renlab2/zhw063",
  "99.MouseBrainPairedTag",
  "allen_brain_atlas/allen_ref_seurat_object",
  "20240625_clean_objects")

regions <- c("AMY", "CPU", "ERC", "HCa", "HCp", "HYP", "NAC",
  "PFC", "VTA_SnR")

# * main
umapList <- lapply(regions, \(r) {
  message("Loading Seurat from ", r)
  seu <- readRDS(file.path(seud, str_glue("{r}.clean.rna.rds")))
  umap_obj <- seu@reductions$umap
  mat <- umap_obj@cell.embeddings |> as.data.frame()
  mat$barcode <- rownames(mat)
  return(mat)
})
umap <- do.call(rbind, umapList)

data.table::fwrite(x = umap,
  file = file.path(outd, "pt.umap.per.region.csv"),
  sep = ",",
  col.names = TRUE,
  row.names = FALSE)

# * test
r <- "AMY"
seu <- readRDS(file.path(seud, str_glue("{r}.clean.rna.rds")))
umap_obj <- seu@reductions$umap
mat <- umap_obj@cell.embeddings |> as.data.frame()
mat$barcode <- rownames(mat)


