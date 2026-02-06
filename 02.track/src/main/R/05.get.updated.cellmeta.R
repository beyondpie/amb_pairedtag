# Updated meta: 20240317
# We added back the some previously filtered cells by lowing our QC thresholds

library(Seurat)

mtx_dir <- file.path(
  "/projects/ps-renlab2/",
  "zhw063", "99.MouseBrainPairedTag",
  "2024Mar_new_mtxs"
)

# * load RNA Seurat
allRNA <- readRDS(
  file.path(mtx_dir, "brain.all.rna.rds"))

cellMeta <- allRNA@meta.data
cellMeta$orig.ident <- colnames(allRNA)
cellMeta <- dplyr::rename(cellMeta, barcode = orig.ident)

data.table::fwrite(x = cellMeta,
  file = file.path(here::here(), "02.track",
    "src/main/resource", "cellMeta.20240321.csv"),
  sep = ",", col.names = T, row.names = F)
