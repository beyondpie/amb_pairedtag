library(Seurat)
library(tidyverse)

# * meta
projd <- here::here()
workd <- file.path(projd, "03.integration")
mrs <- c("AMY", "CPU", "HYP", "HIP", "ERC", "NAC", "VTA", "PFC")
tfneu_region_dir <- file.path(
  workd, "out",
  str_glue("tfneu_vf_region_cca_k5")
)

# * function
# L5 per region mapping to super type
load_sumtf_seu <- function(r) {
  readRDS(file.path(tfneu_region_dir,
    str_glue("sumtfneu_vf_{r}_cca_k5.seu.rds")))
}

load_anchor_seu <- function(r) {
  file.path(tfneu_region_dir, paste0("tf_", r),
    "tf.anchors.with-cca-kac5.rds") |>
    readRDS(file = _)
}

load_query_seu <- function(r) {
  file.path(tfneu_region_dir, paste0("tf_", r),
    "query.with.tf-cca-kac5_on-cl.rds") |>
    readRDS(file = _)
}

# * ggplot
mytheme <- theme(
  panel.border = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  plot.title = element_text(colour = "black", hjust = 0.5,
    size = 15),
  axis.text = element_blank(),
  axis.ticks = element_blank(),
  axis.title = element_text(colour = "black", size = 12),
  axis.line = element_line(colour = "black"),
  legend.position = "right",
  legend.title = element_text(colour = "black", size = 13),
  legend.text = element_text(colour = "black", size = 12)
)


# * main
for (r in mrs) {
  querySeu <- load_query_seu(r)
  tfscoremat <- querySeu@assays$prediction.score.id@data
  tfscoremax <- querySeu$predicted.id.score

  tfSeu <- load_sumtf_seu(r)
  tfSeu$tfscore <- 0.0
  tfSeu$tfscore[!tfSeu$isRef] <- tfscoremax[
    with(tfSeu@meta.data, barcode[!isRef])]

  p_tech <- DimPlot(tfSeu, reduction = "UMAP", raster = TRUE,
    group.by = "tech")

  p_tfscore <- FeaturePlot(tfSeu, reduction = "UMAP", raster = TRUE,
    features = "tfscore", cells = with(tfSeu@meta.data, barcode[!isRef]),
    cols = c("white", "red"))

  p_umap <- p_tech + p_tfscore
  ggsave(
    filename = file.path(tfneu_region_dir,
      str_glue("sumtfscore_vf_{r}_cca_k5.pdf")),
    plot = p_umap, width = 12, height = 6)
}
