library(Seurat)
library(tidyverse)
library(reticulate)

# set python env
use_condaenv(condaenv = "sa2", 
             conda = "/home/szu/miniforge3/bin/conda")
sc <- import("scanpy")
ad <- import("anndata")
pd <- import("pandas")
np <- import("numpy")

# set env
projd <- "/mnt/tscc2/szu/projects/amb_pairedtag"
setwd(file.path(projd, "03.integration", "src/main/python"))

# meta and functions 
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

# main
adata <- ad$read_h5ad("test_integration/adata/adata_merge_al_pt_at.h5ad")

barcodes <- adata$obs_names$to_list()
umap <- adata$obsm['X_umap']
rownames(umap) <- barcodes
colnames(umap) <- c("UMAP1", "UMAP2")
umap <- as.data.frame(umap)

# add batch info
batch <- adata$obs['batch']$batch
umap$batch <- batch
batch_colors <- c(
  "allen" = "red",
  "atac" = "blue",
  "pt" = "darkgreen"
)

# plot

umap_batch <- ggplot(data = umap, aes(x = UMAP1, y = UMAP2, color = batch)) +
  geom_point(size = 1, shape = 19, alpha = 0.5) +
  guides(color = guide_legend(override.aes = list(size=8))) +
  mytheme +
  scale_color_brewer(palette = "Paired")
