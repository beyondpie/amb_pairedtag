library(tidyverse)
library(Seurat)

# meta data
projd <- here::here()
tfd <- file.path(projd, "03.integration", "out")
ptseud <- file.path(projd, "data", "pairedtag_seurat")
commonGenes <- data.table::fread(
  file = file.path(
    projd, "meta",
    "common_gene_markers.csv"
  ),
  sep = ",", header = TRUE, data.table = FALSE
)

## ptcellMeta1 <- file.path(
##   projd, "meta",
##   "pairedtag.meta.csv"
## ) |>
##   data.table::fread(file = _, sep = ",", header = TRUE,
##     data.table = FALSE)
## ptcellMeta2 <- file.path(
##   projd, "meta",
##   "pt.barcode.meta.L5.csv"
## ) |>
##   data.table::fread(file = _, sep = ",", header = TRUE,
##     data.table = FALSE)
## ptcellMeta <- merge(x = ptcellMeta1, y = ptcellMeta2,
##   by = "barcode", all.x = FALSE, all.y = TRUE)
## data.table::fwrite(
##   x = ptcellMeta,
##   file = file.path(projd, "meta", "pairedtag.cell.meta.all.csv"),
##   sep = ",", col.names = TRUE, row.names = FALSE
## )
ptcellMeta <- data.table::fread(
  file = file.path(projd, "meta", "pairedtag.cell.meta.all.csv"),
  sep = ",", header = TRUE, data.table = FALSE
)

# load neuron and nn seurat object
nn_ptseu <- readRDS(file.path(ptseud, "ptRNA.nn.k8.L5ds30.rds"))
neu_ptseu <- readRDS(file.path(ptseud, "ptRNA.neu.k8.L5ds50.rds"))
# merge two seurats
ptseu <- merge(nn_ptseu, neu_ptseu,
  add.cell.ids = NULL)

# load transfer label info
nn_tfsum <- file.path(tfd, "tfnn_k8_cca_k5",
  "sum_tfnn_k8_cca_k5.seurat.rds") |>
  readRDS(file = _)
neu_tfsum <- file.path(
  tfd, "tfneu_k8_cca_k5",
  "sum_tfneu_k8_cca_k5.seurat.rds"
) |>
  readRDS(file = _)

# get L5 to supertype info
tfMeta <- rbind(nn_tfsum@meta.data, neu_tfsum@meta.data)
L5_sp <- tfMeta[!tfMeta$isRef, ] |>
  group_by(L5) |>
  summarise(sp = unique(supertype_id_label)) |>
  as.data.frame() |>
  x => `rownames<-`(x, paste0("L5-", x$L5))
ptcellMeta$sp <- L5_sp[paste0("L5-", ptcellMeta$L1_2_3_4_5), "sp"]
rownames(ptcellMeta) <- ptcellMeta$barcode
ptseu$sp <- ptcellMeta[colnames(ptseu), "sp"]

# perform norm and scale exp in pseudo-bulk level
# (supertype level)
pseudo_ptseu <- AggregateExpression(
  object = ptseu,
  group.by = "sp",
  return.seurat = TRUE
)

saveRDS(pseudo_ptseu,
  file.path(tfd, "pt.pseudo_seurat.supertype.rds"))

# output
