library(tidyverse)
library(ComplexHeatmap)

projd <- here::here()
ptGeneCPMfnm <- file.path(
  projd, "data", "amb_PT_RNA_pseudobulk",
  "pt.RNAseq.pseudobulk.CPM.scbyg.rds")
ptscMetafnm <- file.path(projd, "meta",
  "PairedTagSubclassMetaFromCellMetaChromHMM.csv")
ptscMeta <- data.table::fread(file = ptscMetafnm,
  sep = ",", data.table = F, header = T) |>
  x => `rownames<-`(x, x$PairedTagName) |>
  x => x[x$ATAC > 0, ]

ptscs <- ptscMeta$PairedTagName

m <- readRDS(file.path(projd, "data",
  "amb_PT_RNA_pseudobulk", "panelE.sc2motif.enrich.rds"))
filterNULLfromList <- function(l) {
  wisnull <- vapply(l, is.null, T)
  if (sum(wisnull) == length(l)) {
    message("All are NULLs, return NULL.")
    NULL
  } else {
    l[!wisnull]
  }
}


getTopKGeneExpPerSubclass <- function(geneMat, ntop) {
  q <- 1 - ntop / ncol(geneMat)
  vapply(rownames(geneMat), \(sc) {
    quantile(x = geneMat[sc, ], probs = q)
  }, 0.0) |>
    setNames(object = _, nm = rownames(geneMat))
}

filterByGeneExtPerSubclass <- function(motifMat, geneMat, ct) {
  r <- motifMat
  tfnms <- colnames(motifMat)
  for(sc in rownames(motifMat)) {
    index <- geneMat[sc, tfnms] < ct[sc]
    r[sc, index] <- 0.0
  }
  return(r)
}

getDiffTFPerSubclass <- function(diffTFmat, ptscs,
                                 padj = 0.05,
                                 log2FC = 0.5) {
  lapply(ptscs, \(sc) {
    subset(
      diffTFmat, (subclass == sc) & (p_val_adj <= padj) & (avg_log2FC >= log2FC))
  }) |>
    setNames(object = _, nm = ptscs) |>
    filterNULLfromList(l = _)
}


geneLogCPM <- readRDS(ptGeneCPMfnm) |>
  log1p() |>
  x => x[ptscs, ]
ctoff <- getTopKGeneExpPerSubclass(geneLogCPM, ntop = 10000)

tfmat <- geneLogCPM[rownames(m), colnames(m)]
# FIXME: here why we have consistent values for Tcf4
# tfmat[tfmat > quantile(tfmat, 0.9)] <- quantile(tfmat, 0.9)

tfmat <- scale(tfmat)


quantile(tfmat, na.rm = T)
low.val.col <- -0.5
high.val.col <- 1.5

legend_labels <- c(round(low.val.col, 1),
  round(high.val.col, 1))

col_fun <- circlize::colorRamp2(
  seq(low.val.col, high.val.col, length = 60),
  viridis::viridis(60)
)

hmTFExp <- ComplexHeatmap::Heatmap(
  matrix = tfmat,
  col = col_fun,
  cluster_columns = F,
  cluster_rows = F,
  show_row_names = T,
  row_names_gp = grid::gpar(fontsize = 7),
  column_names_gp = grid::gpar(fontsize = 7),
  show_column_names = T,
  top_annotation = NULL,
  left_annotation = NULL,
  use_raster = T,
  show_heatmap_legend = T,
  na_col = "#440154FF",
  heatmap_legend_param = list(
    title = "scaled logCPM",
    at = c(low.val.col, high.val.col),
    labels = legend_labels,
    direction = "horizontal"
  )
)

hmTFExp
saveRDS(object = tfmat, file = file.path(
  projd, "99.figures", "out", "fig4", "panelE.sc2motif.geneExp.rds"))
