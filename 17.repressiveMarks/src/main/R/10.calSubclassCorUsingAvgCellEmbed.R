library(tidyverse)
library(tmpRpkg)
library(reticulate)
use_python("/Users/szu/ann/bin/python")
pd <- import("pandas", convert = T)


# * meta
projd <- here::here()
workd <- file.path(projd, "17.repressiveMarks")
outd <- file.path(workd, "out")
figd <- file.path(workd, "figure")
genomeRanged <- file.path(outd, "allGenomicRange")

ptscMeta <- tmpRpkg::loadptscMeta() |>
  x => x[x$ATAC > 0, ] |>
  x => `rownames<-`(x, x$pairedTagName)
allenscs <- ptscMeta$AllenIdName


# * functions
loadCellEmbed <- function(mod) {
  pd$read_csv(file.path(genomeRanged,
    str_glue("pt.singlecell.emebd.{mod}.csv")), sep = ",",
  header = 0L, index_col = 0L)
}

loadGroupSimMat <- function(mod, groupBy = "sc", usingEmbed = "SpectralEmbed") {
  fnm <- file.path(genomeRanged,
    str_glue("pt.{groupBy}.{mod}.cor.{usingEmbed}.csv"))
  r <- data.table::fread(file = fnm, sep = ",", header = T, data.table = F) |>
    x => `rownames<-`(x, colnames(x))
  # order the subclasses
  ids <- rownames(r) |>
    strsplit(x = _, split = " ") |>
    x => vapply(x, \(y) {as.integer(y[1])}, 1L)
  r[order(ids, decreasing = F), order(ids, decreasing = F)]
}

plotSimMat <- function(scSimK27me3, scSimK9me3, lowVal, highVal, name) {
  ncolor <- 30
  col1 <- circlize::colorRamp2(seq(lowVal, highVal, length = ncolor),
    viridis::viridis(ncolor))
  legendLab <- c(
    round(lowVal, 1),
    round(highVal, 1)
  )
  col2 <- circlize::colorRamp2(seq(lowVal, highVal, length = ncolor),
    viridis::magma(ncolor))


  hmK27me3 <- ComplexHeatmap::Heatmap(
    mat = scSimK27me3,
    col = col1,
    show_row_names = F,
    show_column_names = F,
    cluster_columns = F,
    cluster_rows = F,
    show_row_dend = F,
    show_column_dend = F,
    use_raster = T,
    row_names_gp = grid::gpar(fontsize = 6),
    column_names_gp = grid::gpar(fontsize = 6),
    heatmap_legend_param = list(
      title = str_glue("H3K27me3 {name}"),
      at = c(lowVal, highVal),
      labels = legendLab,
      direction = "horizontal"
    ),
    column_title = "H3K27me3"
  )
  hmK9me3 <- ComplexHeatmap::Heatmap(
    mat = scSimK9me3,
    col = col2,
    show_row_names = T,
    show_column_names = T,
    cluster_columns = F,
    cluster_rows = F,
    show_row_dend = F,
    show_column_dend = F,
    use_raster = T,
    row_names_gp = grid::gpar(fontsize = 6),
    column_names_gp = grid::gpar(fontsize = 6),
    heatmap_legend_param = list(
      title = str_glue("H3K9me3 {name}"),
      at = c(lowVal, highVal),
      labels = legendLab,
      direction = "horizontal"
    ),
    column_title = "H3K9me3"
  )
  hmK27me3 + hmK9me3
}


# * main
# 1. sc level, embed
embed <- "SpectralEmbed"
scSimK27me3 <- loadGroupSimMat(mod = "H3K27me3",
  groupBy = "sc", usingEmbed = embed) |>
  x => x[allenscs, allenscs]
scSimK9me3 <- loadGroupSimMat(mod = "H3K9me3",
  groupBy = "sc", usingEmbed = embed) |>
  x => x[allenscs, allenscs]
(
  hm <- plotSimMat(scSimK27me3, scSimK9me3,
    lowVal = 0.6, highVal = 1.0,
    name = str_glue("Subclass Pearson Cor Using {embed}"))
)

# 2. class level, embed
embed <- "SpectralEmbed"
scSimK27me3 <- loadGroupSimMat(mod = "H3K27me3",
  groupBy = "class", usingEmbed = embed)
scSimK9me3 <- loadGroupSimMat(mod = "H3K9me3",
  groupBy = "class", usingEmbed = embed)
(
  hm <- plotSimMat(scSimK27me3, scSimK9me3,
    lowVal = 0.6, highVal = 1.0,
    name = str_glue("class Pearson Cor Using {embed}"))
)
