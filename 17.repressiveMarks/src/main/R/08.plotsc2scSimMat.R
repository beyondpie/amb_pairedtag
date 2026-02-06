library(tidyverse)
library(tmpRpkg)
library(ComplexHeatmap)

# * meta
projd <- here::here()
workd <- file.path(projd, "17.repressiveMarks")
outd <- file.path(workd, "out")
figd <- file.path(workd, "figure")
upeakd <- file.path(projd, "data", "pairedtag_peak")
scSimd <- file.path(outd, "allGenomicRange")
mods <- c("H3K27me3", "H3K9me3")
ntop <- 20000
simMethod <- c("RawCor", "ReluCor", "MinMaxCor",
  "MinMaxEucScaleLog", "MinMaxRawEuc",
  "RawSpearmanCor", "Jaccard")

# * functions
loadSimMat <- function(mod, simMethod) {
  fnm <- if(simMethod == "Jaccard") {
    file.path(scSimd,
      str_glue("sc2scSimMat.{mod}.{simMethod}.csv"))
  } else {
    file.path(scSimd,
      str_glue("sc2scSimMat.{mod}.{simMethod}.{ntop}.csv"))
  }
  r <- data.table::fread(file = fnm, sep = ",", header = T, data.table = F)
  m <- as.matrix(r)
  rownames(m) <- colnames(r)
  colnames(m) <- rownames(m)
  return(m)
}

plotSimMat <- function(scSimK27me3, scSimK9me3, lowVal, highVal, name){
  ncolor <- 30
  col <- circlize::colorRamp2(seq(lowVal, highVal, length = ncolor),
    viridis::viridis(ncolor))
  legendLab <- c(
    round(lowVal, 1),
    round(highVal, 1)
  )
  hmK27me3 <- ComplexHeatmap::Heatmap(
    mat = scSimK27me3,
    col = col,
    show_row_names = F,
    show_column_names = T,
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
      direction = "vertical"
    )
  )
  hmK9me3 <- ComplexHeatmap::Heatmap(
    mat = scSimK9me3,
    col = col,
    show_row_names = F,
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
      direction = "vertical"
    )
  )
  hmK27me3 + hmK9me3
}

# * main
for(s in simMethod) {
  simMatK27me3 <- loadSimMat("H3K27me3", s)
  simMatK9me3 <- loadSimMat("H3K9me3", s)
  lowVal <- 0.5
  highVal <- 0.9
  (
    h <- plotSimMat(simMatK27me3, simMatK9me3,
      lowVal = lowVal, highVal = highVal, name = s
      )
  )
}

withr::with_pdf(new = file.path(figd, "sc2scCor.H3K27me3.H3K9me3.pdf"),
  code = {draw(h)}, width = 15, height = 7
  )
