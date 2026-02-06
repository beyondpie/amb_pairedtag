library(tidyverse)
library(tmpRpkg)
library(ComplexHeatmap)
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
getscCorMat <- function(pvMat, method = "pearson") {
  scs <- rownames(pvMat)
  # 2xn matrix
  comscs <- combn(scs, 2)

  r <- matrix(
    data = 1.0, nrow = length(scs), ncol = length(scs),
    dimnames = list(
      scs, scs
    )
  )

  for (j in seq_len(ncol(comscs))) {
    sc1 <- comscs[1, j]
    sc2 <- comscs[2, j]
    cr <- cor(pvMat[sc1, ], pvMat[sc2, ], method= method)
    r[sc1, sc2] <- cr
    r[sc2, sc1] <- r[sc1, sc2]
  }
  return(r)
}

plotSimMat <- function(scSimK27me3, scSimK9me3, lowVal, highVal, name){
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
cntK27me3 <- pd$read_csv(file.path(genomeRanged, "aggK27me3.pd.csv"),
  sep = ",", header = 0L, index_col = 0L)
matK27me3 <- as.matrix(cntK27me3)[allenscs, ]
rpmK27me3 <- matK27me3 * 1e4 / rowSums(matK27me3)

sc2scCorK27me3 <- getscCorMat(rpmK27me3)

sc2scSPCorK27me3 <- getscCorMat(rpmK27me3, method = "spearman")

cntK9me3 <- pd$read_csv(file.path(genomeRanged, "aggK9me3.pd.csv"),
  sep = ",", header = 0L, index_col = 0L)
matK9me3 <- as.matrix(cntK9me3)[allenscs, ]
rpmK9me3 <- matK9me3 * 1e4 / rowSums(matK9me3)

sc2scCorK9me3 <- getscCorMat(rpmK9me3)
sc2scSPCorK9me3 <- getscCorMat(rpmK9me3, method = "spearman")


lowVal <- 0.0
highVal <- 1
(
  h <- plotSimMat(sc2scCorK27me3, sc2scCorK9me3,
    lowVal = lowVal, highVal = highVal, name = "50Kbp genome bin"
    )
)


withr::with_pdf(
 new = file.path(figd, "sc2scCor.K27me3-K9me3.details.pdf") ,
 code = {draw(h)},
 width = 30, height = 16
)



lowVal <- 0.0
highVal <- 1
(
  h2 <- plotSimMat(sc2scSPCorK27me3, sc2scSPCorK9me3,
    lowVal = lowVal, highVal = highVal, name = "50Kbp genome bin"
    )
)


withr::with_pdf(
 new = file.path(figd, "sc2scSPCor.K27me3-K9me3.details.pdf") ,
 code = {draw(h2)},
 width = 30, height = 16
)

# ** try class-level
cntK27me3 <- pd$read_csv(file.path(genomeRanged, "aggClassK27me3.pd.csv"),
  sep = ",", header = 0L, index_col = 0L)
matK27me3 <- as.matrix(cntK27me3)
rpmK27me3 <- matK27me3 * 1e4 / rowSums(matK27me3)
sc2scCorK27me3 <- getscCorMat(rpmK27me3)


cntK9me3 <- pd$read_csv(file.path(genomeRanged, "aggClassK9me3.pd.csv"),
  sep = ",", header = 0L, index_col = 0L)
matK9me3 <- as.matrix(cntK9me3)
rpmK9me3 <- matK9me3 * 1e4 / rowSums(matK9me3)

sc2scCorK9me3 <- getscCorMat(rpmK9me3)


lowVal <- 0.1
highVal <- 0.9
(
  h <- plotSimMat(sc2scCorK27me3, sc2scCorK9me3,
    lowVal = lowVal, highVal = highVal, name = "50Kbp genome bin"
    )
)


withr::with_pdf(
 new = file.path(figd, "class2classCor.K27me3-K9me3.details.pdf") ,
 code = {draw(h)},
 width = 15, height = 8
)

