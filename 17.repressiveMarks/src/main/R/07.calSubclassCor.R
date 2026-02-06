library(tidyverse)
library(tmpRpkg)
library(ComplexHeatmap)

# * meta
projd <- here::here()
workd <- file.path(projd, "17.repressiveMarks")
outd <- file.path(workd, "out")
figd <- file.path(workd, "figure")
upeakd <- file.path(projd, "data", "pairedtag_peak")
sc2upWigValued <- file.path(outd, "allGenomicRange")
mods <- c("H3K27me3", "H3K9me3")
ntop <- 20000

# load pairedtag subclass meta
ptscMeta <- tmpRpkg::loadptscMeta() |>
  x => x[x$ATAC > 0, ] |>
  x => `rownames<-`(x, x$pairedTagName)
ptscs <- ptscMeta$PairedTagName

# * functions
loaduPeak <- function(mod) {
  fnm <- file.path(upeakd, str_glue("{mod}.merged.all.blv2.me.peak"))
  r <- data.table::fread(file = fnm, header = T, sep = "\t", data.table = F)
  pnm <- with(r, paste(chrom, paste(startFrom, endTo, sep = "-"), sep = ":"))
  rownames(r) <- pnm
  return(r)
}

loadscPeakValue <- function(sc, mod) {
  fnm <- file.path(sc2upWigValued,
    str_glue("subclass2unifiedPeakWigValue_{mod}"),
    str_glue("{sc}.{mod}.broadPeak.wigValue.bed"))
  data.table::fread(file = fnm, header = F, sep = "\t",
    data.table = F) |>
    setNames(object = _, nm = c("chrom", "startFrom", "endTo", "score")) |>
    x => `rownames<-`(x, paste(x$chrom,
      paste(x$startFrom, x$endTo, sep = "-"), sep = ":"))
}

getsc2PeakValueMat <- function(scPeakValueList, uPeak) {
  m <- matrix(data = 0.0,
    nrow = length(scPeakValueList), ncol = length(uPeak),
    dimnames = list(
      names(scPeakValueList),
      uPeak
  ))
  for (i in seq_along(scPeakValueList)) {
    sc <- names(scPeakValueList)[i]
    r <- scPeakValueList[[i]]
    m[sc, rownames(r)] <- r$score
  }
  return(m)
}

getscCorMat <- function(pvMat) {
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
    cr <- cor(pvMat[sc1, ], pvMat[sc2, ])
    r[sc1, sc2] <- cr
    r[sc2, sc1] <- r[sc1, sc2]
  }
  return(r)
}


plotscCor <- function(scCorMat, ql = 0.1, qh = 0.7,
                      t = "dist", show_row_names = F, show_column_names = F) {
  lowVal <- quantile(scCorMat, ql)
  highVal <- quantile(scCorMat, qh)
  ncolor <- 60
  legendLab <- c(
    round(lowVal, 1),
    round(highVal, 1)
  )
  ComplexHeatmap::Heatmap(
    mat = scCorMat,
    col = circlize::colorRamp2(seq(lowVal, highVal, length = ncolor),
      viridis::viridis(ncolor)),
    show_row_names = show_row_names,
    show_column_names = show_column_names,
    cluster_columns = F,
    cluster_rows = F,
    use_raster = T,
    row_names_gp = grid::gpar(fontsize = 6),
    column_names_gp = grid::gpar(fontsize = 6),
    heatmap_legend_param = list(
      title = t,
      at = c(lowVal, highVal),
      labels = legendLab,
      direction = "vertical"
    )
  )
}

scaleLogPeakValueMat <- function(m) {
  logm <- log1p(m)
  zeroColIndex <- colSums(m) == 0
  m1 <- m[, !zeroColIndex]
  zeroVarColIndex <- vapply(seq_len(ncol(m1)), \(j) {
    v <- var(m1[, j])
  }, 0.0) |>
    x => x == 0
  m1[, !zeroVarColIndex] |>
    scale(x = _, center = T, scale = T)
}

logPeakValueMat <- function(m) {
  logm <- log1p(m)
  zeroColIndex <- colSums(m) == 0
  m1 <- m[, !zeroColIndex]
  zeroVarColIndex <- vapply(seq_len(ncol(m1)), \(j) {
    v <- var(m1[, j])
  }, 0.0) |>
    x => x == 0
  m1[, !zeroVarColIndex]
}

# * main

# ** H3K27me3
mod <- "H3K27me3"

# load unified peak
upeakH3K27me3 <- loaduPeak(mod)
scPeakValues <- lapply(ptscs, \(sc) {
  loadscPeakValue(sc, mod)
}) |>
  setNames(object = _, nm = ptscs)
mat <- getsc2PeakValueMat(scPeakValues, rownames(upeakH3K27me3))
slmat <- logPeakValueMat(mat)

distMat <- dist(slmat, "euclidean") |>
  as.matrix()
# ndistmat <- (max(distMat) - distMat) / max(distMat)
hpndist <- plotscCor(-distMat / length(ptscs), qh = 0.95, ql = 0.3,
  show_row_names = T, show_column_names = F
  )
hpndist



# ** H3K9me3
mod <- "H3K9me3"
upeak2 <- loaduPeak(mod)
scPeakValues2 <- lapply(ptscs, \(sc) {
  loadscPeakValue(sc, mod)
}) |>
  setNames(object = _, nm = ptscs)
mat2 <- getsc2PeakValueMat(scPeakValues2, rownames(upeak2))
slmat2 <- logPeakValueMat(mat2)
distMat2 <- dist(slmat2, "euclidean") |>
  as.matrix()
hpndist2 <- plotscCor(-distMat2 / length(ptscs), qh = 0.95, ql = 0.3,
  show_row_names = T, show_column_names = T
  )
hpndist2

hpndist + hpndist2
