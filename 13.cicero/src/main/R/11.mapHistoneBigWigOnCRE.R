# Get the RPKM of different histones on CREs
# We extend each CRE up and down stream 500 bps

# We also provide eRNA mapping here, where
# we use raw GRE without any extension.

suppressMessages({
  suppressWarnings({
    library(tidyverse)
    library(tmpRpkg)
    # detach("package:tmpRpkg", unload = T)
    library(GenomicRanges)
    library(S4Vectors)
    library(rtracklayer)
    library(ComplexHeatmap)
    library(reticulate)
    library(future)
    library(future.apply)
  })
})

options(future.globals.maxSize = 50 * 1024^3)
future::plan(multicore, workers = 30)
# * meta
projd <- here::here()
workd <- file.path(projd, "13.cicero")
outd <- file.path(workd, "out")
ptRNAbwd <- file.path(projd, "data", "ptRNAbam", "bigwig")
smartRNAbwd <- file.path(projd, "data", "mba_smartseq_bigwig")
chromAnnotd <- tmpRpkg:::CREAnnotd

ptscMeta <- tmpRpkg::loadptscMeta() |>
  x => x[x$ATAC > 0, ] |>
  x => `rownames<-`(x, x$PairedTagName)
scs <- ptscMeta$PairedTagName

# * functions
getlogRPKM <- function(sc, h, CREGR) {
  bwGR <- tmpRpkg::loadHistoneBigWig(sc, h)
  message(str_glue("get {sc}'s {h} log RPKM."))
  s <- getbwSignal(bwGR = bwGR, peakGR = CREGR)
  data.frame(
    CRE = names(s),
    RPKM = s
  ) |>
    mutate(logRPKM = log1p(RPKM)) |>
    x => `rownames<-`(x, x$CRE)
}

getPairedTageRNAlogRPKM <- function(sc, CREGR) {
  bwGR <- tmpRpkg::loadBigWigGR(
    file.path(ptRNAbwd, str_glue("{sc}.RPKM.bw")))
  message(str_glue("get {sc}'s eRNA RPKM."))
  pteRNAall <- getbwSignal(bwGR = bwGR, peakGR = CREGR)
  data.frame(
    CRE = names(pteRNAall),
    RPKM = pteRNAall
  ) |>
    mutate(logRPKM = log1p(RPKM))
}

getlogRPKMat <- function(sc2signalList, CREs, scs) {
  m <- matrix(
    data = NA, nrow = length(CREs), ncol = length(scs),
    dimnames = list(CREs, scs)
  )
  for (sc in scs) {
    m[, sc] <- sc2signalList[[sc]][CREs, "logRPKM"]
  }
  return(m)
}

# * main
# all CREs will be considered
allCRE <- tmpRpkg::loadAllCREAsGR()
allExtCRE <- tmpRpkg::extendGR(allCRE, up = 500, down = 500)

# 1. histones
histones <- c("H3K27ac", "H3K4me1", "H3K27me3", "H3K9me3")
for (h in histones) {
  y <- future_lapply(scs, \(sc) {
    getlogRPKM(sc, h, allExtCRE)
  }) |>
    setNames(object = _, nm = scs)

  x <- getlogRPKMat(
    y,
    CREs = mcols(allExtCRE)$name, scs
  )

  saveRDS(
    object = x,
    file = file.path(outd, str_glue("{h}.logRPKM.allCRE.pbysc.rds"))
  )
}

# 2. eRNA
sc2pteRNA <- future_lapply(scs, \(sc) {
  getPairedTageRNAlogRPKM(sc, allCRE)
}) |> setNames(object = _, nm = scs)
x <- getlogRPKMat(sc2pteRNA, CREs = mcols(allCRE)$name, scs)
saveRDS(
  object = x,
  file = file.path(outd, "eRNA.PairedTag.logRPKM.allCRE.pbysc.rds")
)
