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

tmpRpkg::setCondaEnv(
  envnm = "sa2stable",
  conda = "/home/szu/mambaforge/bin/conda"
)
pd <- reticulate::import(module = "pandas", convert = F)
options(future.globals.maxSize = 50 * 1024^3)
future::plan(multicore, workers = 30)

# * meta
projd <- here::here()
workd <- file.path(projd, "13.cicero")
figd <- file.path(workd, "figure", "eRNA")
outd <- file.path(workd, "out", "eRNA")
ptRNAbwd <- file.path(projd, "data", "ptRNAbam", "bigwig")
chromAnnotd <- tmpRpkg:::CREAnnotd

ptscMeta <- tmpRpkg::loadptscMeta() |>
  x => x[x$ATAC > 0, ] |>
  x => `rownames<-`(x, x$PairedTagName)
scs <- ptscMeta$PairedTagName

atacPeakd <- file.path(projd, "data", "chromHMM", "subclass_peak")
stateNameOrd <- tmpRpkg:::stateNameOrd
mm10promoter <- tmpRpkg::loadPromoter()
stateName2Color <- tmpRpkg::getChromStateColor()

# * functions
geteRNAlogRPKM <- function(sc, CREGR) {
  bwGR <- tmpRpkg::loadBigWigGR(
    file.path(ptRNAbwd, str_glue("{sc}.RPKM.bw"))
  )
  message(str_glue("get {sc}'s eRNA RPKM."))
  pteRNAall <- getbwSignal(bwGR = bwGR, peakGR = CREGR)
  data.frame(
    CRE = names(pteRNAall),
    pteRNA = pteRNAall
  ) |>
    mutate(logptRPKM = log1p(pteRNA))
}

geteRNAMat <- function(sceRNAs, CREs, scs) {
  m <- matrix(
    data = NA, nrow = length(CREs), ncol = length(scs),
    dimnames = list(CREs, scs)
  )
  for (sc in scs) {
    m[, sc] <- sceRNAs[[sc]][CREs, "logptRPKM"]
  }
  return(m)
}

reduceCREbyeRNA <- function(sc, CRElist, fdPbyS, topk = 2000, eRNAfd = 1.1) {
  message(str_glue("reduce CREs by eRNA: {sc}."))
  m <- CRElist[[sc]]
  f <- fdPbyS[m, sc]
  index <- (f >= eRNAfd)
  if (sum(index) < 1) {
    message(str_glue("{sc} has no CREs under eRNA fd {eRNAfd}."))
    return(NULL)
  }
  mm <- m[index]
  return(mm[seq_len(min(topk, sum(index)))])
}

# * main
# * ATAC ComplexHeatmap
state <- "Chr-A"
# logCPM
xPbyS <- loadATACPseudoBulkMat(tmpRpkg:::snATACPMpbyc, pd = pd) |>
  log1p(x = _)
# pass q-value 0.05
allEnhDARCREs <- readRDS(file.path(
  projd,
  "13.cicero", "out", str_glue("distal.{state}.DE.ATACPeak.rds")
))

## 0. get order of subclass based on their ids.
ordscs <- names(allEnhDARCREs) # already ordered

## 1. get subclass-specific enhancers
# further filtered by log2fd
enhDARCREs <- lapply(ordscs, tmpRpkg::retrieveCRE,
  deATAC = allEnhDARCREs,
  log2fd = 0.2, n = 999999
) |>
  setNames(object = _, nm = ordscs)
DAREnh <- unlist(enhDARCREs) |>
  x => x[!duplicated(x)]

## get eRNA logRPKM signals
DAREnhGR <- transformRegion2GR(DAREnh)

# The matrix loading can be optimized using the latest generated eRNA mat.
sceRNAlogRPKMs <- future_lapply(ordscs, geteRNAlogRPKM, DAREnhGR) |>
  setNames(object = _, nm = ordscs)
mPbyS <- geteRNAMat(sceRNAlogRPKMs, DAREnh, ordscs)
saveRDS(mPbyS, file = file.path(
  outd, "eRNA.logRPKM.filteredDAREhn.pbysc.rds"
))

## sdOfm <- vapply(seq_len(nrow(mPbyS)), \(i) sd(mPbyS[i, ], na.rm = T), 0.1) |>
##   setNames(object = _, nm = allCREGRs)
## save to txt file
data.table::fwrite(
  x = mPbyS,
  file = file.path(outd, "mat", "mat.csv"),
  sep = ",", row.names = F, col.names = F
)
write.table(
  x = rownames(mPbyS),
  file = file.path(outd, "mat", "rownames.txt"), quote = F, row.names = F,
  col.names = F
)
write.table(
  x = colnames(mPbyS),
  file = file.path(outd, "mat", "colnames.txt"), quote = F, row.names = F,
  col.names = F
)

## load fold change from scala script
fdPbyS <- data.table::fread(
  file = file.path(outd, "mat", "fd.csv"),
  header = F, sep = ",", data.table = F
) |>
  setNames(object = _, nm = colnames(mPbyS)) |>
  x => `rownames<-`(x, rownames(mPbyS))

## now get CREs that differential and having eRNA signals
getOrdCRE <- function(topk = 1000, eRNAfd = 1.1) {
  partCRElist <- lapply(ordscs, \(sc) {
    reduceCREbyeRNA(
      sc = sc, enhDARCREs, fdPbyS,
      topk = topk, eRNAfd = eRNAfd
    )
  })
  names(partCRElist) <- ordscs
  partCRElist |>
    unlist() |>
    x => x[!duplicated(x)]
}
partCREs <- getOrdCRE(topk = 1000, eRNAfd = 1.1)
# 044_OB_Dopa_Gaba has no CREs under eRNA fd 1.1.

ordscs <- ordscs[!ordscs %in% c("044_OB_Dopa_Gaba")]
## 2. map pt subclass to atac-seq subclass
allenClMeta <- loadAllenClMeta()
scsATAC <- vapply(allenClMeta$subclass_label, getATACSubclassName, "a")
scspt <- vapply(allenClMeta$subclass_id_label, getPairedTagSubclassName, "b")
scpt2ATAC <- data.frame(
  pt = scspt,
  atac = scsATAC
) |>
  unique() |>
  x => `rownames<-`(x, x$pt)
ordscsATAC <- scpt2ATAC[ordscs, "atac"]

## 3-1. get the corresponding logcpm matrix
xATAC <- xPbyS[partCREs, ordscsATAC]

## 3-2. get histone data
cicero_outd <- file.path(projd, "13.cicero", "out")
H3K27ac_logRPKM_allCREbysc <- readRDS(file.path(
  cicero_outd,
  "H3K27ac.logRPKM.allCRE.pbysc.rds"
))
H3K27me3_logRPKM_allCREbysc <- readRDS(file.path(
  cicero_outd,
  "H3K27me3.logRPKM.allCRE.pbysc.rds"
))
H3K4me1_logRPKM_allCREbysc <- readRDS(file.path(
  cicero_outd,
  "H3K4me1.logRPKM.allCRE.pbysc.rds"
))
H3K9me3_logRPKM_allCREbysc <- readRDS(file.path(
  cicero_outd,
  "H3K9me3.logRPKM.allCRE.pbysc.rds"
))
xK27ac <- H3K27ac_logRPKM_allCREbysc[partCREs, ordscs]
xK27me3 <- H3K27me3_logRPKM_allCREbysc[partCREs, ordscs]
xK4me1 <- H3K4me1_logRPKM_allCREbysc[partCREs, ordscs]
xK9me3 <- H3K9me3_logRPKM_allCREbysc[partCREs, ordscs]

## 4. plot Heatmap
genSimpleHeatmap <- function(x,
                             lowq = 0.005,
                             highq = 0.995,
                             title = "logCPM") {
  lowval <- quantile(x, probs = 0.005)
  highval <- quantile(x, probs = 0.995)
  enhColor <- circlize::colorRamp2(
    seq(lowval, highval, length = 60),
    viridis::viridis(60)
  )
  legendEnh <- list(
    title = title,
    at = c(lowval, highval),
    labels = round(c(lowval, highval), 2),
    direction = "horizontal"
  )

  ComplexHeatmap::Heatmap(
    matrix = x,
    col = enhColor,
    cluster_rows = F,
    cluster_columns = F,
    show_row_names = F,
    show_column_names = F,
    use_raster = T,
    show_heatmap_legend = T,
    heatmap_legend_param = legendEnh,
    width = 1
  )
}
n <- 1000
hmATAC <- genSimpleHeatmap(xATAC, lowq = 0.005, highq = 0.995)
withr::with_pdf(file.path(
  figd,
  str_glue("{state}.distal.ATAC.logcpm.ds{n}.heatmap.pdf")
), {
  print(hmATAC)
})

hmK27ac <- genSimpleHeatmap(xK27ac, lowq = 0.005, highq = 0.995)
withr::with_pdf(file.path(
  figd,
  str_glue("{state}.distal.K27ac.logcpm.ds{n}.heatmap.pdf")
), {
  print(hmK27ac)
})

hmK27me3 <- genSimpleHeatmap(xK27me3, lowq = 0.005, highq = 0.995)
withr::with_pdf(file.path(
  figd,
  str_glue("{state}.distal.K27me3.logcpm.ds{n}.heatmap.pdf")
), {
  print(hmK27me3)
})

hmK4me1 <- genSimpleHeatmap(xK4me1, lowq = 0.005, highq = 0.995)
withr::with_pdf(file.path(
  figd,
  str_glue("{state}.distal.K4me1.logcpm.ds{n}.heatmap.pdf")
), {
  print(hmK4me1)
})

hmK9me3 <- genSimpleHeatmap(xK9me3, lowq = 0.005, highq = 0.995)
withr::with_pdf(file.path(
  figd,
  str_glue("{state}.distal.K9me3.logcpm.ds{n}.heatmap.pdf")
), {
  print(hmK9me3)
})

# * plot eRNA signals
## 2. transform to logRPKM mat
m <- mPbyS[partCREs, ordscs]

## 3. scale across subclasses
m <- t(scale(t(m)))
m[is.na(m)] <- min(m, na.rm = T)

## 4. draw heatmap
lowval <- quantile(m, probs = 0.3)
highval <- quantile(m, probs = 0.9)
eRNAColor <- circlize::colorRamp2(
  seq(lowval, highval, length = 30),
  viridis::viridis(30)
)
legendeRNA <- list(
  title = "scaled logRPKM",
  at = c(lowval, highval),
  labels = round(c(lowval, highval), 2),
  direction = "horizontal"
)

hmeRNA <- ComplexHeatmap::Heatmap(
  matrix = m,
  col = eRNAColor,
  cluster_rows = F,
  cluster_columns = F,
  show_row_names = F,
  show_column_names = F,
  use_raster = T,
  show_heatmap_legend = T,
  heatmap_legend_param = legendeRNA,
  width = 1
)

withr::with_pdf(file.path(
  figd,
  str_glue("{state}.distal.ATAC.PairedTageRNA.logRPKM.ds{n}.heatmap.pdf")
), {
  print(hmeRNA)
})

hmsum <- hmATAC + hmeRNA + hmK27ac + hmK4me1

withr::with_pdf(file.path(
  figd,
  str_glue("EnhDAR.sum.heatmap.pdf")
), {
  draw(hmsum,
    column_title = "DAR Enhancers",
    column_title_gp = gpar(col = "black", fontsize = 16),
    row_title = "Enhancers", row_title_gp = gpar(col = "black", fontsize = 16),
    merge_legend = T
  )
})

saveRDS(object = list(
  ATAC = xATAC, K27ac = xK27ac,
  K4me1 = xK4me1, eRNA = m
), file = file.path(projd,
  "13.cicero", "out", "eRNA", "Fig4.PanelD.ds1000.rds"))
