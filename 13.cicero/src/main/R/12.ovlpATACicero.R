suppressMessages({
  suppressWarnings({
    library(tidyverse)
    library(tmpRpkg)
    # detach("package:tmpRpkg", unload = T)
    library(GenomicRanges)
    library(S4Vectors)
    library(rtracklayer)
    library(ComplexHeatmap)
  })
})

# * meta
projd <- here::here()
workd <- file.path(projd, "13.cicero")
cicero_outd <- file.path(workd, "out")
# outd <- file.path(workd, "out")
# figd <- file.path(workd, "figure", "cicero")
outd <- file.path(projd, "99.figures", "out", "fig4")
figd <- outd

snATACd <- file.path(projd, "data", "snATAC")
ptscMeta <- loadptscMeta() |>
  x => `rownames<-`(x, x$PairedTagName)
sc2nmod <- tmpRpkg::getPairedTagSubclass2Nmod(ptscMeta)
## filter number of DNA cells for subclasses
# ptscs <- sc2nmod$sc[sc2nmod$all >= 1000]
ptscs <- ptscMeta$PairedTagName

atac2ptsc <- data.frame(
  atac = gsub("^\\d+_", "", ptscMeta$ATACName),
  pt = ptscMeta$PairedTagName,
  me = gsub("^\\d+_", "", ptscMeta$DNAmethName)
) |> x => `rownames<-`(x, x$pt)

geneCPM <- getCPMOfptRNA(tmpRpkg:::ptPseudoRNAcountfnm)
# prepare a csv file format for python's pandaframe to use
data.table::fwrite(
  x = geneCPM,
  file = file.path(
    projd, "data", "pairedtag_ann",
    "pt.RNAseq.scbyg.CPM.csv"
  ),
  row.names = T, col.names = T
)
# for later usage.
saveRDS(geneCPM, file.path(projd, "data", "pairedtag_ann",
  "pt.RNAseq.pseudobulk.CPM.scbyg.rds"))

geneScaledLogCPM <- log1p(geneCPM) |>
  scale()
distal2gene <- loadCRE2Gene(tmpRpkg:::ppdcWithCorfnm,
  genes = colnames(geneCPM)
)

eRNAd <- file.path(projd, "13.cicero", "out", "eRNA")
eRNAlogRPKMfnm <- file.path(eRNAd,
  "eRNA.logRPKM.filteredDAREhn.pbysc.rds")
eRNAmat <- readRDS(eRNAlogRPKMfnm)
fdPbyS <- data.table::fread(
  file = file.path(tmpRpkg:::eRNAfdpath, "fd.csv"),
  header = F, sep = ",", data.table = F
) |>
  setNames(object = _, nm = colnames(eRNAmat)) |>
  x => `rownames<-`(x, rownames(eRNAmat))

# * functions
reduceCREbyeRNA <- function(sc, CRElist, fdPbyS,
                            topk = 2000, eRNAfd = 1.1) {
  message(str_glue("reduce CREs by eRNA: {sc}."))
  m <- CRElist[[sc]]
  f <- fdPbyS[m[m %in% rownames(fdPbyS)], sc]
  index <- (f >= eRNAfd)
  if (sum(index) < 1) {
    message(str_glue("{sc} has no CREs under eRNA fd {eRNAfd}."))
    return(NULL)
  }
  mm <- m[index]
  return(mm[seq_len(min(topk, sum(index)))])
}

drawHeatmap <- function(x, lowQ = 0.005, highQ = 0.995, outf) {
  lowval <- quantile(x, probs = lowQ)
  highval <- quantile(x, probs = highQ)
  color <- circlize::colorRamp2(
    seq(lowval, highval, length = 60),
    viridis::viridis(60)
  )
  lgd <- list(
    title = "logCPM",
    at = c(lowval, highval),
    labels = round(c(lowval, highval), 2),
    direction = "horizontal"
  )

  hm <- ComplexHeatmap::Heatmap(
    matrix = x,
    col = color,
    cluster_rows = F,
    cluster_columns = F,
    show_row_names = F,
    show_column_names = F,
    use_raster = T,
    show_heatmap_legend = T,
    heatmap_legend_param = lgd,
    width = 1
  )
  withr::with_pdf(outf, {
    print(hm)
  })
  return(hm)
}

# * main
deATAChrA <- readRDS(tmpRpkg:::snATAChrADEfnm)
deEnhList <- lapply(names(deATAChrA), retrieveCRE,
  deATAC = deATAChrA, log2fd = 0.2, n = 999999
) |>
  setNames(object = _, nm = names(deATAChrA))
allDEnh <- unique(unlist(deEnhList))
ppecs <- loadAllppec(tmpRpkg:::snATACppecd)

allppec <- lapply(ppecs, \(ppec) {
  with(ppec, paste(dchr, paste(dstart, dend, sep = "-"), sep = ":"))
}) |>
  unlist() |>
  unique()
allDEppec <- intersect(allDEnh, allppec)

# 1. check positive proximal-enhancer result
# including ATAC, H3K27ac, H3K4me1, H3K27me3
# - load data
atacPbyS <- loadATACPseudoBulkMat(
  f = tmpRpkg:::snATACPMpbyc,
  env = "sa2stable"
) |>
  log1p()

pteRNA_logRPKM_allCREbysc <- readRDS(file.path(
  cicero_outd,
  "eRNA.PairedTag.logRPKM.allCRE.pbysc.rds"
))
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

# - prepare order of CREs
n <- 1000
partCREList_eRNAfd <- lapply(ptscs, \(sc) {
  
})
ordCREList <- lapply(ptscs, \(sc) {
  enhs <- deEnhList[[sc]]
  index <- enhs %in% allDEppec
  if (sum(index) < 1) {
    message(sc, " has no DEppecs. ")
    return(NULL)
  }
  enhs[index][1:min(sum(index), n)]
}) |>
  setNames(object = _, nm = ptscs) |>
  filterNULLfromList(l = _)

partCREList <- lapply(names(ordCREList), \(sc) {
  reduceCREbyeRNA(
    sc = sc, ordCREList,
    fdPbyS = eRNAmat, topk = n,
    eRNAfd = 1.1)
}) |>
  setNames(object = _, nm = names(ordCREList)) |>
  filterNULLfromList()

allOrdCRE <- unlist(partCREList) |>
  x => x[!duplicated(x)] |>
  x => x[x %in% rownames(eRNAmat)]

ordptsc <- names(partCREList) |>
  x => x[atac2ptsc[x, "atac"] %in% colnames(atacPbyS)]

atacPbyS4hm <- atacPbyS[allOrdCRE, atac2ptsc[ordptsc, "atac"]]

# - plot Heatmap
drawHeatmap(atacPbyS4hm, outf = file.path(figd, "test.ATAC.pdf"))
drawHeatmap(H3K27acEnhByS[allOrdCRE, ordptsc],
  outf = file.path(figd, "test.H3K27ac.pdf")
)

# results for H3K27me3 cannot be shown well
## drawHeatmap(H3K27me3EnhByS[allOrdCRE, atac2ptsc$pt],
##   file.path(figd, "test.H3K27me3.pdf"), lowQ = 0.000, highQ = 0.3)

drawHeatmap(H3K4me1EnhByS[allOrdCRE, ordptsc],
  lowQ = 0.05, highQ = 0.95,
  file.path(figd, "test.H3K4me1.pdf")
)

# 2. add gene expression to ppdc
ordGene <- distal2gene[allOrdCRE, "geneByPeak"] |>
  x => x[x %in% colnames(geneScaledLogCPM)]
ptRNAmat <- geneScaledLogCPM[ordptsc, ordGene] |>
  t()
drawHeatmap(ptRNAmat,
  lowQ = 0.2, highQ = 0.95,
  outf = file.path(figd, "test.RNA.pdf")
)

# 3. add eRNA to ppdc
eRNAmat4hm <- eRNAmat[allOrdCRE, ordptsc] |>
  x => t(scale(t(x)))
eRNAmat4hm[is.na(eRNAmat4hm)] <- min(eRNAmat4hm, na.rm = T)

drawHeatmap(eRNAmat4hm,
  lowQ = 0.3, highQ = 0.9,
  outf = file.path(figd, "test.eRNA.pdf")
)

# save all the heatmaps
saveRDS(object = list(
  ATAC = atacPbyS4hm,
  geneExp = ptRNAmat,
  eRNA = eRNAmat4hm,
  H3K4me1 = H3K4me1EnhByS[allOrdCRE, ordptsc],
  H3K27ac = H3K27acEnhByS[allOrdCRE, ordptsc]
), file = file.path(
  figd, "PanelD.DAREnhCRE.allmat.rds"
))

# 4. add RNA exp from Allen
## save it as csv file for python pandas to use.
allenRNAlogCPM <- readRDS(
  file.path(projd, "data", "allen", "sa2.allen.avg.logCPM.gbysc.rds")
)

data.table::fwrite(
  x = t(allenRNAlogCPM),
  file = file.path(
    projd, "data", "allen",
    "sa2.allen.avg.logCPM.scbyg.csv"
  ),
  sep = ",", row.names = T, col.names = T
)
