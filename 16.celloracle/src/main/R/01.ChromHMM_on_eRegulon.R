# ChromHMM on eRegulon
# 1. ChrA-based enhancers: positive correlation on eRegulon
# 2. subclass-specific TFs:
#    - ChrA-based enhancers + DAR + eRegulon

suppressMessages({
  suppressWarnings({
    library(tidyverse)
    # detach("package:tmpRpkg", unload = T)
    library(tmpRpkg)
    library(GenomicRanges)
    library(S4Vectors)
    library(rtracklayer)
    library(ComplexHeatmap)
    library(reticulate)
  })
})

tmpRpkg::setCondaEnv()
mypy <- import("tmpPypkg")
pd <- import("pandas")
myad <- import("anndata")

## fnm <- "{projd}/data/pairedtag_ann/pt.RNAseq.scbyg.CPM.csv"
## x <- pd$read_csv(fnm, sep = ",", header = 0L, index_col = 0L) |>
##   reticulate::py_to_r(x = _)

# * meta
projd <- tmpRpkg:::projd
CREAnnotd <- file.path(projd, "06.ChromHMM/out",
  "CREAnnot250429")
annotd <- file.path(
  projd, "06.ChromHMM", "out", "CREAnnot250429")
eRegulond <- file.path(
  projd, "16.celloracle", "out", "eRegulon")

workd <- file.path(projd, "16.celloracle")
outd <- file.path(workd, "out")

ptscMeta <- loadptscMeta() |>
  subset(x = _, subset = ATAC > 0) |>
  x => `rownames<-`(x, x$PairedTagName)

mm10promoter <- loadPromoter()
AllenTFOrd <- loadAllenTFOrd()

# * funcitons
load_eRegulon <- function(ptsc) {
  fnm <- file.path(
    eRegulond,
    str_glue("cellOracle.{ptsc}.eRegulon.csv")
  )
  data.table::fread(
    file = fnm, sep = ",", header = T,
    data.table = F
  )
}

retrieveCREv2 <- function(sc, deATAC, log2fd = 0.2,
                          n = 99999999L) {
  m <- deATAC[[sc]]
  index <- m[, 2] >= log2fd
  if (sum(index) < 1) {
    message(str_glue("{sc} has no CREs under log2fd {log2fd}."))
    index <- seq_len(nrow(m))
  }
  r <- m[order(m[index, 2], decreasing = T), ]
  colnames(r) <- c("CRE", "log2fd", "p", "q")
  return(r[seq_len(min(nrow(r), n)), ])
}

getDESum <- function(deList, log2fd = 0.2) {
  deCREList <- lapply(names(deList), retrieveCREv2,
    deATAC = deList, log2fd = log2fd, n = 999999
  ) |>
    setNames(object = _, nm = names(deList))
  deSum <- lapply(names(deCREList), \(ptsc) {
    r <- deCREList[[ptsc]]
    r["subclass"] <- ptsc
    r
  }) |> do.call(rbind, args = _)
  return(deSum)
}

# * main
# * check positive correlations
# ** data loadings
# 1. ChrA-based enhancers
scEnhs <- lapply(ptscMeta$PairedTagName, \(sc) {
  rawAnnot <- readSummitAnnot(
    file.path(CREAnnotd, str_glue("{sc}.CRE.annotBySummit.bed")),
    mm10promoter,
    sc
  )
  annotCREs <- rawAnnot |>
    x => x[x$state == "Chr-A", ] |>
    x => x[x$pd == "distal", ]
}) |>
  setNames(object = _, nm = ptscMeta$PairedTagName)

scEnhdf <- lapply(ptscMeta$PairedTagName, \(sc) {
  data.frame(
    CRE = scEnhs[[sc]]$CRE,
    subclass = sc
  )
}) |> do.call(rbind, args = _)

data.table::fwrite(scEnhdf, file.path(
  projd, "data", "enhancer",
  "enhancer.csv"
), sep = ",", col.names = F, row.names = F)

# - without removing distal
scChrAs <- lapply(ptscMeta$PairedTagName, \(sc) {
  rawAnnot <- readSummitAnnot(
    file.path(CREAnnotd, str_glue("{sc}.CRE.annotBySummit.bed")),
    mm10promoter,
    sc
  )
  annotCREs <- rawAnnot |>
    x => x[x$state == "Chr-A", ]
}) |>
  setNames(object = _, nm = ptscMeta$PairedTagName)

scChrAdf <- lapply(ptscMeta$PairedTagName, \(sc) {
  data.frame(
    CRE = scChrAs[[sc]]$CRE,
    subclass = sc
  )
}) |> do.call(rbind, args = _)

data.table::fwrite(scChrAdf, file.path(
  projd, "data", "enhancer",
  "ChromHMMChrA.CRE.csv"
), sep = ",", col.names = F, row.names = F)

# add Chr-O CREs
scChrOs <- lapply(ptscMeta$PairedTagName, \(sc) {
  rawAnnot <- readSummitAnnot(
    file.path(CREAnnotd, str_glue("{sc}.CRE.annotBySummit.bed")),
    mm10promoter,
    sc
  )
  annotCREs <- rawAnnot |>
    x => x[x$state == "Chr-O", ]
}) |>
  setNames(object = _, nm = ptscMeta$PairedTagName)

scChrOdf <- lapply(ptscMeta$PairedTagName, \(sc) {
  data.frame(
    CRE = scChrOs[[sc]]$CRE,
    subclass = sc
  )
}) |> do.call(rbind, args = _)

data.table::fwrite(scChrOdf, file.path(
  projd, "data", "enhancer",
  "ChromHMMChrO.CRE.csv"
), sep = ",", col.names = F, row.names = F)


# 2. eRegulon loading
sceRegulons <- lapply(ptscMeta$PairedTagName, load_eRegulon) |>
  setNames(object = _, nm = ptscMeta$PairedTagName)

# 3. get enhancer-specific eRegulons
scEnhRegulons <- lapply(ptscMeta$PairedTagName, \(sc) {
  eRegulon <- sceRegulons[[sc]]
  enh <- scEnhs[[sc]]
  r <- subset(eRegulon, subset = CRE %in% enh$CRE)
  if (nrow(r) < 1) {
    message("No eReuglon found for ", sc)
  }
  return(r)
}) |>
  setNames(object = _, nm = ptscMeta$PairedTagName)
saveRDS(
  object = scEnhRegulons,
  file = file.path(outd, "scellOracle.enhancer.subclass-eRegulon.rds")
)

# ** show the data
eRegulonSum <- lapply(ptscMeta$PairedTagName, \(sc) {
  eRegulon <- sceRegulons[[sc]]
  r <- data.frame(
    CRE = eRegulon$CRE,
    motif = eRegulon$motif,
    tf = eRegulon$tf,
    gene = eRegulon$gene,
    coef = eRegulon$coefCellOracle,
    cor = eRegulon$corSubclass,
    conn = eRegulon$connCicero,
    subclass = sc,
    label = "others"
  )
  r$label[r$CRE %in% scEnhs[[sc]]$CRE] <- "enhancer"
  r$label[r$CRE %in% scChrOs[[sc]]$CRE] <- "ChrO"
  return(r)
}) |> do.call(what = rbind, args = _) |>
  x => unique(x)

data.table::fwrite(
  x = eRegulonSum,
  file = file.path(
    outd, "all.eRegulon.csv"
  ),
  sep = ",",
  col.names = T,
  row.names = F
)

table(eRegulonSum$label)
## enhancer   others 
##  2104532  4138122 

# 400 unique TFs
# 2842 unique genes
# 59090 enhancers out of unique 133,045 CREs

# 1. in each subclass, # of positive links
sc2npos <- vapply(
  ptscMeta$PairedTagName, \(sc) {
    y <- subset(eRegulonSum, subclass == sc)
    enhy <- subset(y, label == "enhancer")
    oy <- subset(y, label != "enhancer")
    npos_enh <- ifelse(test = nrow(enhy) > 0,
      yes = sum(enhy$coef > 0) / nrow(enhy), no = 0.0
    )
    npos_o <- ifelse(test = nrow(oy) > 0,
      yes = sum(oy$coef > 0) / nrow(oy), no = 0.0
    )
    c(nrow(enhy), npos_enh * 100, nrow(oy), npos_o * 100)
  },
  FUN.VALUE = c(100, 0.1, 100, 0.01)
)
# no significant diff
# 62.7%
sc2npos <- t(sc2npos)

# 2. in each subclass, positive corrs
# no significant diff
# 0.05
sc2coef <- vapply(ptscMeta$PairedTagName, \(sc) {
  y <- subset(eRegulonSum, subclass == sc)
  enhy <- subset(y, label == "enhancer")
  oy <- subset(y, label != "enhancer")
  c(mean(enhy$coef), mean(oy$coef))
}, FUN.VALUE = c(0.0, 0.0)) |>
  x => t(x)

# * subclass-specific TFs
# ** data loadings
# 1. load DARs of ATAC-seq
deATAC <- readRDS(tmpRpkg:::snATACDEfnm)
deCRESum <- getDESum(deATAC, log2fd = 0.2)
# 767393 unique DAR
data.table::fwrite(
  x = deCRESum,
  file = file.path(
    projd, "data", "enhancer",
    "DAR.log2fd0.2.csv"
  ),
  sep = ",",
  col.names = T
)

# 1*. load DAR Enhancers of ATAC-seq
deATAChrA <- readRDS(tmpRpkg:::snATAChrADEfnm)
deEnhSum <- getDESum(deATAChrA, log2fd = 0.2)
# 470176 DAR enhancers
data.table::fwrite(
  x = deEnhSum,
  file = file.path(
    projd, "data", "enhancer",
    "DAR.enhancer.log2fd0.2.csv"
  ),
  sep = ",",
  col.names = T
)

# 2. load AllenRNA or PairedTag RNA data.
ptRNA_logCPM_scbyg <- mypy$globalvar$load_PairedTag_RNA_logCPM_scbyg()
ptRNAlogCPM <- py_to_r(ptRNA_logCPM_scbyg)
attr(ptRNAlogCPM, "pandas.index") <- NULL

# ** get differential TFs by intersecting
#    DARs, ChrA-based enhancers, TFs

## > dim(deRegulonSum)
## 2933246 x 9
## > dim(eRegulonSum)
## 6242654 x 9

deRegulonSum <- merge(
  x = eRegulonSum,
  y = deCRESum,
  by.x = c("CRE", "subclass"),
  by.y = c("CRE", "subclass")
)
rawColnms <- colnames(deRegulonSum)
newColnms <- rawColnms
newColnms[seq(
  from = length(newColnms) - 3,
  to = length(newColnms),
)] <- c("CRElabel", "DARlog2fd", "DARpval", "DARq")
colnames(deRegulonSum) <- newColnms
data.table::fwrite(
  x = deRegulonSum, file = file.path(outd, "DAReRegulon.csv"),
  sep = ","
)

## table(deRegulonSum$CRElabel)
## enhancer   others
##   788529   2144717
## length(unique(deRegulonSum$CRE))
## [1] 97237
## > length(unique(deEnhRegulonSum$CRE))
## [1] 33295
deEnhRegulonSum <- subset(deRegulonSum, CRElabel == "enhancer")
data.table::fwrite(
  x = deEnhRegulonSum, file = file.path(outd, "DAR.enhancer.eRegulon.csv"),
  sep = ","
)


# 47433 uniq enhancers out of 90631 rows.
deEnhInRegulon <- deEnhRegulonSum[, c("CRE", "subclass")] |> unique()

# dim: 33903 x 2
# 399 unique TFs
deEnhTF <- deEnhRegulonSum[, c("tf", "subclass")] |> unique()
