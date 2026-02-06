library(tidyverse)
library(tmpRpkg)

projd <- here::here()
# test the results of filterByRanges in GenomicRange

DARfnm <- file.path(projd, "16.celloracle/out/DAR.log2fd0.2.csv")
ChrAfnm <- file.path(projd, "data/enhancer/ChromHMMChrA.CRE.csv")

ptscMeta <- tmpRpkg::loadptscMeta()
ptscs <- ptscMeta$PairedTagName

allDAR <- data.table::fread(file = DARfnm,
  header = T, sep = ",", data.table = F)

allChrA <- data.table::fread(file = ChrAfnm,
  header = F, sep = ",", data.table = F)

sc <- ptscs[153]
r <- intersect(with(allDAR, CRE[subclass == sc]),
  allChrA[allChrA$V2 == sc, "V1"])
