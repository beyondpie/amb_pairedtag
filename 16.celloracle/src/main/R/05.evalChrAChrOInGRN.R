library(tidyverse)
library(tmpRpkg)

# * meta
projd <- here::here()
eRegulon <- data.table::fread(
  file = file.path(projd, "16.celloracle", "out", "all.eRegulon.csv"),
  header = T, sep = ",", data.table = F
)
ptscMeta <- loadptscMeta() |>
  subset(x = _, subset = ATAC > 0) |>
  x => `rownames<-`(x, x$PairedTagName)
ptscs <- ptscMeta$PairedTagName


# * load allCREs
allCRE <- data.table::fread(
  file = file.path(projd, "data", "chromHMM",
    "allCRE.amb.PairedTag.annot.tsv"),
  header = T, sep = "\t", data.table = F
)
allCRE$CRE <- with(allCRE, paste(chrom, paste(
  startFrom, endTo,  sep = "-"
) ,sep = ":"))

# * if ChrAs more enriched in DARs
# 0.4352
rDARChrA <- with(allCRE, sum(
  isDAR == "DAR" & chromHMMState == "Chr-A") / sum(chromHMMState == "Chr-A"))
# 0.4295
rDARChrO <- with(allCRE, sum(
  isDAR == "DAR" & chromHMMState == "Chr-O") / sum(chromHMMState == "Chr-O"))

prop.test(
  c(sum(allCRE$isDAR == "DAR" & allCRE$chromHMMState == "Chr-A"),
    sum(allCRE$isDAR == "DAR" & allCRE$chromHMMState == "Chr-O")),
  c(sum(allCRE$chromHMMState == "Chr-A"),
    sum(allCRE$chromHMMState == "Chr-O")))

# * if ChrAs more enriched in CellOracle
r <- lapply(ptscs, \(sc) {
  ## pCRE <- with(eRegulon, CRE[coef > 0 & subclass == sc]) |>
  ##   unique()
  pCRE <- with(eRegulon, CRE[subclass == sc]) |>
    unique()
  CREChrA <- with(allCRE, CRE[subclass == sc & chromHMMState == "Chr-A"])
  CREChrO <- with(allCRE, CRE[subclass == sc & chromHMMState == "Chr-O"])
  rChrA <- length(intersect(CREChrA, pCRE)) / length(CREChrA)
  rChrO <- length(intersect(CREChrO, pCRE)) / length(CREChrO)  
  data.frame(
    sc = sc,
    nChrA = length(CREChrA),
    nChrO = length(CREChrO),
    posChrA = length(intersect(CREChrA, pCRE)),
    posChrO = length(intersect(CREChrO, pCRE)),
    rposChrA = rChrA,
    rposChrO = rChrO
  )
}) |> do.call(rbind, args = _)

# if using all the positive eRegulon
# 0.025 
sum(r$posChrA) / sum(r$nChrA)
# 0.028
sum(r$posChrO) / sum(r$nChrO)

# if using all the eRegulons
# 0.027
# 0.031

# * check eRegulon CRE themselves
eRegulonCREAnnots <- lapply(ptscs, \(sc) {
  e <- with(eRegulon, CRE[subclass == sc]) |>
    unique()
  c <- subset(allCRE, subclass == sc) |>
    x => `rownames<-`(x, x$CRE)
  table(c[e, "chromHMMState"]) |>
    as.data.frame(stringsAsFactors = FALSE) |>
    setNames(object = _, nm = c("state", "count"))
}) |>
  setNames(object = _, nm = ptscs)


