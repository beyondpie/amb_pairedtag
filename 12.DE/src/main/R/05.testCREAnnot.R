library(tidyverse)
library(tmpRpkg)

# * meta
projd <- here::here()
DARfnm <- file.path(projd, "16.celloracle", "out", "DAR.log2fd0.2.csv")
promoterfnm <- file.path("/projects/ps-renlab2/szu",
  "/genome/mouseGenePromoter.GRCm38-M25.promoter.1kbAwayTSS.tsv")
CREAnnotd <- file.path(projd, "/06.ChromHMM/out/CREAnnot250429")
CREd <- file.path(projd, "data", "chromHMM", "CRE")

CAdivfnm <- file.path(projd, "meta", "mba.hba.CAdiver.orthInHg38.bed")
CAcnsrvfnm <- file.path(projd, "meta", "mba.hba.CAcons.orthInHg38.rmdup.bed")

CREHomerAnntfnm <- file.path(projd, "meta", "allCREHomerAnnot.txt")
CREHomerAnntCleanUpfnm <- file.path(projd, "meta", "CREHomerAnnotCleanUp.rds")

# * function
loadCREAnnot <- function(sc) {
  data.table::fread(
    file = file.path(CREAnnotd, str_glue("{sc}.CRE.annotBySummit.bed")),
    header = F, sep = "\t", data.table = F
  )
}

loadCRE <- function(sc, onlyDistal = F) {
  r <- data.table::fread(
    file = file.path(CREd, str_glue("{sc}.CRE.ChromHMM.DAR.distal.annot.tsv")),
    header = T, sep = "\t", data.table = F
  )
  if (onlyDistal) {
    r <- r[r$isDistal == "distal", ]
  }
  return(r)
}

checkDAR <- function(CRE, DARs, sc) {
  dCRE <- subset(CRE, isDAR == "DAR") |>
    x => with(x, paste(chrom,
      paste(startFrom, endTo, sep = "-"), sep = ":"))
  all(dCRE %in% with(DARs, CRE[subclass == sc]))
}

checkAnnot <- function(CRE, aCRE) {
  if (nrow(CRE) != nrow(aCRE)) {
    return(FALSE)
  }
  x <- with(CRE,
    paste(
      paste(chrom, paste(startFrom, endTo, sep = "-"), sep = ":"),
      chromHMMState, sep = "@")
  )
  y <- with(aCRE,
    paste(paste(V1, paste(V2, V3, sep = "-"), sep = ":"),
      V4, sep = "@")
  )
  return(all(x %in% y))
}

checkDistal <- function(CRE, promoters) {
  gHead <- c("chrom", "startFrom", "endTo", "isDistal") 
  gCRE <- tmpRpkg::loadGRfromBed(
    beds = CRE[, gHead], colnms = gHead)
  gPromoter <- tmpRpkg::loadGRfromBed(
    beds = promoters[, 1:4])
  hits <- GenomicRanges::findOverlaps(query = gCRE, subject = gPromoter)
  queryIndex <- S4Vectors::queryHits(hits)
  n <- sum(CRE[queryIndex, "isDistal"] == "proximal") == length(queryIndex)
  r <- CRE[queryIndex, ] |>
    x => x[x$isDistal != "proximal", ]
  s <- tmpRpkg::loadGRfromBed(beds = r[, gHead], colnms = gHead) |>
    GenomicRanges::findOverlaps(query = _, subject = gPromoter) |>
    y => gPromoter[S4Vectors::subjectHits(x = y)]
  print(r)
  print(s)
  if (nrow(r) < 5) {
    return(T)
  }
  return(F)
}


highTEGlutscs <- c(
  "016_CA1-ProS_Glut",
  "017_CA3_Glut",
  "001_CLA-EPd-CTX_Car3_Glut",
  "037_DG_Glut",
  "010_IT_AON-TT-DP_Glut",
  "002_IT_EP-CLA_Glut",
  "007_L2-3_IT_CTX_Glut",
  "008_L2-3_IT_ENT_Glut",
  "009_L2-3_IT_PIR-ENTl_Glut",
  "019_L2-3_IT_PPP_Glut",
  "011_L2_IT_ENT-po_Glut",
  "006_L4-5_IT_CTX_Glut",
  "003_L5-6_IT_TPE-ENT_Glut",
  "022_L5_ET_CTX_Glut",
  "005_L5_IT_CTX_Glut",
  "030_L6_CT_CTX_Glut",
  "004_L6_IT_CTX_Glut",
  "029_L6b_CTX_Glut",
  "014_LA-BLA-BMA-PA_Glut",
  "151_TH_Prkcd_Grin2c_Glut"
)

# * main
ptscMeta <- tmpRpkg::loadptscMeta() |>
  x => x[x$ATAC >0, ]
ptscs <- ptscMeta$PairedTagName
rownames(ptscMeta) <- ptscMeta$ATACName
highTEptscs <- ptscMeta[highTEGlutscs, "PairedTagName"]

DARs <- data.table::fread(
  file = DARfnm, header = T, sep = ",", data.table = F)

promoters <- data.table::fread(
  file = promoterfnm, header = F, sep = "\t", data.table = F
)

CAdiv <- data.table::fread(
  file = CAdivfnm, header = F,
  sep = "\t", data.table = F
) |>
  setNames(object = _, nm = c("chrom", "startFrom", "endTo", "mouseCRE")) |>
  x => `rownames<-`(x, x$mouseCRE)

CAcnsrv <- data.table::fread(file = CAcnsrvfnm, header = F,
  sep = "\t", data.table = F) |>
  setNames(object = _, nm = c("chrom", "startFrom", "endTo", "mouseCRE")) |>
  x => `rownames<-`(x, x$mouseCRE)

allCREHomerAnnot <- data.table::fread(
  file = CREHomerAnntfnm, header = T, sep = "\t", data.table = F)
allCREHomerAnnotCleanUp <- readRDS(CREHomerAnntCleanUpfnm)
rownames(allCREHomerAnnotCleanUp) <- allCREHomerAnnotCleanUp$PeakID

for (sc in ptscs){
  message("checking ", sc)
  CRE <- loadCRE(sc)
  aCRE <- loadCREAnnot(sc)
  a1 <- checkDAR(CRE, DARs, sc)
  a2 <- checkAnnot(CRE, aCRE)
  a3 <- checkDistal(CRE, promoters)
  if(all(c(a1, a2, a3))) {
    message(sc, " pass test.")
  } else {
    message(sc, " does not pass test.")
  }
}


allCRE <- lapply(ptscs, loadCRE) |>
  do.call(rbind, args = _)

allCREnm <- with(allCRE, paste(chrom,
  paste(startFrom, endTo, sep = "-"), sep = ":"))

uCREnm <- unique(allCREnm) # 932029
length(intersect(uCREnm, rownames(CAdiv)))  # 333699
length(intersect(uCREnm, rownames(CAcnsrv))) # 221395

table(allCRE$convervation[
  allCRE$chromHMMState == "Chr-A"]) / sum(allCRE$chromHMMState == "Chr-A")

## CAcnsrv   CAdiv   mouse 
## 3113051 1208244 1238782 
##   CAcnsrv     CAdiv     mouse 
## 0.5598935 0.2173071 0.2227994 

table(allCRE$convervation[
  allCRE$chromHMMState == "Chr-O"]) / sum(allCRE$chromHMMState == "Chr-O")
## CAcnsrv   CAdiv   mouse 
## 3039870 1984280 2156402 
##   CAcnsrv     CAdiv     mouse 
## 0.4233477 0.2763409 0.3003115

# * add CRE conservation on distal elements
allDistalCRE <- lapply(ptscs, \(sc) {loadCRE(sc, onlyDistal = T)}) |>
  do.call(rbind, args = _)

allDistalCREnm <- with(allDistalCRE, paste(chrom,
  paste(startFrom, endTo, sep = "-"), sep = ":"))

uCREnm <- unique(allDistalCREnm) # 874563
length(intersect(uCREnm, rownames(CAdiv)))  # 320900
length(intersect(uCREnm, rownames(CAcnsrv))) # 195181


table(allDistalCRE$convervation[
  allDistalCRE$chromHMMState == "Chr-A"]) /
  sum(allDistalCRE$chromHMMState == "Chr-A")

## CAcnsrv   CAdiv   mouse 
## 3113051 1208244 1238782 
##   CAcnsrv     CAdiv     mouse 
## 0.4699197 0.2700034 0.2600768

table(allDistalCRE$convervation[
  allDistalCRE$chromHMMState == "Chr-O"]) /
  sum(allDistalCRE$chromHMMState == "Chr-O")
## CAcnsrv   CAdiv   mouse 
## 3039870 1984280 2156402 
##   CAcnsrv     CAdiv     mouse 
## 0.3692140 0.3098343 0.3209516



# * check TE enrichment of different chromstates.
allDistalCRE$superFam <- allDistalCREHomerAnnotCleanUp[allDistalCREnm, "AnnoSuperFam1"]
TEcols <- c("LINE", "LTR", "SINE", "Other Repeats")

highTEState <- lapply(highTEptscs, \(sc) {
  t <- subset(allCRE, subclass == sc & superFam %in% TEcols)
  table(t$chromHMMState) |>
    as.data.frame(stringsAsFactors = F) |>
    setNames(object = _, nm = c("state", "counts"))
}) |>
  do.call(rbind, args = _)

## 1 Chr-A   201017
## 2 Chr-B    10479
## 3 Chr-O   351559
## 4 Chr-R    19967
## 5 Hc-H      2230
## 6 Hc-P       386
## 7 ND        4905
highTEStateRatio <- highTEState |>
  group_by(state) |>
  summarise(
    cumCount = sum(counts)
  )

ustates <- c("Chr-A", "Chr-O", "Chr-B", "Chr-R", "Hc-H", "Hc-P", "ND")

TERatioChrA <- lapply(ptscs, \(sc) {
  t <- subset(allCRE, subclass == sc & chromHMMState == "Chr-A")
  data.frame(
    nTE = sum(t$superFam %in% TEcols),
    n = nrow(t),
    r = sum(t$superFam %in% TEcols) / nrow(t),
    sc = sc
  )
}) |> do.call(rbind, args = _)

rChrAHighTE <- TERatioChrA[TERatioChrA$sc %in% highTEptscs, ]
## mean(rChrAHighTE$r)
## [1] 0.1144803


rChrANotHighTE <- TERatioChrA[!TERatioChrA$sc %in% highTEptscs, ]
## > mean(rChrANotHighTE$r)
## [1] 0.06143528


TERatioChrO <- lapply(ptscs, \(sc) {
  t <- subset(allCRE, subclass == sc & chromHMMState == "Chr-O")
  data.frame(
    nTE = sum(t$superFam %in% TEcols),
    n = nrow(t),
    r = sum(t$superFam %in% TEcols) / nrow(t),
    sc = sc
  )
}) |> do.call(rbind, args = _)

rChrOHighTE <- TERatioChrO[TERatioChrO$sc %in% highTEptscs, ]
## mean(rChrOHighTE$r)
## [1] 0.2033028

rChrONotHighTE <- TERatioChrO[!TERatioChrO$sc %in% highTEptscs, ]
## mean(rChrONotHighTE$r)
## [1] 0.1237797

## with(TERatioChrA, sum(nTE) / sum(n))
## [1] 0.09313846
## > 
## with(TERatioChrO, sum(nTE) / sum(n))
## > [1] 0.1421733

# * save data
data.table::fwrite(x = allCRE,
  file = file.path(projd, "data", "chromHMM", "allCRE.amb.PairedTag.annot.tsv"),
  sep = "\t", row.names = F, col.names = T)
