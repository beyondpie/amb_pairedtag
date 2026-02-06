library(Seurat)
library(tidyverse)
s <- str_glue

# * meta
projd <- here::here()
snATACd <- file.path(projd, "data", "snATAC")
snATACtfd <- file.path(snATACd, "snATACTransferLabel")
snATACMeta <- readRDS(
  file.path(snATACd, "mba.whole.cell.meta.v9.7.rds"))
L4topL4 <- snATACMeta[, c("L4", "pL4")] |>
  unique() |>
  x => `rownames<-`(x, x$L4)

allenClMeta <- data.table::fread(
  file = file.path(projd, "meta",
    "AIT21_annotation_freeze_081523.tsv"),
  sep = "\t", header = TRUE, data.table = FALSE
) |>
  x => `rownames<-`(x, paste0("cl-", x[, 1]))

allPeaks <- data.table::fread(
  file = file.path(snATACd, "wmb.snATAC.peak.srt.bed"),
  header = FALSE, data.table = FALSE
) |>
  setNames(object = _, nm = c("chrom", "start", "end", "name"))

# * function
transformAllenLabel2str <- function(allenLabels) {
  vapply(allenLabels, \(a) {
    gsub(" |-|/", "_", a) 
  }, "a")
}

map2snATACscstr <- function(sc) {
  str_replace_all(sc, "^\\d+ ", "") |>
    x => str_replace_all(x, " +",  "_") |>
    x => str_replace_all(x, "/", "-") 
}

sortAllenLabelById <- function(allenLabels) {
  ids <- vapply(allenLabels, \(a) {
    as.integer(str_split_1(a, "_")[1])
  }, 1)
  allenLabels[order(ids)]
}

loadpL4Peak <- function(pL4) {
  data.table::fread(file = file.path(snATACd,
    "sa2.pL4.final.peak.srt", s("{pL4}.bed")),
    header = FALSE, data.table = FALSE
  ) |>
    setNames(object = _, nm = c("chrom", "start", "end", "name"))
}

loadscPeak <- function(sc) {
  data.table::fread(file = file.path(snATACd,
    "sa2.subclassv3.final.peak.srt",
    str_glue("{sc}.bed")
  ), header = FALSE, data.table = FALSE) |>
    setNames(object = _,
      nm = c("chrom", "start", "end", "name"))
}

# * main
# sp 2 sc
sp2sc <- allenClMeta[, c("supertype_id_label",
  "subclass_id_label")] |>
  unique() |>
  setNames(object = _, nm = c("sp", "sc")) |>
  x => `rownames<-`(x, x$sp)
sp2sc$sp <- transformAllenLabel2str(sp2sc$sp)
sp2sc$sc <- transformAllenLabel2str(sp2sc$sc)

# pL4 to sc
pL4tosc <- snATACMeta[, c("pL4", "subclass_id_label_v3")] |>
  unique() |>
  setNames(object = _, nm = c("pL4", "sc")) |>
  x => `rownames<-`(x, x$pL4)
pL4tosc$sc <- transformAllenLabel2str(pL4tosc$sc)

# get groups (subclass + supertype) of pariedtag
ptDNAd <- file.path(projd, "data", "ptDNAbam", "bam")
ptGroups <- list.dirs(ptDNAd, full.names = FALSE, recursive = FALSE)

# map groups to pL4, then get the pL4-level peaks
nn2sp <- data.table::fread(
  file = file.path(snATACtfd, "snATAC_nnL4tosp_240814.csv"),
  header = TRUE, sep = ",", data.table = FALSE
) |> x => `rownames<-`(x, x[, 1])

neu2sc <- unique(snATACMeta[
  snATACMeta$NT_v3 != "NN",
  c("L4", "subclass_id_label_v3")
]) |>
  setNames(object = _, c("L4", "sc"))

nnsp <- unique(transformAllenLabel2str(nn2sp$sp)) |>
  sortAllenLabelById()
neusc <- unique(transformAllenLabel2str(neu2sc$sc)) |>
  sortAllenLabelById()

sum(ptGroups %in% c(nnsp, neusc))

ptGroups[!(ptGroups %in% c(nnsp, neusc))]
##  [1] "077_CEA_BST_Gal_Avp_Gaba"   "092_TMv_PMv_Tbx3_Hist_Gaba"
##  [3] "096_PVHd_Gsc_Gaba"          "103_PVHd_DMH_Lhx6_Gaba"
##  [5] "105_TMd_DMH_Foxd2_Gaba"     "107_DMH_Hmx2_Gaba"
##  [7] "108_ARH_PVp_Tbx3_Gaba"      "122_LHA_MEA_Otp_Glut"
##  [9] "123_DMH_Nkx2_4_Glut"        "125_DMH_Hmx2_Glut"
## [11] "140_PMd_LHA_Foxb1_Glut"     "143_MM_ant_Foxb1_Glut"
## [13] "191_PAG_MRN_Rln3_Gaba"
groups <- ptGroups[ptGroups %in% c(nnsp, neusc)]

# map groups to pL4
nn2sp$sc <- NULL
nn2sp$sp <- transformAllenLabel2str(nn2sp$sp)
colnames(nn2sp) <- c("L4", "group")
neu2sc$sc <- transformAllenLabel2str(neu2sc$sc)
colnames(neu2sc) <- c("L4", "group")
rawL4toGroup <- rbind(nn2sp, neu2sc)
rawL4toGroup$pL4 <- L4topL4[rawL4toGroup$L4, "pL4"]
pL4toGroup <- unique(rawL4toGroup[, c("pL4", "group")])

data.table::fwrite(pL4toGroup,
  file = file.path(snATACtfd, "allpL4toneusc_nnsp.csv"),
  col.names = TRUE, row.names = FALSE, sep = ",")

pL4toGroup.pt <- pL4toGroup[pL4toGroup$group %in% groups, ]

# merge pL4's peaks (double check all of them in our final peak lists)
allpL4Peaks <- lapply(pL4toGroup$pL4, loadpL4Peak)
names(allpL4Peaks) <- pL4toGroup$pL4
# length(unique(do.call(rbind, allpL4Peaks)[, 4]))
peaksPtGroups <- allpL4Peaks[pL4toGroup.pt$pL4] |>
  do.call(what = rbind, args = _) |>
  x => unique(x[,4])

data.table::fwrite(
  data.frame(name = peaksPtGroups),
  file = file.path(snATACd,
    "snATAC.peaks.included.pairedtag.csv"),
  col.names = FALSE, row.names = FALSE)

# get the specificity score (p-value or other measurements)

# get the enrichment score (pseudo-bulk level CPM) 
sc.pt <- pL4tosc[pL4toGroup.pt$pL4, "sc"] |>
  unique()
scLabelPt <- str_replace_all(sc.pt, "^\\d+_", "")

cpm_pbysc <- data.table::fread(
  file.path(snATACd, "cpm_peakBysubclass.csv"),
  sep = ",", header = TRUE, data.table = FALSE)
rownames(cpm_pbysc) <- cpm_pbysc$V1
cpm_pbysc$V1 <- NULL
old_names <- colnames(cpm_pbysc)
colnames(cpm_pbysc) <- str_replace_all(old_names,"-", "_")
cpmPeakByScPt <- cpm_pbysc[peaksPtGroups, scLabelPt]

# get subclass-level peaks
# check original results
pL4tosc <- snATACMeta[, c("pL4", "subclass_id_label_v3")] |>
  unique() |>
  setNames(object = _, c("pL4", "sc")) |>
  x => `rownames<-`(x, x$pL4)
pL4tosc$rawsc <- pL4tosc$sc
ptsc <- transformAllenLabel2str(pL4tosc$sc)
pL4tosc$sc <- ptsc
snATACsc <- map2snATACscstr(pL4tosc$sc)
pL4tosc$snATACsc <- snATACsc

allsc <- unique(pL4tosc$sc)

allsc2peaks <- lapply(allsc, \(sc) {
  pL4s <- pL4tosc[pL4tosc$sc == sc, "pL4"]
  allpL4Peaks[pL4s] |>
    do.call(what = rbind, args = _) |>
    x => unique(x[, 4])
}) |>
  setNames(object = _, nm = allsc)

rawsc2peaks <- lapply(unique(pL4tosc$snATACsc),loadscPeak) |>
  setNames(object = _, nm = allsc)

# TRUE
all(vapply(allsc, \(sc) {
  length(allsc2peaks[[sc]]) == nrow(rawsc2peaks[[sc]])
}, TRUE))

# TRUE
all(vapply(allsc, \(sc) {
   all(allsc2peaks[[sc]] %in% rawsc2peaks[[sc]][, 4])
}, TRUE))

sc2peaks.pt <- rawsc2peaks[sc.pt]

a <- cpmPeakByScPt[sc2peaks.pt[[sc.pt[1]]]$name,
  gsub("^\\d+_", "", sc.pt[1])]

b <- cpmPeakByScPt[
  !rownames(cpmPeakByScPt) %in% sc2peaks.pt[[sc.pt[1]]]$name,
  gsub("^\\d+_", "", sc.pt[1])]

# get CPM for each subclass's peaks

