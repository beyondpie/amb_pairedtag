# set R lib search path
# use this to fix circlize error in my computer
.libPaths(rev(.libPaths()))

library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(ggpubr)



# * functions
readTEsSignal <- function(f) {
  r <- data.table::fread(file = f, header = F,
    sep = "\t", data.table = F)
  r[, 5]
}

readscAvgTESignal <- function(sc, h) {
  d <- file.path(sc2TEd, str_glue("{sc}.{h}"))
  tebg <- mean(readTEsSignal(file.path(d, "tebg.tsv")))
  teChrA <- mean(readTEsSignal(file.path(d, "teChrA.tsv")))
  teChrO <- mean(readTEsSignal(file.path(d, "teChrO.tsv")))
  teCRE <- mean(readTEsSignal(file.path(d, "teCRE.tsv")))
  data.frame(
    teCRE = teCRE,
    teChrA = teChrA,
    teChrO = teChrO,
    tebg = tebg
  )
}

readsc2Subfam <- function(h, withNaN = F, useAll = F) {
  # Note: without _ChrO, that is for ChrA
  ## suffix <- ifelse(useAll, "All",
  ##   ifelse(withNaN, "withNaN", "withoutNaN"))
  suffix <- ifelse(useAll, "All_ChrO",
    ifelse(withNaN, "withNaN", "withoutNaN"))
  d <- file.path(outd, "ovlpTE", str_glue("sc2TEsubfamily_{h}_{suffix}"))
  scs <- data.table::fread(
    file = file.path(d, "scs.txt"),
    header = F,
    data.table = F
  )$V1
  sf <- data.table::fread(
    file = file.path(d, "TEgroup.txt"),
    header = F,
    data.table = F
  )$V1
  m <- data.table::fread(
    file = file.path(d, "mat.csv"),
    header = F,
    sep = ",",
    data.table = F
  ) |> as.matrix()
  rownames(m) <- scs
  colnames(m) <- sf

  scs_ord <- vapply(scs,
    \(sc) as.double(str_split_1(sc, "_")[1]), FUN.VALUE = 0.0) |>
    order(decreasing = F)
  m <- m[scs_ord, ]
  return(m)
}

getPairedTagSubclassName <- function(rawnm) {
  gsub("  ", " ", x = rawnm) |>
    gsub(" ", "_", x = _) |>
    gsub("/", "_", x = _) |>
    gsub("-", "_", x = _)
}

getCEMBATACSubclassName <- function(rawnm) {
  gsub("^\\d+", "", x = rawnm) |>
    gsub("  ", " ", x = _) |>
    gsub(" ", "_", x = _) |>
    gsub("/", "-", x = _) |>
}

readPandasCSV <- function(f) {
  r <- data.table::fread(file = f, sep = ",",
    header = TRUE, data.table = FALSE)
  rownames(r) <- r$V1
  # remove rowname columns
  r[, -1]
}

hc4TE <- function(m, dist = "cor", cmthd = "ward.D2") {
  cordist <- c("pearson", "kendall", "spearman")
  if (any(grepl(dist, cordist))) {
    d <- cor(d, method = dist) |>
      x => 1-x |>
      as.dist()
  } else {
    d <- dist(x = m, method = dist)
  }
  hclust(d, method = cmthd) |>
    x => x$order
}

orderTEByClusterWithFam <- function(sc2te, dist = "cor", cmthd = "ward.D2") {
  fams <- vapply(colnames(sc2te), \(nm) {
    str_split_1(nm, ":")[1]
  }, "a")
  # ? order of fams
  ufams <- unique(fams)
  
  lapply(ufams, \(f) {
    index <- fams == f
    if (sum(index) < 2) {
      return(colnames(sc2te)[which(index)])
    }
    m <- sc2te[, fams == f]
    index <- hc4TE(t(m), dist, cmthd)
    colnames(m)[index]
  }) |> unlist()
}

gen_hm <- function(h, m, colQmin = 0.05, colQmax = 0.99, color = "#D85641") {
  col_fun <- colorRamp2(
    c(quantile(m, colQmin), quantile(m, colQmax)),
    c("white", color))

  p <- Heatmap(m, name = "mat", col = col_fun,
    row_names_gp = gpar(fontsize = 5),
    column_names_gp = gpar(fontsize = 3),
    cluster_rows = F,
    column_km = 9
  )

  withr::with_pdf(
    new = file.path(outd, str_glue("ptscs2TEsubfam_{h}.pdf")),
    code = draw(p),
    width = 40,
    height = 10
  )
}

genhm2 <- function(h, m, colQmin = 0.05,
                   colQmax = 0.99, color = "#D85641", topAnnot) {
  col_fun <- colorRamp2(
    c(quantile(m, colQmin), quantile(m, colQmax)),
    c("white", color))

  p <- Heatmap(m, name = "mat", col = col_fun,
    row_names_gp = gpar(fontsize = 5),
    column_names_gp = gpar(fontsize = 3),
    cluster_rows = F,
    cluster_columns =  F,
    top_annotation = topAnnot
  )

  withr::with_pdf(
    new = file.path(outd, str_glue("ptscs2TEsubfam_{h}.pdf")),
    code = draw(p),
    width = 40,
    height = 10
  )
  p
}

# * main
## ## * Read Kangli's data
## peakAnnotTE <- readRDS(file.path(resourced,
##   "peak.anno.TEs.rds"))
## peakAnnotRename <- readRDS(file.path(resourced,
##   "peak.anno.rename.rds"))
## # 1062 TE subfamliy, family, superfamily
## TEDict <- readRDS(file.path(resourced,
##   "TE_dict.rds"))
## highTEppdc <- readRDS(file.path(resourced,
##   "TE.cCREs_pos.pdc_highTE_vs_other_test.rds"))
root <- here::here()
workd <- file.path(root, "17.repressiveMarks")
figd <- file.path(workd, "fig")
outd <- file.path(workd, "out")
resourced <- file.path(workd, "src/main/resources")
hs <- c("H3K27me3", "H3K9me3")
sc2TEd <- file.path(outd, "ovlpTE/sc2TE")

fnm <- file.path(root, "meta", "mm10TE_class2family2subfamily.csv")
mm10TE_class2fam2sf <- data.table::fread(
  file = fnm, header = F, sep = ",", data.table = F) |>
  setNames(object = _, nm = c("class", "family", "subfamily")) |>
  x => `rownames<-`(x, paste(x$family, x$subfamily, sep = ":"))

TEClassOrd <- c("DNA", "LINE", "LTR", "SINE")
TEClass2Family <- mm10TE_class2fam2sf[, c("class", "family")] |>
  unique()

# order by string
classFamily <- with(TEClass2Family, paste(class, family, sep = ":"))
classFamilyOrd <- order(classFamily) |>
  setNames(object = _, nm = classFamily)

subfam2classFam <- data.frame(
  sf = with(mm10TE_class2fam2sf, paste(family, subfamily, sep = ":")),
  cf = with(mm10TE_class2fam2sf, paste(class, family, sep = ":"))
) |>
  x => `rownames<-`(x, x$sf) |>
  x => mutate(.data = x , ord = classFamilyOrd[x$cf])

# * load subclass hierarchical structure
AllenClMetafnm <- file.path(root, "meta",
  "AIT21_annotation_freeze_081523.tsv")
AllenClMeta <- data.table::fread(
  file = AllenClMetafnm, sep = "\t", header = T,
  data.table = F
) |>
  x => `rownames<-`(x, paste0("cl-", x[, 1]))

sc2class <- AllenClMeta[, c("subclass_id_label", "class_id_label")] |>
  unique()

ptscs_ <- vapply(sc2class[, 1], getPairedTagSubclassName, "a")
sc2class$subclass_id_label <- ptscs_
rownames(sc2class) <- sc2class$subclass_id_label

sc2ntype <- AllenClMeta[,
  c("subclass_id_label", "nt_type_label")] |>
  unique() |>
  mutate(.data = _, ptsc = getPairedTagSubclassName(subclass_id_label)) |>
  x => x[!duplicated(x$subclass_id_label), ] |>
  x => `rownames<-`(x, x$ptsc)

sc2ntype[is.na(sc2ntype$nt_type_label), "nt_type_label"] <- "NN"
sc2ntype[is.na(sc2ntype$nt_type_label), "Glut-GABA"] <- "GABA-Glut"

classColors <- c(
  S4Vectors::Pairs("NN", "#90D5E4"),
  S4Vectors::Pairs("Glut", "#89C75F"),
  S4Vectors::Pairs("GABA", "#F37B7D"),
  S4Vectors::Pairs("Dopa", "#D51F26"),
  S4Vectors::Pairs("Nora", "#D8A767"),
  S4Vectors::Pairs("Sero", "#C06CAB"),
  S4Vectors::Pairs("Chol", "#89288F"),
  S4Vectors::Pairs("GABA-Glut", "#8C4E3C"),
  S4Vectors::Pairs("Glut-GABA", "#8C4E3C"),
  S4Vectors::Pairs("GABA-Chol", "#89288F"),
  S4Vectors::Pairs("GABA-Dopa", "#5B7C3D"),
  S4Vectors::Pairs("GABA-Gly", "#8A9FD1"),
  S4Vectors::Pairs("GABA-Hist", "#000080"),
  S4Vectors::Pairs("Glut-Chol", "#704E89"),
  S4Vectors::Pairs("Glut-Dopa", "#FB8072")
) |> as.data.frame() |>
  setNames(object = _, nm = c("class", "color"))
myClassColor <- classColors$color |>
  setNames(object = _, nm = classColors$class)



# * load subclass region info
atacsc2region <- data.table::fread(
  file = file.path(root, "meta", "atacSubclassRegionInfo.csv"),
  sep = ",",
  header = F,
  data.table = F
) |>
  setNames(object = _, nm = c("sc", "region", "score")) |>
  x => `rownames<-`(x, x$sc)

# from CEMBA ATAC repo
largeRegionColors <- c(
  S4Vectors::Pairs("Telencephalon", "#00688B"),
  S4Vectors::Pairs("Diencephalon", "#F15F30"),
  S4Vectors::Pairs("Midbrain", "#74E44B"),
  S4Vectors::Pairs("Hindbrain", "#788FC8"),
  S4Vectors::Pairs("Cerebellum", "#DEB34C"),
  S4Vectors::Pairs("Non-Telencephalon", "#999999")
)

## * load subclass level TE repressive marks
ptscs <- list.files(sc2TEd, full.names = F, no.. = T, include.dirs = F) |>
  x => str_replace(x, ".H3K27me3|.H3K9me3|.H3K27ac|.H3K4me1", "") |>
  unique() |>
  sort()

# * load # of cells in H3K27ac
sc2hcnt <- readPandasCSV(file.path(root, "meta",
  "amb_PT_cell_count_table.csv"))
rownames(sc2hcnt) <- vapply(rownames(sc2hcnt), getPairedTagSubclassName, "a")

# * load sc2subfam data from H3K27ac
sc2sf_K27ac_All <- readsc2Subfam(h = "H3K27ac", withNaN = F, useAll = T)

# filter by # of cells in H3K27ac
ptscs <- rownames(sc2hcnt)[sc2hcnt >= 100] |>
  x => intersect(x, rownames(sc2sf_K27ac_All))
sc2sf_K27ac_All <- sc2sf_K27ac_All[ptscs, ]

# filter subfamilies with small signals
sumSigOfTE <- colSums(sc2sf_K27ac_All)
plot(sort(sumSigOfTE, decreasing = F))
# (minSigOfTE <- sort(sumSigOfTE, decreasing = F)[150])
(minSigOfTE <- quantile(sumSigOfTE, 0.05))
(maxSigOfTE <- quantile(sumSigOfTE, 0.95))

TEsfs <- which(sumSigOfTE >= minSigOfTE) |>
  x => intersect(x, which(sumSigOfTE <= maxSigOfTE))

# * given a family, perform clustering or hierarchical clustering
# on subfamilies.
# Then sort by class and family order
teord <- orderTEByClusterWithFam(sc2sf_K27ac_All[ptscs, TEsfs],
  dist = "euclidean", cmthd = "ward.D2") |>
  x => x[order(subfam2classFam[x, "ord"])]


TEClassColor <- brewer.pal(name = "Set2", n = 4) |>
  setNames(object = _, nm = TEClassOrd)
TEFamilyUniq <- mm10TE_class2fam2sf[teord, "family"] |>
  unique()

qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors,
  rownames(qual_col_pals)))
set.seed(2025)
TEFamilyColor <- sample(col_vector, size = length(TEFamilyUniq)) |>
  setNames(object = _, nm = TEFamilyUniq)
# * add teord with colors
tehmAnnot <- ComplexHeatmap::HeatmapAnnotation(
  TEClass = mm10TE_class2fam2sf[teord, "class"],
  TEFamily = mm10TE_class2fam2sf[teord, "family"],
  col = list(
    TEClass = TEClassColor,
    TEFamily = TEFamilyColor
  ),
  annotation_name_side = "left"
)

# TODO * get log-level # of TEs in each subfamily

# * perform heatmap on each family (not subfamily)
pK27ac <- genhm2("H3K27ac", sc2sf_K27ac_All[ptscs, teord],
  colQmin = 0.01, colQmax = 0.99, color = "#D85641",
  topAnnot = tehmAnnot)


# * joint heatmap with all the families

# * load sc2subfam data from other histones
sc2sf_K27me3_All <- readsc2Subfam(h = "H3K27me3", withNaN = F, useAll = T)
pK27me3 <- genhm2("H3K27me3", sc2sf_K27me3_All[ptscs, teord],
  colQmin = 0.01, colQmax = 0.99, color = "#E7AA23", topAnnot = NULL)

sc2sf_K4me1_All <- readsc2Subfam(h = "H3K4me1", withNaN = F, useAll = T)
pK4me1 <- genhm2("H3K4me1", sc2sf_K4me1_All[ptscs, teord],
  colQmin = 0.01, colQmax = 0.99, color = "#90C03E", topAnnot = NULL)

sc2sf_K9me3_All <- readsc2Subfam(h = "H3K9me3", withNaN = F, useAll = T)
pK9me3 <- genhm2("H3K9me3", sc2sf_K9me3_All[ptscs, teord],
  colQmin = 0.01, colQmax = 0.99, color = "#5551A2", topAnnot = NULL)

p <- pK27ac %v% pK4me1 %v% pK27me3 %v% pK9me3

withr::with_pdf(
  # for ChrA
  ## new = file.path(outd,
  ##   str_glue("ptscs2TEsubfam_K27ac-K4me1-K27me3-K9me3.pdf")),
  new = file.path(outd,
    str_glue("ptscs2TEsubfam_K27ac-K4me1-K27me3-K9me3_ChrO.pdf")),
  code = draw(p, heatmap_legend_side = "left"),
  width = 40,
  height = 40
)

# * get ChrA / ChrO ovlp with TE ratio per subclasses
allCREfnm <- file.path(root, "data", "CRE",
  "allCRE.amb.PairedTag.annot.tsv")
CREovlpTEfnm <- file.path(workd, "out", "ovlpTE",
  "CRE-TE-TECat.tsv")

allCREs <- data.table::fread(
  file = allCREfnm, header = T, sep = "\t", data.table = F
)

CREovlpTE <- data.table::fread(
  file = CREovlpTEfnm, header = T, sep = "\t",
  data.table = F
)
colnames(CREovlpTE) <- c("chrom", "startFrom", "endTo",
  "TEChrom", "TEStartFrom", "TEEndTo", "TEClass",
  "TESubfamilyHOMER", "TESubfamily"
  )

CREStrOvlpTE <- with(CREovlpTE,
  paste(chrom, paste(startFrom, endTo, sep = "-"), sep = ":")) |>
  unique()

atacscs <- unique(allCREs$subclass) |>
  sort()

CREovlpTEStat <- lapply(atacscs, \(sc) {
  sCREs <- allCREs[allCREs$subclass == sc, ]
  cres <- with(sCREs,
    paste(chrom, paste(startFrom, endTo, sep = "-"), sep = ":"))
  rownames(sCREs) <- cres

  nCRE <- nrow(sCREs)
  nChrACRE <- sum(sCREs$chromHMMState == "Chr-A")
  nChrOCRE <- sum(sCREs$chromHMMState == "Chr-O")
  
  sCREovlpTE <- intersect(CREStrOvlpTE, cres)
  nTECRE <- length(sCREovlpTE)
  if (nTECRE <= 0) {
    nTEChrACRE <- 0
    nTEChrOCRE <- 0
  } else {
    s <- sCREs[sCREovlpTE, , drop = F]
    nTEChrACRE <- sum(s$chromHMMState == "Chr-A")
    nTEChrOCRE <- sum(s$chromHMMState == "Chr-O")
  }
  
  data.frame(
    subclass = sc,
    nCRE = nCRE,
    nChrACRE = nChrACRE,
    nChrOCRE = nChrOCRE,
    nTECRE = nTECRE,
    nTEChrACRE = nTEChrACRE,
    nTEChrOCRE = nTEChrOCRE,
    rTECRE = nTECRE / nCRE,
    rTEChrACRE = nTEChrACRE / nChrACRE,
    rTEChrOCRE = nTEChrOCRE / nChrOCRE
  )
}) |> do.call(what = rbind, args = _) |>
  x => `rownames<-`(x, x$subclass)

# * add ntype
CREovlpTEStat$ntype <- sc2ntype[CREovlpTEStat$subclass, "nt_type_label"]

with(CREovlpTEStat, subclass[order(rTEChrACRE, decreasing = T)]) |>
  head(x = _, n = 20)

(
  histTEChrA <- ggplot(data = CREovlpTEStat, aes(x = rTEChrACRE, fill = ntype)) +
  geom_histogram(bins = 45) +
  scale_fill_manual(values = myClassColor)
)


with(CREovlpTEStat, subclass[order(rTEChrOCRE, decreasing = T)]) |>
  head(x = _, n = 20)

(
  histTEChrO <- ggplot(data = CREovlpTEStat, aes(x = rTEChrOCRE, fill = ntype)) +
  geom_histogram(bins = 45) +
  scale_fill_manual(values = myClassColor)
)


with(CREovlpTEStat, subclass[order(rTECRE, decreasing = T)]) |>
  head(x = _, n = 20)

(
  histTE <- ggplot(data = CREovlpTEStat, aes(x = rTECRE, fill = ntype)) +
  geom_histogram(bins = 45) +
  scale_fill_manual(values = myClassColor)
)

hist <- ggarrange(histTE,
  histTEChrA,
  histTEChrO, nrow = 3, ncol = 1, common.legend =T, legend = "right")
ggsave(filename = file.path(outd, "ratioOfCREOvlpWithTE.pdf"),
  plot = hist, width = 10, height = 18)

withr::with_pdf(
  new = file.path(outd, str_glue("ratioOfCREOvlpWithTE.pdf")),
  code = {hist},
  width = 10,
  height = 30
)

saveRDS(CREovlpTEStat, file = file.path(outd, "ratioOfCREOvlpWithTE.rds"))


# * highTE versus non-highTE

schighTE <- data.table::fread(file.path(root, "meta", "highTESubclass.csv"),
  header = T, sep = ",", data.table = FALSE) |>
  x => x[,2] |>
  unlist()
CREovlpTEStat <- readRDS(file.path(outd, "ratioOfCREOvlpWithTE.rds"))

CREovlpTEStat$highTEGlut <- "highTEGlut"
CREovlpTEStat$highTEGlut[!CREovlpTEStat$subclass %in% schighTE] <- "NotHighTEGlut"

highTEGlutColor <- c("red", "gray") |>
  setNames(object = _, c("highTEGlut", "NotHighTEGlut"))
(
  histTEChrA <- ggplot(data = CREovlpTEStat, aes(x = rTEChrACRE, fill = highTEGlut)) +
  geom_histogram(bins = 45, alpha = 0.5) +
  scale_fill_manual(values = highTEGlutColor)
)

(
  histTEChrO <- ggplot(data = CREovlpTEStat, aes(x = rTEChrOCRE, fill = highTEGlut)) +
  geom_histogram(bins = 45, alpha = 0.5) +
  scale_fill_manual(values = highTEGlutColor)
)

(
  histTE <- ggplot(data = CREovlpTEStat, aes(x = rTECRE, fill = highTEGlut)) +
  geom_histogram(bins = 45, alpha = 0.5) +
  scale_fill_manual(values = highTEGlutColor)
)

hist <- ggarrange(histTE,
  histTEChrA,
  histTEChrO, nrow = 3, ncol = 1, common.legend =T, legend = "right")

hist

# * reCover previous histone modifications on TEs
avgTESignals <- lapply(
  c("H3K27ac", "H3K4me1", "H3K27me3", "H3K9me3"), \(h) {
  lapply(ptscs, \(sc) {
    readscAvgTESignal(sc, h) |>
      x => mutate(.data = x, sc = sc)
  }) |>
    x => do.call(rbind, x) |>
    x => mutate(.data = x, h = h)
}) |>
  x => do.call(rbind, x)


boxplot(avgTESignals[avgTESignals$h == "H3K27ac",
  c("teCRE", "teChrA", "teChrO", "tebg")])

boxplot(avgTESignals[avgTESignals$h == "H3K27me3",
  c("teCRE", "teChrA", "teChrO", "tebg")])


boxplot(avgTESignals[avgTESignals$h == "H3K4me1",
  c("teCRE", "teChrA", "teChrO", "tebg")])

boxplot(avgTESignals[avgTESignals$h == "H3K9me3",
  c("teCRE", "teChrA", "teChrO", "tebg")])

saveRDS(avgTESignals, file = file.path(outd, "avgTESignalsOnHistones.rds"))

# t <- readRDS(file.path(outd, "avgTESignalsOnHistones.rds"))
