# check cell-type specific ChrAs' target genes
# if they are keep cell type specificity

.libPaths(
  c("/opt/homebrew/Caskroom/miniforge/base/lib/R/library",
    "/Users/szu/Library/R/arm64/4.5/library"))
library(tidyverse)
library(reticulate)
library(ComplexHeatmap)

reticulate::use_python(
  python = "/opt/homebrew/Caskroom/miniforge/base/bin/python")
pd <- reticulate::import(module = "pandas", convert = F)
np <- reticulate::import(module = "numpy")

getATACSubclassName <- function(rawnm) {
  gsub("^\\d+_", "", x = rawnm) |>
    gsub("  ", " ", x = _) |>
    gsub(" ", "_", x = _) |>
    gsub("/", "-", x = _)
}

genSimpleHeatmap <- function(x,
                             lowq = 0.005,
                             highq = 0.995,
                             title = "logCPM") {
  lowval <- quantile(x, probs =  lowq)
  highval <- quantile(x, probs = highq)
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
projd <- here::here()
workd <- file.path(projd, "16.celloracle")
outd <- file.path(workd, "out")

# load subclass name meta
ptscMeta <- data.table::fread(
  file = file.path(projd, "meta",
    "PairedTagSubclassMetaFromCellMetaChromHMM.csv"),
  sep = ",",
  header = T,
  data.table = F) |>
  x => x[x$ATAC > 0, ] |>
  x => `rownames<-`(x, x$PairedTagName)

ptscMeta$ATACName2 <- gsub("^\\d+_", "", ptscMeta$ATACName)
atac2ptscName <- data.frame(
  ATAC = ptscMeta$ATACName2,
  PairedTag = ptscMeta$PairedTagName
) |>
  x => `rownames<-`(x, x$ATAC)

# load subclass DAR ChrA CREs
ptCREs <- data.table::fread(
  file.path(projd, "data", "CRE",
    "allCRE.amb.PairedTag.annot.tsv"),
  header = T, sep = "\t", data.table = F)

DARChrAs <- ptCREs[with(ptCREs,
  (chromHMMState == "Chr-A") & (isDistal == "distal") &
    (isDAR == "DAR")), ]

# load pseudobulk-level gene expression 
sc2gCPM <- readRDS(file.path(projd, "data",
  "amb_PT_RNA_pseudobulk",
  "pt.RNAseq.pseudobulk.CPM.scbyg.rds"))


# load previous proximal distal interactions
sc2pdc <- lapply(ptscMeta$PairedTagName, \(sc) {
  atacsc <- ptscMeta[sc, "ATACName2"]
  data.table::fread(
    file = file.path(projd, "data/snATAC/sa2pdc_bedpe",
      str_glue("sa2subclass.{atacsc}.pdc.bedpe")),
    sep = "\t",
    data.table = F
  ) |>
    mutate(.data = _,
      CRE = paste(V4, paste(V5, V6, sep = "-"),
        sep = ":"),
      gene = gsub("\\|.+$", "", V7),
      pearson = V8) |>
    x => x[, c("CRE", "gene", "pearson")]
}) |>
  setNames(object = _, nm = ptscMeta$PairedTagName)


# load snATAC pseuduobulk expressions
sc2ATAC <- pd$read_csv(file.path(projd, "data/snATAC",
  "cpm_peakBysubclass.csv"), index_col = 0L, sep = ",") |>
  reticulate::py_to_r(x = _)
scs <- sc2ATAC$columns$to_numpy()
cres <- sc2ATAC$index$to_numpy()
atacCPM <- sc2ATAC$to_numpy()
colnames(atacCPM) <- scs
rownames(atacCPM) <- cres

newatacCPM <- atacCPM[ , atac2ptscName[
  colnames(atacCPM), "PairedTag"] %in% ptscMeta$PairedTagName] |>
  x => `colnames<-`(x, atac2ptscName[
    colnames(x), "PairedTag"])

# * get subclass: CRE-gene 
c2g <- lapply(ptscMeta$PairedTagName, \(ptsc) {
  m <- sc2pdc[[ptsc]]
  darCREs <- DARChrAs[DARChrAs$subclass == ptsc, 1:3] |>
    x => with(x, paste(chrom,
      paste(startFrom, endTo, sep = "-"), sep = ":"))|>
    x => intersect(x, m$CRE) |>
    x => x[order(as.vector(newatacCPM[x, ptsc]), decreasing = T)]
  m[m$CRE %in% darCREs, ] |>
    x => x[x$gene %in% colnames(sc2gCPM), ] |>
    x => mutate(.data = x,
      gCPM = sc2gCPM[ptsc, x$gene]) |>
    group_by(CRE) |>
    slice_max(pearson) |>
    as.data.frame()
}) |> setNames(object = _, nm = ptscMeta$PairedTagName)

# * prepare heatmap
c2g_ds <- lapply(c2g, \(x) {
  x[1:min(100, nrow(x)), ]
}) |> do.call(what = rbind, args = _) |>
  x => x[!duplicated(x$CRE), ]


ataclogCPM_CRE2sc <- newatacCPM[
  c2g_ds$CRE, ptscMeta$PairedTagName] |>
  x => log1p(x)

# for ATAC heatmap
# load(".RData")
lowval <- quantile(ataclogCPM_CRE2sc, 0.05)
highval <- quantile(ataclogCPM_CRE2sc, 0.99)
atacColor <- circlize::colorRamp2(
  seq(lowval, highval, length = 60),
  viridis::viridis(60)
)
legend <- list(
  title = "ATAC logCPM",
  at = c(lowval, highval),
  labels =  round(c(lowval, highval), 2),
  direction = "horizontal"
)

(
  hmATAC <- Heatmap(
    matrix = ataclogCPM_CRE2sc,
    col = atacColor,
    cluster_rows = F,
    cluster_columns = F,
    show_row_names = F,
    show_column_names = F,
    use_raster = T,
    show_heatmap_legend = T,
    heatmap_legend_param = legend,
    width = 1
  )
)

# for gene heatmap
g2sc_zscore <- log1p(sc2gCPM[
  ptscMeta$PairedTagName, c2g_ds$gene]) |>
  scale() |>
  t()

lowval <- quantile(g2sc_zscore, 0.1)
highval <- quantile(g2sc_zscore, 0.9)
gColor <- circlize::colorRamp2(
  seq(lowval, highval, length = 60),
  viridis::viridis(60)
)
legend <- list(
  title = "Gene zscore",
  at = c(lowval, highval),
  labels =  round(c(lowval, highval), 2),
  direction = "horizontal"
)

(
  hmGene <- Heatmap(
    matrix = g2sc_zscore,
    col = gColor,
    cluster_rows = F,
    cluster_columns = F,
    show_row_names = F,
    show_column_names = F,
    use_raster = T,
    show_heatmap_legend = T,
    heatmap_legend_param = legend,
    width = 1
  )
)

# * get differential genes for heatmap
markerd <- file.path(projd, "16.celloracle", "out",
  "subclass_marker_genes_csvs")
markerscs <- list.files(markerd, full.names = F) |>
  x => gsub("_marker_genes\\.csv", "", x)
allDEGenes <- lapply(markerscs, \(sc) {
  data.table::fread(
    file = file.path(markerd, str_glue("{sc}_marker_genes.csv")),
    header = T, sep = ",", data.table = F
  )
}) |>
  x => setNames(object = x, nm = markerscs)
# filter genes to get subclasses' markers
padj <- 0.01
log2fc <- 1.0

markers <- lapply(markerscs, \(sc) {
  y <- allDEGenes[[sc]] |>
    x => with(x, gene[p_val_adj <= padj & avg_log2FC >= log2fc])
  if (length(y) < 1) {
    NULL
  } else {
    y
  }
}) |>
  x => setNames(object = x, nm = gsub("-", "_", markerscs)) |>
  x => Filter(Negate(is.null), x)

# overlap DEGene with pdc
c2g <- lapply(ptscMeta$PairedTagName, \(ptsc) {
  m <- sc2pdc[[ptsc]]
  darCREs <- DARChrAs[DARChrAs$subclass == ptsc, 1:3] |>
    x => with(x, paste(chrom,
      paste(startFrom, endTo, sep = "-"), sep = ":"))|>
    x => intersect(x, m$CRE) |>
    x => x[order(as.vector(newatacCPM[x, ptsc]), decreasing = T)]
  m[m$CRE %in% darCREs, ] |>
    x => x[x$gene %in% colnames(sc2gCPM), ] |>
    x => x[x$gene %in% markers[[ptsc]], ] |>
    x => mutate(.data = x,
      gCPM = sc2gCPM[ptsc, x$gene]) |>
    group_by(CRE) |>
    slice_max(pearson) |>
    as.data.frame()
}) |> setNames(object = _, nm = ptscMeta$PairedTagName)

c2g2 <- Filter(\(x) nrow(x) > 0, c2g)

# * prepare heatmap
c2g_ds <- lapply(c2g2, \(x) {
  x[1:min(100, nrow(x)), ]
}) |> do.call(what = rbind, args = _) |>
  x => x[!duplicated(x$CRE), ]

ataclogCPM_CRE2sc <- newatacCPM[
  c2g_ds$CRE, names(c2g2)] |>
  x => log1p(x)

hmATAC <- genSimpleHeatmap(ataclogCPM_CRE2sc, 0.05, 0.99)
g2sc_zscore <- log1p(sc2gCPM[
  names(c2g2), c2g_ds$gene]) |>
  scale() |>
  t()
hmGene <- genSimpleHeatmap(g2sc_zscore, 0.1, 0.9)
## legend <- list(
##   title = "Gene zscore",
##   at = c(lowval, highval),
##   labels =  round(c(lowval, highval), 2),
##   direction = "horizontal"
## )
hmATAC + hmGene

saveRDS(
  object = list(
    atachm = ataclogCPM_CRE2sc,
    genehm = g2sc_zscore),
  file = file.path(projd,
    "13.cicero", "out", "DARChrACRE_DEGene.heatmap.rds"))

# now save the CRE-gene-subclass
c2gdf <- lapply(names(c2g2), \(ptsc) {
  c2g2[[ptsc]] |>
    mutate(ptsc = ptsc)
}) |> do.call(what = rbind, args = _)

c2gdsdf <- lapply(names(c2g2), \(ptsc) {
  c2g2[[ptsc]] |>
    mutate(ptsc = ptsc) |>
    x => x[1:min(100, nrow(c2g2[[ptsc]])), ]
}) |>
  x => do.call(what = rbind, args = x) |>
  x => x[!duplicated(x$CRE), ]

data.table::fwrite(
  x = c2gdf, file = file.path(projd,
    "13.cicero", "out",
    "ptscs_DARChrACRE-maxppdc-DEGene.heatmap.csv"),
  sep = ",",
  row.names = F, col.names = T
)

data.table::fwrite(
  x = c2gdsdf, file = file.path(projd,
    "13.cicero", "out",
    "ptscs_DARChrACRE-maxppdc-DEGene_ds100.heatmap.csv"),
  sep = ",",
  row.names = F, col.names = T
)

# * 20260131: generate the gene expressions for Fig4d.

sc2pdc <- lapply(ptscMeta$PairedTagName, \(sc) {
  atacsc <- ptscMeta[sc, "ATACName2"]
  data.table::fread(
    file = file.path(projd, "data/snATAC/sa2pdc_bedpe",
      str_glue("sa2subclass.{atacsc}.pdc.bedpe")),
    sep = "\t",
    data.table = F
  ) |>
    mutate(.data = _,
      CRE = paste(V4, paste(V5, V6, sep = "-"),
        sep = ":"),
      gene = gsub("\\|.+$", "", V7),
      pearson = V8) |>
    x => x[, c("CRE", "gene", "pearson")]
}) |>
  setNames(object = _, nm = ptscMeta$PairedTagName)

a <- readRDS(file.path(projd, "13.cicero",
  "out", "Fig4.PanelD.ds1000.rds"))

# overlap CREs in Fig4d with pdc and markers
nowCREs <- rownames(a$ATAC)
names(allDEGenes) <- gsub("-", "_", names(allDEGenes))
padj <- 0.01
log2fc <- 0.25

markers <- lapply(names(allDEGenes), \(sc) {
  y <- allDEGenes[[sc]] |>
    x => with(x, gene[p_val_adj <= padj & avg_log2FC >= log2fc])
  if (length(y) < 1) {
    NULL
  } else {
    y
  }
}) |>
  x => setNames(object = x, nm = gsub("-", "_", markerscs)) |>
  x => Filter(Negate(is.null), x)


c2g <- lapply(colnames(a$K27ac), \(ptsc) {
  m <- sc2pdc[[ptsc]]
  darCREs <- nowCREs[nowCREs %in% m$CRE]
  
  y <- allDEGenes[[ptsc]] |>
    x => x[x$gene %in% markers[[ptsc]], ] |>
    x => `rownames<-`(x, x$gene)

  m[m$CRE %in% darCREs, ] |>
    x => x[x$gene %in% colnames(sc2gCPM), ] |>
    x => x[x$gene %in% markers[[ptsc]], ] |>
    # x => x[x$gene %in% y$gene, ] |>
    x => mutate(.data = x,
      gCPM = sc2gCPM[ptsc, x$gene],
      logfc = y[x$gene, "avg_log2FC"]
      ) |>
    group_by(CRE) |>
    #slice_max(pearson) |>
    slice_max(logfc) |>
    as.data.frame()
}) |> setNames(object = _, nm = colnames(a$K27ac)) |>
  Filter(\(x) nrow(x) >0, x = _)

c2g_ds <- do.call(rbind, c2g) |>
  x => x[!duplicated(x$CRE), ] |>
  x => `rownames<-`(x, x$CRE)
rCREs <- c2g_ds$CRE |>
  x => rownames(a$ATAC)[rownames(a$ATAC) %in% x]
rptscs <- names(c2g)
c2g_ds <- c2g_ds[rCREs, ]



g2sc_zscore <- log1p(sc2gCPM[
  names(c2g), c2g_ds$gene]) |>
  scale() |>
  t()

hmGene <- genSimpleHeatmap(
  g2sc_zscore, lowq = 0.25, highq = 0.9)

xATAC <- a$ATAC
colnames(xATAC) <- colnames(a$K27ac)
xATAC <- xATAC[rCREs, rptscs]

hmATAC <- genSimpleHeatmap(xATAC, lowq = 0.005, highq = 0.995)
hmK27ac <- genSimpleHeatmap(
  a$K27ac[rCREs, rptscs], lowq = 0.05, highq = 0.99)
hmK4me1 <- genSimpleHeatmap(a$K4me1[rCREs, rptscs], lowq = 0.005,
  highq = 0.995)
hmeRNA <- genSimpleHeatmap(a$eRNA[rCREs, rptscs],
  lowq = 0.3, highq = 0.9)
(hmsum <- hmATAC + hmK27ac + hmK4me1 + hmeRNA + hmGene)



