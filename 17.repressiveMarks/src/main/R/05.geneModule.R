library(tidyverse)
library(tmpRpkg)
library(ComplexHeatmap)
library(cluster)
library(Seurat)
library(topGO)
Sys.setenv("_R_USE_PIPEBIND_" = TRUE)

# * meta
projd <- here::here()
workd <- file.path(projd, "17.repressiveMarks")
rtsd <- file.path(workd, "out")
outd <- file.path(workd, "out")
figd <- file.path(workd, "figure")
mod <- "H3K27me3"
qtop <- 0.05
cumdyCutOff <- 0.9

# load subclass meta
ptscMeta <- tmpRpkg::loadptscMeta() |>
  x => x[x$ATAC > 0, ] |>
  x => `rownames<-`(x, x$pairedTagName)
ptscs <- ptscMeta$PairedTagName
# load Allen's Gene lists
AllenGened <- file.path(projd, "meta", "wmbAllenUpdate")
tfDE <- data.table::fread(file = file.path(AllenGened, "TF534.txt"),
                          head = F, data.table = F)$V1
tfMod <- data.table::fread(file = file.path(AllenGened,
                                            "TF261CoExpModule52.tsv"),
                           head = F, sep = "\t", data.table = F) |>
  setNames(object = _, nm = c("mod", "gene", "score")) |>
  x => `rownames<-`(x, x$gene)
gDE <- data.table::fread(file = file.path(AllenGened, "DEG8460.txt"),
                         head = F, data.table = F)$V1

# * functions
loadRTS <- function(fnm) {
  data.table::fread(
    file = fnm, sep = "\t", header = F, data.table = F
  ) |> setNames(
    object = _,
    nm = c("chrom", "startFrom", "endTo", "gene", "rts", "strand")
  )
}

calculate_wss <- function(data, k) {
  kmeans_result <- kmeans(data, centers = k, nstart = 25, iter.max = 100)
  return(kmeans_result$tot.withinss)
}

calculate_silhouette <- function(data, k) {
  kmeans_result <- kmeans(data, centers = k, nstart = 25, iter.max = 100)
  silhouette_result <- silhouette(kmeans_result$cluster, dist(data))
  return(mean(silhouette_result[, 3]))
}

#' @param gMods vector of strings
#' @return matrix, sc by mods
getAvgExpModule <- function(xscbyg, gMods) {
  umods <- unique(gMods)
  vapply(umods, \(m) {
    xscbyg[, gMods %in% c(m), drop = F] |>
      rowMeans()
  }, FUN.VALUE = rep(0.0, nrow(xscbyg))) |>
    setNames(object = _, nm = umods) |>
    x => `rownames<-`(x, rownames(xscbyg))
}

#' @param t double, threshold to assign mod to the subclass
#' @param r double \in (0, 1), proportion of subclasses
#'   global modues should cover
#' @return vector of strings, module names
getGlobalMod <- function(xscbymod, t, r = 0.7) {
  nsc <- colSums(xscbymod >= t)
  if (sum(nsc >= nrow(xscbymod) * r) >= 1) {
    colnames(xscbymod) |>
      x => x[nsc >= nrow(xscbymod) * r]
  } else {
    NULL
  }
}

#' Assign a top row id for each module based on
#' their max values.
getTopRowIdPerMod <- function(xscbymod) {
  vapply(seq_len(ncol(xscbymod)), \(m) {
    which.max(xscbymod[, m])
  }, FUN.VALUE = 1L) |>
    setNames(object = _, nm = colnames(xscbymod))
}

#' NOTE: global modules are NOT considered here.
assignRowId4Mods <- function(xscbymod) {
  mods <- colnames(xscbymod)
  rids <- getTopRowIdPerMod(xscbymod)
  ursrt <- sort(unique(rids))
  mods <- lapply(ursrt, \(i) {
    index <- rids == i
    if (sum(index) == 1) {
      mods[index]
    } else {
      m <- mods[index]
      x <- xscbymod[i, index]
      m[order(x, decreasing = T)]
    }
  }) |>
    unlist()
}

getOrderOfGene <- function(ordOfMod, geneMod) {
  r <- factor(geneMod, levels = ordOfMod)
  order(r, decreasing = F)
}


# * main
# * load RTS
rts <- loadRTS(file.path(
  rtsd,
  str_glue("gene2RTS.{mod}.top{qtop}.all.bed")
)) |>
  x => x[!duplicated(x$gene), ] |>
  x => `rownames<-`(x, x$gene)
dy <- abs(diff(rts$rts))
cumdy <- cumsum(dy)
ng <- min(which(cumdy >= cumdyCutOff))
minRTS <- rts$rts[ng]
topRTS <- rts[rts$rts >= minRTS, ]
data.table::fwrite(
  x = topRTS, file = file.path(workd, "out", "H3K27me3.topRTS.gene.bed"),
  sep = "\t", col.names = F, row.names = F
)

cpmscbyg <- readRDS(file.path(
  projd, "data", "pairedtag_ann",
  "pt.RNAseq.pseudobulk.CPM.scbyg.rds"
)) |>
  x => x[ptscs, ]

# * load gene expression
logcpmscbyg <- readRDS(file.path(
  projd, "data", "pairedtag_ann",
  "pt.RNAseq.pseudobulk.CPM.scbyg.rds"
)) |>
  x => log1p(x)

geneTopRTS <- intersect(topRTS$gene, colnames(logcpmscbyg))
X <- logcpmscbyg[ptscs, geneTopRTS] |>
  scale(x = _, center = T, scale = T)

# * PCA
# This PCA is to reduce dim for genes (based on their cells features.)
xPCA <- RunPCA(X, npcs = 50, weight.by.var = T, assay = "RNA")

pct <- xPCA@stdev / sum(xPCA@stdev) * 100
cmu <- cumsum(pct)
plot(x = 1:length(cmu), cmu)
plot(x = 1:length(xPCA@stdev), xPCA@stdev)
co1 <- which(cmu > 90 & pct < 5)[1]
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1),
  decreasing = T
)[1] + 1
pcs <- max(co1, co2) # 42

# * k-means to get gene modules
kmeansBIC <- function(fit) {
  m <- ncol(fit$centers) # Number of features
  n <- length(fit$cluster) # Number of data points
  k <- nrow(fit$centers) # Number of clusters
  D <- fit$tot.withinss # Within-cluster sum of squares
  return(D + log(n) * m * k) # BIC calculation
}

k_values <- 10:80
x4km <- xPCA@cell.embeddings[, 1:pcs]
BICs <- sapply(k_values, \(k) {
  r <- kmeans(x4km,
    centers = k,
    nstart = 25, iter.max = 100
  )
  kmeansBIC(r)
})

plot(k_values, BICs,
  type = "b",
  xlab = "Number of Clusters (k)",
  ylab = "BIC",
  main = "BIC"
)

# * choose k based on BIC
k <- 34
kmeans_cols <- kmeans(x4km, centers = k, iter.max = 100, nstart = 25)
# kmeans_cols to data.frame
kmeansCols <- data.frame(
  cluster = paste0("g", kmeans_cols$cluster),
  gene = names(kmeans_cols$cluster)
) |> x => `rownames<-`(x, x$gene)
data.table::fwrite(
  x = kmeansCols,
  file = file.path(outd, str_glue("H3K27me3.RTS.kmeans.{k}.geneModule.csv")),
  sep = ",", row.names = F, col.names = T
)

# * rank modules
avgExpMod <- getAvgExpModule(X, kmeansCols$cluster)
rowIds <- assignRowId4Mods(avgExpMod)
columnOrd <- getOrderOfGene(
  ordOfMod = rowIds,
  geneMod = kmeansCols$cluster
)

Xord <- X[, columnOrd]
saveRDS(object = Xord, 
        file = file.path(outd, str_glue("{mod}.RTS.scbygene.ord.rds")))
# save ordered gene and their modules
g2modOrd <- data.frame(
  gene = colnames(Xord),
  module = kmeansCols[colnames(Xord), "cluster"]
)
data.table::fwrite(
  x = g2modOrd,
  file = file.path(outd, str_glue("{mod}.RTS.geneModule.ord.csv")),
  sep = ",", row.names = F, col.names = F
)

# * draw Heatmap
low.val.col <- quantile(Xord, -1.06, na.rm = T)
high.val.col <- quantile(Xord, 0.96, na.rm = T)
legend_labels <- c(
  round(low.val.col, 1),
  round(high.val.col, 1)
)
colors <- ggsci::pal_igv()(k) |>
  setNames(object = _, nm = paste("g", seq_len(k), sep = ""))
ha <- ComplexHeatmap::HeatmapAnnotation(
  group = kmeansCols[colnames(Xord), "cluster"],
  col = list(group = colors)
)
hmExp <- ComplexHeatmap::Heatmap(
  mat = Xord,
  col = circlize::colorRamp2(
    seq(low.val.col, high.val.col, length = 50),
    viridis::viridis(50)
  ),
  show_row_names = F,
  show_column_names = F,
  cluster_columns = F,
  cluster_rows = F,
  use_raster = T,
  row_names_gp = grid::gpar(fontsize = 6),
  column_names_gp = grid::gpar(fontsize = 4),
  heatmap_legend_param = list(
    title = "scaled logCPM",
    at = c(low.val.col, high.val.col),
    labels = legend_labels,
    direction = "vertical"
  ),
  top_annotation = ha
)
hmExp

hmExp <- ComplexHeatmap::Heatmap(
  mat = Xord,
  col = circlize::colorRamp2(
    seq(low.val.col, high.val.col, length = 50),
    viridis::viridis(50)
  ),
  show_row_names = T,
  show_column_names = T,
  cluster_columns = F,
  cluster_rows = F,
  use_raster = T,
  row_names_gp = grid::gpar(fontsize = 6),
  column_names_gp = grid::gpar(fontsize = 6),
  heatmap_legend_param = list(
    title = "scaled logCPM",
    at = c(low.val.col, high.val.col),
    labels = legend_labels,
    direction = "vertical"
  ),
  top_annotation = ha
)
hmExp

# * what if we focus on TFs
xTFDEOrd <- Xord[, colnames(Xord) %in% tfDE]
low.val.col <- quantile(xTFDEOrd, 0.06, na.rm = T)
high.val.col <- quantile(xTFDEOrd, 0.96, na.rm = T)
legend_labels <- c(
  round(low.val.col, 1),
  round(high.val.col, 1)
)
colors <- ggsci::pal_igv()(k) |>
  setNames(object = _, nm = paste("g", seq_len(k), sep = ""))
ha <- ComplexHeatmap::HeatmapAnnotation(
  group = kmeansCols[colnames(xTFDEOrd), "cluster"],
  col = list(group = colors)
)
hmTFDE <- ComplexHeatmap::Heatmap(
  mat = xTFDEOrd,
  col = circlize::colorRamp2(
    seq(low.val.col, high.val.col, length = 50),
    viridis::viridis(50)
  ),
  show_row_names = T,
  show_column_names = T,
  cluster_columns = F,
  cluster_rows = F,
  use_raster = T,
  row_names_gp = grid::gpar(fontsize = 6),
  column_names_gp = grid::gpar(fontsize = 6),
  heatmap_legend_param = list(
    title = "scaled logCPM",
    at = c(low.val.col, high.val.col),
    labels = legend_labels,
    direction = "vertical"
  ),
  top_annotation = ha
)
hmTFDE
# * how we label the gene markers
# * how we perform GO term analysis
data(geneList)
g1 <- kmeansCols$gene[kmeansCols$cluster == "g1"]
allGenes <- intersect(rts$gene, colnames(logcpmscbyg)) |>
  x => rts[x, "rts"]

sampleGOdata <- new(
  "topGOdata",
  description = "Sample GO data",
  ontology = "BP",
  allGenes = allGenes,
  geneSel = g1,
  annot = annFUN.org,
  mapping = "org.Mm.eg.db",
  ID = "symbol",
  nodeSize = 5
)
