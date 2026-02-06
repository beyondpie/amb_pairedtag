suppressMessages({
  suppressWarnings({
    library(cicero)
    library(tidyverse)
    library(reticulate)
    library(future.apply)
    library(Matrix)
    library(GenomicRanges)
    library(rtracklayer)
    library(ComplexHeatmap)
    library(tmpRpkg)
  })
})
# set _R_USE_PIPEBIND as true if you want to use pipebind
Sys.setenv("_R_USE_PIPEBIND_" = TRUE)

# * meta
projd <- here::here()
workd <- file.path(projd, "13.cicero")
outd <- file.path(workd, "out")
outfigd <- file.path(workd, "figure")
supenhd <- file.path(
  projd, "10.superEnhancer", "out", "ROSE_H3K27ac_bed"
)
H3K27acDEd <- file.path(
  projd, "12.DE", "out", "H3K27ac_sa2LRT_DE"
)
H3K27acExprds <- file.path(
  projd, "data", "pairedtag_ann",
  "pt.H3K27ac.pseudobulk.count.sc.bestSPM.rds"
)
RNAexprds <- file.path(
  projd, "data", "pairedtag_ann",
  "pt.RNAseq.pseudobulk.count.sc.allgenesymbol.rds"
)
H3K27acDEPeakd <- file.path(outd, "H3K27acDEPeak")
H3K27acDECiceroPeak <- file.path(outd, "H3K27acDECiceroPeak")

# * set ggplot themes
theme_my_minimal <- theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.position = "right",
    axis.title = element_text(colour = "black", family = "serif"),
    axis.text = element_text(colour = "black", family = "serif"),
    plot.title = element_text(
      colour = "black", family = "serif",
      hjust = 0.5
    )
  )

# * load cicero correlations
cicero_cor <- data.table::fread(
  file = file.path(outd, "all.H3K27ac.pdc.cor.csv"),
  sep = ",", header = F, data.table = F
)
cicero_shufcor <- data.table::fread(
  file = file.path(outd, "shuf.H3K27ac.pdc.cor.csv"),
  sep = ",", header = F, data.table = F
)
thres_cor <- with(
  cicero_shufcor,
  mean(V3) + 1.0 * sd(V3)
)
cicero_filtercor <- cicero_cor[cicero_cor$V3 >= thres_cor, ]

# * load super enhancers
supenhfnms <- list.files(supenhd, include.dirs = F, no.. = T) |>
  y => y[grepl("^\\d+", y, perl = T)]
index <- vapply(supenhfnms, \(s) {
  as.integer(str_split_1(s, "_")[1])
}, 1)
supenhfnms <- supenhfnms[index <= 400]
groups <- lapply(supenhfnms, \(s) {
  str_split_1(s, "\\.")[1]
}) |> unlist()

supenhdf_List <- lapply(supenhfnms, \(s) {
  data.table::fread(
    file = file.path(supenhd, s), sep = "\t", header = F,
    data.table = F, col.names = c("chr", "start", "end")
  )
}) |>
  setNames(object = _, nm = groups) |>
  filterNULLfromList()

# * load subclass-specifc enhancers
K27acDEfnms <- list.files(H3K27acDEd, include.dirs = F, no.. = T) |>
  y => y[grepl("^\\d+", y, perl = T)]

groups <- gsub("_sc_diffpeaks.txt", "", K27acDEfnms)

K27acDEdf_List <- lapply(K27acDEfnms, \(s) {
  data.table::fread(
    file = file.path(H3K27acDEd, s), sep = ",", header = T,
    data.table = F
  )
}) |> setNames(object = _, nm = groups)

qval <- 0.05
enhdf_List <- lapply(K27acDEdf_List, \(d) {
  index <- d[, 4] <= qval
  if (sum(index) < 1) {
    NULL
  } else {
    d[index, ]
  }
}) |> filterNULLfromList()

# * check the overlap with subclass-specific enhancers
ovlpEnhdfList <- lapply(enhdf_List, \(d) {
  index <- d[, 1] %in% cicero_filtercor[, 1]
  if (sum(index) < 1) {
    NULL
  } else {
    d[index, ]
  }
}) |> filterNULLfromList()

# 49325 out of 60346 enhancers showed cell-type specific pattern
# * TODO check the overlap with superEnhancer

# * plot heatmap
## ** organize order of enhancers
## ovlpEnhdfList is sorted by default based on the names
ord_groups <- names(ovlpEnhdfList)
ord_enhs <- lapply(ord_groups, \(og) {
  ovlpEnhdfList[[og]][, 1]
}) |>
  unlist() |>
  unique()
## subclass by enhancers
raw_enh_count <- readRDS(H3K27acExprds)
raw_enh_count <- as.matrix(raw_enh_count)
raw_enh_cpm <- (raw_enh_count * 10^6) / rowSums(raw_enh_count)
enhlogcpm <- raw_enh_cpm[ord_groups, ord_enhs] |>
  as.matrix() |>
  log1p()

## ** organize the genes based on the enhancers
enh2gene_cor <- cicero_filtercor[
  cicero_filtercor$V1 %in% ord_enhs,
]
enh2gene <- enh2gene_cor |>
  group_by(.data = _, V1) |>
  summarize(.data = _, topgene = V2[which.max(V3)])
rownames(enh2gene) <- enh2gene$V1

## ordeded by ord_enhs
ord_gene <- with(enh2gene, topgene[match(ord_enhs, V1)])
## read gene expression
raw_gene_count <- readRDS(RNAexprds)
raw_gene_count <- as.matrix(raw_gene_count)
raw_gene_cpm <- (raw_gene_count * 10^6) / rowSums(raw_gene_count)
genelogcpm <- raw_gene_cpm[ord_groups, ord_gene] |>
  log1p()

## * save matrix for heatmaps
saveRDS(enhlogcpm,
  file = file.path(outd, "ord.H3K27ac.enh.logcpm.rds")
)
saveRDS(genelogcpm, file = file.path(outd, "ord.gene.logcpm.rds"))

## draw heatmap
enhlogcpm <- readRDS(file.path(outd, "ord.H3K27ac.enh.logcpm.rds"))
genelogcpm <- readRDS(file.path(outd, "ord.gene.logcpm.rds"))

## ** draw heatmap of enhlogcpm
lowval <- quantile(enhlogcpm, probs = 0.005)
highval <- quantile(enhlogcpm, probs = 0.995)
enhColor <- circlize::colorRamp2(
  seq(lowval, highval, length = 60),
  viridis::viridis(60)
)
legendEnh <- list(
  title = "logCPM (H3K27ac)",
  at = c(lowval, highval),
  labels = round(c(lowval, highval), 2),
  direction = "horizontal"
)

hmEnh <- ComplexHeatmap::Heatmap(
  matrix = t(enhlogcpm),
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

## ** draw heatmap of genelogcpm
genelogcpm <- readRDS(file.path(outd, "ord.gene.logcpm.rds"))
geneScaledLogCPM <- scale(genelogcpm)
genelogcpm <- geneScaledLogCPM
lowval <- quantile(genelogcpm, probs = 0.2)
highval <- quantile(genelogcpm, probs = 0.995)
geneColor <- circlize::colorRamp2(
  seq(lowval, highval, length = 60),
  viridis::magma(60)
)
legendGene <- list(
  title = "scaled logCPM (RNA)",
  at = c(lowval, highval),
  labels = round(c(lowval, highval), 2),
  direction = "horizontal"
)
hmRNA <- ComplexHeatmap::Heatmap(
  matrix = t(genelogcpm),
  col = geneColor,
  cluster_rows = F,
  cluster_columns = F,
  show_row_names = F,
  show_column_names = F,
  use_raster = T,
  show_heatmap_legend = T,
  heatmap_legend_param = legendGene,
  width = 1
)

## ** merge hmEnh and hmRNA
hm_Enh_RNA <- ComplexHeatmap::draw(
  hmEnh + hmRNA,
  heatmap_legend_side = "bottom",
  merge_legend = T
)
## save the heatmap to outfigd
withr::with_pdf(file.path(outfigd, "H3K27ac_enh_gene_heatmap.pdf"), {
  print(hm_Enh_RNA)
})

# * output DE and DE-cicero peaks
## * DE peaks
H3K27acDEPeaks <- lapply(K27acDEdf_List, \(t) {
  strsplit(t[, 1], split = ":|-", perl = T) |>
    do.call(what = rbind, args = _) |>
    as.data.frame() |>
    setNames(object = _, nm = c("chr", "startFrom", "endTo")) |>
    mutate(name = t[, 1])
}) |>
  setNames(object = _, nm = names(K27acDEdf_List))

invisible(lapply(names(H3K27acDEPeaks), \(sc) {
  t <- H3K27acDEPeaks[[sc]]
  data.table::fwrite(
    x = t,
    file = file.path(
      H3K27acDEPeakd,
      str_glue("{sc}.K27ac.DE.q0.05.bed")
    ),
    col.names = F, row.names = F, sep = "\t"
  )
}))

## * create the background for DE peaks
allDEPeaks <- do.call(rbind, H3K27acDEPeaks)
allDEPeaks <- allDEPeaks[!duplicated(allDEPeaks[, 4]), ]
H3K27acDEbgPeaks <- lapply(H3K27acDEPeaks, \(t) {
  allDEPeaks[!(allDEPeaks$name %in% t$name), ]
})
invisible(lapply(names(H3K27acDEPeaks), \(sc) {
  t <- H3K27acDEbgPeaks[[sc]]
  data.table::fwrite(
    x = t,
    file = file.path(H3K27acDEPeakd, str_glue("{sc}.K27ac.DE.bg.bed")),
    col.names = F, row.names = F, sep = "\t"
  )
}))

## * DE cicero peaks
ovlpEnhs <- lapply(ovlpEnhdfList, \(t) {
  strsplit(t[, 1], split = ":|-", perl = T) |>
    do.call(what = rbind, args = _) |>
    as.data.frame() |>
    setNames(object = _, nm = c("chr", "startFrom", "endTo")) |>
    mutate(name = t[, 1])
}) |>
  setNames(object = _, nm = names(ovlpEnhdfList))

invisible(lapply(names(ovlpEnhs), \(sc) {
  t <- ovlpEnhs[[sc]]
  data.table::fwrite(
    x = t,
    file = file.path(
      H3K27acDECiceroPeak,
      str_glue("{sc}.K27ac.DE.cicero.bed")
    ),
    col.names = F, row.names = F, sep = "\t"
  )
}))

# * Run homer2 outside of this script
# * Summarize Homer2 results
## prepare motif matrix
field <- "logp"
tf <- \(i) {
  -i * log10(exp(1))
}

# ** load H3K27ac DE bgGenome
homerd_K27acDEbgG <- file.path(
  workd, "out",
  "homer2_H3K27ac_DE_bgGenome"
)
allsubd <- list.files(homerd_K27acDEbgG,
  full.names = F, all.files = F, recursive = F,
  no.. = T, include.dirs = F
)
scs <- vapply(allsubd, \(i) str_split_1(i, "\\.")[1], "001")

motif_K27acDEbgG <- lapply(allsubd, \(i) {
  loadHomer2KnownResults(file.path(homerd_K27acDEbgG, i), rmdupmotif = T)
}) |>
  setNames(object = _, nm = scs)

motifmat_K27acDEbgG <- getMotifMat(
  motif_K27acDEbgG,
  field = field, tf = tf)
saveRDS(
  motifmat_K27acDEbgG,
  file.path(
    outd, "MotifMat_H3K27acDEbgGenome.rds"
  )
)

# ** load H3K27ac DE bgOthers
homerd_K27acDEbgO <- file.path(
  workd, "out",
  "homer2_H3K27ac_DE_bgOthers"
)
allsubd <- list.files(homerd_K27acDEbgO, full.names = F, no.. = T)
scs <- vapply(allsubd, \(i) str_split_1(i, "\\.")[1], "001")
motif_K27acDEbgO <- lapply(allsubd, \(i) {
  loadHomer2KnownResults(file.path(homerd_K27acDEbgO, i), rmdupmotif = T)
}) |>
  setNames(object = _, nm = scs)

## prepare motif matrix
motifmat_K27acDEbgO <- getMotifMat(
  motif_K27acDEbgO,
  field = field, tf = tf)
saveRDS(
  motifmat_K27acDEbgO,
  file.path(
    outd, "MotifMat_H3K27acDEbgOthers.rds"
  )
)

# ** load H3K27ac DE Cicero bgGenome
homerd_K27acDECicerobyG <- file.path(
  workd, "out",
  "homer2_H3K27ac_DEandCicero_bgGenome"
)
allsubd <- list.files(homerd_K27acDECicerobyG, full.names = F, no.. = T)
scs <- vapply(allsubd, \(i) str_split_1(i, "\\.")[1], "001")
motif_K27acDECicerobgG <- lapply(allsubd, \(i) {
  loadHomer2KnownResults(file.path(homerd_K27acDECicerobyG, i), rmdupmotif = T)
}) |>
  setNames(object = _, nm = scs)

## prepare motif matrix
motifmat_K27acDECicerobgG <- getMotifMat(
  motif_K27acDECicerobgG,
  field = field, tf = tf)

saveRDS(
  motifmat_K27acDECicerobgG,
  file.path(
    outd, "MotifMat_H3K27acDECicerobgGenome.rds"
  )
)

# ** load H3K27ac DE Cicero bgOthers
homerd_K27acDECicerobyO <- file.path(
  workd, "out",
  "homer2_H3K27ac_DEandCicero_bgOthers"
)
allsubd <- list.files(homerd_K27acDECicerobyO, full.names = F, no.. = T)
scs <- vapply(allsubd, \(i) str_split_1(i, "\\.")[1], "001")
motif_K27acDECicerobgO <- lapply(allsubd, \(i) {
  loadHomer2KnownResults(file.path(homerd_K27acDECicerobyO, i), rmdupmotif = T)
}) |>
  setNames(object = _, nm = scs)

## prepare motif matrix
motifmat_K27acDECicerobgO <- getMotifMat(
  motif_K27acDECicerobgO,
  field = field, tf = tf)

saveRDS(
  motifmat_K27acDECicerobgO,
  file.path(
    outd, "MotifMat_H3K27acDECicerobgOthers.rds"
  )
)

# ** load H3K27ac super enhancer bg Genome
homerd_K27acSupEnhbyG <- file.path(
  workd, "out",
  "homer2_H3K27ac_SupEnh_bgGenome"
)
allsubd <- list.files(homerd_K27acSupEnhbyG, full.names = F, no.. = T)
scs <- vapply(allsubd, \(i) str_split_1(i, "\\.")[1], "001")
motif_K27acSupEnhbgG <- lapply(allsubd, \(i) {
  loadHomer2KnownResults(file.path(homerd_K27acSupEnhbyG, i), rmdupmotif = T)
}) |>
  setNames(object = _, nm = scs)

## prepare motif matrix
motifmat_K27acSupEnhbyG <- getMotifMat(
  motif_K27acSupEnhbgG,
  field = field, tf = tf)

saveRDS(
  motifmat_K27acSupEnhbyG,
  file.path(
    outd, "MotifMat_H3K27acSupEnhbgGenome.rds"
  )
)

# * heatmap for motif enrichments
## ** visualize H3K27ac DE by 
enhlogcpm <- readRDS(file.path(outd, "ord.H3K27ac.enh.logcpm.rds"))
motifmat_K27acDEbgG <- readRDS(
  file.path(outd, "MotifMat_H3K27acDEbgGenome.rds")
)
motifmat <- motifmat_K27acDEbgG[rownames(enhlogcpm), ]
# p_enrich <- 10e-10
## 64 motifs not enriched in any subclasses
notEnrichedMotifs <- (motifmat < 10) |>
  x => colnames(x)[colSums(x) == nrow(x)]
## 22 motifs showed global enrichments
globalEnrichedMotifs <- (motifmat >= 10) |>
  x => colnames(x)[colSums(x) >= 0.9 * nrow(x)]
motifmat <- motifmat[
  ,
  !(colnames(motifmat) %in% c(notEnrichedMotifs, globalEnrichedMotifs))
]


## ** order the motifs
groupsMat <- groupColsByHierarchical(motifmat, k = 50, sumfn = mean)
ordgroups <- getOrderOfRowsByOrdersOfCols(groupsMat$enrich, topk = 1)
ordmotifs <- getOrderOfElementsByGroupOrder(groupsMat$group, ordgroups)

matmotif <- motifmat[, ordmotifs]
ordmotifs_K27DEbgG <- colnames(matmotif)



## remove (xxx) in the colnames of matmotif
colnames(matmotif) <- gsub("\\(.*\\)", "", colnames(matmotif))
# "OCT:OCT(POU,Homeobox,IR1)"            
# "OCT:OCT(POU,Homeobox)" 
# will be considered as the same motif here.
# now motif dim is: 116 subclasses x 373 motifs
matmotif <- matmotif[ , !duplicated(colnames(matmotif))]
saveRDS(matmotif, file.path(outd, "Heatmap_MotifMat_H3K27acDEbgGenome.rds"))
# draw heatmap of matmotif
# show posisions randomly half of the motifs
mark_pos <- sample(1:ncol(matmotif), ncol(matmotif) / 2)

motif_annot <- columnAnnotation(
  motif = anno_mark(at = mark_pos,
  labels = colnames(matmotif)[mark_pos],
  labels_gp = grid::gpar(fontsize = 7))
)
legend_motif <- list(
  title = "-log10(p-value)",
  at = c(0, 100),
  labels = c(0, 100),
  direction = "horizontal"  
)
hmMotif <- ComplexHeatmap::Heatmap(
  matrix = matmotif,
  col = circlize::colorRamp2(c(0, 100), c("white", "red")),
  cluster_rows = F,
  cluster_columns = F,
  show_row_names = F,
  show_column_names = F,
  use_raster = T,
  show_heatmap_legend = T,
  top_annotation = motif_annot,
  left_annotation = NULL,
  heatmap_legend_param = legend_motif,
  na_col = "white",
  width = 1
  )

show(hmMotif)
# save hmMotif to pdf
withr::with_pdf(file.path(outfigd, "H3K27acDEbgGenome_MotifHeatmap.pdf"), width = 25, height = 15, {
  print(hmMotif)
})

# * repeat the heatmap analysis for H3K27acDEbgOthers
motifmat_K27acDEbgO <- readRDS(
  file.path(outd, "MotifMat_H3K27acDEbgOthers.rds")
)
motifmat <- motifmat_K27acDEbgO[rownames(enhlogcpm), ]
# remove the columns having NA elements
motifmat <- motifmat[ , !apply(motifmat, 2, anyNA)]

# p_enrich <- 10e-10
## 83 motifs not enriched in any subclasses
notEnrichedMotifs <- (motifmat < 10) |>
  x => colnames(x)[colSums(x) == nrow(x)]

## 5 motifs showed global enrichments
globalEnrichedMotifs <- (motifmat >= 10) |>
  x => colnames(x)[colSums(x) >= 0.9 * nrow(x)]
motifmat <- motifmat[
  ,
  !(colnames(motifmat) %in% c(notEnrichedMotifs, globalEnrichedMotifs))
]

groupsMat <- groupColsByHierarchical(motifmat, k = 50, sumfn = mean)
ordgroups <- getOrderOfRowsByOrdersOfCols(groupsMat$enrich, topk = 1)
ordmotifs <- getOrderOfElementsByGroupOrder(groupsMat$group, ordgroups)

# matmotif <- motifmat[, ordmotifs]
matmotif <- motifmat[ , ordmotifs_K27DEbgG[ordmotifs_K27DEbgG %in% colnames(motifmat)]]
## remove (xxx) in the colnames of matmotif
colnames(matmotif) <- gsub("\\(.*\\)", "", colnames(matmotif))
# 116 subclasses x 372 motifs
matmotif <- matmotif[ , !duplicated(colnames(matmotif))]
saveRDS(matmotif, file.path(outd, "Heatmap_MotifMat_H3K27acDEbgOthers.rds"))
# draw heatmap of matmotif
# show posisions randomly half of the motifs
mark_pos <- sample(1:ncol(matmotif), ncol(matmotif) / 2)


motif_annot <- columnAnnotation(
  motif = anno_mark(at = mark_pos,
  labels = colnames(matmotif)[mark_pos],
  labels_gp = grid::gpar(fontsize = 7))
)
legend_motif <- list(
  title = "-log10(p-value)",
  at = c(0, 100),
  labels = c(0, 100),
  direction = "horizontal"  
)
hmMotif <- ComplexHeatmap::Heatmap(
  matrix = matmotif,
  col = circlize::colorRamp2(c(0, 100), c("white", "red")),
  cluster_rows = F,
  cluster_columns = F,
  show_row_names = F,
  show_column_names = F,
  use_raster = T,
  show_heatmap_legend = T,
  top_annotation = motif_annot,
  left_annotation = NULL,
  heatmap_legend_param = legend_motif,
  na_col = "white",
  width = 1
  )

show(hmMotif)
# save hmMotif to pdf
withr::with_pdf(file.path(outfigd, "H3K27acDEbgOthers_MotifHeatmap.pdf"), width = 25, height = 15, {
  print(hmMotif)
})

# * repeat the heatmap analysis for H3K27ac SupEnhbgGenome
# meta
hmd <- file.path(projd, "10.superEnhancer", "out", "heatmap")
supenhd <- file.path(projd, "10.superEnhancer", "out", "ROSE")
scs <- list.files(supenhd, full.names = FALSE, no.. = TRUE, all.files = FALSE)
K27acAllSE <- readRDS(
    file.path(hmd, "signalH3K27ac.allsuperEnhancer.scList.rds"))
allSEmat <- do.call(rbind, K27acAllSE)
rowNormAllSEmat <- (allSEmat * 10^5) / rowSums(allSEmat)
scaleAllSEmat <- scale(rowNormAllSEmat, center = TRUE, scale = TRUE)
sel_scs <- scs[vapply(scs, \(sc) {
    strsplit(sc, "_", fixed = TRUE) |> 
    unlist() |>
    x => as.integer(x[1])
}, 1L) < 500]
str(sel_scs)

motifmat_K27acSupEnhbyG <- readRDS(
  file.path(outd, "MotifMat_H3K27acSupEnhbgGenome.rds")
)

motifmat <- motifmat_K27acSupEnhbyG[sel_scs[sel_scs %in% rownames(motifmat_K27acSupEnhbyG)], ]
# remove the columns having NA elements
motifmat <- motifmat[ , !apply(motifmat, 2, anyNA)]
# 131 not enriched motifs
notEnrichedMotifs <- (motifmat < 10) |>
  x => colnames(x)[colSums(x) == nrow(x)]
## 0 motifs showed global enrichments
globalEnrichedMotifs <- (motifmat >= 10) |>
  x => colnames(x)[colSums(x) >= 0.9 * nrow(x)]
motifmat <- motifmat[
  ,
  !(colnames(motifmat) %in% c(notEnrichedMotifs, globalEnrichedMotifs))
]
groupsMat <- groupColsByHierarchical(motifmat, k = 50, sumfn = mean)
ordgroups <- getOrderOfRowsByOrdersOfCols(groupsMat$enrich, topk = 1)
ordmotifs <- getOrderOfElementsByGroupOrder(groupsMat$group, ordgroups)

matmotif <- motifmat[, ordmotifs]
## remove (xxx) in the colnames of matmotif
colnames(matmotif) <- gsub("\\(.*\\)", "", colnames(matmotif))
# 163 subclasses x 328 motifs
matmotif <- matmotif[ , !duplicated(colnames(matmotif))]

saveRDS(matmotif, file.path(outd, "Heatmap_MotifMat_H3K27acSupEnhbgGenome.rds"))
# draw heatmap of matmotif
# show posisions randomly half of the motifs
mark_pos <- sample(1:ncol(matmotif), ncol(matmotif) / 2)
motif_annot <- columnAnnotation(
  motif = anno_mark(at = mark_pos,
  labels = colnames(matmotif)[mark_pos],
  labels_gp = grid::gpar(fontsize = 7))
)
legend_motif <- list(
  title = "-log10(p-value)",
  at = c(0, 100),
  labels = c(0, 100),
  direction = "horizontal"  
)
hmMotif <- ComplexHeatmap::Heatmap(
  matrix = matmotif,
  col = circlize::colorRamp2(c(0, 100), c("white", "red")),
  cluster_rows = F,
  cluster_columns = F,
  show_row_names = F,
  show_column_names = F,
  use_raster = T,
  show_heatmap_legend = T,
  top_annotation = motif_annot,
  left_annotation = NULL,
  heatmap_legend_param = legend_motif,
  na_col = "white",
  width = 1
  )

show(hmMotif)
# save hmMotif to pdf
withr::with_pdf(file.path(outfigd, "H3K27acSupEnhbgGenome_MotifHeatmap.pdf"), width = 20, height = 15, {
  print(hmMotif)
})
