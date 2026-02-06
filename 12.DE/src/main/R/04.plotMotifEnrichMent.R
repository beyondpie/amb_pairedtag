library(tidyverse)
library(ComplexHeatmap)
library(tmpRpkg)

# * meta
projd <- here::here()
workd <- file.path(projd, "12.DE")
figd <- file.path(workd, "figure")
outd <- file.path(workd, "out")
DARChrAHomerd <- file.path(outd, "DARChrAHomer")
ptGeneCPMfnm <- file.path(
  projd, "data", "pairedtag_ann",
  "pt.RNAseq.pseudobulk.CPM.scbyg.rds")

ptscMeta <- tmpRpkg::loadptscMeta() |>
  x => x[x$ATAC > 0, ]

# * check tf Gene Mat
# 2026-01-01
ptGeneCPMfnm <- file.path(
  projd, "data", "amb_PT_RNA_pseudobulk",
  "pt.RNAseq.pseudobulk.CPM.scbyg.rds")
ptscMetafnm <- file.path(projd, "meta",
  "PairedTagSubclassMetaFromCellMetaChromHMM.csv")
ptscMeta <- data.table::fread(file = ptscMetafnm,
  sep = ",", data.table = F, header = T) |>
  x => `rownames<-`(x, x$PairedTagName) |>
  x => x[x$ATAC > 0, ]

## sc2nmod <- tmpRpkg::getPairedTagSubclass2Nmod(ptscMeta)
## sc2n <- table(seuTF@meta.data$subclass)
## ptscs <- with(sc2nmod, sc[H3K27ac >= 200 & H3K4me1 >= 200])

ptscs <- ptscMeta$PairedTagName

# * function
filterNULLfromList <- function(l) {
  wisnull <- vapply(l, is.null, T)
  if (sum(wisnull) == length(l)) {
    message("All are NULLs, return NULL.")
    NULL
  } else {
    l[!wisnull]
  }
}
extractSigMotif <- function(motifMat, sc, cutoff = 10, ntop = 5) {
  s <- motifMat[sc, ] |>
    setNames(object = _, nm = colnames(motifMat))
  ords <- s[order(s, decreasing = T)] |>
    x => x[seq_len(ntop)]
  x <- ords >= cutoff
  if (sum(x) < 1) {
    message("no significant motifs found: ", sc)
    return(NULL)
  }
  names(ords)[x]
}

getTopKGeneExpPerSubclass <- function(geneMat, ntop) {
  q <- 1 - ntop / ncol(geneMat)
  vapply(rownames(geneMat), \(sc) {
    quantile(x = geneMat[sc, ], probs = q)
  }, 0.0) |>
    setNames(object = _, nm = rownames(geneMat))
}

filterByGeneExtPerSubclass <- function(motifMat, geneMat, ct) {
  r <- motifMat
  tfnms <- colnames(motifMat)
  for(sc in rownames(motifMat)) {
    index <- geneMat[sc, tfnms] < ct[sc]
    r[sc, index] <- 0.0
  }
  return(r)
}

modifyMotifnm <- function(motifnm) {
  gsub("\\(.*", "", motifnm) |>
    gsub("COUP-TFII", "NR2F2", x = _) |>
    gsub("-distal", "", x = _ ) |>
    gsub("\\+.*", "", x = _) |>
    gsub("-AP1", "", x = _) |>
    gsub("NF-E2", "Nfe2", x = _) |>
    gsub("-halfsite", "", x = _) |>
    gsub("n-Myc", "Mycn", x = _) |>
    gsub("c-Myc", "Myc", x = _) |>
    gsub("Nkx2\\.1", "Nkx2-1", x = _) |>
    gsub("Nkx2\\.2", "Nkx2-2", x = _) |>
    gsub("Nkx2\\.5", "Nkx2-5", x = _) |>
    gsub("Nkx3\\.1", "Nkx3-1", x = _) |>
    gsub("Nkx6\\.1", "Nkx6-1", x = _) |>
    gsub("\\+il21", "", x = _) |>
    gsub("HIF-1b", "Arnt", x = _) |>
    gsub("HIF-1a", "Hif1a", x = _) |>
    gsub("\\+1bp", "", x = _) |>
    gsub("OCT:OCT-short", "OCT", x = _) |>
    gsub("OCT:OCT", "OCT", x = _) |>
    gsub("RBPJ:Ebox", "Rbpj", x = _) |>
    gsub("PU\\.1", "Spi1", x = _) |>
    gsub("ZNF143\\|STAF", "Zfp143", x = _) |>
    gsub("NFkB2-p52", "NFkb2", x = _) |>
    gsub("NFkB-p65", "Rela", x = _) |>
    gsub("AP-2alpha", "Tfap2a", x = _)
}

getDiffTFPerSubclass <- function(diffTFmat, ptscs,
                                 padj = 0.05,
                                 log2FC = 0.5) {
  lapply(ptscs, \(sc) {
    subset(
      diffTFmat, (subclass == sc) & (p_val_adj <= padj) & (avg_log2FC >= log2FC))
  }) |>
    setNames(object = _, nm = ptscs) |>
    filterNULLfromList(l = _)
}

# * main
# 1. load motif result on subclass level
motifResults <- lapply(ptscs, \(sc) {
  file.path(DARChrAHomerd, sc) |>
    tmpRpkg::loadHomer2KnownResults(
      resultd = _, rmdupmotif = T)
}) |>
  setNames(object = _, nm = ptscs)

neglog10fn <- \(i) {
  -i * log10(exp(1))
}

motifMat <- tmpRpkg::getMotifMat(
  motifResults, field = "logp", tf = neglog10fn)
saveRDS(motifMat, file.path(outd, "rawMotifEnrichMat.neglog10pval.rds"))

# 2. filter motifMat values
# a. filter motif result by TF expressions
geneLogCPM <- readRDS(ptGeneCPMfnm) |>
  log1p() |>
  x => x[ptscs, ]
ctoff <- getTopKGeneExpPerSubclass(geneLogCPM, ntop = 10000)

# get tfnms from motifs
rawMotifnms <- colnames(motifMat)
mdfnms <- modifyMotifnm(rawMotifnms)

# after modified names, some of them have repeat in colnames of motifMat
motifMat2 <- motifMat[, !duplicated(mdfnms)]
colnames(motifMat2) <- modifyMotifnm(colnames(motifMat2))

motifMat3 <- motifMat2[,
  tolower(colnames(motifMat2)) %in% tolower(colnames(geneLogCPM))] |>
  x => x[ , !duplicated(tolower(colnames(x)))]

# update names as gene names
geneNames <- colnames(geneLogCPM) |>
  g => g[match(x = tolower(colnames(motifMat3)),
    table = tolower(colnames(geneLogCPM)))]
colnames(motifMat3) <- geneNames

# then filter by gene expression
motifMat4 <- filterByGeneExtPerSubclass(
  motifMat = motifMat3, geneMat = geneLogCPM, ct = ctoff)

## scaleMotifMat <- diag(1 / rowSums(motifMat4)) %*% motifMat4 |>
##   x => scale(x)

# b. filter commonly and no-used motifs
commonMotif <- colnames(motifMat4) |>
  x => x[colSums(motifMat4 >= 10) >= (0.7 * nrow(motifMat4))]
noUsedMotif <- colnames(motifMat4) |>
  x => x[colSums(motifMat4 >= 10) <= 0]

motifMat5 <- motifMat4[ ,
  !(colnames(motifMat4) %in% c(commonMotif, noUsedMotif))]

# c. filter by p-value size
# motifMat4[motifMat4 <= 10] <- 0
motifMat5[motifMat5 >= 100] <- 100


# 2. get subclass-level motifs
diffTF <- data.table::fread(
  file = file.path(projd, "16.celloracle", "out", "diff.tf.wilcox.csv"),
  sep = ",",
  header = T,
  data.table = F)

# check
# 2026-01-01
diffTF <- data.table::fread(
  file = file.path(projd, "data", "amb_PT_RNA_pseudobulk", "diff.tf.wilcox.csv"),
  sep = ",",
  header = T,
  data.table = F)

sc2diffTF <- getDiffTFPerSubclass(diffTF, ptscs, padj = 0.05, log2FC = 0.5)

sigMotif <- lapply(ptscs, \(sc) {
  m1 <- extractSigMotif(motifMat5, sc,
    # cutoff: neglog10pval
    cutoff = 10.0, ntop = 10)
  # further filter by diff TFs
  m1[m1 %in% sc2diffTF[[sc]]$tf]
}) |>
  setNames(object = _, nm = ptscs) |>
  tmpRpkg::filterNULLfromList(l = _)
  
m <- motifMat5[names(sigMotif), unique(unlist(sigMotif))]

# 3. plot heatmap of motif enrichments
motifColor = circlize::colorRamp2(
  c(0, 10, 50), c("white", "gray90", "red"))
## pos <- c(1:2,4,6,37:38, 59:60, 63, 68, 70:72, 79:80,82:83,86:87, 125, 170, 195, 198, 202, 214, 314)
## motif_annot <- columnAnnotation(
##   motif = anno_mark(at = pos,
##     labels = colnames(expMotifMat)[pos],
##   labels_gp = grid::gpar(fontsize = 9))
## )
hmMotif <- ComplexHeatmap::Heatmap(
  mat = m,
  ## mat = motifMat,
  col = motifColor,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_column_dend = FALSE,
  row_names_gp = grid::gpar(fontsize = 7),
  column_names_gp = grid::gpar(fontsize = 7),
  # top_annotation = motif_annot,
  left_annotation = NULL,
  use_raster = T,
  show_heatmap_legend = TRUE,
  heatmap_legend_param = list(
    title = latex2exp::TeX(r"(-$\log_{10}$p)"),
    at = c(10, 50),
    labels = c(10, 50),
    direction = "horizontal"
  ),
  na_col = "white" ,
)
hmMotif
saveRDS(object = m, file = file.path(
  projd, "99.figures", "out", "fig4", "panelE.sc2motif.enrich.rds"))

# 4. plot corresponding TF heatmap

# check
# 2026-01-01
m <- readRDS(file.path(projd, "data",
  "amb_PT_RNA_pseudobulk", "panelE.sc2motif.enrich.rds"))
tfmatpre <- geneLogCPM[rownames(m), colnames(m)]
tfmatpre[ , "Tcf4"]

tfmat <- geneLogCPM[rownames(m), colnames(m)]
# FIXME: here why we have consistent values for Tcf4
# tfmat[tfmat > quantile(tfmat, 0.9)] <- quantile(tfmat, 0.9)

tfmat <- scale(tfmat)

quantile(tfmat, na.rm = T)
low.val.col <- -0.5
high.val.col <- 1.5

legend_labels <- c(round(low.val.col, 1),
  round(high.val.col, 1))

col_fun <- circlize::colorRamp2(
  seq(low.val.col, high.val.col, length = 60),
  viridis::viridis(60)
)

hmTFExp <- ComplexHeatmap::Heatmap(
  matrix = tfmat,
  col = col_fun,
  cluster_columns = F,
  cluster_rows = F,
  show_row_names = T,
  row_names_gp = grid::gpar(fontsize = 7),
  column_names_gp = grid::gpar(fontsize = 7),
  show_column_names = T,
  top_annotation = NULL,
  left_annotation = NULL,
  use_raster = T,
  show_heatmap_legend = T,
  na_col = "#440154FF",
  heatmap_legend_param = list(
    title = "scaled logCPM",
    at = c(low.val.col, high.val.col),
    labels = legend_labels,
    direction = "horizontal"
  )
)

hmTFExp
saveRDS(object = tfmat, file = file.path(
  projd, "99.figures", "out", "fig4", "panelE.sc2motif.geneExp.rds"))

withr::with_pdf(
  new = file.path(figd, "motifEnrich.TFExp.seperateHeatmap.pdf"),
  code = {draw(hmMotif + hmTFExp)},
  width = 20, height = 15
)
hmMotif + hmTFExp

# 5. optimize the figure
motifEnrich <- tmpRpkg::to3.matrix(
  m, names = c("subclass", "tf", "neglog10pval"),
  int2str = F, factor2str = T)

tfExp <- tmpRpkg::to3.matrix(
  tfmat, names = c("subclass", "tf", "scaledLogCPM"),
  int2str = F, factor2str = T
)

m4plot <- merge(motifEnrich, tfExp,
  by.x = c("subclass", "tf"),
  by.y = c("subclass", "tf"))

m4plot$subclass <- factor(x = m4plot$subclass,
  levels = rev(rownames(m)))
m4plot$tf <- factor(x = m4plot$tf,
  levels = colnames(m))
m4plot$neglog10pval[m4plot$neglog10pval >= 50] <- 50
m4plot$scaledLogCPM[is.na(m4plot$scaledLogCPM)] <- (-1.0)
m4plot$scaledLogCPM[m4plot$scaledLogCPM < 0.1] <- 0.0
m4plot$scaledLogCPM[m4plot$scaledLogCPM > 1.2] <- 1.2


mytheme <- theme(
  panel.border = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  plot.title = element_text(colour = "black", hjust = 0.5,
    size = 15),
  axis.text = element_blank(),
  axis.ticks = element_blank(),
  axis.title = element_text(colour = "black", size = 12),
  axis.line = element_line(colour = "black"),
  legend.position = "right",
  legend.title = element_text(colour = "black", size = 13),
  legend.text = element_text(colour = "black", size = 12)
)
p <- ggplot(data = m4plot,
  aes(x = tf, y = subclass, fill = neglog10pval)) +
  geom_tile(color = "white", linewidth = 0.15) +
  # scale_fill_distiller(palette = "RdBu", limits = c(0, 50)) +
  scale_fill_gradient(low = "white", high = "red", limits = c(0, 50)) +
  geom_point(aes(x = tf, y = subclass, size = scaledLogCPM),
    alpha = ifelse(m4plot$scaledLogCPM <= 0.8, 0.0, 0.5)) +
  scale_size_continuous(
    breaks = c(0.8, 1.2),
    limits = c(0.8, 1.2),
    range = c(0, 2)
  ) + 
  coord_fixed() +
  mytheme

ggsave(filename = file.path(figd, "DAREnh.MotifEnrich.dotSizedByTFExp.pdf"),
  plot = p, width = 10, height = 20)
