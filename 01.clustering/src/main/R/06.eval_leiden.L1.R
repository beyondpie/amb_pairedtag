library(tidyverse)
library(rlang)
Sys.setenv("_R_USE_PIPEBIND_" = TRUE)

projroot <- here::here()
rdir <- file.path(projroot, "package/R")
import::from(.from = "colors.R", .directory = rdir,
  ArchRPalettes)

projroot <- here::here()
rdir <- file.path(projroot, "package/R")
metadir <- file.path(projroot, "meta")
workdir <- file.path(projroot, "01.clustering")

# * afqc configs
cellmetafnm <- file.path(metadir, "pt.barcode.meta.with.qc.v2.csv")
L1metafnm <- file.path(workdir, "out/afqc/L1",
  "RNA_k8_npc50_k50_L0_0.csv")
silhtmetafnm <- file.path(workdir, "out/afqc/L1",
  "RNA_k8_npc50_k50_L0_0_silht.csv")
nr <- 10
ns <- 4
outdir <- file.path(workdir, "out/afqc/L1/figures")

# * nodlt configs
## cellmetafnm <- file.path(metadir, "pairedtag.meta.v1.csv")
## L1metafnm <- file.path(workdir, "out/nodlt/L1",
##   "RNA_k8_npc50_k50_L0_0.csv")
## silhtmetafnm <- file.path(workdir, "out/nodlt/L1",
##   "RNA_k8_npc50_k50_L0_0_silht.csv")
## nr <- 10
## ns <- 6
## outdir <- file.path(workdir, "out/nodlt/L1/figures")

getSils <- function(r) {
  rownms <- paste(paste0("r", format(r, nsmall = 1)),
    paste0("seed", seq(0, nr - 1)), sep = "_")
  silhtmeta[rownms, paste0("seed", seq(0, ns - 1))]
}
getAvgSil <- function(r) {
  s <- getSils(r)
  ## s_1 <- s > 0
  ## avgs <- rowSums(s * s_1) / rowSums(s_1)
  avgs <- rowMeans(s)
  return(mean(avgs))
}

# load file
cellmeta <- data.table::fread(
  file = cellmetafnm,
  sep = ",",
  header = TRUE,
  data.table = FALSE
)
rownames(cellmeta) <- cellmeta$barcode

L1meta <- data.table::fread(
  file = L1metafnm,
  sep = ",",
  header = TRUE,
  data.table = FALSE)
rownames(L1meta) <- L1meta$barcode

silhtmeta <- data.table::fread(
  file = silhtmetafnm,
  sep = ",",
  header = TRUE,
  data.table = FALSE)
rownames(silhtmeta) <- silhtmeta$reso_seed

# select reso based on silht
rs <- seq(0.1, 2, 0.1)

avgSils <- vapply(seq(0.1, 2, 0.1),
  FUN = getAvgSil, FUN.VALUE = 0.1)

# * for afqc
# r = 0.1
r_maxsil <- rs[which.max(avgSils)]
sils_argmax <- getSils(r_maxsil)
# seed = 1
seed_maxsil <- seq(0, nr - 1)[
  which.max(rowMeans(sils_argmax))]


umap.df <- data.frame(
  barcode = L1meta$barcode,
  UMAP1 = L1meta$umap1,
  UMAP2 = L1meta$umap2,
  BrainRegion = factor(cellmeta[rownames(L1meta), "brainregion"]),
  HistoneModality = factor(cellmeta[rownames(L1meta), "modality"]),
  Sex = factor(cellmeta[rownames(L1meta), "sex"]),
  Rep = factor(
    gsub("Female|Male", "", cellmeta[rownames(L1meta), "rep"])),
  Exp = factor(cellmeta[rownames(L1meta), "exp"]),
  Leiden = factor(L1meta[,
    paste(paste0("r", format(r_maxsil, nsmall = 1)),
      paste0("seed", seed_maxsil), sep = "_")])
)

umap.df.ds <- umap.df[
  sample(seq_len(nrow(umap.df)), size = 100000), ]

# plot umap
dir.create(outdir)

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

(
  pu_sex <- ggplot(
    data = umap.df.ds,
    mapping = aes(x = UMAP1, y = UMAP2,
      color = Sex)) +
    geom_point(size = 1, shape = 19, alpha = 0.5) +
    labs(title = "Adult PairedTag 1.5 million snRNA-seq on Sex",
      x = "UMAP1", y = "UMAP2") +
    guides(color = guide_legend(override.aes = list(size = 8))) +
    mytheme
)

ggsave(filename = file.path(outdir, "RNA_L1_UMAP_sex.pdf"),
  plot = pu_sex, width = 9, height = 7)

(
  pu_histone <- ggplot(
    data = umap.df.ds,
    mapping = aes(x = UMAP1, y = UMAP2,
      color = HistoneModality)) +
    geom_point(size = 1, shape = 19, alpha = 0.5) +
    labs(title = "Adult PairedTag 1.5 million snRNA-seq on Histone Mods",
      x = "UMAP1", y = "UMAP2") +
    guides(color = guide_legend(override.aes = list(size = 8))) +
    mytheme +
    scale_color_brewer(palette="Spectral")
)

ggsave(filename = file.path(outdir, "RNA_L1_UMAP_histone.pdf"),
  plot = pu_histone, width = 9, height = 7)

(
  pu_region <- ggplot(
    data = umap.df.ds,
    mapping = aes(x = UMAP1, y = UMAP2,
      color = BrainRegion)) +
    geom_point(size = 1, shape = 19, alpha = 0.5) +
    labs(title = "Adult PairedTag 1.5 million snRNA-seq on Brain Region",
      x = "UMAP1", y = "UMAP2") +
    guides(color = guide_legend(override.aes = list(size = 8))) +
    mytheme +
    scale_color_brewer(palette = "Set1")
)

ggsave(filename = file.path(outdir, "RNA_L1_UMAP_region.pdf"),
  plot = pu_region, width = 9, height = 7)

(
  pu_rep <- ggplot(
    data = umap.df.ds,
    mapping = aes(x = UMAP1, y = UMAP2,
      color = Rep)) +
    geom_point(size = 1, shape = 19, alpha = 0.5) +
    labs(title = "Adult PairedTag 1.5 million snRNA-seq on Replicate",
      x = "UMAP1", y = "UMAP2") +
    guides(color = guide_legend(override.aes = list(size = 8))) +
    mytheme
)

ggsave(filename = file.path(outdir, "RNA_L1_UMAP_rep.pdf"),
  plot = pu_rep, width = 9, height = 7)

p_sum <- ggpubr::ggarrange(
  pu_sex, pu_rep, pu_histone, pu_region,
  labels = c("A", "B", "C", "D"),
  ncol = 2, nrow = 2
)
ggsave(filename = file.path(outdir, "RNA_L1_UMAP_sum.pdf"),
  plot = p_sum, width = 15, height = 12)
(
  pu_leiden <- ggplot(
    data = umap.df.ds,
    mapping = aes(x = UMAP1, y = UMAP2,
      color = Leiden)) +
    geom_point(size = 0.7, shape = 19, alpha = 0.8) +
    labs(title = paste(
      "Adult PairedTag snRNA-seq on Leiden",
      paste(paste("r:", r_maxsil), paste("seed:", seed_maxsil)),
      sep = "\n"),
    x = "UMAP1", y = "UMAP2") +
    guides(color = guide_legend(override.aes = list(size = 8))) +
    mytheme +
    scale_color_manual(values = unname(ArchRPalettes$bear))
)

ggsave(filename = file.path(outdir, "RNA_L1_UMAP_leiden.pdf"),
  plot = pu_leiden, width = 10, height = 10)
# select reso based on silht

cellmeta <- data.table::fread(
  file = cellmetafnm,
  sep = ",",
  header = TRUE
)
data.table::setkey(cellmeta, "barcode")
data.table::setkey(L1meta, "barcode")
L1meta <- data.table::fread(
  file = L1metafnm,
  sep = ",",
  header = TRUE
)
leiden_col <- "r0.1_seed1"
L1meta <- L1meta[ , c("barcode", "r0.1_seed1", "umap1", "umap2")]
data.table::setkey(L1meta, "barcode")

cellmetaL1 <- merge(cellmeta, L1meta)
cellmetaL1$L1 <- cellmetaL1$r0.1_seed1
cellmetaL1$r0.1_seed1 <- NULL
data.table::fwrite(x = cellmetaL1,
  file = file.path(metadir, "pt.barcode.meta.withL1.csv"),
  sep = ",",
  row.names = FALSE,
  col.names = TRUE
  )

# =============================================
# codes below are mainly used for nodlt
# to double check the qc.
# =============================================

## # * for nodlt
## # r = 0.2, with sil score 0.2
## r_maxsil <- rs[which.max(avgSils)]
## sils_argmax <- getSils(r_maxsil)
## # seed = 9
## seed_maxsil <- seq(0, nr - 1)[
##   which.max(rowMeans(sils_argmax))]

# * check doublet cells in the UMAP.
dlts <- data.table::fread(
  file.path(metadir, "sa2_scrublet_8k_all.csv"),
  header = TRUE,
  sep = ",",
  data.table = FALSE)
rownames(dlts) <- dlts$barcode
dlts$isdlt <- ifelse(test = dlts$dlt_prob >= 0.1, "yes", "no")
umap.dlt <- data.frame(
  barcode = L1meta$barcode,
  UMAP1 = L1meta$umap1,
  UMAP2 = L1meta$umap2,
  dlt_prob = dlts[L1meta$barcode, "dlt_prob"],
  dlt_score = dlts[L1meta$barcode, "dlt_score"],
  is_dlt = factor(dlts[L1meta$barcode, "isdlt"])
)
## ndlt = sum(umap.dlt$is_dlt)
rownames(umap.dlt) <- umap.dlt$barcode

umap.dlt.ds <- umap.dlt[
  sample(seq_len(nrow(umap.dlt)), size = 100000), ]

(
  umap.dlt.prob <- ggplot(
    data = umap.dlt.ds,
    mapping = aes(x = UMAP1, y = UMAP2,
      color = dlt_prob)) +
    geom_point(size = 0.6, shape = 19, alpha = 0.5) +
    guides(color = guide_legend(override.aes = list(size = 8))) +
    ## scale_color_gradient(low = "grey", high = "red")+
    scale_color_gradient2(low = "blue", high = "red",
      mid = "white", space = "Lab", midpoint = 0.05) +
    mytheme
)

# too slow
(
  umap.dlt <- ggplot(
    data = umap.dlt.ds,
    mapping = aes(x = UMAP1, y = UMAP2,
      color = is_dlt)) +
    geom_point(size = 0.6, shape = 19, alpha = 0.5) +
    guides(color = guide_legend(override.aes = list(size = 8))) +
    scale_color_manual(values = c("grey", "red")) +
    ## labs(title = "Doublet removal with SnapATAC2 using 8k features",
    ## subtitle = str_glue("In total, {ndlt} doublets.")) +
    ## scale_color_gradient2(low = "blue", high = "red",
    ##   mid = "white", space = "Lab", midpoint = 0.5) +
    mytheme
)

# ** define doublets per sample
# match the ratio defined in the table

# * select the non-replicable group
# identification of high-expressed genes
# r0.2_seed9
## L1reso <- str_glue("r{r_maxsil}_seed{seed_maxsil}")
(

  pu_cl11 <- ggplot(
    data = umap.df.ds,
    mapping = aes(x = UMAP1, y = UMAP2,
      color = Leiden)) +
    geom_point(size = 0.7, shape = 19, alpha = 0.8) +
    labs(title = paste(
      "Adult PairedTag 1.5 million snRNA-seq on Leiden",
      paste(paste("r:", r_maxsil), paste("seed:", seed_maxsil)),
      sep = "\n"),
    x = "UMAP1", y = "UMAP2") +
    guides(color = guide_legend(override.aes = list(size = 8))) +
    mytheme +
    scale_color_manual(values = c("11" = "red"))
)

## output the cell barcodes
barcodes_cl11 <- with(L1meta, barcode[r0.2_seed9 == 11])
write.table(barcodes_cl11,
  file = file.path(projroot, "meta", "lowrep_L1cl-11.barcode.txt"),
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE)

# get marker genes in barcodes_cl11
# use down-sampled pairedtag Seurat
library(Seurat)
library(SeuratObject)
library(reticulate)
library(Matrix)
reticulate::use_condaenv(
  condaenv = "sa2",
  conda = "/home/szu/miniforge3/bin/mamba",
  required = TRUE
)
sa2 <- reticulate::import(module = "snapatac2")
sa2annfnm <- file.path(projroot, "01.clustering",
  "out/scRNAseq_sa2_all.ann.h5ad")
nds <- 12000

barcodes_cl11 <- read.table(
  file = file.path(projroot, "meta", "lowrep_L1cl-11.barcode.txt"),
  header = FALSE)$V1
outdir <- file.path(projroot, "01.clustering", "out",
  "nodlt", "L1", "pt_seurat")
outs5dir <- file.path(outdir, "pt_RNA_nds12k_L1_s5")
outs5matdir <- file.path(outs5dir, "_mat")
for (d in c(outdir, outs5dir, outs5matdir)) {
  dir.create(d)
}
sa2anndsfnm <- file.path(outdir, "pt_RNA_sa2_nds12k_L1.h5ad")
outs5fnm <- file.path(outs5dir, "pt_RNA_nds12k_L1.seurat.rds")
out.allmarkers.fnm <- file.path(outs5dir, "pt_RNA_nds12k_L1_k8.markers.rds")

sa2ann <- sa2$read(sa2annfnm)
barcodesds <- L1meta |>
  dplyr::group_by(r0.2_seed9) |>
  dplyr::slice_sample(n = nds) |>
  x => x$barcode
# table(L1meta[barcodesds, "r0.2_seed9"])
sa2annds <- sa2ann$subset(obs_indices = barcodesds,
  out = sa2anndsfnm)
sa2ann$close()
sa2annds$close()

sa2annds <- sa2$read(sa2anndsfnm, backed = NULL)
mat <- t(sa2annds$X)
rownames(mat) <- sa2annds$var_names$to_list()
colnames(mat) <- sa2annds$obs_names$to_list()
mat <- MatrixExtra::as.csc.matrix(mat) |>
  BPCells::convert_matrix_type(matrix = _, type = "uint32_t")
BPCells::write_matrix_dir(mat = mat, dir = outs5matdir,
  overwrite = TRUE)

cellmetads <- L1meta[colnames(mat),
  c("barcode", "umap1", "umap2", "r0.2_seed9")]
cellmetads$exp <- sa2annds$obs[colnames(mat), "sample"]
cellmetads$L1 <- cellmetads$r0.2_seed9

d <- BPCells::open_matrix_dir(outs5matdir)
s5 <- Seurat::CreateSeuratObject(counts = d,
  assay = "RNA", meta.data = cellmetads)

s5 <- Seurat::NormalizeData(object = s5,
  normalization.method = "LogNormalize")
s5 <- Seurat::ScaleData(object = s5)
Seurat::Idents(s5) <- s5$L1
s5[["umap"]] <- Seurat::CreateDimReducObject(
  embeddings = as.matrix(s5@meta.data[, c("umap1", "umap2")]),
  key = "UMAP_",
  assay = DefaultAssay(s5)
)
saveRDS(s5, outs5fnm)

allenk8fnm <- file.path(projroot, "03.integration",
  "src/main/resource",
  "AIT21_k8_markers.txt")
k8genes <- data.table::fread(allenk8fnm, header = FALSE,
  data.table = FALSE)$V1

k8genes <- intersect(k8genes, rownames(s5))

all.markers <- Seurat::FindAllMarkers(
  object = s5,
  features = k8genes,
  logfc.threshold = 0.5,
  test.use = "wilcox",
  slot = "data",
  verbose = TRUE
)
saveRDS(all.markers, out.allmarkers.fnm)

# double check marker genes
s5 <- readRDS(outs5fnm)
s5@assays$RNA$data@matrix@matrix@matrix@matrix@matrix@matrix@dir <-
  outs5matdir
s5@assays$RNA$scale.data@matrix@matrix@matrix@matrix@matrix@matrix@matrix@matrix@matrix@dir <-
  outs5matdir

all.markers <- readRDS(out.allmarkers.fnm)
all.markers$cluster <- levels(all.markers$cluster)[all.markers$cluster]
markers.cl11 <- all.markers[
(all.markers$cluster == "11") &
  (all.markers$avg_log2FC > 0.5) & (all.markers$pct.1 > 0.1), ]

markers.cl11 <- markers.cl11[order(markers.cl11$avg_log2FC, decreasing = TRUE), ]
markers.cl11 <- markers.cl11[order(markers.cl11$pct.1, decreasing = TRUE), ]
Seurat::DimPlot(s5, reduction = "umap")
Seurat::FeaturePlot(s5,
  features = markers.cl11[1:20, "gene"],
  reduction = "umap",
  slot = "scale.data")

write.table(markers.cl11,
  file = file.path(outs5dir, "markers.cl11.csv"),
  quote = FALSE, row.names = FALSE, col.names = TRUE,
  sep = ",")

# * double check the doublet removal
barcode2dlt.filtered <- data.table::fread(
  file.path(projroot, "meta",
    "barcode2dlt.filtered.L1.k8.dltscore.csv"),
  sep = ",", header = TRUE
)
data.table::setkey(barcode2dlt.filtered, "barcode")

# umap.df defined early in this script
umap.df$idlt <- ifelse(
  !umap.df$barcode %in% barcode2dlt.filtered$barcode,
  "yes", "no")
ndlt <- sum(umap.df$idlt == "yes")
nall <- nrow(umap.df)
nndlt <- nall - ndlt

umap.df.ds <- umap.df[
  sample(seq_len(nrow(umap.df)), size = 500000), ]
(
  pu_dlt <- ggplot(
    data = umap.df.ds,
    mapping = aes(x = UMAP1, y = UMAP2,
      color = idlt)) +
    geom_point(data = umap.df.ds[umap.df.ds$idlt == "no",],
      size = 1, shape = 19, alpha = 0.5) +
    geom_point(data = umap.df.ds[umap.df.ds$idlt == "yes", ],
      size = 0.8, shape = 19, alpha = 0.8) +
    labs(title =
           str_glue("SnapATAC2 Scrublet with scRNAseq on K8: ",
             "{ndlt} dlts {round(ndlt *100/ nall, 2)}%, ",
             "{nndlt} nndlts among {nall} cells"),
      x = "UMAP1", y = "UMAP2") +
    guides(color = guide_legend(override.aes = list(size = 8))) +
    mytheme +
    scale_color_manual(values = c("yes" = "red", "no" = "gray90"))
  ## scale_color_brewer(palette="Spectral")
)

# double check the score / prob
exp2dlthres <- read.table(
  file.path(projroot, "meta", "exp2dlthres.L1.k8.csv"),
  header = TRUE, sep = ",")
rownames(exp2dlthres) <- exp2dlthres$exp

exp2dltr <- read.table(
  file.path(projroot, "meta", "doublet.rate.est.csv"),
  header = TRUE, sep = ","
)

ae <- "Exp12"
index <- barcode2dlt.filtered$exp == ae
sum(index)
sum(barcode2dlt.filtered$dlt_score[index] <=
      exp2dlthres[ae, "dltscore_thres"])
sum(barcode2dlt.filtered$dlt_prob[index] <= 0.5)

sum(cellmeta$exp  == ae)
1 - sum(index) / sum(cellmeta$exp  == ae)
