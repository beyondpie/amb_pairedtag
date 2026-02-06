# Seurat with BPCells:
# Ref:
# - https://satijalab.org/seurat/articles/seurat5_bpcells_interaction_vignette
library(Seurat)
library(SeuratObject)
library(tidyverse)
# for merge list of Seurat objects
library(scCustomize)
library(BPCells)
library(ggplot2)
# support R new feature in 4.2
projdir <- here::here()
rdir <- file.path(projdir, "package/R")
import::from(.from = "integration.R", .directory = rdir,
  convertAnn2Seurat5,
  downsampleSeurat,
  toSeuratInMemory,
  isOnDiskMat.Seurat,
  calVarOfFea.Seurat,
  getvf.Seurat,
  setVariableFeatures)

import::from(.from = "colors.R", .directory = rdir,
  SnapATACPalette,
  ArchRPalettes)

# * configs
workdir <- file.path(projdir, "03.integration")
pts5fnm <- file.path(projdir, "data",
  "pairedtag_seurat",
  "nds0.25m_balanceregion",
  "pt_RNA_nds0.25m_balanceregion.rds")

allens5dir <- file.path(projdir, "data", "allen_seurat5")
pt2allenRegionfnm <- file.path(workdir, "src/main/resource",
  "pairedtagRegion2AllenFile.csv")

allenk8fnm <- file.path(workdir, "src/main/resource",
  "AIT21_k8_markers.txt")
k8genes <- data.table::fread(allenk8fnm, header = FALSE,
  data.table = FALSE)$V1
allenannotfnm <- file.path(workdir, "src/main/resource",
  "AIT21_annotation_freeze_081523.tsv")

gp <- "all"
nPCA <- 50
tfmethod <- "cca"
kAnchor <- 50
feaGroup <- "k8"

outdir <- file.path(workdir, "out")
dir.create(outdir, showWarnings = FALSE)
suballens5fnm <- file.path(outdir,
  "allen_10xv3_nds0.25m_balance-ptregion.rds")
tfanchorfnm <- file.path(outdir,
  str_glue("{gp}_{feaGroup}_{tfmethod}_{kAnchor}.balance-pt-region.anchors.rds"))

tfLabelfnm <- file.path(outdir,
  str_glue("{gp}_{feaGroup}_{tfmethod}_{kAnchor}.balance-pt-region.tfLabel.cl.rds"))
## pt2allenr <- data.table::fread(file = pt2allenRegionfnm,
##   header = TRUE, sep = ",", data.table = FALSE)
## all10xv3fnms <- with(pt2allenr, seurats5[grepl("10xv3", seurats5)])
## names(all10xv3fnms) <- gsub(".rds", "", basename(all10xv3fnms))
## all10xv3S5s <- lapply(all10xv3fnms, readRDS)
## names(all10xv3S5s) <- names(all10xv3fnms)
## ptregions <- with(pt2allenr, ptr[grepl("10xv3", seurats5)])
## ncells <- vapply(all10xv3S5s, ncol, 1)

## s510xv3 <- scCustomize::Merge_Seurat_List(
##   list_seurat = all10xv3S5s)
## s510xv3$barcode <- gsub("^_", "", colnames(s510xv3))

## # downsample by region
## ndp <- 250000
## nr <- length(unique(ptregions))
## ndpavg <- ceiling(ndp / nr)
## barcodes <- s510xv3$barcode
## regions <- lapply(seq_along(all10xv3S5s), \(i) {
##   rep(ptregions[i], ncells[i]) }) |>
##   do.call(what  = "c", args = _)
## allenb2r <- data.frame(
##   barcode = barcodes,
##   region = regions
## )

## sub_allenb2r <- allenb2r |>
##   dplyr::group_by(region) |>
##   dplyr::slice_sample(n = ndpavg)

## sub_s510xv3 <- subset(s510xv3,
##   subset = barcode %in% sub_allenb2r$barcode)

## sub_s510xv3 <- JoinLayers(sub_s510xv3)
## # about 5 minutes, 30G
## sub_s510xv3[["RNA"]]$counts <- as(
##   object = sub_s510xv3[["RNA"]]$counts, Class = "dgCMatrix")
## sub_s510xv3 <- SeuratObject::RenameCells(
##   sub_s510xv3, new.names = sub_s510xv3$barcode)
## # takes about 10 minutes with 3G on disk.
## # maybe we can save the JoinLayers
## saveRDS(sub_s510xv3, suballens5fnm)

sub_s510xv3 <- readRDS(suballens5fnm)
# * load downsampled sa2 ann with region
# add UMAP and cl1 result to meta
# check Snap25, Gad1, Gad2 and other neuronal markers
ptSeu <- readRDS(pts5fnm)
ptSeu[["RNA"]]$counts <- as(
  object = ptSeu[["RNA"]]$counts, Class = "dgCMatrix")

# * transfer label
ref <- sub_s510xv3
query <- ptSeu

# set NormalizeData before this
# to set data field for the two seurats
ref <- Seurat::NormalizeData(
  ref,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  margin = 1)
query <- Seurat::NormalizeData(
  query,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  margin = 1
)

features <- setVariableFeatures(
  ref = ref,
  query = query,
  features = k8genes,
  onlayer = "counts",
  eps = 0.0001
)

# 740G at most; takes about 30-60 minutes
# sometimes use 10 - 100 CPUs in mediator
# Found 2203029 anchors
# After this function, 220G mem on hold
anchors <- Seurat::FindTransferAnchors(
  reference = ref,
  query = query,
  scale = TRUE,
  reduction = tfmethod,
  normalization.method = "LogNormalize",
  npcs = nPCA,
  dims = seq_len(nPCA),
  l2.norm = TRUE,
  k.anchor = kAnchor,
  features = features,
  verbose = TRUE
)
saveRDS(anchors, tfanchorfnm)

# 320g in mem
tfLabel <- Seurat::TransferData(
  anchorset = anchors,
  refdata = ref$cl,
  query = query,
  dims = seq_len(nPCA),
  weight.reduction = tfmethod,
  l2.norm = FALSE
)
# add transfer label and scores to sa2 S5
# 1.8G in total on disk
saveRDS(tfLabel, tfLabelfnm)

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

outfigdir <- file.path(workdir, "out", "figures")
dir.create(outfigdir)

allenannot <- data.table::fread(
  file = allenannotfnm, sep = "\t", header = TRUE,
  data.table = FALSE)
rownames(allenannot) <- allenannot$cl
anchors <- readRDS(tfanchorfnm)
refquerymeta <- anchors@object.list[[1]]@meta.data
refmeta <- refquerymeta[!is.na(refquerymeta$cl), ]
# 211
length(table(allenannot[refmeta$cl, "subclass_id"]))

ptSeu <- readRDS(tfLabelfnm)
ptSeu$predicted.scid <- allenannot[
  ptSeu$predicted.id, "subclass_id"]
ptSeu$predicted.sc.lb <- allenannot[
  ptSeu$predicted.id, "subclass_label"]
# 195
length(table(ptSeu$predicted.scid))
ptSeu$predicted.clid <- allenannot[
  ptSeu$predicted.id, "class_id"]
ptSeu$predicted.cl.lb <- allenannot[
  ptSeu$predicted.id, "class_label"]
ptSeu$predicted.spid <- allenannot[ptSeu$predicted.id,
  "supertype_id"]
ptSeu$predicted.sp.lb <- allenannot[ptSeu$predicted.id,
  "supertype_label"]

umapdf <- data.frame(
  barcode = ptSeu$barcode,
  UMAP1 = ptSeu$umap1,
  UMAP2 = ptSeu$umap2,
  pred.cl = factor(ptSeu$predicted.cl.lb),
  pred.sp = factor(ptSeu$predicted.sp.lb),
  pred.sc = factor(ptSeu$predicted.sc.lb),
  pred.id = factor(ptSeu$predicted.id),
  pred.id.score = factor(ptSeu$predicted.id.score)
)

umapdf.ds <- umapdf[sample(seq_len(nrow(umapdf)), size = 100000), ]
(
  umap.cl <- ggplot(data = umapdf.ds,
    mapping = aes(x = UMAP1, y = UMAP2, color = pred.cl)) +
    geom_point(size = 0.8, shape = 19, alpha = 0.7) +
    labs(
      title = "Integration with Allen: 0.25 million cells from each dataset",
      x = "UMAP1", y = "UMAP2") +
    guides(color = guide_legend(override.aes = list(size = 6))) +
    mytheme +
    scale_color_manual(values = SnapATACPalette)
)

sc.colors <- grDevices::colorRampPalette(
  colors = c("blue", "purple", "yellow", "red", "orange", "white"),
  bias = 1,
  space = "rgb",
  interpolate = "linear",
  alpha = FALSE)(length(unique(umapdf$pred.sc)))

(
  umap.sc <- ggplot(data = umapdf.ds,
    mapping = aes(x = UMAP1, y = UMAP2, color = pred.sc)) +
    geom_point(size = 0.8, shape = 19, alpha = 0.7) +
    labs(
      title = "Integration with Allen: 0.25 million cells from each dataset",
      x = "UMAP1", y = "UMAP2") +
    guides(color = guide_legend(override.aes = list(size = 6))) +
    mytheme +
    theme(legend.position = "none") +
    scale_color_manual(values = sc.colors)
)

# TODO: plot GABA, Glut, NN

# check the corresponding gene markers
# TODO: co-embedding with the two modalities

# * add test data set for transfer label in tools
allenMeta <- data.table::fread(
  file = file.path(workdir, "src/main/resource",
    "AIT21_annotation_freeze_081523.tsv"), sep = "\t",
  header = TRUE, data.table = FALSE
)
rownames(allenMeta) <- paste0("c", allenMeta$cl)
sub_s510xv3 <- readRDS(suballens5fnm)
sub_s510xv3$cl <- paste0("c", sub_s510xv3$cl)
sub_s510xv3$clidlabel <- allenMeta[sub_s510xv3$cl, "cluster_id_label"]
sub_s510xv3$spidlabel <- allenMeta[sub_s510xv3$cl, "supertype_id_label"]
sub_s510xv3$scidlabel <- allenMeta[sub_s510xv3$cl, "subclass_id_label"]
sub_s510xv3$cidlabel <- allenMeta[sub_s510xv3$cl, "class_id_label"]

barcodes <- sub_s510xv3@meta.data |>
  dplyr::group_by(scidlabel) |>
  dplyr::slice_sample(n = 50) |>
  x => x$barcode

ref <- sub_s510xv3[ , barcodes]
saveRDS(ref, file.path(workdir, "src/test/resource/allen_ref.rds"))

ptSeu <- readRDS(pts5fnm)

ptMeta <- data.table::fread(
  file = file.path(projdir, "meta", "pt.barcode.meta.withL4.csv"),
  sep = ",", header = TRUE, data.table = FALSE
)
rownames(ptMeta) <- ptMeta$barcode

query <- ptSeu[, rownames(ptMeta)] |>
  x => x[, sample(seq_len(ncol(x)), size = ncol(ref))]
query[["RNA"]]$counts <- as(
  object = query[["RNA"]]$counts, Class = "dgCMatrix")

query$L2 <- ptMeta[colnames(query), "L1_2"]
query$L3 <- ptMeta[colnames(query), "L1_2_3"]
query$L4 <- ptMeta[colnames(query), "L1_2_3_4"]
saveRDS(query, file.path(workdir, "src/test/resource",
  "pt_query.rds"))
