#!/usr/bin/env Rscript

Sys.setenv("_R_USE_PIPEBIND_" = TRUE)
library(optparse)
library(tidyverse)
library(Seurat)
library(patchwork)
library(ggrastr)

args <- OptionParser(
  usage = "Summarize Transfer Label from Allen scRNA-seq data",
  description = "
 Generate UMAP (major region; tech) and consensus matrix.
 NOTE:
 - R version >= 4.2
 - Seurat version >= 5.0
 - subclass and supertype in the output of seurat is already based
   on vote strategy; but cl is the original predicted.cl.

 TODO:
 add region to UMAP.
") |>
  add_option(
    opt_str = c("-a", "--acfnm"), type = "character",
    help = "anchor filename"
  ) |>
  add_option(
    opt_str = c("-q", "--qrfnm"), type = "character",
    help = "query with tf filename"
  ) |>
  add_option(
    opt_str = c("--ametafnm"), type = "character",
    help = "allen meta fnm"
  ) |>
  add_option(
    opt_str = c("--pmetafnm"), type = "character",
    help = "paired-tag meta fnm"
  ) |>
  add_option(
    opt_str = c("--acmetafnm", type = "character",
      default = NULL,
      help = "allen single-cell meta fnm default NULL")
  ) |>
  add_option(
    opt_str = c("--outseufnm"), type = "character",
    help = "outfnm of seurat object for both ref and query"
  ) |>
  add_option(
    opt_str = c("--outumaprefix"), type = "character",
    default = "./tf",
    help = "prefix with dir for out umap figure , default is ./tf"
  ) |>
  add_option(
    opt_str = c("--outcnssprefix", type = "character",
      help = "prefix with dir of consensus matrix files and figures")
  ) |>
  add_option(
    opt_str = c("--rd", type = "character",
      default = "cca.l2",
      help = "reduction filed in anchor, default cca.l2")
  ) |>
  add_option(
    opt_str = c("--dmtrc"), type = "character",
    default = "euclidean",
    help = "distance metric for UMAP, default enclidean"
  ) |>
  add_option(
    opt_str = c("--mdist"), type = "double",
    default = 0.3,
    help = "min dist for UMAP, default 0.3"
  ) |>
  add_option(
    opt_str = c("--kumap"), type = "integer",
    default = 30,
    help = "neighbors in umap, default 30"
  ) |>
  add_option(
    opt_str = c("--lightcolor"), type = "character",
    default = "white",
    help = "light color for consensus map"
  ) |>
  add_option(
    opt_str = c("--darkcolor"), type = "character",
    default = "red",
    help = "dark color for consensus map"
  ) |>
  parse_args(object = _)

# * debug parameters
## projdir <- here::here()
## args$acfnm <- file.path(projdir, "03.integration",
##   "out", "tfneu_k8_pca_k50", "tf.anchors.with-pcaproject-kac50.rds")
## args$qrfnm <- file.path(projdir, "03.integration",
##   "out", "tfneu_k8_pca_k50", "query.with.tf-pcaproject-kac50_on-cl.rds")
## args$ametafnm <- file.path(projdir, "meta",
##   "AIT21_annotation_freeze_081523.tsv")
## args$pmetafnm <- file.path(projdir, "meta",
##   "pt.barcode.meta.L5.csv")
## args$outseufnm <- file.path(projdir, "03.integration",
##   "out", "tfneu_k8_pca_k50", "sum_tfneu_k8_pca_k50.seurat.rds")
## args$outumaprefix <- file.path(projdir, "03.integration",
##   "out", "tfneu_k8_pca_k50", "sum_tfneu_k8_pca")
## args$outcnssprefix <- file.path(projdir, "03.integration",
##   "out", "tfneu_k8_pca_k50", "sum_tfneu_k8_pca")
## args$rd <- "pcaproject.l2"
## args$dmtrc <- "euclidean"
## args$mdist <- 0.3
## args$kumap <- 15
## args$lightcolor <- "white"
## args$darkcolor <- "red"
## args$acmetafnm <- file.path(
##   projdir, "data/allen", "allen.10xv3.cell.meta.csv"
## )

outumapfig <- paste(args$outumaprefix,
  args$rd, "UMAParams",
  args$dmtrc, args$mdist, args$kumap, "pdf",
  sep = ".")

# * functions
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

toNamedArray.1dtable <- function(t) {
  tmp <- as.data.frame(table(t),
    stringsAsFactors = FALSE)
  r <- tmp[, 2]
  names(r) <- tmp[, 1]
  return(r)
}

to3.matrix <- function(mat,
                       names = c("row", "column", "score"),
                       int2str = TRUE,
                       factor2str = TRUE) {
  r <- reshape2::melt(mat)
  colnames(r) <- names
  if (factor2str) {
    if (is.factor(r[, 1])) {
      r[, 1] <- as.character(r[, 1])
    }
    if (is.factor(r[, 2])) {
      r[, 2] <- as.character(r[, 2])
    }
  }
  if (int2str) {
    if (is.numeric(r[, 1])) {
      r[, 1] <- as.character(r[, 1])
    }
    if (is.numeric(r[, 2])) {
      r[, 2] <- as.character(r[, 2])
    }
  }
  return(r)
}

getMaxColScoreWithName <- function(mat) {
  maxScores <- apply(mat, 2, max)
  rnames <- rownames(mat)
  maxNames <- apply(mat, 2, \(x) {
    rnames[which.max(x)]
  })
  names(maxScores) <- maxNames
  return(maxScores)
}


# * main
# 1. check arguments
message("Check arguments.")
for (f in with(args, c(acfnm, qrfnm, ametafnm, pmetafnm))) {
  if (!file.exists(f)) {
    stop("Cannot find ", f)
  }
}
for (d in with(args, c(outseufnm, outumapfig, outcnssprefix))) {
  dir.create(dirname(d), showWarnings = FALSE)
}

# 2. prepare meta for plotting
anchor <- readRDS(args$acfnm)
intseu <- anchor@object.list[[1]]

anchorMeta <- intseu@meta.data
anchorMeta$barcode <- gsub(
  "_reference|_query", "", rownames(anchorMeta)
)

anchorMeta$isRef <- FALSE
anchorMeta[grepl("reference", rownames(anchorMeta)), "isRef"] <- TRUE
rownames(anchorMeta) <- anchorMeta$barcode
anchorMeta$tech <- "PairedTagRNA"
anchorMeta$tech[anchorMeta$isRef] <- "AllenRNA"

# 2.1 umap calculation
message("Run umap.")
rd <- anchor@object.list[[1]]@reductions[[args$rd]]@cell.embeddings
rdUMAP <- RunUMAP(
  object = rd,
  umap.method = "uwot",
  n.neighbors = args$kumap,
  metric = args$dmtrc,
  min.dist = args$mdist
)
rownames(rdUMAP@cell.embeddings) <- anchorMeta$barcode
rdumap <- rdUMAP@cell.embeddings

# 2.2 cell meta join
# 2.2.1 allen's meta data
allenMeta <- data.table::fread(
  file = args$ametafnm,
  header = TRUE,
  data.table = FALSE
)
rownames(allenMeta) <- paste0("cl-", allenMeta$cl)
addRefColwithNA <- function(m,
                            coltoAdd = "subclass_id_label",
                            colnew = coltoAdd) {
  if (!"cl" %in% colnames(m)) {
    message("cl column not in allen meta data.")
    message("load allen single-cell meta data: ", args$acmetafnm)
    acmeta <- data.table::fread(
      args$acmetafnm, sep = ",", header = TRUE,
      data.table = FALSE
    )
    rownames(acmeta) <- acmeta$barcode
    m$cl <- NA
    rbarcodes <- with(m, barcode[isRef])
    m[rbarcodes, "cl"] <- acmeta[rbarcodes, "cl"]
  }
  m[[colnew]] <- NA
  clRef <- paste0("cl-", with(m, cl[isRef]))
  m[m$isRef, colnew] <- allenMeta[clRef, coltoAdd]
  return(m)
}
anchorMeta <- addRefColwithNA(anchorMeta, coltoAdd = "subclass_id_label",
  colnew = "subclass_id_label")
anchorMeta <- addRefColwithNA(anchorMeta, coltoAdd = "class_id_label",
  colnew = "class_id_label")
anchorMeta <- addRefColwithNA(anchorMeta, coltoAdd = "F",
  colnew = "FemaleRef")

# 2.2.2 pairedtag's meta data
ptMeta <- data.table::fread(
  file = args$pmetafnm,
  sep = ",",
  header = TRUE,
  data.table = FALSE
)
rownames(ptMeta) <- ptMeta$barcode
addQueryColwithNA <- function(m,
                              coltoAdd = "L1",
                              colnew = coltoAdd) {
  m[[colnew]] <- NA
  barcodeQuery <- with(m, barcode[!isRef])
  m[barcodeQuery, colnew] <- ptMeta[barcodeQuery, coltoAdd]
  return(m)
}
anchorMeta <- addQueryColwithNA(anchorMeta, coltoAdd = "L1")
anchorMeta <- addQueryColwithNA(anchorMeta, coltoAdd = "L1_2",
  colnew = "L2")
anchorMeta <- addQueryColwithNA(anchorMeta, coltoAdd = "L1_2_3",
  colnew = "L3")
anchorMeta <- addQueryColwithNA(anchorMeta, coltoAdd = "L1_2_3_4",
  colnew = "L4")
anchorMeta <- addQueryColwithNA(anchorMeta, coltoAdd = "L1_2_3_4_5",
  colnew = "L5")

# add UMAP to the meta
# all(rownames(anchorMeta) == rownames(rdumap))
anchorMeta$UMAP1 <- rdumap[, 1]
anchorMeta$UMAP2 <- rdumap[, 2]

# 3. plot consensus matrix
q <- readRDS(args$qrfnm)
# get predicted cls per barcode
b2pcl <- data.frame(
  barcode = colnames(q),
  cl = paste0("cl-", q@meta.data$predicted.id)
) |> x => `rownames<-`(x, x$barcode)
# get L5-level specific cl annotation based on vote.
b2pcl$L5 <- paste0("L5-", anchorMeta[b2pcl$barcode, "L5"])

get_mode_cl <- function(x) {
  names(which.max(table(x)))
}

# so use vote strategy to map cl to L5
L5tocl <- b2pcl |> group_by(L5) |>
  summarise(cl = get_mode_cl(cl)) |>
  as.data.frame() |>
  x => `rownames<-`(x, x$L5)

# transform the cluster-level annotation to barcode annotation
b2pcl$clvoted <- L5tocl[b2pcl$L5, "cl"]
rlvls <- c("class_id_label",
  "supertype_id_label", "subclass_id_label")
for (i in rlvls) {
  b2pcl[[i]] <- allenMeta[b2pcl$cl, i]
}

# add the annotations to anchorMeta
anchorMeta[rownames(b2pcl), "cl"] <- gsub("cl-", "", b2pcl$cl)
for (i in rlvls) {
  anchorMeta[rownames(b2pcl), i] <- allenMeta[b2pcl$clvoted, i]
}

# get consensus matrix
getAvgVoteMat <- function(rlvl = "cl") {
  rg <- unique(b2pcl[[rlvl]])
  qg <- unique(b2pcl$L5)
  getAvgVote <- function(g) {
    t <- g[[rlvl]]
    r <- rep(0.0, length(rg)) |>
      setNames(object = _, nm = rg)
    tmp <- toNamedArray.1dtable(t)
    r[names(tmp)] <- (tmp / length(t))
    r
  }
  vapply(qg, \(i) {
    g <- b2pcl[b2pcl$L5 == i, ]
    getAvgVote(g)
  }, FUN.VALUE = rep(0.0, length(rg)))
}

tfmat.cl <- getAvgVoteMat(rlvl = "cl")
tfmat.sp <- getAvgVoteMat(rlvl = "supertype_id_label")
tfmat.sc <- getAvgVoteMat(rlvl = "subclass_id_label")

# draw consensus matrix on cl, supertype and subclass
prepareDotPlot4TransferLabel <- function(tfmat,
                                         refOrder = NULL,
                                         names = c("row", "column", "score"),
                                         ignoreEmptyRef = TRUE) {
  maxScore <- getMaxColScoreWithName(mat = tfmat)
  query2ref <- data.frame(
    query = colnames(tfmat),
    ref = names(maxScore),
    row.names = colnames(tfmat)
  )
  if (is.null(refOrder)) {
    message("refOrder is null, will use default numeric order for it.")
    refOrder <- rownames(tfmat)
    refOrder <- refOrder[order(as.integer(refOrder))]
  } else {
    refOrder <- refOrder[refOrder %in% rownames(tfmat)]
  }
  if (ignoreEmptyRef) {
    # message("remove refs not having query mapped to.")
    tfmat <- tfmat[rownames(tfmat) %in% query2ref$ref, ]
  }
  queryOrder <- query2ref$query[
    order(factor(query2ref$ref, levels = refOrder))]

  meltmat <- to3.matrix(tfmat, names)
  meltmat[, 1] <- factor(meltmat[, 1], levels = refOrder)
  meltmat[, 2] <- factor(meltmat[, 2], levels = queryOrder)
  # reduce size of meltmat
  meltmat <- meltmat[meltmat[, 3] > 0, ]
  return(meltmat)
}

plot_cnss_tile <- function(rlvl = "cl", refOrder, tfmat) {
  m <- prepareDotPlot4TransferLabel(
    tfmat = tfmat,
    ignoreEmptyRef = TRUE,
    refOrder = refOrder,
    names = c("Allen", "PairedTag", "score"))

  lowSimScore <- quantile(m$score, 0.1)
  highSimScore <- max(m$score)

  ggplot(data = m, aes(x = PairedTag, y = Allen)) +
    geom_tile_rast(aes(color = score)) +
    mytheme +
    theme(
      panel.background = element_rect(fill = args$lightcolor),
      axis.text.x = element_blank(),
      axis.text.y = element_blank()) +
    scale_color_gradient(
      low = args$lightcolor, high = args$darkcolor,
      limits = c(lowSimScore, highSimScore), na.value = args$lightcolor) +
    xlab(str_glue("Paired-Tag RNA {nlevels(m$PairedTag)} L5")) +
    ylab(str_glue("Allen 10Xv3 RNA {nlevels(m$Allen)} {rlvl}"))
}

raw_rcls <- rownames(tfmat.cl)
clord <- raw_rcls[order(as.integer(gsub("cl-", "", raw_rcls)))]

p_cnss_cl <- plot_cnss_tile(rlvl = "cl", refOrder = clord,
  tfmat = tfmat.cl)

raw_rscs <- rownames(tfmat.sc)
int_rscs <- map_int(raw_rscs, \(i) {as.integer(
  str_split_1(i, " ")[1])
})
scord <- raw_rscs[order(int_rscs)]
p_cnss_sc <- plot_cnss_tile(rlvl = "subclass", refOrder = scord,
  tfmat = tfmat.sc)

raw_rsps <- rownames(tfmat.sp)
int_rsps <- map_int(raw_rsps, \(i) {as.integer(
  str_split_1(i, " ")[1])
})
spord <- raw_rsps[order(int_rsps)]
p_cnss_sp <- plot_cnss_tile(rlvl = "super type", refOrder = spord,
  tfmat = tfmat.sp)

p_cnss <- p_cnss_cl + p_cnss_sp + p_cnss_sc
message("save consensus matrix pdf.")
ggsave(filename = paste0(args$outcnssprefix, ".consensus.pdf"),
  plot = p_cnss, width = 30, height = 10)
message("save consensus matrix result.")
saveRDS(object = list(cl = tfmat.cl, sp = tfmat.sp, sc = tfmat.sc),
  file = paste0(args$outcnssprefix, ".consensus.mat.rds"))

# 4. draw umap
# create a new seurat for easily plot
cnts <- intseu@assays$RNA@counts
colnames(cnts) <- anchorMeta$barcode

aseu <- CreateSeuratObject(
  counts = cnts,
  project = "TransferLabel",
  assay = "RNA",
  names.delim = "_",
  meta.data = anchorMeta
)
aseu[["UMAP"]] <- rdUMAP

class_ignore <- table(aseu$class_id_label[aseu$isRef]) |>
  x => names(x)[x < 50]
cells <- aseu$barcode |>
  x => x[!(aseu$class_id_label %in% class_ignore)]
p_tech_1 <- DimPlot(aseu, reduction = "UMAP", split.by = "tech",
  raster = TRUE, group.by = "tech")
p_tech_2 <- DimPlot(aseu, reduction = "UMAP", raster = TRUE,
  group.by = "tech")
## p_tech <- p_tech_1 + p_tech_2 +
##   plot_layout(design = "AAB", guides = "collect")
subseu <- aseu[ , cells]
p_class <- DimPlot(subseu,
  reduction = "UMAP",
  ## cells = cells,
  split.by = "tech", group.by = "class_id_label", raster = TRUE,
  label = TRUE, label.box = TRUE, label.size = 3, repel = TRUE) +
  theme(legend.position = "none")

design <- "
AAABB
CCCCC
"
p_umap <- p_tech_1 + p_tech_2 + p_class +
  plot_layout(design = design, guides = "collect")
message("save tech umap to ", outumapfig)
ggsave(filename = outumapfig, plot = p_umap, width = 16, height = 16)
message("save seurat with ref and query to ", args$outseufnm)
saveRDS(aseu, file = args$outseufnm)
message("All done. Good luck!")

