library(tidyverse)
library(Seurat)

# * config
# home_tscc <- "/tscc/nfs/home/szu"
home_mediator <- "/home/szu"

# conda <- file.path(home_tscc, "miniforge3/bin/conda")
conda <- file.path(home_mediator, "miniforge3/bin/conda")

# sys_path <- "/tscc/projects/ps-renlab2"
sys_path <- "/projects/ps-renlab2"
projd <- here::here()
mrs <- c("AMY", "CPU", "HYP", "HIP", "ERC", "NAC", "VTA", "PFC")

# * file meta
allen_seurat_dir <- file.path(
  sys_path, "zhw063/99.MouseBrainPairedTag",
  "allen_brain_atlas/allen_ref_seurat_object"
)

allen_k8 <- data.table::fread(
  file = file.path(projd, "meta", "AIT21_k8_markers.txt"),
  header = FALSE, data.table = FALSE
)$V1

pt_ann_dir <- file.path(
  projd, "data", "pairedtag_ann")


# * load allen neuronal seurat objects
loadAllenNeuSeurat <- function(mr) {
  seu <- readRDS(file.path(
    allen_seurat_dir,
    str_glue("20240522.{mr}.clean.allen.rna.rds")
  ))
  mc <- seu@meta.data$manualcheck
  table(mc) |>
    x => paste(names(x), x, sep = ":", collapse = " ") |>
    x => message(x)
  cell_class <- vapply(
    seu@meta.data$subclass_label,
    \(nm) {
      str_split_1(nm, pattern = " ") |>
        tail(x = _, n = 1)
    }, "NN")
  table(cell_class) |>
    x => paste(names(x), x, sep = ":", collapse = " ") |>
    x => message(x)
  k8 <- intersect(allen_k8, rownames(seu))
  seu <- seu[k8, (mc != "discard") & (!grepl("LQ|NN", cell_class))]
  message(str_glue("{mr}: {ncol(seu)} neuronal cells."))
  seu[["RNA"]]$data <- NULL
  seu[["RNA"]]$scale.data <- NULL
  
  outf <- file.path(
    projd, "data", "allen_seurat",
    str_glue("allen.10xv3.{mr}.neu.clean.all.rds")
  )
  message("saved to: ", outf)
  saveRDS(seu, outf)
  return(seu)
}

mrs <- c("AMY", "CPU", "HYP", "HIP", "ERC", "NAC", "VTA", "PFC")
allenNeuSeus <- lapply(mrs, loadAllenNeuSeurat) |>
  x => setNames(object = x, nm = mrs)

# * load pairedtag neuronal cells
# seurat objects are prepared in python script

# * downsample neuronal cells for each major regions
projd <- here::here()
mrs <- c("AMY", "CPU", "HYP", "HIP", "ERC", "NAC", "VTA", "PFC")
## ptMeta <- data.table::fread(
##   file = file.path(projd, "meta",
##   "pairedtag.cell.meta.all.with.init.tf.csv"),
##   sep = ",", header = TRUE, data.table = FALSE) |>
##   x => `rownames<-`(x, x$barcode)

allen_dir <- file.path(projd, "data", "allen_seurat")
pt_dir <- file.path(projd, "data", "pairedtag_seurat",
  "neu_seu_region")

allenSeus <- lapply(mrs, \(r) {
  readRDS(file.path(allen_dir, str_glue("allen.10xv3.{r}.neu.clean.all.rds")))
}) |>
  setNames(object = _, nm = mrs)

ptSeus <- lapply(mrs, \(r){
  readRDS(file.path(pt_dir, str_glue("pt.neu.{r}.all.seu.rds")))
}) |>
  setNames(object = _, nm = mrs)

# check data
## r <- "AMY"
## allenseu <- allenSeus[[r]]
## ptseu <- ptSeus[[r]]

check_dp_balance <- function(r, n_allen = 100, n_pt = 50) {
  alm <- allenSeus[[r]]@meta.data
  ds_allen <- alm |>
    group_by(cl) |>
    slice_sample(n = n_allen)
  ncl <- length(unique(alm$cl))
  message("Brain Region: ", r)
  message(str_glue("allen: {ncl} cls, {nrow(alm)} cells."))
  message(str_glue("after downsampling: {nrow(ds_allen)} cells."))
  message("----------------------")
  
  ptm <- ptSeus[[r]]@meta.data
  nL5 <- length(unique(ptm$L1_2_3_4_5))
  ds_pt <- ptm |>
    group_by(L1_2_3_4_5) |>
    slice_sample(n = n_pt)
  message(str_glue("pt: {nL5} L5s, {nrow(ptm)} cells."))
  message(str_glue("after downsampling: {nrow(ds_pt)} cells."))
  message("======================")
}

na <- 100
np <- 50
for (r in mrs) {
  check_dp_balance(r, n_allen = na, n_pt = np)
}

num4ds <- data.frame(
  region = mrs,
  nallen = c(250, 1500, 150, 1000, 700, 450, 200, 800),
  npt = c(60, 20, 300, 30, 50, 20, 10, 650)
) |>
  x => `rownames<-`(x, x$region)

for (r in mrs) {
  check_dp_balance(r,
    n_allen = num4ds[r, "nallen"],
    n_pt = num4ds[r, "npt"]
  )
}

data.table::fwrite(
  x = num4ds, file = file.path(projd, "03.integration/out",
    "int.region.specific.downsample.nums.csv"),
  sep = ",", row.names = FALSE, col.names = TRUE)

# * perform downsample on the two datasets.
downsample_seu <- function(r, isallen = TRUE) {
  seu <- if (isallen) {
    allenSeus[[r]]
  } else {
    ptSeus[[r]]
  }
  outfnm <- if(isallen) {
    file.path(projd, "data", "allen_seurat",
      str_glue("allen.10xv3.{r}.neu.ds.seu.rds"))
  } else {
    file.path(projd, "data", "pairedtag_seurat",
      "neu_seu_region", str_glue("pt.neu.{r}.ds.seu.rds"))
  }
  nds <- ifelse(isallen, num4ds[r, "nallen"], num4ds[r, "npt"])
  bygroup <- ifelse(isallen, "cl", "L1_2_3_4_5")
  message(
    str_glue("nds {nds} under {bygroup}."))
  seu_meta <- seu@meta.data
  seu_meta$barcode <- colnames(seu)
  barcodes <- seu_meta |>
    group_by(across(bygroup)) |>
    slice_sample(n = nds) |>
    x => x$barcode
  sub_seu <- subset(seu, cells = barcodes)
  message(str_glue("output file: {outfnm}."))
  saveRDS(sub_seu, file = outfnm)
}

for (r in mrs) {
  message("downsample for allen.")
  downsample_seu(r, isallen = TRUE)
  message("downsample for pairedtag.")
  downsample_seu(r, isallen = FALSE)
}

# * check VTA regions
aclmeta <- data.table::fread(
  file = file.join(projd, "meta", ""))
seu <- readRDS(file.path(
    allen_seurat_dir,
    "20240522.VTA.clean.allen.rna.rds"
))
mc <- seu@meta.data$manualcheck
table(mc) |>
  x => paste(names(x), x, sep = ":", collapse = " ") |>
  x => message(x)
cell_class <- vapply(
  seu@meta.data$subclass_label,
  \(nm) {
    str_split_1(nm, pattern = " ") |>
      tail(x = _, n = 1)
  }, "NN")
table(cell_class) |>
  x => paste(names(x), x, sep = ":", collapse = " ") |>
  x => message(x)

# * check pairted meta
ptmeta <- data.table::fread(
  file.path(projd, "meta",  "pairedtag.cell.meta.all.csv"),
  sep = ",", header = TRUE, data.table = FALSE)

# * Update AMY part
seu <- readRDS(file.path(allen_seurat_dir,
  "new.AMY.allen.rna.rds"))
seu <- subset(seu, method == "10Xv3")
mc <- seu@meta.data$manualcheck
table(mc) |>
  x => paste(names(x), x, sep = ":", collapse = " ") |>
  x => message(x)
cell_class <- vapply(
  seu@meta.data$subclass_label,
  \(nm) {
    str_split_1(nm, pattern = " ") |>
      tail(x = _, n = 1)
  }, "NN")
table(cell_class) |>
  x => paste(names(x), x, sep = ":", collapse = " ") |>
  x => message(x)
k8 <- intersect(allen_k8, rownames(seu))
seu <- seu[k8, (mc != "discard") & (!grepl("LQ|NN", cell_class))]
message(str_glue("AMY: {ncol(seu)} neuronal cells."))
seu[["RNA"]]$data <- NULL
seu[["RNA"]]$scale.data <- NULL
outf <- file.path(projd, "data",
  "allen_seurat", "allen.10xv3.AMY.neu.clean.all.rds")
saveRDS(seu, outf)

# now check how to set downsampleing number
pt_dir <- file.path(projd, "data", "pairedtag_seurat",
  "neu_seu_region")
ptSeu <- readRDS(file.path(pt_dir, "pt.neu.AMY.all.seu.rds"))

check_dp_balance <- function(n_allen = 100, n_pt = 50) {
  alm <- seu@meta.data
  ds_allen <- alm |>
    group_by(cl) |>
    slice_sample(n = n_allen)
  ncl <- length(unique(alm$cl))
  message(str_glue("allen: {ncl} cls, {nrow(alm)} cells."))
  message(str_glue("after downsampling: {nrow(ds_allen)} cells."))
  ptm <- ptSeu@meta.data
  nL5 <- length(unique(ptm$L1_2_3_4_5))
  ds_pt <- ptm |>
    group_by(L1_2_3_4_5) |>
    slice_sample(n = n_pt)
  message(str_glue("pt: {nL5} L5s, {nrow(ptm)} cells."))
  message(str_glue("after downsampling: {nrow(ds_pt)} cells."))
  message("======================")
}

check_dp_balance(n_allen = 400, n_pt = 130)
## allen: 958 cls, 106408 cells.
## after downsampling: 76255 cells.
## pt: 1744 L5s, 201552 cells.
## after downsampling: 75503 cells.


downsample_seu <- function(isallen = TRUE) {
  s <- if (isallen) {
    seu
  } else {
    ptSeu
  }
  outfnm <- if(isallen) {
    file.path(projd, "data", "allen_seurat",
      "allen.10xv3.AMY.neu.ds.seu.rds")
  } else {
    file.path(projd, "data", "pairedtag_seurat",
      "neu_seu_region", "pt.neu.AMY.ds.seu.rds")
  }
  nds <- ifelse(isallen, 400, 130)
  bygroup <- ifelse(isallen, "cl", "L1_2_3_4_5")
  message(
    str_glue("nds {nds} under {bygroup}."))
  seu_meta <- s@meta.data
  seu_meta$barcode <- colnames(s)
  barcodes <- seu_meta |>
    group_by(across(bygroup)) |>
    slice_sample(n = nds) |>
    x => x$barcode
  sub_seu <- subset(s, cells = barcodes)
  message(str_glue("output file: {outfnm}."))
  saveRDS(sub_seu, file = outfnm)
}

downsample_seu(isallen = TRUE)
downsample_seu(isallen = FALSE)

# check AMY
a <- readRDS(file.path(projd, "data", "allen_seurat",
  "allen.10xv3.AMY.neu.ds.seu.rds"))


