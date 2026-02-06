# SnapATAc2 AnnData to Seurat
# - Support downsampling on regions

library(BPCells)
library(Seurat)
library(stringr)
library(purrr)
library(reticulate)
library(Matrix)
library(MatrixExtra)
reticulate::use_condaenv(
  condaenv = "sa2",
  conda = "/home/szu/miniforge3/bin/mamba",
  required = TRUE
)
sa2 <- reticulate::import(module = "snapatac2")

# * configs
projdir <- here::here()
workdir <- file.path(projdir, "03.integration")
outdir <- file.path(projdir, "data/pairedtag_seurat")
dir.create(outdir)

ptmetafnm <- file.path(projdir, "meta",
  "pairedtag.meta.v1.csv")
sa2annfnm <- file.path(projdir, "01.clustering",
  "out/scRNAseq_sa2_all.ann.h5ad")
L1meta <- data.table::fread(
  file = file.path(projdir, "01.clustering", "out/nodlt/L1",
    "RNA_k8_npc50_k50_L0_0.csv"), sep = ",", header = TRUE)
# use r0.2_seed9
data.table::setkey(L1meta, barcode)


# so 0.5 million cells for integration analysis in total
nds <- 250000
sa2annds_fnm <- file.path(outdir,
  "pt_RNA_sa2_nds0.25m_balanceregion.h5ad")
if (file.exists(sa2annds_fnm)) {
  message(sa2annds_fnm, " exists and will be removed.")
  file.remove(sa2annds_fnm)
}

# load pt meta
ptmeta <- data.table::fread(file = ptmetafnm,
  header = TRUE, sep = ",")
data.table::setkey(ptmeta, barcode)

# * load SnapATAC2 AnnData using reticulate
sa2ann <- sa2$read(filename = sa2annfnm, backed = "r")

# * downsampling and keep balance of the region
# for integration
barcodes <- sa2ann$obs_names
annmeta <- ptmeta[barcodes, ]
# treat HCa and HCp as HIP
annmeta$brainregion[
  annmeta$brainregion %in% c("HCa", "HCp")] <- "HIP"
nr <- length(unique(annmeta$brainregion))
ndsper <- ceiling(nds / nr)

subannmeta <- annmeta |>
  dplyr::group_by(brainregion) |>
  dplyr::slice_sample(n = ndsper)
subannmeta <- as.data.frame(subannmeta)
rownames(subannmeta) <- subannmeta$barcode

sa2annds <- sa2ann$subset(
  obs_indices = subannmeta$barcode,
  out = sa2annds_fnm
)
sa2ann$close()
sa2annds$close()

# * save to Seurat v5
outs5dir <- file.path(outdir, "nds0.25m_balanceregion")
outs5matdir <- file.path(outs5dir, "_mat")
dir.create(outs5dir)
dir.create(outs5matdir)
outs5fnm <- file.path(outs5dir,
  "pt_RNA_nds0.25m_balanceregion.rds")
# read into memory so that metadata is readable by R.
sa2annds <- sa2$read(sa2annds_fnm, backed = NULL)
mat <- t(sa2annds$X)
rownames(mat) <- sa2annds$var_names$to_list()
colnames(mat) <- sa2annds$obs_names$to_list()
mat <- MatrixExtra::as.csc.matrix(mat)
mat <- BPCells::convert_matrix_type(mat, type = "uint32_t")
BPCells::write_matrix_dir(mat = mat, dir = outs5matdir,
  overwrite = TRUE)
cellmeta <- subannmeta[sa2annds$obs_names$to_list(), ]
cellmeta$umap1 <- L1meta[cellmeta$barcode, "umap1"]$umap1
cellmeta$umap2 <- L1meta[cellmeta$barcode, "umap2"]$umap2
cellmeta$L1 <- L1meta[cellmeta$barcode, "r0.2_seed9"]$r0.2_seed9

d <- BPCells::open_matrix_dir(outs5matdir)
s5 <- Seurat::CreateSeuratObject(counts = d,
  assay = "RNA", meta.data = cellmeta)
saveRDS(s5, file = outs5fnm)

# * prepare anndata as test sample
library(dplyr)
nds <- 5000
outsa2annfnm <- file.path(workdir, "src/test/resource",
  "sa2ann_test.h5ad")
outsa2meta <- file.path(workdir, "src/test/resource",
  "sa2ann_cellmeta.csv")

sa2ann <- sa2$read(filename = sa2annfnm, backed = "r")
barcodes <- sa2ann$obs_names
annmeta <- ptmeta[barcodes, ]

# treat HCa and HCp as HIP
annmeta$brainregion[
  annmeta$brainregion %in% c("HCa", "HCp")] <- "HIP"
nr <- length(unique(annmeta$brainregion))
ndsper <- ceiling(nds / nr)

dsannmeta <- annmeta |> group_by(brainregion) |> slice_sample(n = ndsper)
dsa2ann <- sa2ann$subset(obs_indices = dsannmeta$barcode,
  out = outsa2annfnm)

data.table::fwrite(dsannmeta,
  file = outsa2meta, col.names = TRUE, row.names = FALSE,
  sep = ",")
