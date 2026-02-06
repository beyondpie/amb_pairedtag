## ** DE analysis on these TFs.
## 1. get seurat object of sampled barcode 2 these TFs
## 2. Wilcox to extract subclass-specific TFs

suppressMessages({
  suppressWarnings({
    library(tidyverse)
    # detach("package:tmpRpkg", unload = T)
    library(tmpRpkg)
    library(GenomicRanges)
    library(S4Vectors)
    library(rtracklayer)
    library(ComplexHeatmap)
    library(reticulate)
    library(Seurat)
    library(Matrix)
  })
})

tmpRpkg::setCondaEnv()
mypy <- import("tmpPypkg")
myad <- import("anndata")

# * meta
projd <- tmpRpkg:::projd
workd <- file.path(projd, "16.celloracle")
outd <- file.path(workd, "out")
deEnheRegulon <- data.table::fread(
  file = file.path(outd, "DAR.enhancer.eRegulon.csv"),
  sep = ",",
  header = T,
  data.table = F
)

ptscMeta <- loadptscMeta() |>
  subset(x = _, subset = ATAC > 0) |>
  x => `rownames<-`(x, x$AllenIdName)

# * function
ann2seu <- function(anntf) {
  
  barcodes <- py_to_r(anntf$obs_names$to_list())
  genes <- py_to_r(anntf$var_names$to_list())
  X <- py_to_r(anntf$X)
  rownames(X) <- barcodes
  colnames(X) <- genes

  # hack obs
  obs <- py_to_r(anntf$obs)
  attr(obs, "pandas.index") <- NULL
  rownames(obs) <- barcodes
  oldsc <- obs$subclass
  newsc <- levels(oldsc)[oldsc]
  obs$subclass <- ptscMeta[levels(oldsc)[oldsc], "PairedTagName"]
  seu <- CreateSeuratObject(
    # gene by cell for Seurat
    counts = t(X),
    project = "Ann2Seu",
    assay = "RNA",
    meta.data = obs
  )
  
  Idents(seu) <- "subclass"
  return(seu)
}

runDE <- function(sc, seutf) {
  tfMarkers <- FindMarkers(
    object = seutf,
    ident.1 = sc,
    test.use = "wilcox", max.cells.per.ident = 20000,
    only.pos = T)
  tfMarkers$tf <- rownames(tfMarkers)
  tfMarkers$subclass <- sc
  return(tfMarkers)
}

# * main
ptannk8 <- myad$read_h5ad(
  filename = file.path(projd, "data", "pairedtag_ann",
    "pt.RNA.ann.k8.ds.raw.h5ad"), backed = NULL)

## k8genes <- py_to_r(ptannk8$var_names$to_list())
## tfs <- unique(deEnheRegulon$tf)
## all(tfs %in% k8genes) # yes
## anntf <- ptannk8[ , tfs]

# * create Seurat object
seuk8 <- ann2seu(ptannk8)
seuk8 <- NormalizeData(seuk8,
  scale.factor = 10000,
  normalization.method = "LogNormalize")
seuk8 <- ScaleData(seuk8, do.scale = T,
  do.center = T, scale.max = 10)
seutf <- seuk8[tfs, ]
saveRDS(seutf, file.path(projd, "16.celloracle", "out",
  "PairedTagRNA.TF.seu.rds"))

# * perform diff test
tfMarkers <- lapply(ptscMeta$PairedTagName, runDE,
  seutf = seutf) |>
  setNames(object = _, nm = ptscMeta$PairedTagName)
tfMarkers <- filterNULLfromList(tfMarkers)
tfMarker <- do.call(rbind, tfMarkers)
data.table::fwrite(x = tfMarker,
  file = file.path(outd, "diff.tf.wilcox.csv"),
  col.names = T, row.names = F)
