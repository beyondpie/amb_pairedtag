## * Master TF
## Diff TF + eRegulon + ChromHMM-ChrA + DECRE

suppressMessages({
  suppressWarnings({
    library(tidyverse)
    # detach("package:tmpRpkg", unload = T)
    library(tmpRpkg)
    library(ComplexHeatmap)
    library(reticulate)
    library(Seurat)
  })
})
tmpRpkg::setCondaEnv()
mypy <- import("tmpPypkg")

# * meta
projd <- tmpRpkg:::projd
workd <- file.path(projd, "16.celloracle")
figd <- file.path(projd, "99.figures/out/fig4")
outd <- figd
ptscMetafnm <- file.path(
  projd, "meta", "PairedTagSubclassMetaFromCellMetaChromHMM.csv")
ptscMeta <- data.table::fread(
  file = ptscMetafnm,
  sep = ",", header = TRUE, data.table = FALSE
) |>
  x => `rownames<-`(x, x$PairedTagName) |>
  x => x[x$ATAC >0, ]
sc2nmod <- tmpRpkg::getPairedTagSubclass2Nmod(ptscMeta)
sc2n <- table(seuTF@meta.data$subclass)
ptscs <- with(sc2nmod, sc[H3K27ac >= 200 & H3K4me1 >= 200])

# * functions

# * main
# 1. get master TF
# - load diff TF
diffTF <- data.table::fread(
  file = file.path(projd, "16.celloracle", "out", "diff.tf.wilcox.csv"),
  sep = ",",
  header = T,
  data.table = F)

# - load diff enhancers at eRegulon
deEnheRegulon <- data.table::fread(
  file = file.path(projd, "16.celloracle", "out", "DAR.enhancer.eRegulon.csv"),
  sep = ",",
  header = T,
  data.table = F
)

# - intersection and check results
scSpecTF <- merge(diffTF, deEnheRegulon,
  by.x = c("subclass", "tf"),
  by.y = c("subclass", "tf")
  )

# - remove low-number subclass
# 2. heatmap of master TF
# - load PairedTag RNA
ptRNA_logCPM_scbyg <- mypy$globalvar$load_PairedTag_RNA_logCPM_scbyg()
ptRNAlogCPM <- py_to_r(ptRNA_logCPM_scbyg)
attr(ptRNAlogCPM, "pandas.index") <- NULL

sc2diffTFmat <- ptRNAlogCPM[unique(scSpecTF$subclass), unique(scSpecTF$tf)] |>
  as.matrix() |>
  scale()

saveRDS(
  object = sc2diffTFmat,
  file = file.path(projd, "16.celloracle", "out",
    "subclassSpecEnhTF.pairedTagRNA.scaledlogCPM.rds")
)

# - ord by subclass and TFs
seuTF <- readRDS(file.path(projd, "16.celloracle",
  "out",
  "PairedTagRNA.TF.seu.rds"))

TFOrd <- lapply(ptscs, \(sc) {
  m <- scSpecTF[scSpecTF$subclass == sc, ]
  tfs <- m$tf
  fc <- m[, "avg_log2FC"]
  index <- order(fc, decreasing = T)
  tfs[index[seq_len(min(length(index), 10))]]
}) |>
  unlist() |>
  unique()

sc2tfmat <- ptRNAlogCPM[ptscs, TFOrd] |>
  as.matrix() |>
  scale()

saveRDS(
  object = sc2tfmat,
  file = file.path(outd,
    "panelE.subclassSpecEnhTF.pairedTagRNA.scaledlogCPM.rds")
)

