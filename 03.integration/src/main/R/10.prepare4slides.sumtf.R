library(Seurat)
library(tidyverse)

projd <- here::here()
workd <- file.path(projd, "03.integration")
mrs <- c("AMY", "CPU", "HYP", "HIP", "ERC", "NAC", "VTA", "PFC")

# * load meta
ptMeta_tfg <- file.path(projd, "meta",
  "pairedtag.cell.meta.all.with.init.tf.csv") |>
  data.table::fread(file = _, sep = ",", header = TRUE, data.table = FALSE) |>
  x => `rownames<-`(x, x$barcode)
ptMeta_tfr <- file.path(projd, "meta",
  "pairedtag.cell.meta.all.with.tfv2.csv") |>
  data.table::fread(file = _, sep = ",", header = TRUE, data.table = FALSE) |>
  x => `rownames<-`(x, x$barcode)
ptMeta_tfg$isNeuFromL1 <- ptMeta_tfr[rownames(ptMeta_tfg), "isNeuFromL1"]

# * load Allen's seurat data
allenSeus <- lapply(mrs, \(r) {
  readRDS(file.path(projd, "data", "allen_seurat",
    paste("allen.10xv3", r, "neu.ds.seu.rds", sep = ".")))
})
names(allenSeus) <- mrs

allenNeuMeta <- lapply(mrs, \(r) {
  m <- allenSeus[[r]]@meta.data
  rownames(m) <- colnames(allenSeus[[r]])
  return(m)
}) |> do.call(rbind, args = _)

n_neu_sp <- allenNeuMeta$supertype_label |>
  unique() |>
  length()
n_neu_sc <- allenNeuMeta$subclass_label |>
  unique() |>
  length()
allenNeuMeta$isNeu <- TRUE

allenNNSeu <- readRDS(
  file.path(projd, "data", "allen_seurat",
    "allen.10xv3.pt.regions.nn.imn.k8.cl.ds1000.rds")
)
allenNNMeta <- allenNNSeu@meta.data
rownames(allenNNMeta) <- colnames(allenNNSeu)
allenNNMeta$isNeu <- FALSE


allenMeta <- rbind(
  allenNeuMeta[ , c("supertype_id_label", "subclass_id_label")],
  allenNNMeta[, c("supertype_id_label", "subclass_id_label")] )

n_allen_sp <- allenMeta$supertype_id_label |>
  unique() |> length()
n_allen_sc <- allenMeta$subclass_id_label |>
  unique() |> length()

# * global transfer label
n_sc <- length(unique(ptMeta_tfg$subclass_id_label))
n_nn_sc <- ptMeta_tfg[!ptMeta_tfg$isNeuFromL1, ] |>
  x => unique(x$subclass_id_label) |>
  x => length(x)
n_neu_sc <- n_sc - n_nn_sc 

n_sp <- length(unique(ptMeta_tfg$supertype_id_label))
n_nn_sp <- ptMeta_tfg[!ptMeta_tfg$isNeuFromL1, ] |>
  x => unique(x$supertype_id_label) |>
  x => length(x)
n_neu_sp <- n_sp - n_nn_sp

n_tfg <- c(n_sp, n_nn_sp, n_neu_sp, n_sc, n_nn_sc, n_neu_sc)

# * region transfer labels
n_sc <- length(unique(ptMeta_tfr$subclass_id_label))
n_nn_sc <- ptMeta_tfr[!ptMeta_tfr$isNeuFromL1, ] |>
  x => unique(x$subclass_id_label) |>
  x => length(x)
n_neu_sc <- n_sc - n_nn_sc 

n_sp <- length(unique(ptMeta_tfr$supertype_id_label))
n_nn_sp <- ptMeta_tfr[!ptMeta_tfr$isNeuFromL1, ] |>
  x => unique(x$supertype_id_label) |>
  x => length(x)
n_neu_sp <- n_sp - n_nn_sp

n_tfr <- c(n_sp, n_nn_sp, n_neu_sp, n_sc, n_nn_sc, n_neu_sc)
