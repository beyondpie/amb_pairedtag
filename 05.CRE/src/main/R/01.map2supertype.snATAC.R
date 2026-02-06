library(Seurat)
library(tidyverse)

# * meta data
projd <- here::here()
workd <- file.path(projd, "05.CRE")
snATACd <- file.path(projd, "data", "snATAC")
snATACtfd <- file.path(snATACd, "snATACTransferLabel")

allenClMeta <- data.table::fread(
  file = file.path(projd, "meta", "AIT21_annotation_freeze_081523.tsv"),
  sep = "\t", header = TRUE, data.table = FALSE
) |>
  x => `rownames<-`(x, x[, 1])

allencl2sp <- unique(allenClMeta[, c("cl", "supertype_id_label")])
rownames(allencl2sp) <- paste0("cl-", allencl2sp$cl)

snATACMeta <- readRDS(file.path(snATACd, "mba.whole.cell.meta.v9.7.rds"))
cl2L4mat <- readRDS(file.path(snATACtfd, "sa2.all.cl2L4.mat.rds"))

# * function
getTopSuperType <- function(l4) {
  cl <- rownames(cl2L4mat)[which.max(cl2L4mat[, l4])]
  allencl2sp[paste0("cl-", cl), "supertype_id_label"]
}

# * Main
# * get nn L4s
L4s <- unique(snATACMeta$L4)
nnL4s <- unique(snATACMeta$L4[snATACMeta$NT_v3 == "NN"])
nnL4tosp <- data.frame(
  L4 = nnL4s,
  sp = vapply(nnL4s, getTopSuperType, "001")
)


nnL4toscFromMeta <- unique(snATACMeta[
  snATACMeta$L4 %in% nnL4s, c("L4", "subclass_id_label_v3")])
rownames(nnL4toscFromMeta) <- nnL4toscFromMeta$L4

nnL4tosp$sc <- nnL4toscFromMeta[nnL4tosp$L4, "subclass_id_label_v3"]
# only 17-5-1-1 is different with subclass_id_label_v3
# https://docs.google.com/spreadsheets/d/1PRaBDlMuFdMAUrcahdyVrIuk0Pwj3Yx2_1DG3NHFMOg/edit?gid=1426252395#gid=1426252395
# from here: we follow subclass id label and correct tf label
# based on the annotation of 17-5-2-1, we label it as VLMC NN_1
nnL4tosp["17-5-1-1", "sp"] <- "1187 VLMC NN_1"
data.table::fwrite(
  x = nnL4tosp,
  file = file.path(snATACtfd, "snATAC_nnL4tosp_240814.csv"),
  col.names = TRUE,
  row.names = FALSE
)


