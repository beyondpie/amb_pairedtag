library(ggplot2)
library(tidyverse)

projroot <- here::here()
workdir <- file.path(projroot, "03.integration")
cellMeta <- data.table::fread(
  file = file.path(projroot, "meta", "pt.barcode.meta.withL4.csv"),
  sep = ",", header = T)
c11barcodes <- data.table::fread(
  file = file.path(projroot, "meta", "lowrep_L1cl-11.barcode.txt"),
  header = F, data.table = FALSE
) |>
  x => x$V1
allenMeta <- data.table::fread(
  file = file.path(projroot, "meta",
    "AIT21_annotation_freeze_081523.tsv"),
  header = TRUE, sep = "\t", data.table = FALSE
)
rownames(allenMeta) <- paste0("c", allenMeta$cl)

# * plot doublet rate per experiment
# * check cluster 11 annotation with integration
tfLabels <- readRDS(file.path(workdir, "out",
  "all_k8_cca_50.balance-pt-region.tfLabel.cl.rds"))
b2cl <- data.frame(
  barcode = rownames(tfLabels@meta.data),
  cl = paste0("c", tfLabels@meta.data$predicted.id)
) |>
  x => `rownames<-`(x, x$barcode)

c112cl <- b2cl[intersect(c11barcodes, b2cl$barcode), "cl"]
c112class <- allenMeta[c112cl, "class_id_label"]
c112subclass <- allenMeta[c112cl, "subclass_id_label"]
a <- table(c112subclass) |>
  as.data.frame(stringsAsFactors = F) |>
  x => x[order(x$Freq, decreasing = T), ]
  
