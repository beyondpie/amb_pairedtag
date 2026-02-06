library(tidyverse)

projd <- here::here()

cellMeta <- file.path(projd, "meta",
  "pairedtag.cell.meta.all.240626.csv") |>
  data.table::fread(file = _, head = TRUE, sep = ",",
    data.table = FALSE)
rownames(cellMeta) <- cellMeta$barcode

# * compared with the global annotations
cellMetaGlobal <- file.path(projd, "meta",
  "pairedtag.cell.meta.all.with.init.tf.csv") |>
  data.table::fread(file = _, head = TRUE, sep = ",",
    data.table = FALSE) |>
  x => `rownames<-`(x, x$barcode)
cellMetaGlobal <- cellMetaGlobal[rownames(cellMeta), ]

all(cellMetaGlobal$nCount_RNA == cellMeta$mapping.nRNA)
all(cellMetaGlobal$sublib == cellMeta$exp.sublib)
all(cellMetaGlobal$modularity == cellMeta$exp.modularity)
all(cellMetaGlobal$L1UMAP1 == cellMeta$cluster.l1.umap.x)
all(cellMetaGlobal$supertype_id_label == cellMeta$gannot.sp)

# * comapred with region annotations
cellMetaRegion <- file.path(projd, "meta",
  "pairedtag.cell.meta.all.with.tfv3.csv") |>
  data.table::fread(file = _, head = TRUE, sep = ",",
    data.table = FALSE) |>
  x => `rownames<-`(x, x$barcode) |>
  x => x[rownames(cellMeta), ]

barcodes <- with(cellMeta, barcode[!is.na(rannot.sp)])
all(cellMetaRegion[barcodes, "supertype_id_label"] ==
      cellMeta[barcodes, "rannot.sp"])

sum(cellMetaRegion$supertype_id_label == cellMeta$annot.sp,
  na.rm = TRUE)

barcodes <- with(cellMeta, barcode[annotQuality == "Good"])
all(cellMetaRegion[barcodes, "supertype_id_label"] ==
      cellMeta[barcodes, "annot.sp"])

all(cellMetaRegion$isNeuFromL1 == cellMeta$isNeuL1)

# * check LQ
barcodes <- with(cellMetaRegion, barcode[grepl("LQ", L5r)])
table(cellMeta[barcodes, "annotQuality"])

cellCheck <- file.path(projd, "03.integration", "out",
  "tfneu_vf_region_cca_k5",
  "amb_PT_allcellmetadata.ZW.20240614.csv"
  ) |>
  data.table::fread(file = _, header = TRUE, sep = ",",
    data.table = FALSE) |>
  x => `rownames<-`(x, x$barcode) 
barcodes <- with(cellCheck, V1[umap_filter != "keep"])
table(cellMeta[barcodes, "annotQuality"])

# * check manual HQ
ERC_enrich <- file.path(projd, "03.integration",
  "out", "tfneu_vf_region_cca_k5", "un_coembed",
  "ERC.enrich.L5r.with.annot.txt") |>
  data.table::fread(file = _, header = FALSE,
    sep = ",") |>
  setNames(object = _, nm = c("L5r", "n", "sp"))

CPU_enrich <- file.path(projd, "03.integration",
  "out", "tfneu_vf_region_cca_k5", "un_coembed",
  "CPU.enrich.L5r.with.annot.txt") |>
  data.table::fread(file = _, header = FALSE,
    sep = ",") |>
  setNames(object = _, nm = c("L5r", "n", "sp"))

HQ_enrich <- rbind(ERC_enrich, CPU_enrich)
rownames(HQ_enrich) <- HQ_enrich$L5r

barcodes <- with(cellMeta, barcode[l5r %in% HQ_enrich$L5r])
table(cellMeta[barcodes, "annotQuality"])

sum(is.na(cellMeta[cellMeta$annotQuality != "Good", "annot.sp"]))
sum(is.na(cellMeta[cellMeta$annotQuality != "Good", "l5r"]))
sum(is.na(cellMeta[cellMeta$annotQuality == "Good", "annot.sp"]))
length(unique(cellMeta[cellMeta$annotQuality == "Good", "annot.sp"]))
length(unique(cellMeta[cellMeta$annotQuality == "Good", "annot.sc"]))
length(unique(cellMeta[cellMeta$annotQuality == "Good", "annot.c"]))

length(unique(cellMeta[(cellMeta$annotQuality == "Good") & (!cellMeta$isNeuL1),
  "annot.sp"]))

a <- table(cellMeta[barcodes, c("l5r", "annot.sp")]) |>
  as.data.frame(stringsAsFactors = FALSE) |>
  reshape2::melt() |>
  subset(x = _, subset = value > 0)

rownames(a) <- a$L5r
all(a[rownames(HQ_enrich), "sp"] == HQ_enrich$sp)


