library(tidyverse)

projd <- here::here()
workd <- file.path(projd, "16.celloracle")
outd <- file.path(workd, "out")

markerd <- file.path(outd, "subclass_marker_genes_csvs")

# * load meta data
# load subclass-level marker gene analysis
markerscs <- list.files(markerd, full.names = F) |>
  x => gsub("_marker_genes\\.csv", "", x)

allDEGenes <- lapply(markerscs, \(sc) {
  data.table::fread(
    file = file.path(markerd, str_glue("{sc}_marker_genes.csv")),
    header = T, sep = ",", data.table = F
  )
}) |>
  x => setNames(object = x, nm = markerscs)

# filter genes to get subclasses' markers
padj <- 0.01
log2fc <- 1

markers <- lapply(markerscs, \(sc) {
  y <- allDEGenes[[sc]] |>
    x => with(x, gene[p_val_adj <= padj & avg_log2FC >= log2fc])
  if (length(y) < 1) {
    NULL
  } else {
    y
  }
}) |>
  x => setNames(object = x, nm = gsub("-", "_", markerscs)) |>
  x => Filter(Negate(is.null), x)

# all eRegulon
eRegulons <- data.table::fread(
  file = file.path(outd, "all.eRegulon.csv"),
  header = T, data.table = F, sep = ","
)

eRegscs <- unique(eRegulons$subclass)

ptscs <- intersect(eRegscs, names(markers))

# for each subclass, check DE overlap with
# ChrA/ChrO-marked TF-downstream genes

r <- lapply(ptscs, \(sc) {
  x <- eRegulons[(eRegulons$subclass == sc) & (eRegulons$coef > 0.0),
    c("tf", "gene", "label")]
  tfChrA <- with(x, tf[label == "enhancer"]) |>
    unique()
  tfChrO <- with(x, tf[label == "ChrO"]) |>
    unique() |>
    setdiff(x = _, y = tfChrA)
  
  gChrA <- with(x, gene[tf %in% tfChrA]) |> unique()
  gChrO <- with(x, gene[tf %in% tfChrO]) |> unique()

  g2ChrA <- with(x, gene[label == "enhancer"]) |>
    unique()
  g2ChrO <- with(x, gene[label == "ChrO"]) |>
    unique() 
  g3ChrO <- with(x, gene[label == "ChrO"]) |>
    unique()  |>
    setdiff(x = _, y = g2ChrA)
  
  data.frame(
    nMarkerChrA = length(intersect(gChrA, markers[[sc]])),
    nMarkerChrO = length(intersect(gChrO, markers[[sc]])),
    nMarker2ChrA = length(intersect(g2ChrA, markers[[sc]])),
    nMarker2ChrO = length(intersect(g2ChrO, markers[[sc]])),
    nMarker3ChrO = length(intersect(g3ChrO, markers[[sc]])),
    sc = sc
  )
}) |> do.call(what = rbind, args = _)

saveRDS(object = r,
  file = file.path(outd, "ChrA-ChrO-TF-gene_intersect_markers.rds"))

# * check if ChrA enriched TF have more centralirities than ChrO marked.
# 2026-01-01
# all eRegulon
eRegulons <- data.table::fread(
  file = file.path(outd, "all.eRegulon.csv"),
  header = T, data.table = F, sep = ","
)

eRegscs <- unique(eRegulons$subclass)

ptscs <- intersect(eRegscs, names(markers))
tfscores <- read.csv(
  file = file.path(workd, "out", "TFCentraScore.txt"),
  header = T,
  row.names = 1, sep = "\t", as.is = T, check.names = F
)

ptscNameMeta <- data.table::fread(
  file = file.path(projd, "meta",
    "PairedTagSubclassMetaFromCellMetaChromHMM.csv"),
  sep = ",", header = T, data.table = F
) |>
  mutate(.data = _, ATACName2 = gsub("^\\d+_", "", ATACName)) |>
  x => `rownames<-`(x, x$ATACName2) |>
  x => x[x$ATAC > 0, ]

# TRUE
all(ptscNameMeta$ATACName2 %in% colnames(tfscores))

tfs <- tfscores[, ptscNameMeta$ATACName2] |>
  setNames(object = _, nm = ptscNameMeta$PairedTagName)
# TRUE
all(eRegscs %in% colnames(tfs))

# assign TF is ChrA or ChrO
nAlltf2sc <- matrix(0.0,
  nrow = length(rownames(tfs)),
  ncol = length(colnames(tfs)),
  dimnames = list(
    rownames(tfs),
    colnames(tfs)
  ))
ChrAtf2sc <- nAlltf2sc
ChrOtf2sc <- nAlltf2sc

for(sc in colnames(tfs)) {
  x <- eRegulons[(eRegulons$subclass == sc) & (eRegulons$coef > 0.0),
    c("gene", "label")]
  gChrA <- with(x, gene[label == "enhancer"]) |>
    unique()
  
  gChrO <- with(x, gene[label == "ChrO"]) |>
    unique() |>
    setdiff(x = _, y = gChrA)
}

alltf <- rownames(tfs)
tfs2sc <- lapply(colnames(tfs), \(sc) {
  x <- eRegulons[(eRegulons$subclass == sc) & (eRegulons$coef > 0.0),
    c("tf", "label")]
  gChrA <- with(x, tf[label == "enhancer"]) |>
    unique()
  
  gChrO <- with(x, tf[label == "ChrO"]) |>
    unique() |>
    setdiff(x = _, y = gChrA)

  avgChrA <- tfs[alltf %in% gChrA, sc] |>
    mean()
  avgChrO <- tfs[alltf %in% gChrO, sc] |>
    mean()

  data.frame(
    tfsChrA = avgChrA,
    tfsChrO = avgChrO,
    sc = sc
  )
}) |>
  x => do.call(what = rbind , args = x)
saveRDS(object = tfs2sc,
  file = file.path(outd, "tfcentralirity_ChrAvsChrO.rds"))





