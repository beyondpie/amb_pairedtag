library(tidyverse)
library(ggpubr)

# * config
projroot <- here::here()
clrec <- data.table::fread(
  file = file.path(projroot, "meta", "pt.scRNAseq.sa2clustering.record.csv"),
  sep = ",", header = T, data.table = F
)

cellMeta <- data.table::fread(
  file = file.path(projroot, "meta", "pt.barcode.meta.withL2.csv"),
  sep = ",", header = T, data.table = F
)
rownames(cellMeta) <- cellMeta$barcode

L2cls <- paste0("c", unique(cellMeta$L1_2))

L3resdir <- file.path(
  projroot, "01.clustering",
  "out/afqc/L3"
)
resprefix <- "RNA_vg2000_npc30_k30_L1_2_"

L3leidens <- lapply(L2cls, \(x) {
  data.table::fread(
    file = file.path(
      L3resdir,
      paste0(resprefix, x, ".leiden.csv")
    ),
    sep = ",", header = T
  )
}) |> do.call(what = rbind, args = _)
rownames(L3leidens) <- L3leidens$barcode

L3silhts <- lapply(L2cls, \(x) {
  data.table::fread(file = file.path(
    L3resdir,
    paste0(resprefix, x, ".silht.csv")
  ), sep = ",", header = T)
})
names(L3silhts) <- L2cls

L3cellmetas <- lapply(L2cls, \(x) {
  data.table::fread(file = file.path(
    L3resdir,
    paste0(resprefix, x, ".cellmeta.leiden.csv")
  ), sep = ",", header = T)
})
names(L3cellmetas) <- L2cls

# * get best reso and its seed for each cluster
seeds <- seq_len(10)-1
rs <- seq(0.1,2,0.1)
getAvgSil <- function(r, m) {
  r <- as.numeric(r)
  cols <- paste0(
    stringr::str_glue("r{format(r, nsmall = 1)}_"), "seed", seeds
  )
  mean(colMeans(m[, cols]))
}

rrssncl <- lapply(L2cls, \(x) {
  m <- L3silhts[[x]] |> as.data.frame()
  max_r <- rs[which.max(map_vec(rs, getAvgSil, m = m))]
  max_r <- as.numeric(max_r)
  cols <- paste0(
    str_glue("r{format(max_r, nsmall=1)}_"),
    "seed", seeds
  )
  mm <- m[, cols]
  max_seed <- seeds[which.max(colMeans(mm))]
  t <- str_glue("r{max_r}_seed{max_seed}")
  ncl <- length(table(L3cellmetas[[x]][, ..t]))
  data.frame(r = max_r, seed = max_seed, ncl = ncl)
}) |> do.call(rbind, args = _)
rownames(rrssncl) <- L2cls

L3clrec <- data.frame(
  clustering_level = "L2",
  cluster = as.integer(gsub("c", "", L2cls)),
  resolution = rrssncl$r,
  seed = rrssncl$seed,
  ncl = rrssncl$ncl,
  use_auto = TRUE
) |> x => x[order(L3clrec$cluster), ]
rownames(L3clrec) <- paste0("c", L3clrec$cluster)


newclrec <- rbind(clrec, L3clrec)
data.table::fwrite(newclrec,
  file = file.path(projroot, "meta",
    "pt.scRNAseq.sa2clustering.record.csv"),
  sep = ",",
  row.names = F, col.names = T
)

# update cellmeta and prepare L4 level clustering
allL3UMAP <- lapply(L3cellmetas, \(x) {
  as.data.frame(x) |>
    x => x[, c("barcode", "umap1", "umap2")] |>
    rename(.data = _, L3_UMAP1 = umap1, L3_UMAP2 = umap2)
}) |>
    do.call(what = rbind, args = _)
rownames(allL3UMAP) <- allL3UMAP$barcode
allL3 <- lapply(names(L3cellmetas), \(nm) {
  x <- as.data.frame(L3cellmetas[[nm]])
  r <- L3clrec[nm, "resolution"]
  s <- L3clrec[nm, "seed"]
  t <- str_glue("r{r}_seed{s}")
  data.frame(
    barcode = x$barcode,
    L3 = x[[t]]
  )
}) |> do.call(what = rbind, args = _)
rownames(allL3) <- allL3$barcode
cellMeta$L3_UMAP1 <- allL3UMAP[cellMeta$barcode, "L3_UMAP1"]
cellMeta$L3_UMAP2 <- allL3UMAP[cellMeta$barcode, "L3_UMAP2"]
cellMeta$L3 <- allL3[cellMeta$barcode, "L3"]
cellMeta$L1_2_3 <- with(cellMeta,
  L1 * 10000 + L2 * 100 + L3)
data.table::fwrite(
  cellMeta,
  file = file.path(projroot, "meta", "pt.barcode.meta.withL3.csv"),
  sep = ",",
  row.names = FALSE, col.names = TRUE
)


