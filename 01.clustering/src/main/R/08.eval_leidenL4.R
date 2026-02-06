library(tidyverse)
library(ggpubr)
# inhibit scientific notation
# for cluster label.
options(scipen=999)

# * config
nmin_cluster <- 50
projroot <- here::here()
clrec <- data.table::fread(
  file = file.path(projroot, "meta",
    "pt.scRNAseq.sa2clustering.record.csv"),
  sep = ",", header = T, data.table = F
)
cellMeta <- data.table::fread(
  file = file.path(
    projroot, "meta",
    "pt.barcode.meta.withL3.csv"
  ),
  sep = ",", header = T, data.table = F
) |>
  x => `rownames<-`(x, x$barcode)

allL3cls <- paste0("c", unique(cellMeta$L1_2_3))
L3size <- table(cellMeta$L1_2_3) |>
  as.data.frame(
    stringsAsFactors = FALSE
  ) |>
  x => `rownames<-`(x, x$Var) |>
  setNames(nm = c("L3", "size"))

L3cls <- paste0("c",
  with(L3size, L3[size >= nmin_cluster]))
L4resdir <- file.path(
  projroot, "01.clustering", "out/afqc/L4"
)
resprefix <- "RNA_vg2000_npc30_k30_L1_2_3_"

L4leidens <- lapply(L3cls, \(x){
  data.table::fread(
    file = file.path(
      L4resdir,
      paste0(resprefix, x, ".leiden.csv")
    ),
    sep = ",", header = T
  )
}) |> do.call(what = rbind, args= _) |>
  x => `rownames<-`(x, x$barcode)

L4silhts <- lapply(L3cls, \(x) {
  data.table::fread(file = file.path(
    L4resdir,
    paste0(resprefix, x, ".silht.csv")
  ), sep = ",", header = T)
}) |> setNames(nm = L3cls)

L4cellmetas <- lapply(L3cls, \(x) {
  data.table::fread(file = file.path(
    L4resdir,
    paste0(resprefix, x, ".cellmeta.leiden.csv")
  ), sep = ",", header = T)
})
names(L4cellmetas) <- L3cls

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

rrssncl <- lapply(L3cls, \(x) {
  m <- L4silhts[[x]] |> as.data.frame()
  max_r <- rs[which.max(map_vec(rs, getAvgSil, m = m))]
  max_r <- as.numeric(max_r)
  cols <- paste0(
    str_glue("r{format(max_r, nsmall=1)}_"),
    "seed", seeds
  )
  mm <- m[, cols]
  max_seed <- seeds[which.max(colMeans(mm))]
  t <- str_glue(
    "r{format(max_r, nsmall = 1)}_seed{max_seed}")
  ncl <- length(table(L4cellmetas[[x]][, ..t]))
  data.frame(r = max_r, seed = max_seed, ncl = ncl)
}) |>
  do.call(rbind, args = _) |>
  x => `rownames<-`(x, L3cls)

L4clrec <- data.frame(
  clustering_level = "L3",
  cluster = as.integer(gsub("c", "", L3cls)),
  resolution = rrssncl$r,
  seed = rrssncl$seed,
  ncl = rrssncl$ncl,
  use_auto = TRUE
) |>
  x => x[order(x$cluster), ] |>
  x => `rownames<-`(x, paste0("c", x$cluster))



newclrec <- rbind(clrec, L4clrec)
data.table::fwrite(newclrec,
  file = file.path(projroot, "meta",
    "pt.scRNAseq.sa2clustering.record.csv"),
  sep = ",",
  row.names = F, col.names = T
)

# update cellmeta and prepare L4 level clustering
allL4UMAP <- lapply(L4cellmetas, \(x) {
  as.data.frame(x) |>
    x => x[, c("barcode", "umap1", "umap2")] |>
    rename(.data = _, L4_UMAP1 = umap1, L4_UMAP2 = umap2)
}) |>
  do.call(what = rbind, args = _) |>
  x => `rownames<-`(x, x$barcode)

allL4 <- lapply(names(L4cellmetas), \(nm) {
  x <- as.data.frame(L4cellmetas[[nm]])
  r <- as.numeric(L4clrec[nm, "resolution"])
  s <- L4clrec[nm, "seed"]
  t <- str_glue("r{format(r, nsmall = 1)}_seed{s}")
  data.frame(
    barcode = x$barcode,
    L4 = x[[t]]
  )
}) |>
  do.call(what = rbind, args = _) |>
  x => `rownames<-`(x, x$barcode)

cellMeta$L4_UMAP1 <- NA
cellMeta$L4_UMAP2 <- NA
cellMeta$L4 <- 0
bs <- intersect(cellMeta$barcode, rownames(allL4UMAP))

cellMeta[bs, "L4_UMAP1"] <- allL4UMAP[bs, "L4_UMAP1"]
cellMeta[bs, "L4_UMAP2"] <- allL4UMAP[bs, "L4_UMAP2"]
cellMeta[bs, "L4"] <- allL4[bs, "L4"]
cellMeta$L1_2_3_4 <- with(cellMeta,
  L1 * 1000000 + L2 * 10000 + L3 * 100 + L4)
  
data.table::fwrite(
  cellMeta,
  file = file.path(projroot, "meta", "pt.barcode.meta.withL4.csv"),
  sep = ",",
  row.names = FALSE, col.names = TRUE
)


