library(tidyverse)
library(rlang)
library(ggpubr)

projroot <- here::here()
# * configs
clrec <- data.table::fread(
  file = file.path(projroot, "meta", "pt.scRNAseq.sa2clustering.record.csv"),
  sep = ",", header = T, data.table = F
)
rownames(clrec) <- with(clrec, paste(
  clustering_level, paste0("c", cluster),
  sep = ":"
))

L2resdir <- file.path(
  projroot, "01.clustering",
  "out/afqc/L2"
)
resprefix <- "RNA_vg3000_npc30_k50_L1_"

L1cls <- paste0("c", seq_len(13) - 1)
# load all the leiden results
L1leidens <- lapply(L1cls, \(x) {
  data.table::fread(
    file = file.path(L2resdir, paste0(resprefix, x, ".leiden.csv")),
    sep = ",",
    header = TRUE
  )
}) |> do.call(what = rbind, args = _)
rownames(L1leidens) <- L1leidens$barcode

# load all the leiden results
L1silhts <- lapply(L1cls, \(x) {
  data.table::fread(
    file = file.path(L2resdir, paste0(resprefix, x, ".silht.csv")),
    sep = ",",
    header = TRUE
  )
})
names(L1silhts) <- L1cls

# leiden different reso with best seed
L1cellmetas <- lapply(L1cls, \(x) {
  data.table::fread(
    file = file.path(L2resdir, paste0(resprefix, x, ".cellmeta.leiden.csv")),
    sep = ",", header = TRUE
  )
})
names(L1cellmetas) <- L1cls

# * double check clustering meta on L2
## check ncl match
all(vapply(L1cls, \(x) {
  rec <- clrec[str_glue("L1:{x}"), , drop = TRUE]
  r <- TRUE
  if (rec$use_auto) {
    rrec <- rec$resolution
    srec <- rec$seed
    nrec <- rec$ncl
    t <- str_glue("r{rrec}_seed{srec}")
    # use .. since it's data.table
    ncl <- length(table(L1cellmetas[[x]][, ..t]))
    r <- ifelse(ncl == nrec, TRUE, FALSE)
  }
  r
}, TRUE))
# check the best silht
seeds <- seq_len(10) - 1
rs <- seq(0.1, 2, 0.1)
getAvgSil <- function(r, m) {
  r <- as.numeric(r)
  cols <- paste0(
    stringr::str_glue("r{format(r, nsmall = 1)}_"), "seed", seeds
  )
  mean(colMeans(m[, cols]))
}

all(
  a <- vapply(L1cls, \(x) {
    rec <- clrec[str_glue("L1:{x}"), , drop = TRUE]
    r <- TRUE
    if (rec$use_auto) {
      rrec <- rec$resolution
      srec <- rec$seed
      nrec <- rec$ncl
      t <- str_glue("r{rrec}_seed{srec}")
      m <- L1silhts[[x]] |> as.data.frame()
      max_r <- rs[which.max(map_vec(rs, getAvgSil, m = m))]
      is_max_r <- ifelse((max_r - rrec) < 0.001, TRUE, FALSE)
      cols <- paste0(str_glue("r{format(rrec, nsmall = 1)}_"), "seed", seeds)
      mm <- m[, cols]
      is_max_seed <- ifelse(seeds[which.max(colMeans(mm))] == srec,
        TRUE, FALSE
      )
      r <- ifelse(is_max_r & is_max_seed, TRUE, FALSE)
    }
    r
  }, TRUE)
)

# * organize cell meta with L2 clustering and umap
allL2UMAP <- lapply(L1cellmetas, \(x) {
  as.data.frame(x) |>
    x => x[, c("barcode", "umap1", "umap2")] |>
    rename(.data = _, L2_UMAP1 = umap1, L2_UMAP2 = umap2)
}) |>
  do.call(what = rbind, args = _)
rownames(allL2UMAP) <- allL2UMAP$barcode

allL2 <- lapply(names(L1cellmetas), \(nm) {
  x <- as.data.frame(L1cellmetas[[nm]])
  rec <- clrec[str_glue("L1:{nm}"), , drop = TRUE]
  r <- rec$resolution
  s <- rec$seed
  t <- str_glue("r{r}_seed{s}")
  data.frame(
    barcode = x$barcode,
    L2 = x[[t]]
  )
}) |> do.call(what = rbind, args = _)
rownames(allL2) <- allL2$barcode

curCellMeta <- data.table::fread(
  file = file.path(projroot, "meta", "pt.barcode.meta.withL1.csv"),
  sep = ",",
  header = TRUE, data.table = FALSE
)

rownames(curCellMeta) <- curCellMeta$barcode
curCellMeta <- rename(curCellMeta, L1_UMAP1 = umap1, L1_UMAP2 = umap2)
curCellMeta$L2_UMAP1 <- allL2UMAP[curCellMeta$barcode, "L2_UMAP1"]
curCellMeta$L2_UMAP2 <- allL2UMAP[curCellMeta$barcode, "L2_UMAP2"]
curCellMeta$L2 <- allL2[curCellMeta$barcode, "L2"]

all(vapply(L1cls, \(x) {
  rec <- clrec[str_glue("L1:{x}"), , drop = TRUE]
  ifelse(rec$ncl == length(
    table(with(curCellMeta, L2[L1 == as.integer(gsub("c", "", x))]))
  ),
  TRUE, FALSE
  )
}, TRUE))

data.table::fwrite(
  curCellMeta,
  file = file.path(projroot, "meta", "pt.barcode.meta.withL2.csv"),
  sep = ",",
  row.names = FALSE, col.names = TRUE
)

# * draw different factors on the L2 UMAPs.
outfig_dir <- file.path(
  projroot, "01.clustering", "out",
  "afqc", "L2", "figures"
)
mytheme <- theme(
  panel.border = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  plot.title = element_text(
    colour = "black", hjust = 0.5,
    size = 15
  ),
  axis.text = element_blank(),
  axis.ticks = element_blank(),
  axis.title = element_blank(),
  axis.line = element_line(colour = "black"),
  legend.position = "right",
  legend.title = element_text(colour = "black", size = 10),
  legend.text = element_text(colour = "black", size = 10),
)

cellMeta <- data.table::fread(
  file = file.path(projroot, "meta", "pt.barcode.meta.withL2.csv"),
  sep = ",", header = T, data.table = F
)
nds <- 10000

all_p_list <- map(seq_len(13) - 1, \(L1cl) {
  umapdf <- cellMeta[
    cellMeta$L1 == L1cl,
    c(
      "L2_UMAP1", "L2_UMAP2", "brainregion", "modality", "exp",
      "sex", "rep", "L2"
    )
  ]
  if (nrow(umapdf) > nds) {
    umapdf <- umapdf[sample(seq_len(nrow(umapdf)), nds), ]
  }
  umapdf$L2 <- factor(umapdf$L2)
  p_list <- map(c("L2", "modality", "brainregion", "exp", "sex", "rep"), \(colnm) {
    ggplot(data = umapdf, mapping = aes(x = L2_UMAP1, y = L2_UMAP2)) +
      geom_point(aes(color = .data[[colnm]]), size = 0.8, shape = 19, alpha = 0.5) +
      mytheme +
      guides(color = guide_legend(override.aes = list(size = 4)))
  })
  nr <- 3
  p <- ggarrange(plotlist = p_list, ncol = length(p_list) / nr, nrow = nr)
  ggsave(
    filename = file.path(outfig_dir, str_glue("L2.c{L1cl}.pdf")),
    plot = p,
    width = 12,
    height = 15
  )
})
# works
# ggarrange(plotlist = all_p_list, ncol = 1, nrow = 13)

# * prepare meta file for further L3 clustering
cellMeta <- data.table::fread(
  file = file.path(projroot, "meta", "pt.barcode.meta.withL2.csv"),
  sep = ",", header = T, data.table = F
)
cellMeta$L1_2 <- with(cellMeta, L1 * 100 + L2)
data.table::fwrite(
  cellMeta,
  file = file.path(projroot, "meta", "pt.barcode.meta.withL2.csv"),
  sep = ",",
  row.names = FALSE, col.names = TRUE
)
