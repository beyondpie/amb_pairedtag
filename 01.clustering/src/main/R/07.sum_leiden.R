library(tidyverse)
library(ggpubr)
Sys.setenv("_R_USE_PIPEBIND_" = TRUE)

# * snakemake configs
leidens <- snakemake@input[['leidens']]
silhts <- snakemake@input[['silhts']]
cellmetafnm <- snakemake@input[['cellmeta']]
nsample_umap <- as.integer(snakemake@params['nsample_umap'])

out_leiden <- snakemake@output[['leiden']]
out_silht <- snakemake@output[['silht']]
out_umap <- snakemake@output[['umap']]
out_cellmeta <- snakemake@output[["cellmeta"]]

# DEBUG
## input_dir <- file.path(
##   here::here(), "01.clustering", "out/test")
## prefix <- "RNA_k8_npc50_k50_L1"
## resos <- c(0.1, 0.2, 0.3)

## leidens <- file.path(input_dir,
##   paste(prefix, paste0("r", resos, ".leiden.csv"),
##     sep = "_"))
## silhts <- file.path(input_dir,
##   paste(prefix, paste0("r", resos, ".silht.csv"),
##   sep = "_"))
## cellmetafnm <- file.path(input_dir, "RNA_k8_npc50_k50_L1_0.csv")
## nsample_umap <- 50000

## out_dir <- input_dir
## out_leiden <- file.path(out_dir, "RNA_vg_npc50_k50_L1_leiden.csv")
## out_silht <- file.path(out_dir, "RNA_vg_npc50_k50_L1_silht.csv")
## out_umap <- file.path(out_dir, "RNA_vg_npc50_k50_L1_umap.pdf")
## out_cellmeta <- file.path(out_dir,
##   "RNA_vg_npc50_k50_L1_cellmeta.leiden.csv")

# * merge leiden and silht

res_leiden <- lapply(leidens, \(x) {
  r <- data.table::fread(file = x, sep = ",", header = TRUE)
  data.table::setkey(r, "barcode")
}) |> reduce(.f = merge, by = "barcode") |>
  as.data.frame()
rownames(res_leiden) <- res_leiden$barcode

data.table::fwrite(res_leiden,
  file = out_leiden, sep = ",", row.names = F, col.names = T)

res_silht <- lapply(silhts, \(x) {
  r <- data.table::fread(file = x, sep = ",", header = TRUE)
}) |> do.call(what = cbind, args = _) |>
  as.data.frame()

data.table::fwrite(res_silht,
  file = out_silht, sep = ",", row.names = F, col.names = T)

# get meta data for leiden and silht
getSeeds <- function(colnms) {
  colnms|>
    x => x[grep("seed", x)] |>
    map_chr(\(x) str_split_1(x, "_")[2]) |>
    gsub("seed", "", x = _) |>
    as.integer() |> unique()
}
getResoCols <- function(r, .colnms = colnames(res_silht)) {
  # BUG: sometimes (not reproducible)
  # r = 1 will be treated as 1 or 1L
  r <- format(as.numeric(r), nsmall = 1)
  ptn <- str_glue("r{format(r)}_")
  .colnms[grep(ptn, .colnms)]
}

# for each resolution, let's fix the seed
# which is the best silht
# in the meanwhile, let's find the best reso
# based on avg of silht across all the seeds
getAvgSils <- function(r) {
  s <- res_silht[, getResoCols(r, colnames(res_silht))]
  mean(rowMeans(s))
}
resos <- colnames(res_leiden) |>
  x => x[grep("seed", x)] |>
  map_chr(\(x) str_split_1(x, "_")[1]) |>
  gsub("r", "", x = _) |>
  as.numeric() |> unique()

avgSils <- map_vec(resos, getAvgSils)
names(avgSils) <- paste0("r", format(resos, nsmall = 1))

maxAvgSil <- max(avgSils)
rr <- resos[which.max(avgSils)]
seeds <- getSeeds(getResoCols(rr))
maxSil <- max(colMeans(res_silht[, getResoCols(rr)]))
ss <- seeds[which.max(
  colMeans(res_silht[, getResoCols(rr)]))]

message(str_glue("maxAvgSil: {round(maxAvgSil, 3)}"))
message(str_glue("r chosen: {format(rr, nsmall=1)}"))
message(str_glue("seed to use: {ss}"))
message(str_glue("maxSil: {round(maxSil, 3)}"))

# now for each reso, find the best seeds
ss_each <- map_int(resos, \(r){
  colnms <- getResoCols(r)
  seeds[which.max(colMeans(res_silht[ , colnms]))]
})
names(ss_each) <- format(resos, nsmall=1)

resoCols <- paste(paste0("r", format(resos, nsmall=1)),
  paste0("seed", ss_each), sep = "_")

# * load cell meta
cellmeta <- data.table::fread(
  file = cellmetafnm, header = T, data.table = F)
rownames(cellmeta) <- cellmeta$barcode
cellmeta <- cbind(cellmeta, res_leiden[rownames(cellmeta), resoCols])
data.table::fwrite(cellmeta, file = out_cellmeta,
  sep = ",", row.names = F, col.names = T)

if (nrow(cellmeta) > nsample_umap) {
  cellmeta <- cellmeta[
    sample(seq_len(nrow(cellmeta)), size = nsample_umap, replace = FALSE), ]
}

# * draw uma
mytheme <- theme(
  panel.border = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  plot.title = element_text(colour = "black", hjust = 0.5,
    size = 12),
  axis.text = element_blank(),
  axis.ticks = element_blank(),
  axis.title = element_blank(),
  axis.line = element_line(colour = "black"),
  legend.position = "none"
)

umap_plots <- lapply(resos, \(r) {
  rf <- format(r, nsmall=1)
  sil <- round(avgSils[paste0("r",rf)], 2)
  rcolnm <- paste(paste0("r", rf),
    paste0("seed", ss_each[rf]), sep = "_")
  d <- cellmeta[, c("barcode", "umap1", "umap2", rcolnm)]
  d[[rcolnm]] <- factor(d[[rcolnm]])
  d$cluster <- d[[rcolnm]]
  ncl <- length(unique(d$cluster))
 p <- ggplot(data = d, mapping = aes(x = umap1, y = umap2)) +
   geom_point(aes(colour = cluster),
     size = 0.7, shape = 19, alpha = 0.7) +
    mytheme
  if (r == rr) {
    p <- p + theme(
     plot.title = element_text(
        colour = "red", face = "bold",
        hjust = 0.5, size = 12)
    )
  }
  p <- p + labs(title = str_glue("{rcolnm} avgSil: {sil} #cl: {ncl}"))
  return(p)
})


umap_all <- ggpubr::ggarrange(
  plotlist = umap_plots,
  ncol = 5,
  nrow = 4
)

ggsave(filename = out_umap, plot = umap_all,
  width = 15, height = 12)
