library(tidyverse)
library(ggrastr)
library(tmpRpkg)
library(ComplexHeatmap)

# * meta
projd <- here::here()
workd <- file.path(projd, "03.integration")
figd <- file.path(workd, "figure")

# * functions
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
  axis.title = element_text(colour = "black", size = 12),
  axis.line = element_line(colour = "black"),
  legend.position = "right",
  legend.title = element_text(colour = "black", size = 13),
  legend.text = element_text(colour = "black", size = 12)
)

to3.matrix <- function(mat,
                       names = c("row", "column", "score"),
                       int2str = TRUE,
                       factor2str = TRUE) {
  r <- reshape2::melt(mat)
  colnames(r) <- names
  if (factor2str) {
    if (is.factor(r[, 1])) {
      r[, 1] <- as.character(r[, 1])
    }
    if (is.factor(r[, 2])) {
      r[, 2] <- as.character(r[, 2])
    }
  }
  if (int2str) {
    if (is.numeric(r[, 1])) {
      r[, 1] <- as.character(r[, 1])
    }
    if (is.numeric(r[, 2])) {
      r[, 2] <- as.character(r[, 2])
    }
  }
  return(r)
}

readcnss <- function(cnssd) {
  x <- data.table::fread(
    file = file.path(cnssd, "mat.csv"),
    header = F, sep = ",", data.table = F
  )
  rnms <- data.table::fread(
    file = file.path(cnssd, "rownames.txt"),
    header = F, data.table = F
  )[, 1]
  cnms <- data.table::fread(
    file = file.path(cnssd, "colnames.txt"),
    header = F, data.table = F, sep = "@"
  )[, 1]
  x |>
    as.matrix() |>
    x => `rownames<-`(x, rnms) |>
    x => `colnames<-`(x, cnms)
}

plot_cnss_tile <- function(m) {
  lowSimScore <- quantile(m$score, 0.05)
  highSimScore <- quantile(m$score, 0.95)

  ggplot(data = m, aes(x = row, y = column)) +
    geom_tile_rast(aes(fill = score), color = "white") +
    mytheme +
    theme(
      panel.background = element_rect(fill = "white"),
      axis.text.x = element_blank(),
      axis.text.y = element_blank()
    ) +
    scale_color_gradient(
      low = "white", high = "red",
      limits = c(lowSimScore, highSimScore), na.value = "white"
    )
}



# * main
for (clevel in c("L3", "L4")) {
  for (group in c("nn", "neu", "all")) {
    cnssd <- file.path(workd, "out",
      str_glue("consensusMat_{group}_{clevel}"))
    cnss <- readcnss(cnssd)
    ## cnss3c <- to3.matrix(
    ##   mat = cnss, names = c("row", "column", "score"),
    ##   int2str = T, factor2str = T
    ## )

    ## cnss3c$row <- factor(cnss3c$row, levels = rownames(cnss))
    ## cnss3c$column <- factor(cnss3c$column, levels = colnames(cnss))

    ## p <- plot_cnss_tile(cnss3c) +
    ##   xlab(str_glue("Paired-Tag RNA {group} {clevel}")) +
    ##   ylab(str_glue("Allen 10Xv3 RNA subclass(neu) + supertype(nn)"))

    ## ggsave(filename = file.path(
    ##   figd,
    ##   str_glue("consensusMat_{group}_{clevel}.pdf")
    ## ), plot = p, width = 30, height = 10)

    probColfn <- circlize::colorRamp2(c(0.1, 0.7), c("white", "red"))
    ph <- ComplexHeatmap::Heatmap(
      matrix = t(cnss),
      name = "ConsensusMat",
      col = probColfn,
      show_row_names = F,
      show_column_names = F,
      cluster_rows = F,
      cluster_columns = F,
      row_names_gp = gpar(fontsize = 15),
      column_names_gp = gpar(fontsize = 15),
      use_raster = T
    )
    withr::with_pdf(
      new = file.path(figd, str_glue("consensusMat_{group}_{clevel}.chp.pdf")),
      code = {
        print(ph)
      },
      width = 8, height = 4
    )
  }
}

# get the cell numbers during integration
library(Seurat)
allend <- file.path(projd, "data", "allen_seurat")
ptd <- file.path(projd, "data", "pairedtag_seurat")
## global integration (both neu and nn)

### nn
# 35098 cells
allen_nn <- readRDS(file.path(allend,
  "allen.10xv3.pt.regions.nn.imn.k8.cl.ds1000.rds"))
# 40469 cells
pt_nn <- readRDS(file.path(ptd, "ptRNA.nn.k8.L5ds30.rds"))

### neu
# 139912 cells
allen_neu <- readRDS(file.path(allend,
  "allen.10xv3.pt.regions.neu.imn.k8.cl.ds100.rds"))
# 145254 cells
pt_neu <- readRDS(file.path(ptd, "ptRNA.neu.k8.L5ds50.rds"))

## by region (neu)
regions <- c(
 "AMY", "CPU", "HYP", "HIP", "ERC", "NAC", "VTA", "PFC" 
)
## AMY has 76255 allen cells.
## CPU has 15995 allen cells.
## HYP has 115263 allen cells.
## HIP has 36498 allen cells.
## ERC has 47294 allen cells.
## NAC has 15470 allen cells.
## VTA has 5238 allen cells.
## PFC has 91285 allen cells.
for (r in regions) {
  allen_r_neu <- readRDS(file.path(allend,
    str_glue("allen.10xv3.{r}.neu.ds.seu.rds")))
  message(r," has ", dim(allen_r_neu)[2], " allen cells.")
}

## AMY has 75503 pt cells.
## CPU has 17206 pt cells.
## HYP has 109394 pt cells.
## HIP has 36405 pt cells.
## ERC has 48582 pt cells.
## NAC has 16780 pt cells.
## VTA has 7291 pt cells.
## PFC has 90762 pt cells.
for (r in regions) {
  pt_r_neu <- readRDS(file.path(ptd, "neu_seu_region",
    str_glue("pt.neu.{r}.ds.seu.rds")))
  message(r," has ", dim(pt_r_neu)[2], " pt cells.")
}
