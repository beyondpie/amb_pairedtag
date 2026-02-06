library(tidyverse)
library(tmpRpkg)

# * meta
projd <- here::here()
workd <- file.path(projd, "17.repressiveMarks")
outd <- file.path(workd, "out")
figd <- file.path(workd, "figure")
fromd <- file.path(outd, "allGenomicRange")

ptscMeta <- tmpRpkg::loadptscMeta() |>
  x => x[x$ATAC > 0, ] |>
  x => `rownames<-`(x, x$pairedTagName)
allenscs <- ptscMeta$AllenIdName

# * functions

loadGroupSimMat <- function(group = "class",
                            mod = "H3K27me3", embed = "Umap") {
  # load all the data
  wthgfnm <- file.path(fromd,
    str_glue("pt.{group}.{mod}.withingroup.dist.{embed}.csv"))
  wthgSim <- data.table::fread(
    file = wthgfnm, header = F, sep = ",", data.table = F) |>
    setNames(object = _, nm = c("subclass", "dist")) |>
    x => `rownames<-`(x, x$subclass)
  btgfnm <- file.path(fromd,
    str_glue("pt.{group}.{mod}.betweengroup.dist.{embed}.csv"))
  btgSim <- data.table::fread(
    file = btgfnm, header = F, sep = ",", data.table = F) |>
    setNames(object = _, nm = c("sc1", "sc2", "dist"))

  # prepare the mat
  groups <- unique(c(btgSim$sc1, btgSim$sc2))

  # sort groups
  ids <- strsplit(groups, split = " ") |>
    x => vapply(x, \(i) as.integer(i[1]), 0L)
  groups <- groups[order(ids)]

  m <- matrix(data = NA,
    nrow = length(groups), ncol = length(groups),
    dimnames = list(
      groups, groups
    )
  )
  # sign in
  for (i in seq_len(nrow(wthgSim))) {
    sc <- wthgSim$subclass[i]
    if (sc %in% rownames(m)) {
      m[sc, sc] <- wthgSim$dist[i]
    }
  }

  for (i in seq_len(nrow(btgSim))) {
    sc1 <- btgSim[i, "sc1"]
    sc2 <- btgSim[i, "sc2"]
    if ((sc1 %in% rownames(m)) && (sc2 %in% rownames(m))) {
      dist <- btgSim[i, "dist"]
      m[sc1, sc2] <- dist
      m[sc2, sc1] <- dist
    }
  }
  return(m)
}

getAvgDistWithinGroup <- function(simMat) {
  a <- vapply(seq_len(nrow(simMat)), \(i) {
    simMat[i, i]
  }, 0.0)
  mean(a)
}

plotSimMat <- function(K27me3Mat, K9me3Mat, lowVal, highVal, name) {
  ncolor <- 30
  col1 <- circlize::colorRamp2(seq(lowVal, highVal, length = ncolor),
    viridis::viridis(ncolor))
  legendLab <- c(
    round(lowVal, 1),
    round(highVal, 1)
  )
  ## col2 <- circlize::colorRamp2(seq(lowVal, highVal, length = ncolor),
  ##   viridis::magma(ncolor))
  col2 <- circlize::colorRamp2(seq(lowVal, highVal, length = ncolor),
    viridis::viridis(ncolor))


  hmK27me3 <- ComplexHeatmap::Heatmap(
    mat = K27me3Mat,
    col = col1,
    show_row_names = F,
    show_column_names = F,
    cluster_columns = F,
    cluster_rows = F,
    show_row_dend = F,
    show_column_dend = F,
    use_raster = T,
    row_names_gp = grid::gpar(fontsize = 6),
    column_names_gp = grid::gpar(fontsize = 6),
    heatmap_legend_param = list(
      title = str_glue("H3K27me3 {name}"),
      at = c(lowVal, highVal),
      labels = legendLab,
      direction = "horizontal"
    ),
    column_title = "H3K27me3"
  )
  hmK9me3 <- ComplexHeatmap::Heatmap(
    mat = K9me3Mat,
    col = col2,
    show_row_names = T,
    show_column_names = T,
    cluster_columns = F,
    cluster_rows = F,
    show_row_dend = F,
    show_column_dend = F,
    use_raster = T,
    row_names_gp = grid::gpar(fontsize = 6),
    column_names_gp = grid::gpar(fontsize = 6),
    heatmap_legend_param = list(
      title = str_glue("H3K9me3 {name}"),
      at = c(lowVal, highVal),
      labels = legendLab,
      direction = "horizontal"
    ),
    column_title = "H3K9me3"
  )
  hmK27me3 + hmK9me3
}
# * main

# class-level similarty
group <- "class"
embed <- "Umap"
K27me3Mat <- loadGroupSimMat(group, mod = "H3K27me3", embed)
K9me3Mat <- loadGroupSimMat(group, mod = "H3K9me3", embed)

avgDistWithinGroup <- min(
  getAvgDistWithinGroup(K27me3Mat),
  getAvgDistWithinGroup(K9me3Mat))

maxDist <- max(max(K27me3Mat), max(K9me3Mat))

(
  p <- plotSimMat(K27me3Mat, K9me3Mat,
    lowVal = avgDistWithinGroup,
    highVal = maxDist,
    name = str_glue("Avg Euclidean distance of {group} with {embed} embed.")
  )
)

# save data for plotting
saveRDS(
  object = list(K27me3 = K27me3Mat, K9me3Mat = K9me3Mat),
  file = file.path(outd, str_glue("pt.{group}.{embed}.dist.rds"))
)


group <- "class"
embed <- "SpectralEmbed"
K27me3Mat <- loadGroupSimMat(group, mod = "H3K27me3", embed)
K9me3Mat <- loadGroupSimMat(group, mod = "H3K9me3", embed)

avgDistWithinGroup <- min(
  getAvgDistWithinGroup(K27me3Mat),
  getAvgDistWithinGroup(K9me3Mat))

maxDist <- max(max(K27me3Mat), max(K9me3Mat))

(
  p <- plotSimMat(K27me3Mat, K9me3Mat,
    lowVal = avgDistWithinGroup,
    highVal = maxDist,
    name = str_glue("Avg Euclidean distance of {group} with {embed} embed.")
  )
)

# save data for plotting
saveRDS(
  object = list(K27me3 = K27me3Mat, K9me3Mat = K9me3Mat),
  file = file.path(outd, str_glue("pt.{group}.{embed}.dist.rds"))
)



# subclass-level similarty
group <- "sc"
embed <- "Umap"
K27me3Mat <- loadGroupSimMat(group, mod = "H3K27me3", embed)
K9me3Mat <- loadGroupSimMat(group, mod = "H3K9me3", embed)

avgDistWithinGroup <- min(
  getAvgDistWithinGroup(K27me3Mat),
  getAvgDistWithinGroup(K9me3Mat))

neuscs <- rownames(K27me3Mat) |>
  x => x[!grepl("NN", x)]
maxDist <- max(
  max(K27me3Mat[neuscs, neuscs]), max(K9me3Mat[neuscs, neuscs]))

(
  p <- plotSimMat(K27me3Mat, K9me3Mat,
    lowVal = avgDistWithinGroup,
    highVal = maxDist,
    name = str_glue("Avg Euclidean distance of {group} with {embed} embed.")
  )
)

# save data for plotting
saveRDS(
  object = list(K27me3 = K27me3Mat, K9me3Mat = K9me3Mat),
  file = file.path(outd, str_glue("pt.{group}.{embed}.dist.rds"))
  )


group <- "sc"
embed <- "SpectralEmbed"
K27me3Mat <- loadGroupSimMat(group, mod = "H3K27me3", embed)
K9me3Mat <- loadGroupSimMat(group, mod = "H3K9me3", embed)

avgDistWithinGroup <- min(
  getAvgDistWithinGroup(K27me3Mat),
  getAvgDistWithinGroup(K9me3Mat))

maxDist <- max(max(K27me3Mat), max(K9me3Mat))

(
  p <- plotSimMat(K27me3Mat, K9me3Mat,
    lowVal = avgDistWithinGroup,
    highVal = maxDist,
    name = str_glue("Avg Euclidean distance of {group} with {embed} embed.")
  )
)

# save data for plotting
saveRDS(
  object = list(K27me3 = K27me3Mat, K9me3Mat = K9me3Mat),
  file = file.path(outd, str_glue("pt.{group}.{embed}.dist.rds"))
)
