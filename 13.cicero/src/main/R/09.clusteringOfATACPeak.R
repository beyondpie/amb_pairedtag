suppressMessages({
  suppressWarnings({
    library(tidyverse)
    library(uwot)
    library(RcppHNSW)
    library(rnndescent)
    library(RcppAnnoy)
    library(igraph)
    library(tmpRpkg)
  })
})

Sys.setenv("_R_USE_PIPEBIND_" = TRUE)

# * meta
projd <- here::here()
workd <- file.path(projd, "13.cicero")

args <- commandArgs(trailingOnly = T)
p <- as.integer(args[1])
k <- as.integer(args[2])
reso <- as.numeric(args[3])
addweight <- args[4]
weight <- if (addweight == "weight") {
  NA
} else {
  NULL
}
sc <- args[5]
embedf <- args[6]
suffix <- str_glue("p{p}-k{k}-r{reso}-w{addweight}")

# * functions
loadSpectralEmbed <- function(fnm) {
  r <- data.table::fread(fnm, sep = ",", header = F, data.table = F) |>
    x => `rownames<-`(x, x[, 1])
  d <- r[, 3:ncol(r)] |> as.matrix()
  list(raw = r, mat = d)
}

# * main
outd <- file.path(workd, "out", "spectral", sc)
# embedfnm <- file.path(outd, "spectralEmbed.k20.csv")
embedObj <- loadSpectralEmbed(fnm = embedf)

# umap
umap <- uwot::umap(
  X = embedObj$mat[, seq_len(p)],
  n_neighbors = k, min_dist = 0.01, nn_method = "hnsw"
)
colnames(umap) <- c("UMAP1", "UMAP2")

# leiden algorithm
knn <- AnnoyNN(embedObj$mat[, seq_len(p)], k = k)

edge_list <- lapply(seq_len(nrow(umap)), \(i) {
  data.frame(
    from = i,
    to = knn$knn[i, ],
    weight = 1 / (1 + knn$dist[i, ])
  )
}) |> do.call(rbind, args = _)

g <- igraph::graph_from_data_frame(edge_list, directed = F)
leiden_clusters <- igraph::cluster_leiden(
  graph = g, objective_function = "modularity",
  resolution = reso, weights = weight, n_iterations = 20
)

cluster_assign <- leiden_clusters$membership

# show UMAP and cluster_assign
r <- data.frame(
  UMAP1 = umap[, 1],
  UMAP2 = umap[, 2],
  cluster = factor(cluster_assign)
)

p <- ggplot(data = r, mapping = aes(x = UMAP1, y = UMAP2, color = cluster)) +
  geom_point(size = 0.1) +
  ggtitle(str_glue("UMAP of ATAC peaks of {sc}")) +
  theme(
    axis.text.x = element_text(size = 14, colour = "black"),
    axis.text.y = element_text(size = 14, colour = "black"),
    plot.title = element_text(size = 15, colour = "black", hjust = 0.5)
  ) +
  scale_color_manual(values = SnapATACPalette) + 
  guides(colour = guide_legend(override.aes = list(size = 5))) +
  theme_minimal(
    base_size = 15
  )

withr::with_pdf(
  new = file.path(outd, str_glue("{sc}.ATAC.peak.UMAP.{suffix}.pdf")),
  code = {
    print(p)
  },
  width = 7, height = 6
)

# output the UMAP and leiden results
rr <- data.frame(
  CRE = embedObj$raw[, 1],
  state = embedObj$raw[, 2],
  UMAP1 = r$UMAP1,
  UMAP2 = r$UMAP2,
  cluster = cluster_assign
)

data.table::fwrite(
  x = rr,
  file = file.path(outd, str_glue("{sc}.ATAC.peak.UMAP.Leiden.{suffix}.csv")),
  sep = ",", col.names = T, row.names = F
)
