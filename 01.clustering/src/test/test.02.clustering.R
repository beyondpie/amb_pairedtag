
library(tidyverse)
library(future)
library(Seurat)

# 10G
options(future.globals.maxSize= 10 * 1024^3)
old_path <- Sys.getenv("PATH")
Sys.setenv(PATH = paste(
  "/home/szu/miniforge3/envs/r/bin", old_path, sep = ":"))

# * configs
npcs <- 50
# based on PCA result in first round
npcs_knn <- 50
knn_k <- 50
knn_method <- "annoy"
reso <- 0.8

# * meta
projdir <- here::here()
workdir <- file.path(projdir, "01.clustering")
outdir <- file.path(workdir, "out")
figdir <- file.path(outdir, "figures")
dir.create(path = figdir)

# * main
s5 <- readRDS(file.path(outdir, "scRNAseq_k8_all.seurat.rds"))
# about 10-20 minutes。
s5 <- Seurat::NormalizeData(object = s5,
  normalization.method = "LogNormalize",
  scale.factor = 10000)
# about 10-20 minutes
# scale data is not sparse, but take lots of memory
s5 <- Seurat::ScaleData(object = s5,
  features = rownames(s5),
  do.scale = TRUE,
  do.center = TRUE)
# about one hour
# and we can save it indepedently
# take about 300G RAM
# No NA values in PCA
s5 <- Seurat::RunPCA(object = s5,
  npcs = 50,
  seed.use = 42,
  features = rownames(s5))
# this takes lots of time and space (91G)
# saveRDS(s5, file.path(outdir, "scRNAseq_k8_all.pca.seurat.rds"))

saveRDS(s5@reductions, file.path(outdir, "scRNAseq_k8_pca.rds"))

# * check PCA components
withr::with_pdf(new = file.path(figdir, "pca.std.scRNAseq_k8.all.pdf"),
  code = {
    plot(x = seq_len(npcs), y = s5@reductions$pca@stdev)
  })

# * run KNN
# about 30-60 minutes

# Error 1  <= use graph.name to fix this error
## rlang::last_trace()
## <error/rlang_error>
## Error in `[[<-`:
## ! `i` must be one of "nn", not "RNA_nn".
## ---
## Backtrace:
##     ▆
##  1. ├─Seurat::FindNeighbors(...)
##  2. └─Seurat:::FindNeighbors.Seurat(...)
##  3.   ├─methods (local) `[[<-`(`*tmp*`, graph.name[[ii]], value = `<named list>`)
##  4.   └─SeuratObject (local) `[[<-`(i = graph.name[[ii]])
## Run rlang::last_trace(drop = FALSE) to see 9 hidden frames.

# Error 2 <= use compute.SNN to skip this
# Cannot compute.SNN, will have std::bad_alloc error
pcaRes <- readRDS(file.path(outdir, "scRNAseq_k8_pca.rds"))
s5@reductions <- pcaRes
s5 <- Seurat::FindNeighbors(
  object = s5,
  k.param = knn_k,
  nn.method = knn_method,
  dims = 1:npcs_knn,
  reduction = "pca",
  compute.SNN = FALSE,
  graph.name = "nn"
)

saveRDS(s5@graphs,
  file = file.path(outdir,
    str_glue("scRNAseq_k8_pca.{knn_method}-pc{npcs_knn}-k{knn_k}.graphs.rds")))

# * run clustering
# It takes some time
# results have NA warnings on reticulate values load
# but all the cluster ids are not NA.
s5 <- Seurat::FindClusters(
  object = s5,
  resolution = reso,
  method = "igraph",
  algorithm = 4,
  n.start = 10,
  n.iter = 10,
  random.seed = 0,
  group.singletons = FALSE,
  verbose = TRUE,
  graph.name = "nn"
)

saveRDS(s5@meta.data,
  file.path(outdir, str_glue(
    "scRNAseq_k8_pca.{knn_method}-pc{npcs_knn}-k{knn_k}.r-{reso}.leiden.rds")
  ))




