suppressMessages({
  suppressWarnings({
    library(cicero)
    library(tidyverse)
    library(reticulate)
    library(tmpRpkg)
    library(future.apply)
    library(Matrix)
  })
})

# this script should accept arguments outside.
# we can run it among NN and Neu seperately
# or we can run it for all the cells
# or we can run it for different classes.

# * meta
setCondaEnv(envnm = "sa2stable")
future::plan(strategy = multisession, workers = 20)
options(future.globals.maxSize = 40 * 1024^3)
ad <- reticulate::import(module = "anndata", convert = FALSE)
projd <- here::here()
outd <- file.path(projd, "13.cicero", "out")
set.seed(2024)
mouse_mm10_genome <- read.table(
  file = file.path(projd, "meta", "mm10.chrom.sizes.lite"),
  header = F, sep = "\t"
)
ncore4LSI_UMAP <- 10

# * setup program input and output
annfnm <- file.path(projd, "12.DE", "out", "ann.H3K27ac.pmat.h5ad")
cdsbfUMAPfnm <- file.path(outd, "cds.all.H3K27ac.bfUMAP.rds")
cdsUMAPfnm <- file.path(outd, "cds.all.H3K27ac.UMAP.rds")
cicerocdsfnm <- file.path(outd, "cicero.cds.all.H3K27ac.rds")
distparamsfnm <- file.path(outd, "distanceParameters.all.H3K27ac.rds")
ciceromodelfnm <- file.path(outd, "cicero.model.all.H3K27ac.rds")
connsfnm <- file.path(outd, "conns.all.H3K27ac.csv")

# * load whole matrix
ann <- ad$read_h5ad(annfnm, backed = "r")

# * to cicero cds
mat <- getFullAnnX(ann)
cds <- getBinarizeMonocle3Obj(mat, isRowBarcode = T)

# tmp to save the data
saveRDS(object = cds, file = cdsbfUMAPfnm)

# * perform UMAP
# TODO: perform low-dimensional reduction in SnapATAC2
# 3 - 4 hours needed for this step
cds <- monocle3::estimate_size_factors(cds)
cds <- monocle3::preprocess_cds(cds, method = "LSI", num_dim = 50)
cds <- monocle3::reduce_dimension(cds,
  max_components = 2,
  reduction_method = "UMAP",
  preprocess_method = "LSI", cores = ncore4LSI_UMAP
)
saveRDS(object = cds, file = cdsUMAPfnm)
message(Sys.time())

# * get cicero cds
## > + Overlap QC metrics:
## Cells per bin: 50
## Maximum shared cells bin-bin: 44
## Mean shared cells bin-bin: 0.00516665008017892
## Median shared cells bin-bin: 0
cds_cicero <- cicero::make_cicero_cds(
  cds = cds,
  reduced_coordinates = SingleCellExperiment::reducedDims(cds)$UMAP,
  k = 50, return_agg_info = F
)
saveRDS(object = cds, file = cicerocdsfnm)
message(Sys.time())

# * run cicero with explict functions in steps.
# will use future for parallization
distance_parameters <- p_estimate_distance_parameter(
  cds,
  window = 5e+05, maxit = 100,
  s = 0.75,
  sample_num = 100,
  distance_constraint = 250000,
  genomic_coords = mouse_mm10_genome,
)
# this takes about 5 hours
message(Sys.time())
# distance_parameter
# [1] 2.273677e-06
distance_parameter <- mean(distance_parameters)
message("Mean of distance parameters: ", distance_parameter)
saveRDS(distance_parameters, distparamsfnm)

# will use future for parallization
# takes about 5 hours under 10 cores
message(Sys.time())
cicero_cds <- readRDS(cicerocdsfnm)
dps <- readRDS(distparamsfnm)
distance_parameter <- mean(dps)

cicero_model <- p_generate_cicero_models(
  cicero_cds,
  distance_parameter,
  s = 0.75,
  window = 5e+05,
  genomic_coords = mouse_mm10_genome
)
message(Sys.time())
saveRDS(cicero_model, ciceromodelfnm)

conns <- assemble_connections(
  cicero_model_list = cicero_model, silent = F
)

conns$Peak2 <- levels(conns$Peak2)[conns$Peak2]
conns <- conns[!is.na(conns$coaccess), ]

data.table::fwrite(
  x = conns,
  file = connsfnm,
  sep = ",",
  col.names = T, row.names = F
)
message("Done. Good luck!")
