suppressPackageStartupMessages({
  library(tidyverse)
  library(GenomicRanges)
  library(reticulate)
  library(tmpRpkg)
  library(ComplexHeatmap)
  library(future)
  library(furrr)
})

# * meta
plan(strategy = "multicore", workers = 5)

use_condaenv(
  condaenv = "base",
  conda = "/home/szu/mambaforge/bin/conda"
)
pybw <- import("pyBigWig")
projd <- "/projects/ps-renlab2/szu/projects/amb_pairedtag"
# h <- "H3K27ac"
# supenhd <- file.path(projd, "10.superEnhancer", "out", "ROSE")

h <- "H3K27me3"
supenhd <- file.path(projd, "10.superEnhancer", "out", "ROSE_H3K27me3")

bwd <- file.path(projd, "data", "ptDNAbam", "bigwig")
scs <- list.files(supenhd,
  full.names = FALSE, no.. = TRUE, all.files = FALSE
)
outd <- file.path(projd, "10.superEnhancer", "out", "heatmap")

# * functions
loadSuperEnhancer <- function(sc) {
  r <- data.table::fread(
    file = file.path(
      supenhd, sc,
      str_glue("{sc}_AllStitched.table.txt")
    ),
    head = TRUE, sep = "\t", skip = "#", data.table = FALSE
  )
  return(r)
}

loadH3K27acBigWig <- function(sc, h = "H3K27ac") {
  r <- pybw$open(
    file.path(bwd, str_glue("{sc}.{h}.e100.bs100.sm300.bw"))
  )
  return(r)
}

# use it for eRNA
isWithinGene <- function(queryGR, up2TSS = 2000) {}

mapBigWigSignalOnCRE <- function(bwhandle, chrom, startFrom, endTo) {
  r <- bwhandle$stats(chrom, startFrom, endTo,
    type = "mean", exact = TRUE
  )
  return(r)
}

mapBigWigSignalOnEnhancers <- function(bwhandle, enhs) {
  vapply(seq_len(nrow(enhs)), \(i) {
    mapBigWigSignalOnCRE(
      bwhandle, enhs[i, "CHROM"], enhs[i, "START"], enhs[i, "STOP"]
    )
  }, FUN.VALUE = 0.0)
}

# * main
# ** load the super enhancers and bigwigs
supenhl <- lapply(scs, loadSuperEnhancer) |>
  setNames(object = _, nm = scs)

## supenhGRl <- lapply(scs, \(sc) {
##   d <- supenhl[[sc]]
##   loadGRfromBed(beds = d[, 2:ncol(d)], header = TRUE)
## }) |> setNames(object = _, nm = scs)

bwl <- lapply(scs, loadH3K27acBigWig, h = h) |>
  setNames(object = _, nm = scs)

# * map H3K27ac bigwig signals on these enhancers
signalH3K27acl <- lapply(scs, \(sc) {
  mapBigWigSignalOnEnhancers(bwl[[sc]], supenhl[[sc]])
}) |> setNames(object = _, nm = scs)

saveRDS(
  signalH3K27acl,
  file.path(outd, str_glue("signal{h}.scEnhancer.scList.rds"))
)

flatEnhs <- do.call(rbind, supenhl)
flatEnhs$sc <- do.call("c", lapply(scs, \(sc) {
  rep(sc, nrow(supenhl[[sc]]))
}))
flatSuperEnhs <- flatEnhs[flatEnhs$isSuper > 0, ]
rownames(flatSuperEnhs) <- with(
  flatSuperEnhs,
  paste(sc,
    paste(CHROM, paste(START, STOP, sep = "-"), sep = ":"),
    sep = "@"
  )
)

signalH3K27ac_allSupEnh <- furrr::future_map(scs, \(sc) {
  r <- mapBigWigSignalOnEnhancers(bwl[[sc]], flatSuperEnhs)
  message("Done: get signals from the bigwig of ", sc)
  return(r)
}, .options = furrr_options(seed = TRUE)) |>
  setNames(object = _, nm = scs)

saveRDS(
  signalH3K27ac_allSupEnh,
  file.path(outd, str_glue("signal{h}.allsuperEnhancer.scList.rds"))
)
