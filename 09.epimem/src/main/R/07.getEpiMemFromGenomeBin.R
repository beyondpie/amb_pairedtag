suppressPackageStartupMessages({
  library(tidyverse)
  library(furrr)
  library(reticulate)
})
Sys.setenv("_R_USE_PIPEBIND_" = TRUE)

# * meta
projd <- here::here()
binSize <- 5000
bwd <- file.path(projd, "data", "ptDNAbam", "bigwig")
outd <- file.path(projd, "09.epimem", "out", "genomeBin")
args <- commandArgs(trailingOnly = TRUE)
sc <- args[1]

# hs <- c("H3K4me1", "H3K27me3")
# suffix <- c("e100.bs100.sm300", "e100.bs100.sm300")

# 2024-11-18 add H3K27ac signals
# hs <- c("H3K27ac")
# suffix <- c("e100.bs100.sm300")

# 2024-11-22 add H3K9me3 signals
hs <- c("H3K9me3")
argsLen <- length(args)
if (argsLen < 2) {
  hs <- c("H3K9me3")
} else {
  hs <- str_split_1(args[2], ",")
}

suffix <- vapply(hs, \(h) {
  if (h == "H3K9me3") {
    "e100.bs100.sm1000"
  } else {
    "e100.bs100.sm300"
  }
}, "e100")

outsuffix <- paste(paste(hs, collapse = "_"), binSize, sep = ".")

# * set python env
use_condaenv(
  condaenv = "base",
  conda = "/home/szu/mambaforge/bin/conda")
pybw <- import("pyBigWig")

# * functions

#' @param chromSizes data.frame
#' Here is an example
#' chrom_sizes <- data.frame(
#'   chrom = c("chr1", "chr2", "chr3"),
#'   size = c(248956422, 242193529, 198295559)
#' )
#' @return List of data.frame
#' the data.frame has two columns "startFrom" and "endTo"
partitionGenome <- function(chromSizes, binSize = 5000) {
  lapply(seq_len(nrow(chromSizes)), \(i) {
    chr <- chromSizes$chrom[i]
    chrSize <- chromSizes$size[i]
    starts <- as.integer(seq(1, chrSize, by = binSize))
    ends <- as.integer(c(starts[-1] - 1, chrSize))
    data.frame(
      startFrom = starts,
      endTo = ends
    )
  }) |> setNames(object = _, nm = chromSizes$chrom)
}


getEpiSignal <- function(bwhs, chr, genomeBin) {
  message(str_glue("{Sys.time()} {chr} getEpiSignal start."))
  s <- vapply(bwhs, \(bwh) {
    vapply(seq_len(nrow(genomeBin)), \(i){
      bwh$stats(chr,
        as.integer(genomeBin$startFrom[i]),
        as.integer(genomeBin$endTo[i]),
        type = "mean", exact = TRUE
      )
    }, 0.0)
  }, rep(0.0, nrow(genomeBin))) |>
    as.data.frame() |>
    setNames(object = _, nm = names(bwhs))
  s$region <- paste(chr,
    with(genomeBin, paste(startFrom, endTo, sep = "-")),
    sep = ":"
  )
  message(str_glue("{Sys.time()} {chr} getEpiSignal finishes."))
  return(s[, c(c("region"), names(bwhs))])
}


#' @param bwhs List of bigwig handles
#' @return data.frame
#' columns are region, epi names in turn
#' region is chr1:1-20 style
runEpiMap <- function(sc, bwhs, genomeBinList) {
  r <- lapply(seq_along(genomeBinList), \(i) {
    chr <- names(genomeBinList)[i]
    getEpiSignal(bwhs, chr, genomeBinList[[i]])
  })
  return(r)
}

# * main
chromSizes <- data.table::fread(
  file = file.path(projd, "meta", "mm10.chrom.sizes.lite"),
  header = FALSE, sep = "\t", data.table = FALSE
) |> setNames(object = _, nm = c("chrom", "size"))
genomeBinList <- partitionGenome(chromSizes, binSize)

## for (sc in c(
##   "001_CLA_EPd_CTX_Car3_Glut",
##   "002_IT_EP_CLA_Glut",
##   "318_Astro_NT_NN", "326_OPC_NN", "327_Oligo_NN",
##   "334_Microglia_NN"
## )) {
message(Sys.time(), " Run runEpiMap for ", sc, ".")
f <- file.path(outd, str_glue("{sc}.{outsuffix}.csv"))
if (file.exists(f)) {
  message(f, " exist, and skip.")
  quit(save = "no", status = 0)
}

bwhs <- lapply(seq_along(hs), \(i) {
  f <- file.path(bwd, str_glue("{sc}.{hs[i]}.{suffix[i]}.bw"))
  if (!file.exists(f)) {
    message(f, " does not exist, and skip.")
    quite(save = "no", status = 0)
  }
  pybw$open(file.path(bwd, str_glue("{sc}.{hs[i]}.{suffix[i]}.bw")))
}) |> setNames(object = _, nm = hs)

r <- runEpiMap(sc, bwhs, genomeBinList)
write.table(paste(colnames(r[[1]]), collapse = ","),
  file = f,
  append = TRUE, quote = FALSE, row.names = FALSE,
  col.names = FALSE
)

for (df in r) {
  lines <- vapply(seq_len(nrow(df)), \(i) {
    paste(df[i, ], collapse = ",")
  }, "chr1:1-2,0.2,0.2")
  write.table(lines,
    file = f, append = TRUE, quote = FALSE,
    row.names = FALSE, col.names = FALSE
  )
}
message(Sys.time(), " Finish runEpiMap for ", sc, ". Good Luck!")
