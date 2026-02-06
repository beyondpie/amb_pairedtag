suppressPackageStartupMessages({
  library(tidyverse)
  library(reticulate)
  library(GenomicRanges)
})

# * set python env
use_condaenv(condaenv = "base",
  conda = "/home/szu/mambaforge/bin/conda")
pybw <- import("pyBigWig")

# * meta
projd <- here::here()
bwd <- file.path(projd, "data", "ptDNAbam", "bigwig")
peakd <- file.path(projd, "data", "pariedtag_peak", "subclass_peak")
epimemoutd <- file.path(projd, "09.epimem", "out")
repressMark <- "H3K27me3"
activeMark <- "H3K4me1"
sc <- "327_Oligo_NN"
# sc <- "326_OPC_NN"

# * functions
loadMappedSignals <- function(sc, from_m, to_m) {
  fnm <- file.path(epimemoutd,
    str_glue("{from_m}_on_{to_m}"),
    str_glue("{sc}_${from_m}-on-${to_m}_e2500.tsv"))
  r <- data.table::fread(
    file = fnm, header = TRUE, sep = "\t")
  return(r)
}

# * main
# given a subclass, load peak-related signals
# H3K27me3 and h3K4me1
r2a <- loadMappedSignals(
  sc, repressMark, active_Mark)
r2r <- loadMappedSignals(sc, repressMark, repressMark)
a2r <- loadMappedSignals(sc, activeMark, repressMark)
a2a <- loadMappedSignals(sc, activeMark, activeMark)

# plot the signals to have a look firstly.

# mark the regions that have both the signals


# scatter plot and show the regions we selected
