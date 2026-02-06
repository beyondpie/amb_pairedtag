library(tidyverse)
library(reticulate)
library(hdf5r)
library(devtools)
# load my tmp package under project: tmpRpkg
load_all()

# * meta
projd <- tmpRpkg::projd
rscd <- file.path(projd, "data", "snATAC")
cpmOutd <- hiscntdr

atacPeakSrt <- loadATACPeakSrt()
peakInPtSrt <- loadATACPeakInPtSrt()
allenClMeta <- loadAllenClMeta()
sp2sc <- load_sp2sc()
sc2scil <- load_sc2scil()

# * main
# -- ATAC-seq data
cpm_pbysc_atac <- load_snATACPMpbyc()

# -- Histone modifications
hm <- "H3K27ac"
out_cpmhm <- file.path(cpmOutd,
  str_glue("cpm.scbyp.{hm}.h5"))
cnt <- readRDS(
  file.path(histcntdr, str_glue("featureCount.{hm}.mat.rds"))
) |>
  mapfc2fcsc(fc = _, sp2scMap = sp2sc) |>
  x => `rownames<-`(x, atacPeakSrt)
sc.final <- intersect(
  colnames(cpm_pbysc_atac), colnames(cnt)) |>
   x => sortAllenLabelById(x)
rcnt <- cnt[peakInPtSrt, sc.final]
cpm_pbyc <- getCPMOfHistones(rcnt)
saveCPM2h5(
  mat = t(cpm_pbyc),
  peaks = peakInPtSrt,
  clusters = sc.final,
  outf = out_cpmhm)
write.table(peakInPtSrt,
  file = file.path(cpmOutd, str_glue("peak.sc.{hm}.txt")),
  quote = FALSE, row.names = FALSE, col.names = FALSE
)
write.table(sc.final,
  file = file.path(cpmOutd, str_glue("subclass.{hm}.txt")),
  quote = FALSE, row.names = FALSE, col.names = FALSE
)

# test h5ad file.
## f <- file.path(projd, "data/ptHistoneCounts",
##   "ATACPeak_rds", "cpm.scbyp.H3K27ac.h5")
## fh5 <- H5File$new(f, mode="r")
## mat <- fh5[["X/mat"]][, ]
## scs <- fh5[["X/rownames"]][]
## peaks <- fh5[["X/colnames"]][]
