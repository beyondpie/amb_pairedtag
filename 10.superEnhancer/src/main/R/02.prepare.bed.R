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
projd <- here::here()
# h <- "H3K27me3"
h <- "H3K27ac"
supenhd <- file.path(projd, "10.superEnhancer", "out", str_glue("ROSE_{h}"))
outd <- file.path(
  projd, "10.superEnhancer", "out",
  str_glue("ROSE_{h}_bed")
)

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

# * main
scs <- list.files(supenhd,
  full.names = FALSE, no.. = TRUE, all.files = FALSE
)

supenhl <- lapply(scs, loadSuperEnhancer) |>
  setNames(object = _, nm = scs)

n_all <- vapply(supenhl, nrow, 9)
n_sup <- vapply(supenhl, \(i) sum(i$isSuper > 0), 9)

flatEnhs <- do.call(rbind, supenhl)
flatEnhs$sc <- do.call("c", lapply(scs, \(sc) {
  rep(sc, nrow(supenhl[[sc]]))
}))
flatSuperEnhs <- flatEnhs[flatEnhs$isSuper > 0, ]
data.table::fwrite(flatSuperEnhs,
  file = file.path(outd, str_glue("all.superEnhancers.{h}.tsv")),
  sep = "\t", col.names = T, row.names = F
)

# all peaks used for background
# for H3K27me3
## allPeaks <- data.table::fread(
##   file = file.path(projd, "data", "pairedtag_peak", "merge_peak",
##     str_glue("{h}.merged.all.blv2.me.peak")),
##   sep = "\t", head = TRUE, data.table = FALSE
## )
# for H3K27ac
allPeaks <- data.table::fread(
  file = file.path(
    projd, "data", "pairedtag_peak", "merge_peak",
    str_glue("{h}.merged.all.blv2.me.bestspm.peak")
  ),
  sep = "\t", head = TRUE, data.table = FALSE
)
data.table::fwrite(
  x = with(allPeaks, data.frame(
    chrom = chrom,
    startFrom = startFrom,
    endTo = endTo
  )),
  file = file.path(outd, str_glue("all.peaks.{h}.bed")),
  sep = "\t", col.names = F, row.names = F
)


# all super enhancers
data.table::fwrite(
  x = with(
    flatSuperEnhs,
    data.frame(
      chrom = CHROM,
      startFrom = START,
      endTo = STOP
    )
  ),
  file = file.path(outd, str_glue("all.superEnhancers.{h}.bed")),
  sep = "\t", col.names = F, row.names = F
)

# save super enhancers per subclass
for (sc in scs) {
  r <- supenhl[[sc]]
  r <- r[r$isSuper > 0, ]
  data.table::fwrite(
    x = with(r, data.frame(
      chrom = CHROM,
      startFrom = START,
      endTo = STOP
    )),
    file = file.path(outd, str_glue("{sc}.superEhancers.{h}.bed")),
    sep = "\t", col.names = F, row.names = F
  )
}

# load peaks called per subclass
# h <- "H3K27ac"
h <- "H3K27me3"
peakd <- file.path(projd, "data", "pairedtag_peak", "subclass_peak")
suffix <- if (h == "H3K27ac") {
  "BestSPM"
} else {
  "bedtoolmerge"
}
for (sc in scs) {
  raw <- data.table::fread(
    file = file.path(peakd, str_glue("{sc}-{h}.{suffix}.peak")),
    header = T, sep = "\t", data.table = FALSE
  )
  data.table::fwrite(
    x = with(raw, data.frame(
      chrom = chrom,
      startFrom = startFrom,
      endTo = endTo
    )),
    file = file.path(
      projd, "data", "pairedtag_peak", "subclass_peak_bed",
      str_glue("{sc}-{h}.{suffix}.bed")
    ),
    col.names = F, row.names = F, sep = "\t"
  )
}

# * prepare enhancers (H3K27ac peaks overlapped with super enhancers)
# for motif analysis

# ** load super enhancers
supEnhBedir <- file.path(
  projd, "10.superEnhancer",
  "out", "ROSE_H3K27ac_bed"
)
scs <- list.files(supEnhBedir) |>
  strsplit(x = _, split = "\\.", perl = T) |>
  lapply(X = _, FUN = \(i) {
    i[1]
  }) |>
  unlist() |>
  x => x[grepl("_", x)]

id_scs <- strsplit(scs, split = "_") |>
  vapply(X = _, FUN = \(i) {
    as.integer(i[1])
  }, 1)
scs <- scs[id_scs < 500]

loadSupEnhBed <- function(sc) {
  suffix <- "superEhancers.H3K27ac.bed"
  data.table::fread(
    file = file.path(supEnhBedir, str_glue("{sc}.{suffix}")),
    header = F, sep = "\t", data.table = F
  )
}
supEnhs <- lapply(scs, loadSupEnhBed) |>
  setNames(object = _, nm = scs)

# ** load enhancers
ptpeakd <- file.path(
  projd, "data", "pairedtag_peak",
  "subclass_peak_bed"
)
loadEnhBed <- function(sc) {
  suffix <- "H3K27ac.BestSPM.bed"
  data.table::fread(
    file = file.path(ptpeakd, str_glue("{sc}-{suffix}")),
    header = F, sep = "\t", data.table = F
  )
}
# scs is defined above in superEnh loading
enhs <- lapply(scs, loadEnhBed) |>
  setNames(object = _, nm = scs)

# ** ovlp between enh and supenhd
ovlpEnhWithSupEnh <- function(enh, supenh) {
  r1 <- loadGRfromBed(beds = enh)
  r2 <- loadGRfromBed(beds = supenh)
  r <- FindOvlpRegionInB(r1, r2) |>
    transformGR2SimpleBed(gr = _)
  return(r)
}

ovlp_enhs <- lapply(scs, \(sc) {
  r <- enhs[[sc]]
  b <- supEnhs[[sc]]
  ovlpEnhWithSupEnh(r, b)
})|>
  setNames(object = _, nm = scs)

# ** output results
outd <- file.path(
  projd, "13.cicero", "out", "H3K27acOvlpSupEnhPeak"
)
invisible({
  lapply(names(ovlp_enhs), \(sc) {
    data.table::fwrite(
      x = ovlp_enhs[[sc]],
      file = file.path(
        outd,
        str_glue("{sc}.H3K27ac.ovlpsupenh.peak.bed")
      ),
      sep = "\t", col.names = F, row.names = F
    )
  })
})


