## if (!require("BiocManager", quietly = TRUE)) {
##   install.packages("BiocManager")
## }
## BiocManager::install("GenomicRanges")
library(GenomicRanges)
library(tidyverse)
Sys.setenv("_R_USE_PIPEBIND_" = TRUE)


# meta
projd <- here::here()
satd <- file.path(projd, "04.peakcalling", "out", "saturation")

loadNarrowPeak <- function(npfile) {
  r <- data.table::fread(npfile,
    sep = "\t",
    header = FALSE, data.table = FALSE
  ) |>
    x => x[!grepl("chrUn|chrM|random", x[,1]), ]
  # temp ignore summit in narrowPeak
  colnames(r)[1:9] <- c(
      "chr", "left", "right", "name", "-10log10q", "strand",
      "fldchng", "-log10p", "-log10q"
  )
  return(r)
}

np2gr <- function(np) {
  allGRs <- GenomicRanges::GRanges(
    seqnames = np[, 1],
    ranges = IRanges::IRanges(
      start = np[, 2],
      end = np[, 3]
    )
  )
  chrs <- allGRs@seqnames
  allGRs[
    !grepl("random|chrM|chrUn", chrs),
  ]
}

capGR <- function(qGR, sGR) {
  hits <- GenomicRanges::findOverlaps(qGR, sGR)
  sindex <- S4Vectors::subjectHits(hits)
  length(unique(sindex)) * 100 / length(sGR)
}

runSatAnalysis <- function(mod = "H3K27ac",
                           celltype = "Sst_Chodl_Gaba",
                           modPattern = "narrowPeak") {
  peakFiles <- list.files(
    path = file.path(satd, celltype),
    pattern = modPattern,
    full.names = FALSE, include.dirs = FALSE, no.. = TRUE
  ) |>
    x => x[grepl(mod, x)]
  nds <- vapply(peakFiles, \(f){
    str_split_1(f, "_|\\.")[5] |>
      gsub(pattern = "ds", replacement = "", x = _) |>
      as.integer()
  }, 1L) |>
    sort()

  peaks <- lapply(nds, \(n) {
    loadNarrowPeak(file.path(
      satd, celltype,
      str_glue("{celltype}_{mod}_ds{n}_peaks.{modPattern}")
    ))
  }) |>
    x => setNames(x, paste0("ds", nds))

  gRs <- lapply(peaks, \(np) {
    np2gr(np)
  }) |>
    x => setNames(x, paste0("ds", nds))

  numps <- vapply(gRs, length, 1L)
  caps <- vapply(gRs, capGR, 90.0, sGR = gRs[[length(gRs)]])
  return(
    data.frame(
      npeak = numps,
      ratio = caps
    )
  )
}

# * H3K27ac
sst_H3K27ac <- runSatAnalysis(
  mod = "H3K27ac",
  celltype = "Sst_Chodl_Gaba",
  modPattern = "narrowPeak"
)
opc_H3K27ac <- runSatAnalysis(
  mod = "H3K27ac",
  celltype = "OPC_NN_1",
  modPattern = "narrowPeak"
)

# * H3K27me3
sst_H3K27me3 <- runSatAnalysis(
  mod = "H3K27me3",
  celltype = "Sst_Chodl_Gaba",
  modPattern = "broadPeak"
)
opc_H3K27me3 <- runSatAnalysis(
  mod = "H3K27me3",
  celltype = "OPC_NN_1",
  modPattern = "broadPeak"
)

# H3K9me3
sst_H3K9me3 <- runSatAnalysis(
  mod = "H3K9me3",
  celltype = "Sst_Chodl_Gaba",
  modPattern = "broadPeak"
)
opc_H3K9me3 <- runSatAnalysis(
  mod = "H3K9me3",
  celltype = "OPC_NN_1",
  modPattern = "broadPeak"
)

# H3K4me1
sst_H3K4me1 <- runSatAnalysis(
  mod = "H3K4me1",
  celltype = "Sst_Chodl_Gaba",
  modPattern = "broadPeak"
)
opc_H3K4me1 <- runSatAnalysis(
  mod = "H3K4me1",
  celltype = "OPC_NN_1",
  modPattern = "broadPeak"
)

# * plot
nReads <- c(
  1e5, 2e5, 5e5, 7.5e5, 1e6, 1.25e6, 1.5e6, 1.75e6, 2e6, 2.5e6, 4e6, 
  5e6, 1e7, 1.5e7, 2e7) |>
  as.integer()

