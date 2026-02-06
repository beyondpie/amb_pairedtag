suppressPackageStartupMessages({
  library(tidyverse)
  library(GenomicRanges)
  library(tmpRpkg)
})
Sys.setenv("_R_USE_PIPEBIND_" = TRUE)

# * meta
projd <- here::here()
# 5kb is good enough
binSize <- 5000
epiMemGenomeBind <- file.path(projd, "09.epimem", "out", "genomeBin")
peakd <- file.path(projd, "data", "pairedtag_peak", "subclass_peak")
# blacklist
blfnm <- file.path(projd, "meta", "mm10-blacklist.v2.bed")
blGR <- loadGRfromBed(
  bedFile = blfnm, header = FALSE,
  colnms = c("chr", "start", "end", "note")
)

# * set arguments
args <- commandArgs(trailingOnly = TRUE)
sc <- args[1]
## scs <- c("001_CLA_EPd_CTX_Car3_Glut",
##   "002_IT_EP_CLA_Glut",
##   "318_Astro_NT_NN", "326_OPC_NN", "327_Oligo_NN",
##   "334_Microglia_NN")
# sc <- "326_OPC_NN"

f <- file.path(
  epiMemGenomeBind,
  str_glue("{sc}.Epimem.H3K4me1.peak.tsv")
)

if (file.exists(f)) {
  message(f, " exist and skip it.")
  quit(save = "no", status = 0)
}

# * functions
transformRegion2GR <- function(region) {
  r <- vapply(region, \(r) {
    str_split_1(r, ":|-")
  }, c("chr", "1", "2")) |>
    t() |>
    as.data.frame() |>
    setNames(object = _, nm = c("chr", "startFrom", "endTo")) |>
    mutate(
      startFrom = as.integer(startFrom),
      endTo = as.integer(endTo),
      name = region
    ) |>
    loadGRfromBed(beds = _, header = FALSE)
  return(r)
}

loadH3K4me1Peak <- function(sc) {
  data.table::fread(
    file = file.path(
      peakd,
      str_glue("{sc}-H3K4me1.bedtoolmerge.peak")
    ),
    sep = "\t", header = TRUE, data.table = FALSE
  ) |>
    x => `rownames<-`(x, x$name)
}

getHighEpiMemGenomeBin <- function(hemGB, foldsd = 1.0, lowQtl = 0.005) {
  # remove sex chromosome
  hemGB <- hemGB[grepl("chr[1-9]+", hemGB$region), ]
  # remove blacklist region
  gr <- transformRegion2GR(hemGB$region)
  regionInbl <- FindOvlpRegionInB(gr, blGR)
  hemGB <- hemGB[!(hemGB$region %in% regionInbl$name), ]
  # remove low epi-signal regions
  lowK4 <- quantile(with(hemGB, H3K4me1[H3K4me1 > 0]), lowQtl)
  lowK27me3 <- quantile(with(hemGB, H3K27me3[H3K27me3 > 0]), lowQtl)
  hemGB <- hemGB[(hemGB$H3K4me1 > lowK4) & (hemGB$H3K27me3 > lowK27me3), ]
  # scale data
  hemGB$scaled_H3K4me1 <- scale(hemGB$H3K4me1, center = TRUE, scale = TRUE)[, 1]
  hemGB$scaled_H3K27me3 <- scale(hemGB$H3K27me3, center = TRUE, scale = TRUE)[, 1]
  # filter using signals
  cutoffK4 <- mean(hemGB$H3K4me1) + foldsd * sd(hemGB$H3K4me1)
  cutoffK27me3 <- mean(hemGB$H3K27me3) + foldsd * sd(hemGB$H3K27me3)
  hemGB <- hemGB[(hemGB$H3K4me1 >= cutoffK4) & (hemGB$H3K27me3 >= cutoffK27me3), ]
  return(hemGB)
}

# * main
# ** load highly epiMem genome bins
epiMemf <- file.path(
  epiMemGenomeBind,
  str_glue("{sc}.H3K4me1_H3K27me3.5000.csv")
)
epiMem <- data.table::fread(
  file = epiMemf, header = TRUE, data.table = FALSE, sep = ","
) |>
  x => `rownames<-`(x, x$region)
hemGB <- getHighEpiMemGenomeBin(epiMem, foldsd = 1.0, lowQtl = 0.005)
hemGR <- transformRegion2GR(hemGB$region)

# ** intersect with H3K4me1 peaks
pK4 <- loadH3K4me1Peak(sc)
pK4GR <- loadGRfromBed(beds = pK4[, 1:4])

# pK4_hem <- FindOvlpRegionInB(pK4GR, hemGR)
# pK4hemRegion <- strGRSeq(pK4_hem)
# pK4hem <- pK4[pK4_hem$name, ]

ovlpGR <- findOverlaps2(pK4GR, hemGR)
bestOvlp <- ovlpGR |>
  x => data.frame(
    pname = x$queryGR$name,
    gname = x$subjectGR$name,
    width = x$ovlpWidth
  ) |>
  group_by(pname) |>
  slice_max(width) |>
  slice_head(n = 1) |>
  as.data.frame() |>
  x => `rownames<-`(x, x$pname)

pK4hem <- pK4[bestOvlp$pname, ]
pK4hem$width <- bestOvlp$width
pK4hem <- cbind(pK4hem, hemGB[bestOvlp$gname, ])
pK4hem$harmonicMean <- with(
  pK4hem,
  2 / (1 / scaled_H3K4me1 + 1 / scaled_H3K27me3)
)

# ** output
data.table::fwrite(pK4hem, file = f, sep = "\t")
