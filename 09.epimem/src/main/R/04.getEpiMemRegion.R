suppressPackageStartupMessages({
  library(tidyverse)
  library(purrr)
  library(S4Vectors)
  library(GenomicRanges)
  # library(BiocManager)
  # BiocManager::install("trackViewer")
  library(trackViewer)
  library(tmpRpkg)
})
# For reload updated package
# detach("package:tmpRpkg", unload = TRUE)
Sys.setenv("_R_USE_PIPEBIND_" = TRUE)

# * meta
projd <- "/projects/ps-renlab2/szu/projects/amb_pairedtag/"
epimemd <- file.path(projd, "09.epimem")
epimem_outd <- file.path(epimemd, "out")
outd <- file.path(epimem_outd, "actRegionWithInactSignal")
extsize <- 2500
fold <- 2
upQ <- 0.995
# sc <- "326_OPC_NN"
# inactive <- c("H3K27me3", "H3K9me3")
# active <- c("H3K27ac", "H3K4me1", "ATAC")
# suffix <- "intersect_all_inact2act"

# * set arguments
args <- commandArgs(trailingOnly = TRUE)
sc <- args[1]
inactive <- str_split_1(args[2], ",")
active <- str_split_1(args[3], ",")
suffix <- args[4]
# - mutual: then use all the combinations
# - sequence: then treat the combination in sequential
#   from inactive and active, and they must have the same length.
combine <- args[5]
# TODO: filter chrX, chrY

message("Subclass: ", sc)
message("Inactive: ", args[2])
message("Active: ", args[3])
message("Suffix: ", suffix)
message("Combine strategy: ", combine)

# * main
if (combine == "mutual") {
  ina2a <- lapply(inactive, \(a) {
    lapply(active, \(b) {
      c(a, b)
    }) |> x => do.call(rbind, args = x)
  }) |>
    x => do.call(rbind, args = x) |>
  as.data.frame() |>
  setNames(object = _, nm = c("inactive", "active"))
} else {
  ina2a <- data.frame(
    inactive = inactive,
    active = active
  )
}
rownames(ina2a) <- with(ina2a, paste(inactive, active, sep = "_"))

# * load signal results
loadSignals <- function(sc, from_m, on_m, extsize) {
  fnm <- file.path(
    epimem_outd, str_glue("{from_m}_on_{on_m}"),
    # NOTE: use "on" now instead of "to" before
    str_glue("{sc}_{from_m}-on-{on_m}_e{extsize}.tsv")
  )
  r <- data.table::fread(
    file = fnm, sep = "\t", header = TRUE, data.table = FALSE
  )
  return(r)
}

signals <- lapply(seq_len(nrow(ina2a)), \(i) {
  loadSignals(
    sc = sc,
    from_m = ina2a[i, 1], on_m = ina2a[i, 2], extsize = extsize
  )
}) |> setNames(object = _, nm = rownames(ina2a))

# * filter active regions with inactive signals.
getOutlier <- function(score, upQ = 0.995, fold = 2) {
  q <- quantile(score, probs = upQ)
  message(str_glue("Max: {max(score)}; Q{upQ}: {q}."))
  v <- score[score <= q]
  m <- mean(v)
  std <- sd(v)
  thres <- m + fold * std
  message(str_glue("mean {m}; std {std}; threshold {thres}."))
  message(str_glue("{sum(score >= thres)} outliers."))
  return(thres)
}

filterSignals <- lapply(signals, \(s) {
  scores <- s[, dim(s)[2]]
  o <- getOutlier(scores, upQ = upQ, fold = fold)
  s[scores >= o, ]
}) |> setNames(object = _, nm = names(signals))

filterGRs <- lapply(filterSignals, \(gr) loadGRfromBed(
  beds = gr, header = TRUE)) |>
  setNames(object = _, nm = names(signals))

# * intersect regions
extendGRs <- lapply(seq_along(filterGRs), \(i) {
  gr <- filterGRs[[i]]
  if (grepl("H3K27ac|ATAC", names(filterGRs[i]))) {
    return(resizeIgnoreStrand(gr, fix = "center", extsize = extsize))
  }
  return(resizeIgnoreStrand(gr, fix = "twoside", extsize = extsize))
}) |>
  setNames(object = _, nm = names(filterGRs))

irs <- purrr::reduce(extendGRs, filterByGR)

# * save results
nm <- names(filterGRs)[1]
raw_irs <- filterGRs[[1]]
raw_irs <- raw_irs[raw_irs$name %in% irs$name]
data.table::fwrite(
  x = transformGenomicRange2DataFrame(raw_irs),
  file = file.path(outd,
    str_glue("{sc}.{nm}.{suffix}.peak")),
  sep = "\t", row.names = FALSE, col.names = TRUE
)
