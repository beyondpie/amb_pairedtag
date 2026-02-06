suppressMessages({
  suppressWarnings({
    library(tidyverse)
    library(future.apply)
    library(Matrix)
    library(GenomicRanges)
    library(tmpRpkg)
  })
})

# * meta
projd <- here::here()
workd <- file.path(projd, "06.ChromHMM")
chmmd <- file.path(workd, "out", "model_bypeak_b200_s18")
ordStates <- expand.grid(c("P", "D"), 1:18, stringsAsFactors = F) |>
  x => paste(x$Var2, x$Var1, sep = "_")
chrom2bin2pd <- readRDS(file.path(
  workd, "out",
  "mm10.chroms.200bin.id2pd.rds"
))
chrs <- paste0("chr", 1:19)

binaryDatad <- file.path(workd, "out", "bBed_bypeak_b200")
outd <- file.path(workd, "out", "m-bpeak-s18_pd-obs")

plan(multicore, workers = 10)
options(future.globals.maxSize = 40 * 1024^3)

# * main
## * get all the subclasses
scs <- list.files(
  path = file.path(binaryDatad),
  full.names = F, no.. = T
) |>
  y => gsub("_chr\\d+_binary\\.txt", "", y) |>
  x => unique(x)

# args <- commandArgs(trailingOnly = T)
# sc <- args[1]
invisible(
  future.apply::future_lapply(scs, \(sc) {
    message("Perform remapping ChromHMM States for ", sc)

    # * read binarized file
    binaryData.sc <- loadChromHMMBinarizedData(
      prefix = str_glue("{binaryDatad}/{sc}"),
      chrs = chrs
    )

    # * read segmentation
    segment <- readChromHMMDenseBed(
      file.path(chmmd, str_glue("{sc}_18_dense.bed"))
    ) |> x => x[order(x$chr, x$start), ]

    # * get bin2state per chrom
    bin2StateList <- lapply(chrs, \(chr) {
      chrom2Bin <- chrom2bin2pd[[chr]]
      mapBin2State.Chrom(
        segment[segment$chr == chr, ], chrom2Bin
      )
    }) |> setNames(object = _, nm = chrs)

    # * merge bin2pd, state, binary data
    bin2all <- lapply(chrs, \(chr) {
      bin2state <- bin2StateList[[chr]]
      bin2pd <- chrom2bin2pd[[chr]]
      bin2data <- binaryData.sc[[chr]]
      dplyr::inner_join(bin2state, bin2pd, by = c("id")) |>
        dplyr::inner_join(x = _, bin2data, by = c("id")) |>
        x => x[, c(
          "chr", "start", "end", "id", "relaTSS", "state",
          "ATAC", "H3K27ac", "H3K27me3", "H3K4me1", "H3K9me3"
        )]
    }) |>
      setNames(object = _, nm = chrs)

    # * save resuls as txt
    invisible(lapply(chrs, \(chr) {
      r <- bin2all[[chr]]
      data.table::fwrite(
        x = r,
        file = file.path(outd, str_glue("{sc}_{chr}_s18-200bin-pd-obs.csv")),
        col.names = T,
        row.names = F,
        sep = ","
      )
    }))
  }) ## end of future.apply
)


## ## * get emission count per chrom
## state2emitCountList <- lapply(chrs, \(chr) {
##   getEmissionCount.Chrom(
##     bin2PromoterDistalList[[chr]],
##     bin2StateList[[chr]],
##     binaryData.sc[[chr]]
##   )
## }) |> setNames(object = _, nm = chrs)

## ## * get avg emission
## state2avgEmit <- calAvgEmissionMat(state2emitCountList) |>
##   x => x[ordStates, ]
