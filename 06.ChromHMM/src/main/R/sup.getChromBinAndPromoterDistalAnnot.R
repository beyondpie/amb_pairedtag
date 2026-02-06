suppressMessages({
  suppressWarnings({
    library(tidyverse)
    library(Matrix)
    library(GenomicRanges)
    library(tmpRpkg)
  })
})

# * meta
projd <- here::here()
workd <- file.path(projd, "06.ChromHMM")
mm10chromSizefnm <- file.path(
  projd,
  "meta", "mm10.chrom.sizes.lite"
)
chromSizes <- data.table::fread(
  mm10chromSizefnm,
  header = F, data.table = F
) |>
  setNames(object = _, nm = c("chr", "size")) |>
  x => `rownames<-`(x, x$chr)
## chromSizes <- chromSizes[
##   !(rownames(chromSizes) %in% c("chrX", "chrY")),
## ]
chrs <- chromSizes$chr

homerd <- file.path("/projects/ps-renlab2/szu/softwares/homer")
homermm10Annotd <- file.path(homerd, "data/genomes", "mm10")
homermm10TSSfnm <- file.path(homermm10Annotd, "mm10.tss")
homerTSS <- readHomerTSS(homermm10TSSfnm, bufferSize = 2000) |>
  x => x[x$chr %in% chromSizes$chr, ]
# order homerTSS by chrom and start
homerTSS <- homerTSS[order(homerTSS$chr, homerTSS$start), ]

# * get chrom2bin
chrom2Bins <- lapply(seq_len(nrow(chromSizes)), \(i) {
  chrom <- chromSizes[i, 1]
  chromSize <- chromSizes[i, 2]
  r <- partitionChrom(
    chrName = chrom,
    chromSize = chromSize,
    binSize = 200,
    excludeLastPartialBin = T
  )
  return(r)
}) |>
  setNames(object = _, nm = chromSizes[, 1])
saveRDS(
  chrom2Bins,
  file.path(workd, "out", "mm10.chroms.200bin.rds")
)

# * get chromBin2PromoterDistal
bin2PromoterDistalList <- lapply(chrs, \(chr) {
  chrom2Bin <- chrom2Bins[[chr]]
  mapBin2PromoterDistal.Chrom(chrom2Bin, homerTSS)
}) |> setNames(object = _, nm = chrs)

saveRDS(bin2PromoterDistalList, file.path(
  workd, "out",
  "mm10.chroms.200bin.promoter-distal.rds"
))

# * merge chrom2Bin and Bin2PromoterDistal
chrom2Bin2PQ <- lapply(chrs, \(chr) {
  chrom2Bin <- chrom2Bins[[chr]]
  bin2PQ <- bin2PromoterDistalList[[chr]]
  chrom2Bin2PQ <- dplyr::inner_join(
    x = chrom2Bin, y = bin2PQ, by = c("id")
  ) |>
    x => x[, c("chr", "start", "end", "id", "relaTSS")] |>
    x => `rownames<-`(x, paste0("id", x$id))
}) |>
  setNames(object = _, nm = chrs)

saveRDS(
  chrom2Bin2PQ,
  file.path(workd, "out", "mm10.chroms.200bin.id2pd.rds")
)
