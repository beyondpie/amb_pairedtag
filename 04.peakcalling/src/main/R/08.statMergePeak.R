library(tidyverse)
projd <- here::here()

mods <- c("H3K27ac", "H3K27me3", "H3K9me3", "H3K4me1")
bedtoolSuffix <- "merged.all.blv2.bedtools"
suffix <- "merged.all.blv2"


getTotalSizeOfPeaks <- function(peaks) {
  sum(peaks$end - peaks$start)
}

mouseChr2Size <- data.table::fread(
  file = file.path(projd, "meta", "mm10.chrom.sizes.lite"),
  header = FALSE
) |> setNames(object = _, nm = c("chr", "size"))

mm10size <- sum(mouseChr2Size$size)

# main
mergePeakd <- file.path(projd,
  "04.peakcalling", "out", "mergePeak")

mod2peaks_bdt <- lapply(
  mods, \(m, head) {
    f <- file.path(mergePeakd, str_glue("{m}.{bedtoolSuffix}.peak"))
    data.table::fread(
      file = f, header = head, sep = "\t",
      data.table = FALSE
    ) |>
      setNames(object = _, nm = c("chr", "start", "end"))
  }, head = FALSE
) |>
  setNames(object = _, nm = mods)

# get number of peaks
mod2n <- vapply(mod2peaks_bdt, nrow, 1)
## mod2n
##  H3K27ac H3K27me3  H3K9me3  H3K4me1 
##   717962   552774   473445   401383 

# get total size
mod2size <- vapply(mods, \(m) {
  getTotalSizeOfPeaks(mod2peaks_bdt[[m]])
}, 1)
##   H3K27ac   H3K27me3    H3K9me3    H3K4me1 
## 667255795 1181067759  650265990 1105272744

mod2prcnt <- mod2size * 100 / mm10size
##  H3K27ac H3K27me3  H3K9me3  H3K4me1 
## 24.48177 43.33365 23.85841 40.55271 


# get peak size
sizeOfmodPeaks <- lapply(mods, \(m) {
  quantile(with(mod2peaks_bdt[[m]], end - start))
}) |> setNames(object = _, nm = mods)

## str(sizeOfmodPeaks)
## List of 4
##  $ H3K27ac : Named num [1:5] 200 248 468 963 70836
##   ..- attr(*, "names")= chr [1:5] "0%" "25%" "50%" "75%"
##  $ H3K27me3: Named num [1:5] 200 384 944 2164 350424
##   ..- attr(*, "names")= chr [1:5] "0%" "25%" "50%" "75%"
##  $ H3K9me3 : Named num [1:5] 200 279 697 1524 164907
##   ..- attr(*, "names")= chr [1:5] "0%" "25%" "50%" "75%"
##  $ H3K4me1 : Named num [1:5] 200 466 1143 2696 323045
##   ..- attr(*, "names")= chr [1:5] "0%" "25%" "50%" "75%"

allMergedPeaks <- data.table::fread(
  file.path(mergePeakd, "all.blv2.merged.bedtools.peak"),
  head = FALSE, data.table = FALSE
) |> setNames(object = _, nm = c("chr", "start", "end"))

sum(allMergedPeaks$end - allMergedPeaks$start)
sum(allMergedPeaks$end - allMergedPeaks$start) * 100/ mm10size
# 67.22258


