library(tidyverse)
library(tmpRpkg)

# * meta
projd <- here::here()
workd <- file.path(projd, "17.repressiveMarks")
figd <- file.path(workd, "figure")
outd <- file.path(workd, "out")
highTEscfnm <- file.path(projd, "meta", "highTESubclass.csv")
highTEptscs <- data.table::fread(file = highTEscfnm,
  header = T, sep = ",", data.table = F) |>
  x => x$PairedTagName

# * functions
loadPeakOvlpTE <- function(mod = "H3K27me3") {
  data.table::fread(
    file = file.path(outd, str_glue("pt.{mod}.peakOvlpTE.ratio.csv")),
    header = T, sep = ",", data.table = F
  ) |>
    x => `rownames<-`(x, x$subclass)
}

getOvlpRatio <- function(ovlpMat, scs, TEs) {
  m <- ovlpMat[scs, TEs, drop = F]
  rowSums(m) * 100 / ovlpMat[scs, "total"]
}

# * main
TEs <- c("LINE", "SINE", "LTR", "Satellite", "Simple_repeat")
matK27me3OvlpTE <- loadPeakOvlpTE("H3K27me3")
highTEK27me3Ratio <- getOvlpRatio(
  matK27me3OvlpTE, intersect(rownames(matK27me3OvlpTE),highTEptscs), c("LINE"))
othersK27me3Ratio <- getOvlpRatio(
  matK27me3OvlpTE, setdiff(rownames(matK27me3OvlpTE), highTEptscs), c("LINE")
)

## quantile(highTEK27me3Ratio)
##       0%      25%      50%      75%     100% 
## 16.17102 23.20091 28.97135 31.56001 40.33832 
## > quantile(othersK27me3Ratio)
##       0%      25%      50%      75%     100% 
## 10.35276 18.67904 20.15930 23.11226 41.19977 

matH3K9me3OvlpTE <- loadPeakOvlpTE("H3K9me3")
highTEH3K9me3Ratio <- getOvlpRatio(
  matH3K9me3OvlpTE, intersect(rownames(matH3K9me3OvlpTE),highTEptscs), c("LINE"))
othersH3K9me3Ratio <- getOvlpRatio(
  matH3K9me3OvlpTE, setdiff(rownames(matH3K9me3OvlpTE), highTEptscs), c("LINE")
)

## quantile(highTEH3K9me3Ratio)
##       0%      25%      50%      75%     100% 
## 44.00193 47.21110 48.37051 49.26789 52.61352 
## > quantile(othersH3K9me3Ratio)
##       0%      25%      50%      75%     100% 
## 26.99309 45.49490 47.46946 50.14706 57.56020 




matK27me3OvlpTE <- loadPeakOvlpTE("H3K27me3")
highTEK27me3Ratio <- getOvlpRatio(
  matK27me3OvlpTE, intersect(rownames(matK27me3OvlpTE),highTEptscs), c("SINE"))
othersK27me3Ratio <- getOvlpRatio(
  matK27me3OvlpTE, setdiff(rownames(matK27me3OvlpTE), highTEptscs), c("SINE")
)

quantile(highTEK27me3Ratio)
##       0%      25%      50%      75%     100% 
## 23.78029 33.69406 35.57700 37.05390 39.42485 

quantile(othersK27me3Ratio)
##       0%      25%      50%      75%     100% 
## 15.17929 24.14251 27.74711 31.20280 40.39401 

matH3K9me3OvlpTE <- loadPeakOvlpTE("H3K9me3")
highTEH3K9me3Ratio <- getOvlpRatio(
  matH3K9me3OvlpTE, intersect(rownames(matH3K9me3OvlpTE),highTEptscs), c("SINE"))
othersH3K9me3Ratio <- getOvlpRatio(
  matH3K9me3OvlpTE, setdiff(rownames(matH3K9me3OvlpTE), highTEptscs), c("SINE")
)

quantile(highTEH3K9me3Ratio)
##       0%      25%      50%      75%     100% 
## 10.49632 12.03153 13.44290 17.38789 22.37284 
quantile(othersH3K9me3Ratio)
##       0%      25%      50%      75%     100% 
## 10.49632 12.03153 13.44290 17.38789 22.37284 



matK27me3OvlpTE <- loadPeakOvlpTE("H3K27me3")
highTEK27me3Ratio <- getOvlpRatio(
  matK27me3OvlpTE, intersect(rownames(matK27me3OvlpTE),highTEptscs), c("LTR"))
othersK27me3Ratio <- getOvlpRatio(
  matK27me3OvlpTE, setdiff(rownames(matK27me3OvlpTE), highTEptscs), c("LTR")
)

quantile(highTEK27me3Ratio)
##       0%      25%      50%      75%     100% 
## 18.94380 28.43762 32.07180 34.65341 37.14411 

quantile(othersK27me3Ratio)
##       0%      25%      50%      75%     100% 
## 13.57362 20.54970 23.82160 26.98153 39.13591 

matH3K9me3OvlpTE <- loadPeakOvlpTE("H3K9me3")
highTEH3K9me3Ratio <- getOvlpRatio(
  matH3K9me3OvlpTE, intersect(rownames(matH3K9me3OvlpTE),highTEptscs), c("LTR"))
othersH3K9me3Ratio <- getOvlpRatio(
  matH3K9me3OvlpTE, setdiff(rownames(matH3K9me3OvlpTE), highTEptscs), c("LTR")
)

quantile(highTEH3K9me3Ratio)
##       0%      25%      50%      75%     100% 
## 22.99633 30.43500 32.63811 34.61432 37.76144 
quantile(othersH3K9me3Ratio)
##       0%      25%      50%      75%     100% 
## 19.59043 27.57943 31.76435 34.82191 42.43885 

# * read CREovlpTE
