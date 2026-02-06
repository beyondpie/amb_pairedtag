library(tidyverse)
library(ComplexHeatmap)

# * meta
projd <- here::here()
rscd <- file.path(projd, "data", "snATAC")
hiscntd <- file.path(projd, "data", "ptHistoneCounts","ATACPeak")
histcntdr <- file.path(projd, "data", "ptHistoneCounts","ATACPeak_rds")
histones <- c("H3K27ac", "H3K4me1", "H3K27me3", "H3K9me3")

atacProj <- "/projects/ps-renlab2/szu/projects/CEMBA2"
atacCPMfnm <- file.path(atacProj,
  "18.snap2_peakcalling", "out/scfilter",
  "cpm_peakBysubclass.csv")

atacPeakSrt <- data.table::fread(
  file.path(projd, "data", "snATAC", "wmb.snATAC.peak.srt.bed"),
  header = FALSE, sep = "\t", data.table = FALSE
)$V4


# * functions
loadStatPeakNMF <- function(fnm) {
  statPeak <- data.table::fread(
    file = fnm,
    header = FALSE, sep = "\t", data.table = FALSE
  )
  colnames(statPeak) <- c(
    "peak", "index", "class0", "contributes",
    "featureScore", "selt_fs_list",
    "selt_med_list")
  statPeak$moduleN <- paste(
    "m", as.integer(statPeak$class0 + 1), sep = "")
  rownames(statPeak) <- statPeak$peak
  return(statPeak)
}

orderStatPeakByModuleRank <- function(statPeak, mod.ordered){
  index.list <- lapply(mod.ordered, function(m) {
    index.raw <- which(statPeak$moduleN %in% m)
    index <- index.raw[
      order(statPeak$featureScore[index.raw], decreasing = TRUE)]
    return(index)
  })
  index.ordered <- unlist(index.list)
  r <- statPeak[index.ordered, ]
  return(r)
}

setHeatmap <- function(mat,
                       low.quantile = 0.01,
                       show.legend = TRUE,
                       legend.direction = "horizontal",
                       show.row.names = TRUE) {
  highval <- quantile(mat, 1 - low.quantile)
  lowval <- quantile(mat, low.quantile)
  colfn <- circlize::colorRamp2(seq(lowval, highval, length = 60),
    viridis::viridis(60))
  ComplexHeatmap::Heatmap(
    matrix = mat,
    col = colfn,
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    show_row_names = show.row.names,
    show_column_names = FALSE,
    use_raster = TRUE,
    top_annotation = NULL,
    left_annotation = NULL,
    show_heatmap_legend = show.legend,
    heatmap_legend_param = list(
      title = latex2exp::TeX(r"($\log$(CPM+1))"),
      at = c(lowval, highval),
      labels = c(round(lowval, 1), round(highval, 1)),
      direction = legend.direction)
  )
}

combineHeatmap <- function(logCPM, mCG, low.quantile = 0.01) {
  hmATAC <- setHeatmap(mat = logCPM)
  hmCG <- setHeatmap(mat = mCG)
  ComplexHeatmap::draw(
    hmATAC + hmCG,
    heatmap_legend_side = "bottom",
    merge_legend = TRUE)
}

readRawFeatureCount <- function(fnm) {
  data.table::fread(file = fnm,
    skip = 1, header = TRUE, sep = "\t", data.table = FALSE)
}

readFeatureCount <- function(fnm) {
  data.table::fread(file = fnm,
    skip = 1, header = TRUE, sep = "\t",
    drop = 1:6, data.table = FALSE)[, 1]
}

readFeatureCountMergeSex <- function(h, g) {
  malefnm <- file.path(hiscntd,
    str_glue("{g}_{h}_Male.featureCount.txt"))
  femalefnm <- file.path(hiscntd,
    str_glue("{g}_{h}_Female.featureCount.txt"))
  fe <- file.exists(femalefnm)
  me <- file.exists(malefnm)
  if(fe & me) {
    m <- readFeatureCount(malefnm)
    f <- readFeatureCount(femalefnm)
    return(m + f)
  }
  if (fe) {
    return(readFeatureCount(femalefnm))
  }
  return(readFeatureCount(malefnm))
}

getGroups4Histone <- function(fnms, h) {
  suf <- "featureCount.txt"
  a <- str_replace(fnms[grepl(h, fnms)],
    str_glue("_{h}_Male.{suf}|_{h}_Female.{suf}"), "")
  return(unique(a))
}

sortGroups <- function(groups) {
  ids <- vapply(groups, \(g) {
    as.integer(str_split_1(g, pattern = "_")[1])
  }, 0)
  groups[order(ids)]
}

downsample.index.1 <- function(index.all,
                               mod.index,
                               mod.ordered = NULL,
                               score.index = NULL,
                               size = 200,
                               seed = 2022) {
  set.seed(seed = seed)
  mods <- if(!is.null(mod.ordered)) {
    mod.ordered
  } else {
    unique(mod.index)
  }
  index.list <- lapply(mods, function(m) {
    index <- mod.index %in% m
    if(sum(index) < 1) {
      warning("Module ", m, " has no features.")
      return(NULL)
    }
    cols <- which(index)
    cols.sampled <- if(!is.null(score.index)) {
      cols[order(score.index[cols], decreasing = TRUE)][
        1:min(length(cols), size)]
    } else {
      sample(cols,size = min(length(cols), size),
        replace = FALSE)
    }
    return(cols.sampled)
  })
  index.list[sapply(index.list, is.null)] <- NULL
  return(unlist(index.list))
}

downsampleStatPeak.1 <- function(statPeak,
                                 mod.ordered = NULL,
                                 size = 200,
                                 seed = 2022) {
  index.sampled <- downsample.index.1(
    index.all = seq_len(nrow(statPeak)),
    mod.index = statPeak$moduleN,
    mod.ordered = mod.ordered,
    score.index = statPeak$featureScore,
    size = size,
    seed = seed
  )
  statPeak.sampled <- statPeak[index.sampled, ]
  return(statPeak.sampled)
}

fastread.csv <- function(fnm, header = TRUE, rowname = TRUE) {
  r <- data.table::fread(file = fnm, sep = ",",
    header = header, data.table = FALSE)
  if (rowname) {
    rownames(r) <- r$V1
    r$V1 <- NULL
  }
  return(r)
}
read.sa2CPM.pbysc <- function(cap = 0.9999,
                              runlog1p = TRUE,
                              scbyp = TRUE) {
  r <- fastread.csv(
    atacCPMfnm,
    header = TRUE, rowname = TRUE)
  r <- as.matrix(r)
  upValue <- quantile(r, cap)
  message("quantile of ", cap, " : ", upValue)
  r[r > upValue] <- upValue
  if (runlog1p) {
    message("log1p transformation.")
    r <- log1p(r)
  }
  if (scbyp) {
    message("output to subclass by peak matrix.")
    r <- t(r)
  }
  return(r)
}

# * load atac and mCG signals

# * get NMF results for all the peaks
allStatPeak <- loadStatPeakNMF(
  fnm = file.path(rscd, "snATACAllPeakNMF",
    "nmfPmat.allpeak.r150.n0.statW"))
modOrd <- data.table::fread(
  file = file.path(rscd, "snATACAllPeakNMF",
    "sa2.allpeak.nmf.module.order.txt"),
  head = FALSE, data.table = FALSE
)$V1
statPeakOrd <- orderStatPeakByModuleRank(allStatPeak, modOrd)
data.table::fwrite(x = statPeakOrd,
  file = file.path(rscd, "snATACAllPeakNMF",
    "sa2.allpeak.nmf.peak.order.csv"),
  col.names = TRUE, row.names = FALSE, sep = ",")

# * get histone CPM
hiscntfnms <- list.files(hiscntd) |>
  x => x[grepl(pattern = "featureCount.txt$", x = x)]

tmp <- readRawFeatureCount(file.path(hiscntd, hiscntfnms[1]))

gH3K27ac <- getGroups4Histone(hiscntfnms, h = "H3K27ac")
gH3K27me3 <- getGroups4Histone(hiscntfnms, h = "H3K27me3")
gH3K4me1 <- getGroups4Histone(hiscntfnms, h = "H3K4me1")
gH3K9me3 <- getGroups4Histone(hiscntfnms, h = "H3K9me3")

# 173 groups
groups <- intersect(gH3K27ac, gH3K27me3) |>
  x => intersect(x = x, y = gH3K4me1) |>
  x => intersect(x = x, y = gH3K9me3)

saveFeatureCount <- function(h) {
  r <- vapply(groups, \(g) {
    message("read featureCounts for ", g)
    readFeatureCountMergeSex(h = h, g = g)
  }, FUN.VALUE = rep(0, nrow(tmp)))
  rownames(r) <- tmp$Geneid
  saveRDS(r, file = file.path(histcntdr,
    str_glue("featureCount.{h}.mat.rds")))
  return(r)
}

fcH3K27ac <- saveFeatureCount(h = "H3K27ac")
fcH3K27me3 <- saveFeatureCount(h = "H3K27me3")
fcH3K4me1 <- saveFeatureCount(h = "H3K4me1")
fcH3K9me3 <- saveFeatureCount(h = "H3K9me3")

loadFeatureCount <- function(h) {
  readRDS(file.path(histcntdr, str_glue("featureCount.{h}.mat.rds")))
}

fcH3K27ac <- loadFeatureCount("H3K27ac")
fcH3K27me3 <- loadFeatureCount(h = "H3K27me3")
fcH3K4me1 <- loadFeatureCount(h = "H3K4me1")
fcH3K9me3 <- loadFeatureCount(h = "H3K9me3")




# * plot histone heatmap
nds <- 10000
nmod <- length(unique(allStatPeak$moduleN))
scale.factor <- 1e6
cap <- 0.99
gOrd <- sortGroups(colnames(fcH3K27ac))

allStatPeak.ds <- downsampleStatPeak.1(
  statPeak = statPeakOrd,
  mod.ordered = NULL,
  size = floor(nds / nmod),
  seed = 2024
)
peaks.ds <- allStatPeak.ds$peak

drawHisHM <- function(cnt, h, low.quantile = 0.01) {
  r <- t(t(cnt * scale.factor) / colSums(cnt))
  upValue <- quantile(r, cap)
  message("quantile of ", cap, " : ", upValue)
  r[r > upValue] <- upValue
  logCPM <- r |>
    log1p() |>
    # scale() |>
    t()
  
  # oldPeaks <- colnames(logCPM)
  colnames(logCPM) <- atacPeakSrt
  logCPM.plot <- logCPM[gOrd, peaks.ds]
  pHM <- setHeatmap(mat = logCPM.plot,
    low.quantile = low.quantile,
    show.legend = TRUE,
    legend.direction = "horizontal")

  withr::with_pdf(
    new = file.path(projd, "figure",
      str_glue("{h}.onATAC.ATACnmf.all.ds.heatmap.pdf")),
    code = {
      print(pHM)
    },
    height = 20,
    width = 20
  )
  return(pHM)
}

drawHisHM(cnt = fcH3K27ac, "H3K27ac")
drawHisHM(cnt = fcH3K27me3, "H3K27me3")
drawHisHM(cnt = fcH3K4me1, "H3K4me1")
drawHisHM(cnt = fcH3K9me3, "H3K9me3")

# * test ATAC
## ggOrd <- str_replace(gOrd, "^..._", "")
## r <- read.sa2CPM.pbysc(cap = 0.9999)
## oldrnm <- rownames(r)
## rownames(r) <- str_replace_all(oldrnm, "-", "_")
## logCPM <- r
## logCPM.plot <- logCPM[ggOrd[ggOrd %in% rownames(logCPM)], peaks.ds]

## pHM <- setHeatmap(mat = logCPM.plot,
##   low.quantile = 0.01,
##   show.legend = TRUE,
##   legend.direction = "horizontal")

## withr::with_pdf(
##   new = file.path(projd, "figure",
##     str_glue("ATAC.onATAC.ATACnmf.all.ds.heatmap.pdf")),
##   code = {
##     print(pHM)
##   },
##   height = 20,
##   width = 20
## )

