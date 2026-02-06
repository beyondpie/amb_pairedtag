library(tidyverse)
library(tmpRpkg)
library(matrixStats)
library(paletteer)

# * meta
projd <- tmpRpkg:::projd
workd <- file.path(projd, "06.ChromHMM")
# figd <- file.path(workd, "figure")
figd <- file.path(projd, "99.figures/out/fig4")
outd <- file.path(workd, "out", "phastConsMatSum")
if (!dir.exists(outd)) {
  dir.create(outd)
}
chromHMMStates <- tmpRpkg:::stateNameOrd
ptscMeta <- tmpRpkg::loadptscMeta()
ptscs <- ptscMeta[ptscMeta$ATAC > 0, "PairedTagName"]

# * functions
loadPhastConsMat <- function(fnm) {
  r <- data.table::fread(fnm, sep = "\t", skip = 1, data.table = F,
    header = F)
  x <- as.matrix(r[, 7:ncol(r)])
  rownames(x) <- r[,4]
  colnames(x) <- NULL
  return(x)
}

loadPhastConsMats4Subclass <- function(sc) {
  lapply(chromHMMStates, \(s) {
    fnm <- file.path(scCREPhastConsMatd,
      str_glue("{sc}.{s}.mat.gz"))
    if (file.exists(fnm)) {
      loadPhastConsMat(fnm)
    } else {
      NULL
    }
  }) |>
    setNames(object = _, nm = chromHMMStates) |>
    tmpRpkg::filterNULLfromList()
}

loadAllCREPhastConsMat4sc <- function(sc) {
  loadPhastConsMat(
    fnm = file.path(
      scCREPhastConsMatd,
      str_glue("{sc}.autosome.CRE.mat.gz")
    ))
}

getMeansOfPhastCons4State <- function(state, sc2matList) {
  r <- lapply(sc2matList, \(xs) {
    if (state %in% names(xs)) {
      colMeans(xs[[state]])
    } else {
      NULL
    }
  }) |>
    tmpRpkg::filterNULLfromList()
  if (!is.null(r)) {
    return(do.call(what = rbind, r))
  }
  return(NULL)
}


# * main
# 1. load PhastCons Matrix
scCREPhastConsMatd <- file.path(workd, "out",
  "subclassCREStateComputeMatrix")
## bgCREPhastConsMatfnm <- file.path(workd, "out",
##   "mm10.backgroundRanges.computeMatrix.mat.gz")
bgCREPhastConsMatfnm <- file.path(workd, "out",
  "mm10.random2.bg.computeMatrix.mat.gz")

sc2mats <- lapply(ptscs, loadPhastConsMats4Subclass) |>
  setNames(object = _, nm = ptscs)
saveRDS(sc2mats,
  file = file.path(outd, "sc2phastConsMat2.rds"))
bgmat <- loadPhastConsMat(bgCREPhastConsMatfnm)
saveRDS(bgmat, file = file.path(outd, "bgPhastConsMat2.rds"))
bgmean <- colMeans(bgmat)

sc2allCREmat <- lapply(ptscs, loadAllCREPhastConsMat4sc) |>
  setNames(object = _, nm = ptscs)

sc2allCREmeanOfscore <- lapply(ptscs, \(sc) {
  colMeans(sc2allCREmat[[sc]])
})

allCREmeanOfscroe <- do.call(rbind, sc2allCREmeanOfscore) |>
  colMeans(x = _)

# 2. plot as mean on one figure
state2meanOfscore <- lapply(chromHMMStates,
  getMeansOfPhastCons4State,
  sc2matList = sc2mats) |>
  setNames(object = _, nm = chromHMMStates) |>
  tmpRpkg::filterNULLfromList()
state2meanOfscore[["Shuffle"]] <- bgmean
state2meanOfscore[["AllCRE"]] <- allCREmeanOfscroe

## organize data for ggplot
meanOfscore <- lapply(names(state2meanOfscore), \(state) {
  m <- state2meanOfscore[[state]]
  y <- if (!is.null(dim(m))) {
    colMeans(m)
  } else {
    m
  }
  data.frame(
    x = seq(1, 500, 10),
    y = y,
    label = state
  )
}) |> do.call(rbind, args = _)
saveRDS(
  object = meanOfscore,
  file = file.path(projd, "99.figures/out/fig4",
    "panelB.meanOfPhastConsScore2.rds")
)


p <- ggplot(
  data = meanOfscore[meanOfscore$label %in% c(
    "Chr-A", "Chr-O", "AllCRE", "Shuffle"), ],
  mapping = aes(x = x, y = y, color = label)
) +
  geom_line() +
  geom_point() +
  scale_color_brewer(palette = "Set1") +
  ggtitle("PhastCons Score")
ggsave(plot = p,
  filename = file.path(figd, "panelB.meanOfPhastCons2.pdf"))








