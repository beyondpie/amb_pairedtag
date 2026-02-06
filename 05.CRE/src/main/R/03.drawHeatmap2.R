library(tidyverse)
library(ComplexHeatmap)
library(devtools)
load_all(file.path(here::here(), "tmpRpkg"))

# * function
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

orderStatPeakByModuleRank <- function(statPeak, mod.ordered) {
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

# * meta
atacPeakSrt <- loadATACPeakSrt()
peakInPtSrt <- loadATACPeakInPtSrt()
allenClMeta <- loadAllenClMeta()
sc2scil <- load_sc2scil()
sp2sc <- load_sp2sc()

# === MAIN ===
# * redraw using the results from H3K27ac NMF
H3K27acNMFMeta <- readRDS(
  file.path(projd, "05.CRE", "out", "nmf", "H3K27ac.r150.meta.rds")
)
scs <- rownames(H3K27acNMFMeta$logCPM4plot)
peaks <- colnames(H3K27acNMFMeta$logCPM4plot)

cpm_scbyp_ATAC <- t(cpm_pbysc_ATAC) |>
  x => x[scs, peaks] |>
  x => log1p(x)

pATACHM <- setHeatmap(
  mat = cpm_scbyp_ATAC,
  ha_col = NULL,
  ha_row = NULL,
  showRowNames = TRUE,
  lowqtl = 0.01,
  highqtl = 0.99,
  legend_show = TRUE,
  legend_direct = "horizontal")

withr::with_pdf(
  new = file.path(projd, "05.CRE", "out/nmf",
    "snATAC.heatmap.onH3K27ac.pdf"),
  code = {
    print(pATACHM)
  },
  height = 28,
  width = 20
)

# test scale version
# row-levle scaled
cpm_scbyp_ATAC <-
  cpm_pbysc_ATAC[peaks, scs] |>
  log1p() |>
  scale() |>
  t()

pATACHM <- setHeatmap(
  mat = cpm_scbyp_ATAC,
  ha_col = NULL,
  ha_row = NULL,
  showRowNames = TRUE,
  lowqtl = 0.01,
  highqtl = 0.99,
  legend_show = TRUE,
  legend_direct = "horizontal")

withr::with_pdf(
  new = file.path(projd, "05.CRE", "out/nmf",
    "snATAC.heatmap.onH3K27ac.localScaleByRow.pdf"),
  code = {
    print(pATACHM)
  },
  height = 28,
  width = 20
)

cpm_scbyp_ATAC <-
  cpm_pbysc_ATAC |>
  log1p() |>
  scale() |>
  t() |>
  x => x[scs, peaks]

pATACHM <- setHeatmap(
  mat = cpm_scbyp_ATAC,
  ha_col = NULL,
  ha_row = NULL,
  showRowNames = TRUE,
  lowqtl = 0.01,
  highqtl = 0.99,
  legend_show = TRUE,
  legend_direct = "horizontal")

withr::with_pdf(
  new = file.path(projd, "05.CRE", "out/nmf",
    "snATAC.heatmap.onH3K27ac.globalScaleByRow.pdf"),
  code = {
    print(pATACHM)
  },
  height = 28,
  width = 20
)

# column scale 
cpm_scbyp_ATAC <-
  cpm_pbysc_ATAC[peaks, scs] |>
  log1p() |>
  t() |>
  scale() |>

pATACHM <- setHeatmap(
  mat = cpm_scbyp_ATAC,
  ha_col = NULL,
  ha_row = NULL,
  showRowNames = TRUE,
  lowqtl = 0.01,
  highqtl = 0.99,
  legend_show = TRUE,
  legend_direct = "horizontal")

withr::with_pdf(
  new = file.path(projd, "05.CRE", "out/nmf",
    "snATAC.heatmap.onH3K27ac.localScaleByCol.pdf"),
  code = {
    print(pATACHM)
  },
  height = 28,
  width = 20
)


# redraw H3K4me1
hm <- "H3K4me1"
fc <- readRDS(file.path(hiscntdr,
  str_glue("featureCount.{hm}.mat.rds")))
rownames(fc) <- atacPeakSrt
fc_sc <- mapfc2fcsc(fc, sp2sc)

cpm_pbysc_hm <- getCPMOfHistones(fc_sc,
  scale.factor = 1e6)

cpm_scbyp_hm <- t(cpm_pbysc_hm) |>
  x => x[scs[scs %in% rownames(x)], peaks] |>
  x => log1p(x)

pHM <- setHeatmap(
  mat = cpm_scbyp_hm,
  ha_col = NULL,
  ha_row = NULL,
  showRowNames = TRUE,
  lowqtl = 0.01,
  highqtl = 0.99,
  legend_show = TRUE,
  legend_direct = "horizontal")

withr::with_pdf(
  new = file.path(projd, "05.CRE", "out/nmf",
    str_glue("snATAC.heatmap.{hm}.pdf")),
  code = {
    print(pHM)
  },
  height = 28,
  width = 20
)

# scale by row
hm <- "H3K4me1"
fc <- readRDS(file.path(hiscntdr,
  str_glue("featureCount.{hm}.mat.rds")))
rownames(fc) <- atacPeakSrt
fc_sc <- mapfc2fcsc(fc, sp2sc)

cpm_pbysc_hm <- getCPMOfHistones(fc_sc,
  scale.factor = 1e6)

cpm_scbyp_hm <-
  log1p(cpm_pbysc_hm) |>
  scale() |>
  t() |>
  x => x[scs, peaks]

pHM <- setHeatmap(
  mat = cpm_scbyp_hm,
  ha_col = NULL,
  ha_row = NULL,
  showRowNames = TRUE,
  lowqtl = 0.01,
  highqtl = 0.99,
  legend_show = TRUE,
  legend_direct = "horizontal"
  )

withr::with_pdf(
  new = file.path(projd, "05.CRE", "out/nmf",
    str_glue("snATAC.heatmap.{hm}.globalScaleByRow.pdf")),
  code = {
    print(pHM)
  },
  height = 28,
  width = 20
)


cpm_scbyp_hm <- cpm_pbysc_hm[peaks, scs] |>
  log1p() |>
  scale() |>
  t()

pHM <- setHeatmap(
  mat = cpm_scbyp_hm,
  ha_col = NULL,
  ha_row = NULL,
  showRowNames = TRUE,
  lowqtl = 0.01,
  highqtl = 0.99,
  legend_show = TRUE,
  legend_direct = "horizontal"
  )

withr::with_pdf(
  new = file.path(projd, "05.CRE", "out/nmf",
    str_glue("snATAC.heatmap.{hm}.localScaleByRow.pdf")),
  code = {
    print(pHM)
  },
  height = 28,
  width = 20
)

# === Draw heatmap on the H3K27ac-identified cell-type specific peaks ====


# * unify snATAC, H3K27ac, H3K4me1 CPM matrix
# i.e., save row and colnames
# - load data
# -- ATAC-seq data
cpm_pbysc_ATAC <- load_snATACPMpbyc()

# -- H3K27ac
fcH3K27ac <- readRDS(file.path(histcntdr, "featureCount.H3K27ac.mat.rds"))
fcH3K27ac_sc <- mapfc2fcsc(fcH3K27ac, sp2sc)
rownames(fcH3K27ac_sc) <- atacPeakSrt

# remove subclasses with overall low CPM on H3K27ac
rm_scs <- c(
  # "041_OB_in_Frmd7_Gaba",
  "042_OB_out_Frmd7_Gaba",
  "044_OB_Dopa_Gaba",
  "065_IA_Mgp_Gaba",
  "066_NDB_SI_ant_Prdm12_Gab",
  "070_LSX_Prdm12_Slit2_Gaba",
  "083_CEA_BST_Rai14_Pdyn_Crh_Gaba",
  "084_BST_SI_AAA_Six3_Slc22a3_Gaba",
  "091_ARH_PVi_Six6_Dopa_Gaba",
  "109_LGv_ZI_Otx2_Gaba",
  "137_PH_an_Pitx2_Glut",
  "138_PH_Pitx2_Glut",
  "141_PH_SUM_Foxa1_Glut"
)

sc.final <- intersect(colnames(cpm_pbysc), colnames(fcH3K27ac_sc)) |>
  x => setdiff(x, rm_scs) |>
  x => sortAllenLabelById(x)

# -- H3K4me1
fcH3K4me1 <- readRDS(file.path(histcntdr, "featureCount.H3K4me1.mat.rds"))
fcH3K4me1_sc <- mapfc2fcsc(fcH3K4me1, sp2sc)
rownames(fcH3K4me1_sc) <- atacPeakSrt

# unify data
cpmATAC_pbysc <- cpm_pbysc[peakInPtSrt, sc.final]
cpmATAC_pbysc <- as.matrix(cpmATAC_pbysc)
fcH3K27ac_pbysc <- fcH3K27ac_sc[peakInPtSrt, sc.final]
cpmH3K27ac_pbysc <- getCPMOfHistones(fcH3K27ac_pbysc)
fcH3K4me1_pbysc <- fcH3K4me1_sc[peakInPtSrt, sc.final]
cpmH3K4me1_pbysc <- getCPMOfHistones(fcH3K4me1_pbysc)

# get NMF result and downsamle
allStatPeak <- loadStatPeakNMF(
  fnm = file.path(snATACd, "snATACAllPeakNMF",
    "nmfPmat.allpeak.r150.n0.statW"))
modOrd <- data.table::fread(
  file = file.path(snATACd, "snATACAllPeakNMF",
    "sa2.allpeak.nmf.module.order.txt"),
  head = FALSE, data.table = FALSE
)$V1
statPeakOrd <- orderStatPeakByModuleRank(allStatPeak, modOrd)
statPeakOrd <- statPeakOrd[statPeakOrd$peak %in% peakInPtSrt, ]


nds <- 10000
nmod <- length(unique(statPeakOrd$moduleN))
scale.factor <- 1e6
cap <- 0.99
gOrd <- colnames(cpmH3K27ac_pbysc)

allStatPeak.ds <- downsampleStatPeak.1(
  statPeak = statPeakOrd,
  mod.ordered = NULL,
  size = floor(nds / nmod),
  seed = 2024
)
peaks.ds <- allStatPeak.ds$peak

# r is cpm
drawHisHM <- function(r, h, low.quantile = 0.01, cap = 0.999) {
  # r <- t(t(cnt * scale.factor) / colSums(cnt))
  upValue <- quantile(r, cap)
  message("quantile of ", cap, " : ", upValue)
  r[r > upValue] <- upValue
  logCPM <- r |>
    log1p() |>
    # scale() |>
    t()
  
  logCPM.plot <- logCPM[, peaks.ds]
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

drawHisHM(r = cpmATAC_pbysc, h = "snATAC", low.quantile = 0.01, cap = 0.999)
drawHisHM(r = cpmH3K27ac_pbysc, h = "H3K27ac", low.quantile = 0.05, cap = 0.999)
drawHisHM(r = cpmH3K4me1_pbysc, h = "H3K4me1", low.quantile = 0.05, cap = 0.999)

drawHisHMscale <- function(r, h, low.quantile = 0.01, cap = 0.999) {
  # r <- t(t(cnt * scale.factor) / colSums(cnt))
  upValue <- quantile(r, cap)
  message("quantile of ", cap, " : ", upValue)
  r[r > upValue] <- upValue
  logCPM <- r |>
    log1p() |>
    scale() |>
    t()
  
  logCPM.plot <- logCPM[, peaks.ds]
  pHM <- setHeatmap(mat = logCPM.plot,
    low.quantile = low.quantile,
    show.legend = TRUE,
    legend.direction = "horizontal")

  withr::with_pdf(
    new = file.path(projd, "figure",
      str_glue("{h}.onATAC.ATACnmf.all.ds.scale.heatmap.pdf")),
    code = {
      print(pHM)
    },
    height = 20,
    width = 20
  )
  return(pHM)
}
drawHisHMscale(r = cpmH3K27ac_pbysc, h = "H3K27ac", low.quantile = 0.01, cap = 0.99)
drawHisHMscale(r = cpmH3K4me1_pbysc, h = "H3K4me1", low.quantile = 0.01, cap = 0.99)

