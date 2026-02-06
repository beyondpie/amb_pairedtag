library(tidyverse)
library(readxl)
library(tmpRpkg)

Sys.setenv("_R_USE_PIPEBIND_" = TRUE)

# * meta
projd <- here::here()
statfnm <- file.path(projd, "meta",
  "allen_enhancer_chromatinstate_for_Songpeng.xlsx")

# * function
statROC <- function(tbl, scMap,
                    manualCheckFilter = NULL,
                    conservFilter = "conserved",
                    predictState = "ChrA") {
  r <- if (is.null(manualCheckFilter)) {
    tbl
  } else {
    tbl[tbl[["manual check"]] %in% manualCheckFilter, ]
  }
  r <- if(conservFilter %in% c("conserved", "not_conserved")) {
    r[r$conservation %in% conservFilter, ]
  } else {
    r
  }
  sc4check <- unique(r$max_subclass_primary)
  roc <- lapply(sc4check, \(a) {
    rows <- r$max_subclass_primary == a
    tcol <- scMap[a, "pt"]
    ncols <- setdiff(scMap$pt, tcol)
    # one column
    tr <- r[rows, tcol]
    # multiple columns
    nr <- r[rows, ncols, drop = F]
    data.frame(
      numtr = nrow(tr),
      numtp = sum(tr == predictState),
      numnr = nrow(nr) * ncol(nr),
      numfn = sum(nr == predictState)
    )
  }) |>
    do.call(rbind, args = _) |>
    x => `rownames<-`(x, sc4check) |>
  setNames(object = _, nm = c("nt", "ntp", "nn", "nfp"))
}

calcROC <- function(roc) {
  data.frame(
    precision = sum(roc$ntp) / (sum(roc$ntp) + sum(roc$nfp)),
    recall = sum(roc$ntp) / sum(roc$nt),
    sensitiviy = sum(roc$ntp) / sum(roc$nt),
    specificy = 1 - sum(roc$nfp) / sum(roc$nn)
  )
}



# * main

# * map subclasses between Allen's and ours
scMap <- data.frame(
  allen = c("L6_IT_Car3", "L5_ET", "Vip",
    "Lamp5", "Lamp5_Lhx6", "Pvalb", "Sst", "Astro", "Oligo", "Endo",
    "L2-3_IT", "Sncg"),
  pt = c("004_L6IT", "022_L5ET", "046_Vip",
    "049_Lamp5", "050_Lamp5_Lhx6", "052_Pvalb", "053_Sst",
    "319_Astro", "327_Oligo", "333_Endo",
    "007_L23IT", "047_Sncg")
) |>
  x => `rownames<-`(x, x$allen)

# * read statfnm content
s <- read_excel(path = statfnm, sheet = "Sheet1")

roc_consrv_ChrA <- statROC(s, scMap, manualCheckFilter = NULL,
  conservFilter = "conserved", predictState = "ChrA")
calroc_consrv_ChrA <- calcROC(roc_consrv_ChrA)

## > calroc_consrv_ChrA
##   precision    recall sensitiviy specificy
##  0.5434783 0.8333333  0.8333333 0.9363636

roc_consrv_ChrO <- statROC(s, scMap, manualCheckFilter = NULL,
  conservFilter = "conserved", predictState = "ChrO")

calroc_consrv_ChrO <- calcROC(roc_consrv_ChrO)

## > calroc_consrv_ChrO
##    precision    recall sensitiviy specificy
##  0.08888889 0.1333333  0.1333333 0.8757576

roc_all_ChrA <- statROC(s, scMap, manualCheckFilter = NULL,
  conservFilter = "none", predictState = "ChrA")
calroc_all_ChrA <- calcROC(roc_all_ChrA)

## > calroc_all_ChrA
##   precision    recall sensitiviy specificy
##  0.5176471 0.6984127  0.6984127 0.9408369

roc_all_ChrO <- statROC(s, scMap, manualCheckFilter = NULL,
  conservFilter = "none", predictState = "ChrO")
calroc_all_ChrO <- calcROC(roc_all_ChrO)

## > calroc_all_ChrO
##   precision    recall sensitiviy specificy
##  0.1782178 0.2857143  0.2857143 0.8802309
